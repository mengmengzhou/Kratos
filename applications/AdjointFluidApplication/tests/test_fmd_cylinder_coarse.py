import os
from KratosMultiphysics import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_MainKratos
import json
from shutil import copyfile

class ControlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

class TestCase(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def removeFile(self, file_path):
        if os.path.isfile(file_path):
            os.remove(file_path)

    def readNodalCoordinates(self, node_id, model_part_file_name):
        with open(model_part_file_name + '.mdpa', 'r') as model_part_file:
            lines = model_part_file.readlines()
        lines = lines[lines.index('Begin Nodes\n'):lines.index('End Nodes\n')]
        line = lines[node_id] # assumes consecutive node numbering starting with 1
        components = line.split()
        if int(components[0]) != node_id:
            raise RuntimeError('Error parsing file ' + model_part_file_name)
        return [float(components[i]) for i in range(1,4)]

    def writeNodalCoordinates(self, node_id, coords,model_part_file_name):
        with open(model_part_file_name + '.mdpa', 'r') as model_part_file:
            lines = model_part_file.readlines()
        node_lines = lines[lines.index('Begin Nodes\n'):lines.index('End Nodes\n')]
        old_line = node_lines[node_id] # assumes consecutive node numbering starting with 1
        components = old_line.split()
        if int(components[0]) != node_id:
            raise RuntimeError('Error parsing file ' + model_part_file_name)
        new_line = '{:5d}'.format(node_id) + ' ' \
             + '{:19.10f}'.format(coords[0]) + ' ' \
             + '{:19.10f}'.format(coords[1]) + ' ' \
             + '{:19.10f}'.format(coords[2]) + '\n'
        lines[lines.index(old_line)] = new_line
        with open(model_part_file_name + '.mdpa', 'w') as model_part_file:
            model_part_file.writelines(lines)

    def getTimeAveragedFlowMalDistribution(self, filename, model_part_file_name):
        with open(filename, 'r') as input_file:
            lines = input_file.readlines()
        lines = lines[1:]

        if len(lines) > 1:
            Dt = (float(lines[1].split()[0]) - float(lines[0].split()[0]))
            T  = (float(lines[-1].split()[0]) - float(lines[0].split()[0])) + Dt
        else:
            Dt = 1.0
            T = 1.0

        avg_result = 0.0
        
        _remove_iterations = 0
        _itr_counter = 0
        for line in lines:
            _itr_counter += 1
            if _itr_counter < _remove_iterations:
                continue
            avg_result = avg_result + float(line.strip().split()[1])

        return avg_result * (Dt / T)

    def ComputeFiniteDifferenceSensitivity( \
                                self, node_ids, step_size, model_part_file_name, file_name):
        sensitivity = []
        # unperturbed flow mal distribution
        copyfile("%s.mdpa.org" % model_part_file_name, "%s.mdpa" % model_part_file_name)
        self.solve(model_part_file_name)
        flow_mal_distribution_0 = self.getTimeAveragedFlowMalDistribution(file_name, model_part_file_name)
        count = 0
        for node_id in node_ids:
            node_sensitivity = []
            copyfile("%s.mdpa.org" % model_part_file_name, "%s.mdpa" % model_part_file_name)
            coord = self.readNodalCoordinates(node_id, model_part_file_name)

            # X + h
            # perturbed_coord = [coord[0] + step_size, coord[1], coord[2]]
            # self.writeNodalCoordinates(node_id, perturbed_coord, model_part_file_name)
            # self.solve(model_part_file_name)
            # flow_mal_distribution = self.getTimeAveragedFlowMalDistribution(file_name)
            # node_sensitivity.append((flow_mal_distribution - flow_mal_distribution_0) / step_size)


            # Y + h
            perturbed_coord = [coord[0], coord[1]+step_size, coord[2]]
            self.writeNodalCoordinates(node_id, perturbed_coord, model_part_file_name)
            self.solve(model_part_file_name)
            flow_mal_distribution = self.getTimeAveragedFlowMalDistribution(file_name, model_part_file_name)
            node_sensitivity.append((flow_mal_distribution - flow_mal_distribution_0) / step_size)

            # Z + h
            perturbed_coord = [coord[0], coord[1], coord[2] + step_size]
            self.writeNodalCoordinates(node_id, perturbed_coord, model_part_file_name)
            self.solve(model_part_file_name)
            flow_mal_distribution = self.getTimeAveragedFlowMalDistribution(file_name, model_part_file_name)
            node_sensitivity.append((flow_mal_distribution - flow_mal_distribution_0) / step_size)
            sensitivity.append(node_sensitivity)
            self.removeFile("%s.mdpa" % model_part_file_name)


            count = count + 1
        return sensitivity

    def createTest(self, parameter_file_name):
        with open(parameter_file_name + '_parameters.json', 'r') as parameter_file:
            project_parameters = Parameters(parameter_file.read())
            parameter_file.close()
        test = test_MainKratos.MainKratos(project_parameters)
        return test

    def solve(self, parameter_file_name):
        test = self.createTest(parameter_file_name)
        test.Solve()

    def test_cylinder(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            model_part_file_name = 'test_fmd_cylinder_coarse/test_fmd_cylinder_coarse'

            # solve fluid
            copyfile("%s.mdpa.org" % model_part_file_name, "%s.mdpa" % model_part_file_name)
            self.solve(model_part_file_name)
            self.removeFile("%s.mdpa" % model_part_file_name)

            # solve adjoint
            copyfile("%s.mdpa.org" % model_part_file_name, "%s_adjoint.mdpa" % model_part_file_name)
            test = self.createTest('test_fmd_cylinder_coarse/test_fmd_cylinder_coarse_adjoint')
            test.Solve()
            self.removeFile("%s_adjoint.mdpa" % model_part_file_name)

            node_ids = [8]

            _adjoint_sensitivity = []
            for i in range(0, len(node_ids)):
                node_id = node_ids[i]
                node_sensitivity = []
                # node_sensitivity.append( \
                #         test.main_model_part.GetNode(node_id).GetSolutionStepValue(\
                #                             SHAPE_SENSITIVITY_X))
                node_sensitivity.append( \
                        test.main_model_part.GetNode(node_id).GetSolutionStepValue(\
                                            SHAPE_SENSITIVITY_Y))
                node_sensitivity.append( \
                        test.main_model_part.GetNode(node_id).GetSolutionStepValue(\
                                            SHAPE_SENSITIVITY_Z))
                _adjoint_sensitivity.append(node_sensitivity)

            # calculate sensitivity by finite difference
            step_size = 0.000000001
            _finite_difference_sensitivity = self.ComputeFiniteDifferenceSensitivity(\
                        node_ids, \
                        step_size, \
                        './test_fmd_cylinder_coarse/test_fmd_cylinder_coarse', \
                        './test_fmd_cylinder_coarse/flow_mal_distribution.data')

            print("Sensitivities calculated from finite difference method: ")
            print(_finite_difference_sensitivity)
            print("Sensitivities calculated from adjoint method: ")
            print(_adjoint_sensitivity)

            for i in range(0, len(node_ids)):
                self.assertAlmostEqual(_adjoint_sensitivity[i][0], _finite_difference_sensitivity[i][0], 5)
                self.assertAlmostEqual(_adjoint_sensitivity[i][1], _finite_difference_sensitivity[i][1], 5)

            self.removeFile("./test_fmd_cylinder_coarse/test_fmd_cylinder_coarse_0.h5")
            self.removeFile("./test_fmd_cylinder_coarse/flow_mal_distribution.data")
            self.removeFile("./test_fmd_cylinder_coarse/test_fmd_cylinder_coarse.time")
            self.removeFile("./test_fmd_cylinder_coarse/test_fmd_cylinder_coarse_adjoint.time")

    def tearDown(self):
        pass

if __name__ == '__main__':
    KratosUnittest.main()
