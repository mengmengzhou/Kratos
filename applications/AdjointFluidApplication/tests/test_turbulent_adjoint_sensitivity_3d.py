import os
from KratosMultiphysics import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_MainKratos
import os

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

    def getTimeAveragedFlowMalDistribution(self, filename):
        with open(filename, 'r') as input_file:
            lines = input_file.readlines()
        lines = lines[1:]

        if len(lines) > 1:
            Dt = -(float(lines[1].split()[0]) - float(lines[0].split()[0]))
            T  = -(float(lines[-1].split()[0]) - float(lines[0].split()[0])) + Dt
        else:
            Dt = 1.0
            T = 1.0

        avg_result = 0.0
        for line in lines:
            avg_result = avg_result + float(line.strip().split()[1])

        return avg_result/T

    def ComputeFiniteDifferenceSensitivity( \
                                self, node_ids, step_size, model_part_file_name, file_name):
        sensitivity = []
        # unperturbed drag
        print('Solving for unperturbed flow mal dsitribution...')
        self.solve(model_part_file_name)
        flow_mal_distribution_0 = self.getTimeAveragedFlowMalDistribution(file_name)
        for node_id in node_ids:
            node_sensitivity = []
            coord = self.readNodalCoordinates(node_id, model_part_file_name)
            # Y + h
            perturbed_coord = [coord[0], coord[1]+step_size, coord[2]]
            self.writeNodalCoordinates(node_id, perturbed_coord, model_part_file_name)
            print('Solving for perturbation in Y direction flow mal dsitribution...')
            self.solve(model_part_file_name)
            flow_mal_distribution = self.getTimeAveragedFlowMalDistribution(file_name)
            node_sensitivity.append((flow_mal_distribution - flow_mal_distribution_0) / step_size)
            # Z + h
            perturbed_coord = [coord[0], coord[1], coord[2] + step_size]
            self.writeNodalCoordinates(node_id, perturbed_coord, model_part_file_name)
            print('Solving for perturbation in Z direction flow mal dsitribution...')
            self.solve(model_part_file_name)
            flow_mal_distribution = self.getTimeAveragedFlowMalDistribution(file_name)
            node_sensitivity.append((flow_mal_distribution - flow_mal_distribution_0) / step_size)
            sensitivity.append(node_sensitivity)
            # return mdpa file to unperturbed state
            self.writeNodalCoordinates(node_id, coord, model_part_file_name)
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

    def test_Nozzle_3D(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # solve fluid
            print('Solving primal problem...')
            self.solve('test_turbulent_sensitivity_3d/nozzle_test_3d')
            # solve adjoint
            print('Solving adjoint problem...')
            test = self.createTest('test_turbulent_sensitivity_3d/nozzle_test_3d_adjoint')
            test.Solve()

            node_id = 74

            Sensitivity = [[]]
            Sensitivity[0].append(test.main_model_part.GetNode(node_id).GetSolutionStepValue(SHAPE_SENSITIVITY_Y))
            Sensitivity[0].append(test.main_model_part.GetNode(node_id).GetSolutionStepValue(SHAPE_SENSITIVITY_Z))

            # calculate sensitivity by finite difference
            step_size = 0.00000001
            FDSensitivity = self.ComputeFiniteDifferenceSensitivity(\
                        [node_id], \
                        step_size, \
                        './test_turbulent_sensitivity_3d/nozzle_test_3d', \
                        './test_turbulent_sensitivity_3d/flow_mal_distribution.data')
            self.assertAlmostEqual(Sensitivity[0][0], FDSensitivity[0][0], 5)
            self.assertAlmostEqual(Sensitivity[0][1], FDSensitivity[0][1], 5)
            self.removeFile("./test_turbulent_sensitivity_3d/nozzle_test_3d_0.h5")
            self.removeFile("./test_turbulent_sensitivity_3d/flow_mal_distribution.data")
            self.removeFile("./test_turbulent_sensitivity_3d/nozzle_test_3d.time")
            self.removeFile("./test_turbulent_sensitivity_3d/nozzle_test_3d_adjoint.time")

    def tearDown(self):
        pass

if __name__ == '__main__':
    KratosUnittest.main()
