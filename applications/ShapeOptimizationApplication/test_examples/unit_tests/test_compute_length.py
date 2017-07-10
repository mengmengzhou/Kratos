# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division 

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.ShapeOptimizationApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# Additional imports
import json as json

# ======================================================================================================================================
# Parameters
# ======================================================================================================================================

# Input parameters
fem_input_filename = "test_1"
cad_geometry_input_filename = "test_1_geometry.json" 
cad_integration_input_filename = "Benchmark_halbkreis_32x16_multipatch_integration_data.json" 

# Output parameters
cad_geometry_output_filename = cad_geometry_input_filename
cad_geometry_output_filename = cad_geometry_output_filename.replace(".json","_updated.json")

# ======================================================================================================================================
# Reading
# ======================================================================================================================================

# Read the FE model
fe_model_part = ModelPart("name_of_empty_mdpa")
fe_model_part.AddNodalSolutionStepVariable(SHAPE_CHANGE_ABSOLUTE)
model_part_io = ModelPartIO(fem_input_filename)
model_part_io.ReadModelPart(fe_model_part)

# Read CAD data
cad_geometry = {}
with open(cad_geometry_input_filename) as cad_data1:
    cad_geometry = json.load(cad_data1)

cad_integration_data = {}
with open(cad_integration_input_filename) as cad_data2:
    cad_integration_data = json.load(cad_data2)    

# ======================================================================================================================================
# Mapping
# ======================================================================================================================================    
# Create CAD-mapper
linear_solver = SuperLUSolver()
# DiagPrecond = DiagonalPreconditioner()
# linear_solver =  BICGSTABSolver(1e-9, 5000, DiagPrecond)
# linear_solver = AMGCLSolver(AMGCLSmoother.GAUSS_SEIDEL, AMGCLIterativeSolverType.BICGSTAB, 1e-9, 300, 2, 10)
mapper = CADMapper(fe_model_part,cad_geometry,cad_integration_data,linear_solver)

[u_a, v_a] = [0,0]
[u_b, v_b] = [1,0]
length = mapper.compute_real_length(0, u_a, v_a, u_b, v_b)
correct_length = 1
relative_error = abs(100*(length-correct_length)/correct_length)
print("length:", length, "\trelative error: {0:0.2f}%".format(relative_error))


[u_a, v_a] = [0,0]
[u_b, v_b] = [0.5,0]
length = mapper.compute_real_length(0, u_a, v_a, u_b, v_b)
correct_length = 0.5
print("length:", length, "\trelative error: {0:0.2f}%".format(relative_error))

[u_a, v_a] = [0,0]
[u_b, v_b] = [1,1]
length = mapper.compute_real_length(0, u_a, v_a, u_b, v_b)
correct_length = 1.414213562
print("length:", length, "\trelative error: {0:0.2f}%".format(relative_error))
