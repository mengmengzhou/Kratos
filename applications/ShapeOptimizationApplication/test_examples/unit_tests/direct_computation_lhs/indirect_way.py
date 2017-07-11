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
fem_input_filename = "optimized_shell_3418_elements"
cad_geometry_input_filename = "Benchmark_halbkreis_32x16_multipatch_geometry.json" 
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
# Testing nearest neighbour + Newton-Raphson
# ======================================================================================================================================    

# Create CAD-mapper
linear_solver = SuperLUSolver()
mapper = CADMapper(fe_model_part,cad_geometry,cad_integration_data,linear_solver)

# Compute nearest points
u_res = 5
v_res = 5
mapper.compute_nearest_points(u_res,v_res)

#########################################################################
file_to_write = "{0}x{1}.txt".format(u_res, v_res)
mapper.print_nearest_points(file_to_write)
#########################################################################

# Compute a matrix
mapper.compute_a_matrix()

# Set shape update to map
for node in fe_model_part.Nodes:
    shape_update = Vector(3)
    shape_update[0] = node.GetSolutionStepValue(SHAPE_CHANGE_ABSOLUTE_X)
    shape_update[1] = node.GetSolutionStepValue(SHAPE_CHANGE_ABSOLUTE_Y)
    shape_update[2] = node.GetSolutionStepValue(SHAPE_CHANGE_ABSOLUTE_Z)
    node.SetValue(SHAPE_CHANGE_ABSOLUTE,shape_update)

# Perform mapping
mapper.map_to_cad_space_2("indirect_way.txt")
