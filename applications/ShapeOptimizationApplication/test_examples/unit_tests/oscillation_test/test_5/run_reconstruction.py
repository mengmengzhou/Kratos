# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division 
from math import sin, cos, pi, fabs
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
# surface definition
# ======================================================================================================================================    
def surf_1(u,v):
    x_1 = sin(pi/4) *  u - cos(pi/4) * v
    # if x_1>0:
    #     return [u,v, -x_1**2]
    # else:
    #     return [u,v, +x_1**2]
    return [u,v,x_1**2]    
# ======================================================================================================================================
# Mapping
# ======================================================================================================================================    

# Create CAD-mapper
linear_solver = SuperLUSolver()
mapper = CADMapper(fe_model_part,cad_geometry,cad_integration_data,linear_solver)

# Set nearest point + shape update to map
# mapper.set_point(id, u, v, updated_x, updated_y, updated_z)
# for u in range(6):
#     for v in range(6):
#         [x,y,z] = surf_1(u, v)
#         mapper.set_point(0, u, v, x, y, z)

u=5
v=0
[x,y,z] = surf_1(u,u)
mapper.set_point(0, u, v, x, y, z)
u=0
v=5
[x,y,z] = surf_1(u,u)
mapper.set_point(0, u, v, x, y, z)
for u in range(6):
    [x,y,z] = surf_1(u,u)
    mapper.set_point(0, u, u, x, y, z)

# u=0
# v=0
# [x,y,z] = surf_1(u,u)
# mapper.set_point(0, u, v, x, y, z)
# u=5
# v=5
# [x,y,z] = surf_1(u,u)
# mapper.set_point(0, u, v, x, y, z)
# for u in range(6):
#     [x,y,z] = surf_1(u,5-u)
#     mapper.set_point(0, u, 5-u, x, y, z)

# # Perform mapping
# mapper.external_map_to_cad_space()

# ======================================================================================================================================
# Writing results
# ======================================================================================================================================

# Output some surface nodes of updated cad geometry
file_to_write = "surface_nodes_of_updated_cad_geometry.txt"
u_resolution = 30
v_resolution = 30
mapper.output_surface_points(file_to_write, u_resolution, v_resolution, -1)

# Output control point update in gid-format 
mapper.output_control_point_displacements()

# Output json file with updated geometry
with open(cad_geometry_output_filename, 'w') as fp:
    json.dump(cad_geometry, fp)
