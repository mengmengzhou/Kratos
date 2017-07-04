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
# Mapping
# ======================================================================================================================================    

# Create CAD-mapper
linear_solver = SuperLUSolver()
# DiagPrecond = DiagonalPreconditioner()
# linear_solver =  BICGSTABSolver(1e-9, 5000, DiagPrecond)
# linear_solver = AMGCLSolver(AMGCLSmoother.GAUSS_SEIDEL, AMGCLIterativeSolverType.BICGSTAB, 1e-9, 300, 2, 10)
mapper = CADMapper(fe_model_part,cad_geometry,cad_integration_data,linear_solver)

# Output some surface nodes of cad geometry
file_to_write = "surface_nodes_of_cad_geometry.txt"
u_resolution = 100
v_resolution = 100
mapper.output_surface_points(file_to_write, u_resolution, v_resolution, -1)

# Compute mapping matrix
u_resolution = 500
v_resolution = 500
mapper.compute_mapping_matrix(u_resolution,v_resolution)

# Set shape update to map
for node in fe_model_part.Nodes:
    shape_update = Vector(3)
    shape_update[0] = node.GetSolutionStepValue(SHAPE_CHANGE_ABSOLUTE_X)
    shape_update[1] = node.GetSolutionStepValue(SHAPE_CHANGE_ABSOLUTE_Y)
    shape_update[2] = node.GetSolutionStepValue(SHAPE_CHANGE_ABSOLUTE_Z)
    node.SetValue(SHAPE_CHANGE_ABSOLUTE,shape_update)

# Apply boundary conditions
penalty_factor_displacement_coupling = 1e4
penalty_factor_rotation_coupling = 1e4
penalty_factor_dirichlet_condition = 1e5
edges_with_specific_dirichlet_conditions = []
# edges_with_specific_dirichlet_conditions = [ [1001,[False,True,False]], [1003,[False,True,False]] ]
edges_with_enforced_tangent_continuity = []
# edges_with_enforced_tangent_continuity = [ [1004, 1e3] ]
mapper.apply_boundary_conditions( penalty_factor_displacement_coupling, 
                                  penalty_factor_rotation_coupling, 
                                  penalty_factor_dirichlet_condition,
                                  edges_with_specific_dirichlet_conditions,
                                  edges_with_enforced_tangent_continuity )

# Perform mapping
mapper.map_to_cad_space()

# Output some surface nodes of updated cad geometry
file_to_write = "surface_nodes_of_updated_cad_geometry.txt"
u_resolution = 250
v_resolution = 250
mapper.output_surface_points(file_to_write, u_resolution, v_resolution, -1)

# Output control point update in gid-format 
mapper.output_control_point_displacements()

# Output json file with updated geometry
with open(cad_geometry_output_filename, 'w') as fp:
    json.dump(cad_geometry, fp)