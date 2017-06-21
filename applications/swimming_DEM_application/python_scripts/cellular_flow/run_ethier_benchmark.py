import KratosSwimmingDEM as script
import sys
import ethier_benchmark_algorithm
import ProjectParameters as pp
import DEM_explicit_solver_var as DEM_parameters
varying_parameters = dict()

irregular_mesh_sizes = {0.2}
regular_mesh_n_points = {}
derivatives_types = range(1, 3)
combinations_that_failed = []
for size in irregular_mesh_sizes.union(regular_mesh_n_points):
    varying_parameters['size_parameter'] = size
    for derivatives_type in derivatives_types:
        varying_parameters['material_acceleration_calculation_type'] = derivatives_type
        varying_parameters['laplacian_calculation_type'] = derivatives_type
        with ethier_benchmark_algorithm.Algorithm() as algorithm:
            test = script.Solution(algorithm, varying_parameters)
            try:
                test.Run()
            except:
                error = sys.exc_info()
                combinations_that_failed.append({'size':size, 'type':derivatives_type, 'error':error})

print()
print('****************************************')

if len(combinations_that_failed):
    print('The following combinations produced an error:')
    print()
    for combination in combinations_that_failed:
        print(combination)
else:
    print('All combinations run without errors')
print('****************************************')
