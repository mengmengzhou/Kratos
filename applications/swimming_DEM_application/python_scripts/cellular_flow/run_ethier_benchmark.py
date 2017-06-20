import KratosSwimmingDEM as script

import ethier_benchmark_algorithm
varying_parameters = dict()

for size in {0.2}:
    varying_parameters['size_parameter'] = size
    for derivatives_type in range(1, 8):
        varying_parameters['material_acceleration_calculation_type'] = derivatives_type
        varying_parameters['laplacian_calculation_type'] = derivatives_type
        test = script.Solution(ethier_benchmark_algorithm, varying_parameters)
        try:
            test.Run()
        except:
            pass
