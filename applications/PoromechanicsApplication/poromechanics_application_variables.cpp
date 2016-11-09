//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#include "poromechanics_application_variables.h"

namespace Kratos
{
//Create Variables 
KRATOS_CREATE_VARIABLE( double, NEWMARK_COEFFICIENT_U )
KRATOS_CREATE_VARIABLE( double, NEWMARK_COEFFICIENT_P )

KRATOS_CREATE_VARIABLE( double, DT_WATER_PRESSURE )
KRATOS_CREATE_VARIABLE( double, NORMAL_FLUID_FLUX )

KRATOS_CREATE_VARIABLE( double, DENSITY_SOLID )
KRATOS_CREATE_VARIABLE( double, BULK_MODULUS_SOLID )
KRATOS_CREATE_VARIABLE( double, BULK_MODULUS_FLUID )
KRATOS_CREATE_VARIABLE( double, PERMEABILITY_XX )
KRATOS_CREATE_VARIABLE( double, PERMEABILITY_YY )
KRATOS_CREATE_VARIABLE( double, PERMEABILITY_ZZ )
KRATOS_CREATE_VARIABLE( double, PERMEABILITY_XY )
KRATOS_CREATE_VARIABLE( double, PERMEABILITY_YZ )
KRATOS_CREATE_VARIABLE( double, PERMEABILITY_ZX )

KRATOS_CREATE_VARIABLE( double, MINIMUM_JOINT_WIDTH )
KRATOS_CREATE_VARIABLE( double, TRANSVERSAL_PERMEABILITY )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( FLUID_FLUX_VECTOR )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( LOCAL_FLUID_FLUX_VECTOR )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( LOCAL_STRESS_VECTOR )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( LOCAL_RELATIVE_DISPLACEMENT_VECTOR )
KRATOS_CREATE_VARIABLE( Matrix, PERMEABILITY_MATRIX )
KRATOS_CREATE_VARIABLE( Matrix, LOCAL_PERMEABILITY_MATRIX )

KRATOS_CREATE_VARIABLE( double, CRITICAL_DISPLACEMENT )

KRATOS_CREATE_VARIABLE(bool, IS_CONVERGED)

KRATOS_CREATE_VARIABLE( Matrix, TOTAL_STRESS_TENSOR )

KRATOS_CREATE_VARIABLE( double, STATE_VARIABLE )
KRATOS_CREATE_VARIABLE( double, ARC_LENGTH_LAMBDA )
KRATOS_CREATE_VARIABLE( double, ARC_LENGTH_RADIUS_FACTOR )

KRATOS_CREATE_VARIABLE( double, TIME_UNIT_CONVERTER )

KRATOS_CREATE_VARIABLE( double, LOCAL_EQUIVALENT_STRAIN )
KRATOS_CREATE_VARIABLE( double, NONLOCAL_EQUIVALENT_STRAIN )
}
