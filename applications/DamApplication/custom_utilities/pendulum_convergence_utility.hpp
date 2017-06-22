//
//   Project Name:        KratosDamApplication    $
//   Last modified by:    $Author: Lorenzo Gracia $
//   Date:                $Date:        June 2017 $
//   Revision:            $Revision:          0.0 $
//

#if !defined(KRATOS_PENDULUM_CONVERGENCE_UTILITY )
#define  KRATOS_PENDULUM_CONVERGENCE_UTILITY

#include <cmath>

#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"
#include "utilities/openmp_utils.h"
#include "utilities/math_utils.h"

#include "dam_application_variables.h"

namespace Kratos
{
    
class PendulumConvergenceUtility
{
    
public:
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    PendulumConvergenceUtility(ModelPart& model_part
                                ) : mr_model_part(model_part)
    {

    }

    ///------------------------------------------------------------------------------------
    
    /// Destructor
    ~PendulumConvergenceUtility() {}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Vector CheckConvergence(double tolerance, double reference)
    {
        KRATOS_TRY;
        Vector Result = ZeroVector(2);   
        Result.resize(2); 
        double value_to_check = mr_model_part.GetNode(173).GetSolutionStepValue(DISPLACEMENT_Y);
        double dif = value_to_check-reference;

        Result[1] = dif;

        if(fabs(dif) < tolerance)
            Result[0] = 1.0;
    
        return Result;

        KRATOS_CATCH("");
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    bool CheckGlobalConvergence(unsigned int number_of_pendulums, double tolerance, double reference)
    {
        KRATOS_TRY;
        bool Result = false;   
        Vector ConvergencePendulumsVector = ZeroVector (number_of_pendulums);
        ConvergencePendulumsVector.resize(number_of_pendulums);

        // Here a for to check all control points
        for(unsigned int i = 0; i<number_of_pendulums; i++)
        {
            double value_to_check = mr_model_part.GetNode(173).GetSolutionStepValue(DISPLACEMENT_Y);
            double dif = value_to_check-reference;
            if(fabs(dif) < tolerance)
                ConvergencePendulumsVector[i]= 1.0;
        }
        
        for(unsigned int j = 0; j<number_of_pendulums; j++)
        {
            if (ConvergencePendulumsVector[j] > 0.5)
                Result = true;
            else 
                break;
        }
        return Result;

        KRATOS_CATCH("");
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
protected:

    /// Member Variables
    ModelPart& mr_model_part;


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

};//Class

} /* namespace Kratos.*/

#endif /* KRATOS_PENDULUM_CONVERGENCE_UTILITY defined */
