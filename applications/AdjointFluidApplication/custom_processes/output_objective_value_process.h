//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//                   Suneth Warnakulasuriya
//

#if !defined(KRATOS_OUTPUT_OBJECTIVE_VALUE_PROCESS_H_INCLUDED )
#define KRATOS_OUTPUT_OBJECTIVE_VALUE_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <iomanip>

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "includes/model_part.h"

// Application includes
#include "custom_utilities/objective_function.h"
#include "custom_utilities/flow_mal_distribution_objective_function.h"

namespace Kratos
{

///@addtogroup AdjointFluidApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Process to output primal solution (e.g., velocity and pressure) needed for adjoint solution.
/**
 * Requires mesh topology and node numbering to remain constant.
 */
template<unsigned int TDim>
class OutputObjectiveValueProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(OutputObjectiveValueProcess);

    typedef ModelPart::NodeType NodeType;

    typedef hsize_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    OutputObjectiveValueProcess(ModelPart& rModelPart, Parameters& rParameters)
    : Process(), mrModelPart(rModelPart)
    {
        KRATOS_TRY

        if (rParameters["objective_type"].GetString().compare("flow_mal_distribution") == 0)
        {
                mObjectiveFunction = static_cast<ObjectiveFunction*>(new FlowMalDistributionObjectiveFunction<TDim>(rParameters));
        }
        else
        {
            KRATOS_THROW_ERROR(std::invalid_argument, "objective type not supported ", rParameters["objective_type"].GetString());
        }

        KRATOS_CATCH("")
    }

    /// Destructor.
    virtual ~OutputObjectiveValueProcess() {
        delete mObjectiveFunction;
    }

    ///@}
    ///@name Operators
    ///@{

    /// This operator simply calls the Execute method.
    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    virtual void Execute() {}

    virtual void ExecuteInitialize()
    {
        KRATOS_TRY

        mObjectiveFunction->Initialize(mrModelPart);

        KRATOS_CATCH("")
    }

    virtual void ExecuteBeforeSolutionLoop()
    {
    }

    virtual void ExecuteInitializeSolutionStep()
    {
        KRATOS_TRY

        mObjectiveFunction->InitializeSolutionStep(mrModelPart);

        KRATOS_CATCH("")
    }

    virtual void ExecuteFinalizeSolutionStep()
    {
        KRATOS_TRY
        
        mObjectiveFunction->FinalizeSolutionStep(mrModelPart);
        
        KRATOS_CATCH("")
    }

    virtual void ExecuteBeforeOutputStep()
    {
    }

    virtual void ExecuteAfterOutputStep()
    {
    }

    virtual void ExecuteFinalize()
    {
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "OutputObjectiveValueProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << this->Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    ObjectiveFunction* mObjectiveFunction;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
}; /* Class OutputObjectiveValueProcess */

///@}

///@name Input and output
///@{

/// input stream function
template<unsigned int TDim>
inline std::istream& operator >> (std::istream& rIStream,
                                  OutputObjectiveValueProcess<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template<unsigned int TDim>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const OutputObjectiveValueProcess<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Adjoint Fluid Application group

}  /* namespace Kratos */

#endif /* KRATOS_OUTPUT_OBJECTIVE_VALUE_H_INCLUDED defined */
