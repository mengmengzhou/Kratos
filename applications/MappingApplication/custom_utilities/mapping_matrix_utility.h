//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_MAPPING_MATRIX_UTILITY_H_INCLUDED )
#define  KRATOS_MAPPING_MATRIX_UTILITY_H_INCLUDED

// System includes

// External includes


// Project includes
#include "includes/define.h"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
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

/// Short class definition.
/** Detail class definition.
*/
class MappingMatrixUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MappingMatrixUtility
    KRATOS_CLASS_POINTER_DEFINITION(MappingMatrixUtility);

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MappingMatrixUtility() { }

    /// Destructor.
    virtual ~MappingMatrixUtility() { }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual void AssembleRHSOrigin(const Variable<double>& rVariable,
                                   Kratos::Flags MappingOptions,
                                   std::unordered_map<int, Node<3>*>& rGlobalVectorIndices)
    {
        for (auto const &entry : rGlobalVectorIndices)
        {
            // mQ_o[entry.first] = entry.second.FastGetCurrentSolutionStepValue(rVaiable);
        }
        // Do Global Assembly
    }

    virtual void AssembleRHSDestination() // Needed for Conservative Mapping
    {


    }
    
    virtual void UpdateOrigin() // Needed for Conservative Mapping
    {


    }

    virtual void UpdateDestination(const Variable<double>& rVariable,
                                   Kratos::Flags MappingOptions,
                                   std::unordered_map<int, Node<3>*>& rGlobalVectorIndices)
    {
        for (auto const &entry : rGlobalVectorIndices)
        {
            // entry.second.FastGetCurrentSolutionStepValue(rVariable) = mQ_d[entry.first];
        }
        // Do Global Assembly
    }



    virtual void Multiply()
    {
        // Do Multiplication based on Spaces stuff
        // mQ_d = mM_do * mQ_o // in case of Mortar one has to do some more operation on mQ_d

    }

    virtual void MultiplyConservative()
    {
        // Do Multiplication based on Spaces stuff
        // mQ_o = mM_do^T * mQ_d // in case of Mortar one has to do some more operation on mQ_d

    }


    // Functions only needed to distributed matrices
    virtual void ConstructEpetraGraph()
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual void AssembleMatrixGlobal()
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }


    ///@}
    ///@name Access
    ///@{

    // TODO Use operator overloading for this?
    virtual double GetEntry(const int RowIndex, const int ColumnIndex)
    {   
        // retreive the value from the matrix
        double value = 5.0;
        return value;
    }

    virtual void SetEntry(const int RowIndex, const int ColumnIndex, const double Value)
    {   
        // set the value in the matrix
        return;
    }

    virtual void GetMatrix()
    {
        // interface function if one really wants to have the matrix
    }


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "MappingMatrixUtility" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MappingMatrixUtility";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


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

    // TSystemMatrixType mM_do;
    // TSystemVectorType mQ_o;
    // TSystemVectorType mQ_d;


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
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    MappingMatrixUtility& operator=(MappingMatrixUtility const& rOther);

    //   /// Copy constructor.
    //   MappingMatrixUtility(MappingMatrixUtility const& rOther){}


    ///@}

}; // Class MappingMatrixUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MappingMatrixUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MappingMatrixUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MAPPING_MATRIX_UTILITY_H_INCLUDED  defined