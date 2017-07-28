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

#if !defined(KRATOS_MORTAR_MAPPER_H_INCLUDED )
#define  KRATOS_MORTAR_MAPPER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "mapper_matrix_based.h"


namespace Kratos
{

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

class MortarMapper : public MapperMatrixBased
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of MortarMapper
    KRATOS_CLASS_POINTER_DEFINITION(MortarMapper);

    ///@}
    ///@name Life Cycle
    ///@{

    MortarMapper(ModelPart& i_model_part_origin, ModelPart& i_model_part_destination,
                          Parameters& rJsonParameters, BuilderAndSolver BuilderAndSolverMdd) : MapperMatrixBased(
                                  i_model_part_origin, i_model_part_destination, rJsonParameters)
    {
        // @Jordi the "Condition_Gauss_Point" can be of the type that Riccardo suggested. Then we can get directly the ShapeFunctionValues that we need
        mpMapperCommunicator->InitializeOrigin(MapperUtilities::Condition_Gauss_Point);
        mpMapperCommunicator->InitializeDestination(MapperUtilities::Condition_Gauss_Point);
        mpMapperCommunicator->Initialize();
        // these three steps correspond to "MapperCommunicator::BuildInterface()"

        // //mpMapperCommunicator->GetBuilderAndMultiplier()->BuildLHS(scheme, modelpart, Mdo); // could be moved to baseclass...?
        // FillMappingMatrix();
        // mpMapperCommunicator->GetBuilderAndSolver()->BuildLHS(scheme, modelpart, mM_dd); // this would also initialize this "BuilderAndSolver" => same as below, should be member of this class?
    }

    /// Destructor.
    virtual ~MortarMapper() { }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void UpdateInterface(Kratos::Flags MappingOptions, double SearchRadius) override
    {

    }

    /* This function maps from Origin to Destination */
    void Map(const Variable<double>& rOriginVariable,
             const Variable<double>& rDestinationVariable,
             Kratos::Flags MappingOptions) override
    {   
        // InterpolateToDestinationMesh(mQ_tmp) // here the Multiplication is done
        
        // // @Jordi this BuilderAndSolver is only needed for Mortar. Can it be a member of this class then?
        // mpMapperCommunicator->GetBuilderAndSolver()->Solve(scheme, modelpart_destination, mM_dd, mQ_d, mQ_tmp);

        // SetNodalValues();
    }

    /* This function maps from Origin to Destination */
    void Map(const Variable< array_1d<double, 3> >& rOriginVariable,
             const Variable< array_1d<double, 3> >& rDestinationVariable,
             Kratos::Flags MappingOptions) override
    {
        
    }

    /* This function maps from Destination to Origin */
    void InverseMap(const Variable<double>& rOriginVariable,
                    const Variable<double>& rDestinationVariable,
                    Kratos::Flags MappingOptions) override
    {

    }

    /* This function maps from Destination to Origin */
    void InverseMap(const Variable< array_1d<double, 3> >& rOriginVariable,
                    const Variable< array_1d<double, 3> >& rDestinationVariable,
                    Kratos::Flags MappingOptions) override
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
    virtual std::string Info() const override
    {
        return "MortarMapper";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MortarMapper";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
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

    void FillMappingMatrix() override
    {
        // // @Jordi I hope this is possible. And I hope we only need to access data in this row.
        // // For NN and NE this is the case, for Mortar I am not sure atm.
        // // local information: row for destination Node


        // mM_do->Reset();


        // // Ask the communicator for the information needed for filling the matrix
        // // NN: Neighbor Indices
        // // NE: Neighbor Indices and corresponding shape function values
        // // Mortar:

        // // @Jordi Here the integration would be done. I would like to follow what Vicente does for Contact
        // // to not reinvent the wheel. He has a utility ("exact_mortar_segmentation_utility.h") which seems to do the
        // // integration, but I want to talk to him when he comes next week.

        // // Unitl Vicente comes next week I will also try to catch up again on the Mortar Theory, esp the triangulaton.

        // // The data that is necessary will be querried from the communicator. Aditya is currently trying to extend the 
        // // communicator for sending arbitrary data, this he needs for his MPC stuff. We could then also use this
        // // Or the global pointers, lets see how it goes

        
        // mpCommunicator->AssembleMappingMatrix(mM_do);
    }

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
    
    // // @Jordi same question as for the base class
    // TSystemMatrixType mM_dd;
    // TSystemVectorType mQ_tmp;

    // strategy

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
    MortarMapper& operator=(MortarMapper const& rOther);

    /// Copy constructor.
    //MortarMapper(MortarMapper const& rOther);

    ///@}

}; // Class MortarMapper

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function


inline std::istream & operator >>(std::istream& rIStream,
                                  MortarMapper& rThis)
{
    return rIStream;
}

/// output stream function

inline std::ostream & operator <<(std::ostream& rOStream,
                                  const MortarMapper& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_MORTAR_MAPPER_H_INCLUDED  defined