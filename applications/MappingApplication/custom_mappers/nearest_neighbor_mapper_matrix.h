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

#if !defined(KRATOS_NEAREST_NEIGHBOR_MAPPER_MATRIX_H_INCLUDED )
#define  KRATOS_NEAREST_NEIGHBOR_MAPPER_MATRIX_H_INCLUDED

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

class NearestNeighborMapperMatrix : public MapperMatrixBased
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of NearestNeighborMapperMatrix
    KRATOS_CLASS_POINTER_DEFINITION(NearestNeighborMapperMatrix);

    ///@}
    ///@name Life Cycle
    ///@{

    NearestNeighborMapperMatrix(ModelPart& i_model_part_origin, ModelPart& i_model_part_destination,
                          Parameters& rJsonParameters) : MapperMatrixBased(
                                  i_model_part_origin, i_model_part_destination, rJsonParameters)
    {
        mpMapperCommunicator->InitializeOrigin(MapperUtilities::Node_Coords);
        mpMapperCommunicator->InitializeDestination(MapperUtilities::Node_Coords);
        mpMapperCommunicator->Initialize();

        FillMappingMatrix();

        //mpMapperCommunicator->GetBuilderAndMultiplier()->BuildLHS(scheme, modelpart, Mdo); // could be moved to baseclass...?
    }

    /// Destructor.
    virtual ~NearestNeighborMapperMatrix() { }

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
        InterpolateToDestinationMesh(q_d)
        SetNodalValues();
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
        return "NearestNeighborMapperMatrix";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NearestNeighborMapperMatrix";
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

    void FillMappingMatrix() override
    {
        // local information: row for destination Node


        mM_do->Reset();


        // Ask the communicator for the information needed for filling the matrix
        // NN: Neighbor Indices
        // NE: Neighbor Indices and corresponding shape function values
        // Mortar:

        
        std::vector<int> neighbor_IDs = mpCommunicator->GetNeighborIDs();
        int i = 0;

        for (auto &local_node : rModelPart.GetCommunicator().LocalMesh().Nodes())
        {
            mM_do(local_node.Id() , neighbor_IDs[i]) = 1.0f; // Fill the local part of the matrix
            ++i;
        }

        mpCommunicator->AssembleMappingMatrix(mM_do);
    }

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
    NearestNeighborMapperMatrix& operator=(NearestNeighborMapperMatrix const& rOther);

    /// Copy constructor.
    //NearestNeighborMapperMatrix(NearestNeighborMapperMatrix const& rOther);

    ///@}

}; // Class NearestNeighborMapperMatrix

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function


inline std::istream & operator >>(std::istream& rIStream,
                                  NearestNeighborMapperMatrix& rThis)
{
    return rIStream;
}

/// output stream function

inline std::ostream & operator <<(std::ostream& rOStream,
                                  const NearestNeighborMapperMatrix& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_NEAREST_NEIGHBOR_MAPPER_MATRIX_H_INCLUDED  defined