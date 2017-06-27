//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//                   Suneth Warnakulasuriya
//

#if !defined(KRATOS_FLOW_MAL_DISTRIBUTION_OBJECTIVE_FUNCTION)
#define KRATOS_FLOW_MAL_DISTRIBUTION_OBJECTIVE_FUNCTION

// System includes
#include <vector>
#include <string>
#include <fstream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/process_info.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"
#include "utilities/normal_calculation_utils.h"

// Application includes
#include "custom_utilities/objective_function.h"

namespace Kratos
{
///@addtogroup AdjointFluidApplication
///@{

///@name Kratos Classes
///@{

/// An objective function for flow mal distribution.
template <unsigned int TDim>
class FlowMalDistributionObjectiveFunction : public ObjectiveFunction
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(FlowMalDistributionObjectiveFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    FlowMalDistributionObjectiveFunction(Parameters& rParameters)
    {
        KRATOS_TRY

        Parameters DefaultParamsSurfaceNormal(R"(
        {
            "objective_type": "flow_mal_distribution",
            "surface_model_part_name": "PLEASE_SPECIFY_MODEL_PART",
            "flow_mal_distribution_direction": "surface_normal",
            "output_file": "flow_mal_distribution"
        })");

        Parameters DefaultParamsCustomDirection(R"(
        {
            "objective_type": "flow_mal_distribution",
            "surface_model_part_name": "PLEASE_SPECIFY_MODEL_PART",
            "flow_mal_distribution_direction": [1.0, 0.0, 0.0],
            "output_file": "flow_mal_distribution"
        })");

        mSurfaceModelPartName = rParameters["surface_model_part_name"].GetString();
        mOutputFilename = rParameters["output_file"].GetString();

        if (rParameters["flow_mal_distribution_direction"].IsArray() == false)
        {
            rParameters.ValidateAndAssignDefaults(DefaultParamsSurfaceNormal);
            if (rParameters["flow_mal_distribution_direction"].GetString().compare("surface_normal")==0)
                mMalDistributionDirection = MAL_DISTRIBUTION_DIRECTION_SURFACE_NORMAL;
            else
                KRATOS_THROW_ERROR(std::runtime_error,
                               "flow_mal_distribution_direction only accepts as text inputs \"surface_normal\" or \"absolute_value\":",
                               rParameters.PrettyPrintJsonString())
        }
        else 
        {
            rParameters.ValidateAndAssignDefaults(DefaultParamsCustomDirection);
            if (rParameters["flow_mal_distribution_direction"].size() != 3)
            {
                KRATOS_THROW_ERROR(std::runtime_error,
                                   "flow_mal_distribution_direction vector does not have size 3:",
                                   rParameters.PrettyPrintJsonString())
            }
            else
            {
                mMalDistributionDirection = MAL_DISTRIBUTION_DIRECTION_CUSTOM_DIRECTION;
                for (unsigned int d = 0; d < TDim; ++d)
                    mDirection[d] = rParameters["flow_mal_distribution_direction"][d].GetDouble();

                if (std::abs(norm_2(mDirection) - 1.0) > 1e-3)
                {
                    const double magnitude = norm_2(mDirection);
                    if (magnitude == 0.0)
                        KRATOS_THROW_ERROR(std::runtime_error,
                                        "flow_mal_distribution_direction is not properly defined.",
                                        "")

                    std::cout
                        << "WARNING: non unit vector detected in \"flow_mal_distribution_direction\": "
                        << rParameters.PrettyPrintJsonString() << std::endl;
                    std::cout << "normalizing \"flow_mal_distribution_direction\"..." << std::endl;

                    for (unsigned int d = 0; d < TDim; d++)
                        mDirection[d] /= magnitude;
                }
            }
        }        

        KRATOS_CATCH("")
    }

    /// Destructor.
    virtual ~FlowMalDistributionObjectiveFunction()
    {
        if (mOutputFileOpenend)
            mOutputFileStream.close();
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    virtual void Initialize(ModelPart& rModelPart)
    {
        KRATOS_TRY

        if (rModelPart.HasSubModelPart(mSurfaceModelPartName) == false)
        {
            KRATOS_THROW_ERROR(
                std::runtime_error,
                "invalid parameters \"surface_model_part_name\": ",
                mSurfaceModelPartName)
        }


        // initialize the variables to zero.
        #pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodesBegin, NodesEnd);

            for (auto it = NodesBegin; it != NodesEnd; ++it)
                it->Set(STRUCTURE, false);
        }
        ModelPart& rSurfaceModelPart = rModelPart.GetSubModelPart(mSurfaceModelPartName);
        // go through all the conditions to identify repeating nodes
        for(auto it = rSurfaceModelPart.ConditionsBegin();
                it != rSurfaceModelPart.ConditionsEnd();
                it++)
        {
            Element::GeometryType& rNodes = it->GetGeometry();
            for(unsigned int in = 0; in<rNodes.size(); in++)
                mNodeSharingCounter[rNodes[in].Id()]++;
        }

        // mark objective surface
        for (auto it = rSurfaceModelPart.NodesBegin();
             it != rSurfaceModelPart.NodesEnd();
             ++it)
            it->Set(STRUCTURE, true);
        
        mObjectiveValue = Calculate(rModelPart);

        mOutputFileStream.open(mOutputFilename + ".data");
        mOutputFileStream<<"#time       FLOW_MAL_DISTRIBUTION"<<std::endl;
        mOutputFileOpenend = true;

        NormalCalculationUtils normal_calculation_obj = NormalCalculationUtils();
        const unsigned int domain_size = rModelPart.GetProcessInfo().GetValue(DOMAIN_SIZE);
        normal_calculation_obj.CalculateOnSimplex(rSurfaceModelPart, domain_size);
        CalculateAverageObjectiveParameters(rModelPart);        

        KRATOS_CATCH("")
    }

    virtual void InitializeSolutionStep(ModelPart& rModelPart)
    {
        KRATOS_TRY

        // allocate auxiliary memory. this is done here instead of Initialize()
        // in case of restart.

        KRATOS_CATCH("")
    }

    virtual void CalculateAdjointVelocityContribution(const Element& rElem,
                                                      const Matrix& rAdjointMatrix,
                                                      Vector& rRHSContribution,
                                                      ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        //Question: about the nodes being accounted by this element and the adjacent element as well??

        if (rRHSContribution.size() != rAdjointMatrix.size1())
            rRHSContribution.resize(rAdjointMatrix.size1(), false);

        const unsigned int NumNodes = rElem.GetGeometry().PointsNumber();

        double inv_coeff = 1.0/(mAverageSurfaceArea*(mN-1)*mObjectiveValue);

        if (mMalDistributionDirection == MAL_DISTRIBUTION_DIRECTION_SURFACE_NORMAL)
        {
            unsigned int LocalIndex = 0;
            for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
            {
                if (rElem.GetGeometry()[iNode].Is(STRUCTURE))
                {
                    const array_1d<double,3>& normal = rElem.GetGeometry()[iNode].FastGetSolutionStepValue(NORMAL, 0);
                    const array_1d<double,3>& velocity = rElem.GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, 0);

                    double inv_numberOfAdjacentElements = 1.0/mNodeSharingCounter[rElem.GetGeometry()[iNode].Id()];

                    for (unsigned int d = 0; d < TDim; d++)
                        rRHSContribution[LocalIndex++] = inv_coeff * inv_numberOfAdjacentElements *
                                                        normal[d] * (
                                                                normal[d]*
                                                                velocity[d] / 
                                                                mAverageSurfaceArea - mAverageVelocity
                                                                );
                }
                else
                {
                    for (unsigned int d = 0; d < TDim; d++)
                        rRHSContribution[LocalIndex++] = 0.0;
                }

                //Adding pressure DOF                        
                rRHSContribution[LocalIndex++] = 0.0;
            }
        } 
        else if (mMalDistributionDirection == MAL_DISTRIBUTION_DIRECTION_CUSTOM_DIRECTION)
        {
            unsigned int LocalIndex = 0;
            for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
            {
                if (rElem.GetGeometry()[iNode].Is(STRUCTURE))
                {
                    const array_1d<double,3>& normal = rElem.GetGeometry()[iNode].FastGetSolutionStepValue(NORMAL, 0);
                    const array_1d<double,3>& velocity = rElem.GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, 0);
                    double magnitude = norm_2(normal);

                    double inv_numberOfAdjacentElements = 1.0/mNodeSharingCounter[rElem.GetGeometry()[iNode].Id()];

                    for (unsigned int d = 0; d < TDim; d++)
                        rRHSContribution[LocalIndex++] = inv_coeff * inv_numberOfAdjacentElements *
                                                        mDirection[d]*
                                                        magnitude*(
                                                                        mDirection[d] *
                                                                        magnitude *
                                                                        velocity[d] / 
                                                                        mAverageSurfaceArea - mAverageVelocity
                                                                    );
                }
                else
                {
                    for (unsigned int d = 0; d < TDim; d++)
                        rRHSContribution[LocalIndex++] = 0.0;
                }

                //Adding pressure DOF                        
                rRHSContribution[LocalIndex++] = 0.0;
            }
        }

        KRATOS_CATCH("")
    }

    virtual void CalculateAdjointAccelerationContribution(const Element& rElem,
                                                          const Matrix& rAdjointMatrix,
                                                          Vector& rRHSContribution,
                                                          ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rRHSContribution.size() != rAdjointMatrix.size1())
            rRHSContribution.resize(rAdjointMatrix.size1(), false);

        for (unsigned int k = 0; k < rRHSContribution.size(); ++k)
            rRHSContribution[k] = 0.0;

        KRATOS_CATCH("")
    }

    virtual void CalculateSensitivityContribution(const Element& rElem,
                                                  const Matrix& rDerivativesMatrix,
                                                  Vector& rRHSContribution,
                                                  ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rRHSContribution.size() != rDerivativesMatrix.size1())
            rRHSContribution.resize(rDerivativesMatrix.size1(), false);

        for (unsigned int k = 0; k < rRHSContribution.size(); ++k)
            rRHSContribution[k] = 0.0;

        KRATOS_CATCH("")
    }

    virtual double Calculate(ModelPart& rModelPart)
    {
        double result = 0.0;

        CalculateAverageObjectiveParameters(rModelPart);

        ModelPart& rSurfaceModelPart = rModelPart.GetSubModelPart(mSurfaceModelPartName);

        if (mMalDistributionDirection == MAL_DISTRIBUTION_DIRECTION_SURFACE_NORMAL)
        {
            #pragma omp parallel reduction(+:result)
            {
                ModelPart::NodeIterator NodesBegin;
                ModelPart::NodeIterator NodesEnd;
                OpenMPUtils::PartitionedIterators(rSurfaceModelPart.Nodes(), NodesBegin, NodesEnd);

                for (auto it = NodesBegin; it != NodesEnd; ++it)
                {
                    const array_1d<double,3>& normal = it->FastGetSolutionStepValue(NORMAL, 0);
                    const array_1d<double,3>& velocity = it->FastGetSolutionStepValue(VELOCITY, 0);

                    result += pow(
                                    (
                                        normal[0]*velocity[0] + 
                                        normal[1]*velocity[1] + 
                                        normal[2]*velocity[2]
                                    ) / mAverageSurfaceArea - mAverageVelocity,
                                    2
                                );
                }
            }

        } 
        else if (mMalDistributionDirection == MAL_DISTRIBUTION_DIRECTION_CUSTOM_DIRECTION)
        {
            #pragma omp parallel reduction(+:result)
            {
                ModelPart::NodeIterator NodesBegin;
                ModelPart::NodeIterator NodesEnd;
                OpenMPUtils::PartitionedIterators(rSurfaceModelPart.Nodes(), NodesBegin, NodesEnd);

                for (auto it = NodesBegin; it != NodesEnd; ++it)
                {
                    const array_1d<double,3>& normal = it->FastGetSolutionStepValue(NORMAL, 0);
                    const array_1d<double,3>& velocity = it->FastGetSolutionStepValue(VELOCITY, 0);
                    
                    double magnitude = norm_2(normal);

                    result += pow(
                                    (
                                        mDirection[0]*velocity[0] + 
                                        mDirection[1]*velocity[1] + 
                                        mDirection[2]*velocity[2]
                                    )*magnitude / mAverageSurfaceArea - mAverageVelocity,
                                    2
                                );
                }
            }

        }

        return sqrt(result/(mN-1));
    }

    virtual void FinalizeSolutionStep(ModelPart& rModelPart)
    {
        ModelPart& rSurfaceModelPart = rModelPart.GetSubModelPart(mSurfaceModelPartName);
        NormalCalculationUtils normal_calculation_obj = NormalCalculationUtils();
        const unsigned int domain_size = rModelPart.GetProcessInfo().GetValue(DOMAIN_SIZE);
        normal_calculation_obj.CalculateOnSimplex(rSurfaceModelPart, domain_size);
        CalculateAverageObjectiveParameters(rModelPart);

        ProcessInfo& rProcessInfo = rModelPart.GetProcessInfo();
        mObjectiveValue = Calculate(rModelPart);
        mOutputFileStream.precision(5);
        mOutputFileStream<<std::scientific<<rProcessInfo[TIME]<<" ";
        mOutputFileStream.precision(15);
        mOutputFileStream<<std::scientific<<mObjectiveValue<<std::endl;        
    }

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::string mSurfaceModelPartName;
    std::string mOutputFilename;
    std::ofstream mOutputFileStream;
    bool mOutputFileOpenend = false;
    array_1d<double, TDim> mDirection;
    std::map<unsigned long,unsigned int> mNodeSharingCounter;
    int mMalDistributionDirection;
    double mAverageSurfaceArea;
    double mN;
    double mAverageVelocity;
    double mObjectiveValue;

    const int MAL_DISTRIBUTION_DIRECTION_SURFACE_NORMAL = 1;
    const int MAL_DISTRIBUTION_DIRECTION_CUSTOM_DIRECTION = 2;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{   
    void CalculateAverageObjectiveParameters(ModelPart& rModelPart)
    {
        ModelPart& rSurfaceModelPart = rModelPart.GetSubModelPart(mSurfaceModelPartName);

        mAverageSurfaceArea = 0.0;
        mAverageVelocity = 0.0;
        mN = 0;

        if (mMalDistributionDirection == MAL_DISTRIBUTION_DIRECTION_SURFACE_NORMAL)
        {
            for (auto it = rSurfaceModelPart.NodesBegin();
                it != rSurfaceModelPart.NodesEnd();
                ++it)
            {
                const array_1d<double,3>& normal = it->FastGetSolutionStepValue(NORMAL);
                double magnitude = norm_2(normal);
                const array_1d<double,3>& velocity = it->FastGetSolutionStepValue(VELOCITY, 0);

                mAverageVelocity += (
                                        normal[0]*velocity[0] + 
                                        normal[1]*velocity[1] + 
                                        normal[2]*velocity[2]
                                    );
                mAverageSurfaceArea += magnitude;
                mN++;
            }
        }
        else if (mMalDistributionDirection == MAL_DISTRIBUTION_DIRECTION_CUSTOM_DIRECTION)
        {
            for (auto it = rSurfaceModelPart.NodesBegin();
                it != rSurfaceModelPart.NodesEnd();
                ++it)
            {
                const array_1d<double,3>& normal = it->FastGetSolutionStepValue(NORMAL);
                const array_1d<double,3>& velocity = it->FastGetSolutionStepValue(VELOCITY, 0);

                double magnitude = norm_2(normal);

                mAverageVelocity += (
                                        mDirection[0]*velocity[0] + 
                                        mDirection[1]*velocity[1] + 
                                        mDirection[2]*velocity[2]
                                    )*magnitude;

                mAverageSurfaceArea += magnitude;
                mN++;
            }
        }

        //Getting the average velocity by dividing flow rate by total objective surface area
        mAverageVelocity /= mAverageSurfaceArea;
        mAverageSurfaceArea /= mN;
    }

    /*void CalculateAreaGradients(const Element& rElem,
                                    std::vector<array_1d<double, 3>&>& rAreaGradients,
                                    ProcessInfo& rProcessInfo)
    {
        
        const unsigned int NumNodes = rElem.GetGeometry().PointsNumber();

        double inv_coeff = 1.0/(mAverageSurfaceArea*(mN-1)*mObjectiveValue);

        Geometry<Node<3> >& pGeometry = rElem.GetGeometry();

        array_1d<double, 3> v1, v2, p;

        unsigned int LocalIndex = 0;
        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            int index_1 = iNode - 1;
            int index_2 = iNode + 1;

            if (iNode < 0)
                index_1 = NumNodes - 1;
            if (iNode > NumNodes - 1)
                index_2 = 0

            p[0] = pGeometry[iNode].X();
            p[1] = pGeometry[iNode].Y();
            p[2] = pGeometry[iNode].Z();
(
            v1[0] = pGeometry[index_1].X() - pGeometry[iNode].X();
            v1[1] = pGeometry[index_1].Y() - pGeometry[iNode].Y();
            v1[2] = pGeometry[index_1].Z() - pGeometry[iNode].Z();

            v2[0] = pGeometry[index_2].X() - pGeometry[iNode].X();
            v2[1] = pGeometry[index_2].Y() - pGeometry[iNode].Y();
            v2[2] = pGeometry[index_2].Z() - pGeometry[iNode].Z();

            double mag_v1 = norm_2(v1);
            double mag_v2 = norm_2(v2);

            double mag_v1_2 = pow(mag_v1, 2);
            double mag_v2_2 = pow(mag_v2, 2);
            double v1_dot_v2 = dot_product(v1, v2);

            double cos_alpha = v1_dot_v2/(mag_v1*mag_v2);
            double sin_alpha = sqrt(1-pow(cos_alpha,2));

            double element_area = 0.5*mag_v1*mag_v2*sin_alpha;

            array_1d<double, 3> gradient_of_v1_dot_v2;

            gradient_of_v1_dot_v2 = -v1 -v2;

            array_1d<double, 3> temp = v1/mag_v1_2 + v2/mag_v2_2;

            rAreaGradients[iNode] = -2*element_area*temp 
                        -(v1_dot_v2*0.5/elemen_area)*(
                            gradient_of_v1_dot_v2 + v1_dot_v2*temp
                            );
            rAreaGradients[iNode] *= 0.5;

        }
    }

    void double dot_product(const array_1d<double, 3>& a, 
                            const array_1d<double, 3>& b)
    {
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }*/
    ///@}
};

///@} // Kratos Classes

///@} // Adjoint Fluid Application group

} /* namespace Kratos.*/

#endif /* KRATOS_FLOW_MAL_DISTRIBUTION_OBJECTIVE_FUNCTION defined */
