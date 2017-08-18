/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

/* *********************************************************
*
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: 29 Mar 2017 $
*   Revision:            $Revision: 1.13 $
*
* ***********************************************************/

#if !defined(KRATOS_VARIABLE_FAST_TRANSFER_UTILITY_INCLUDED )
#define  KRATOS_VARIABLE_FAST_TRANSFER_UTILITY_INCLUDED

//System includes
#ifdef _OPENMP
#include <omp.h>
#endif

//External includes
#include "boost/smart_ptr.hpp"
#include "boost/timer.hpp"
#include "boost/progress.hpp"

//Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "includes/element.h"
#include "geometries/geometry.h"
#include "spaces/ublas_space.h"
#include "structural_application.h"

namespace Kratos
{
class VariableFastTransferUtility
{
public:
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;
    typedef LinearSolver<SparseSpaceType, DenseSpaceType> LinearSolverType;
    typedef Element::GeometryType GeometryType;
    typedef GeometryType::PointType NodeType;
    typedef NodeType::PointType PointType;

    /**
     * Constructor.
     */
    VariableFastTransferUtility(ModelPart::ElementsContainerType& pElements, LinearSolverType::Pointer pLinearSolver,
        const double& Dx, const double& Dy, const double& Dz)
    : mpElements(pElements), mpLinearSolver(pLinearSolver), mEchoLevel(0), mDx(Dx), mDy(Dy), mDz(Dz)
    {
        // extract the active nodes
        for(ModelPart::ElementsContainerType::ptr_iterator it = pElements.ptr_begin();
                it != pElements.ptr_end(); ++it)
        {
            if( (*it)->GetValue(IS_INACTIVE) == false || (*it)->Is(ACTIVE) )
            {
                for(std::size_t i = 0; i < (*it)->GetGeometry().size(); ++i)
                {
                    mActiveNodes.insert( (*it)->GetGeometry()[i].Id() );
                }
            }
        }

        // assign each node an id. That id is the row of this node in the global L2 projection matrix
        std::size_t cnt = 0;
        for(std::set<std::size_t>::iterator it = mActiveNodes.begin(); it != mActiveNodes.end(); ++it )
        {
            mNodeRowId[*it] = cnt++;
        }

        std::cout << "VariableFastTransferUtility created" << std::endl;
    }

    /**
     * Destructor.
     */
    virtual ~VariableFastTransferUtility() {}

    void SetEchoLevel(int Level)
    {
        mEchoLevel = Level;
    }

    int GetEchoLevel()
    {
        return mEchoLevel;
    }

protected:

    LinearSolverType::Pointer mpLinearSolver;

    //**********AUXILIARY FUNCTION**************************************************************
    //******************************************************************************************
    /// Find an element in pMasterElements contains rSourcePoint and assign it to pTargetElement. The rLocalTargetPoint is the local point in pTargetElement of rSourcePoint
    /// REMARK: we should disable the move mesh flag if we want to search in the reference configuration
    bool SearchPartnerWithBin( const PointType& rSourcePoint, ModelPart::ElementsContainerType& pMasterElements,
            Element::Pointer& pTargetElement, PointType& rLocalTargetPoint ) const
    {
        ModelPart::ElementsContainerType pMasterElementsCandidates;

        // get the containing elements from the bin
        int ix = (mDx != 0.0) ? (int) floor(rSourcePoint.X()) / mDx : 0;
        int iy = (mDy != 0.0) ? (int) floor(rSourcePoint.Y()) / mDy : 0;
        int iz = (mDz != 0.0) ? (int) floor(rSourcePoint.Z()) / mDz : 0;

        SpatialKey key(ix, iy, iz);
        std::map<SpatialKey, std::set<std::size_t> >::const_iterator it_bin_elements = mBinElements.find(key);

        if(it_bin_elements != mBinElements.end())
        {
            for(std::set<std::size_t>::const_iterator it = it_bin_elements->second.begin(); it != it_bin_elements->second.end(); ++it )
            {
                Element::GeometryType& r_geom = pMasterElements[*it].GetGeometry();

                bool is_inside = r_geom.IsInside( rSourcePoint, rLocalTargetPoint );
                if( is_inside )
                {
                    pTargetElement = pMasterElements(*it);
                    return true;
                }
            }
        }

        std::cout << " !!!! WARNING: NO ELEMENT FOUND TO CONTAIN " << rSourcePoint << " !!!! " << std::endl;
        return false;
    }

    //**********AUXILIARY FUNCTION**************************************************************
    //******************************************************************************************
    /// Construct the L2 projection matrix
    void ConstructLHSMatrix(SparseSpaceType::MatrixType& rA)
    {
        // set up system
        if(rA.size1() != mActiveNodes.size() || rA.size2() != mActiveNodes.size())
            rA.resize(mActiveNodes.size(), mActiveNodes.size());
        noalias(rA) = ZeroMatrix(mActiveNodes.size(), mActiveNodes.size());

        int number_of_threads = 1;
#ifdef _OPENMP
        number_of_threads = omp_get_max_threads();
#endif
        vector<unsigned int> element_partition;
        CreatePartition(number_of_threads, mpElements.size(), element_partition);
//        boost::progress_display show_progress( mpElements.size() );

        // create the structure for M a priori
        ConstructMatrixStructure(rA, mpElements, mNodeRowId);

#ifdef _OPENMP
        //create the array of lock
        std::vector< omp_lock_t > lock_array(rA.size1());
        unsigned int system_size = rA.size1();
        for(unsigned int i = 0; i < system_size; ++i)
            omp_init_lock(&lock_array[i]);
#endif

#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for(int k = 0; k < number_of_threads; ++k)
        {
            ModelPart::ElementsContainerType::ptr_iterator it_begin = mpElements.ptr_begin() + element_partition[k];
            ModelPart::ElementsContainerType::ptr_iterator it_end = mpElements.ptr_begin() + element_partition[k+1];

            for( ModelPart::ElementsContainerType::ptr_iterator it = it_begin; it != it_end; ++it )
            {
                if( (*it)->GetValue(IS_INACTIVE) == true && !(*it)->Is(ACTIVE) )
                    continue;

                unsigned int dim = (*it)->GetGeometry().WorkingSpaceDimension();

                const GeometryType::IntegrationPointsArrayType& integration_points
                    = (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());

                GeometryType::JacobiansType J(integration_points.size());
                J = (*it)->GetGeometry().Jacobian(J, (*it)->GetIntegrationMethod());

                const Matrix& Ncontainer = (*it)->GetGeometry().ShapeFunctionsValues((*it)->GetIntegrationMethod());

                double DetJ;
                for(unsigned int point = 0; point < integration_points.size(); ++point)
                {
                    DetJ = MathUtils<double>::Det(J[point]);

                    double dV = DetJ*integration_points[point].Weight();

                    for(unsigned int prim = 0 ; prim < (*it)->GetGeometry().size(); ++prim)
                    {
                        unsigned int row = mNodeRowId[(*it)->GetGeometry()[prim].Id()];
#ifdef _OPENMP
                        omp_set_lock(&lock_array[row]);
#endif
                        for(unsigned int sec = 0 ; sec < (*it)->GetGeometry().size(); ++sec)
                        {
                            unsigned int col = mNodeRowId[(*it)->GetGeometry()[sec].Id()];
                            rA(row, col) += Ncontainer(point, prim)*Ncontainer(point, sec) * dV;
                        }
#ifdef _OPENMP
                        omp_unset_lock(&lock_array[row]);
#endif
                    }
                }
//                ++show_progress;
            }
        }
    }

    //**********AUXILIARY FUNCTION**************************************************************
    //******************************************************************************************
    void ConstructMatrixStructure(SparseSpaceType::MatrixType& A,
            ModelPart::ElementsContainerType& pElements, std::map<std::size_t, std::size_t>& NodeRowId) const
    {
        std::size_t equation_size = A.size1();
        std::vector<std::vector<std::size_t> > indices(equation_size);

        Element::EquationIdVectorType ids;
        for(ModelPart::ElementsContainerType::iterator i_element = pElements.begin();
                i_element != pElements.end() ; ++i_element)
        {
            if( !(i_element)->GetValue( IS_INACTIVE ) || (i_element)->Is(ACTIVE) )
            {
                ids.resize((i_element)->GetGeometry().size());
                for(unsigned int i = 0; i < (i_element)->GetGeometry().size();  ++i)
                {
                    ids[i] = NodeRowId[(i_element)->GetGeometry()[i].Id()];
                }

                for(std::size_t i = 0 ; i < ids.size() ; ++i)
                {
                    if(ids[i] < equation_size)
                    {
                        std::vector<std::size_t>& row_indices = indices[ids[i]];
                        for(std::size_t j = 0 ; j < ids.size() ; j++)
                        {
                            if(ids[j] < equation_size)
                            {
                                AddUnique(row_indices,ids[j]);
                            }
                        }
                    }
                }
            }
        }

        //allocating the memory needed
        int data_size = 0;
        for(std::size_t i = 0 ; i < indices.size() ; i++)
        {
            data_size += indices[i].size();
        }
        A.reserve(data_size, false);

        //filling with zero the matrix (creating the structure)
#ifndef _OPENMP
        for(std::size_t i = 0 ; i < indices.size() ; i++)
        {
            std::vector<std::size_t>& row_indices = indices[i];
            std::sort(row_indices.begin(), row_indices.end());

            for(std::vector<std::size_t>::iterator it = row_indices.begin(); it != row_indices.end() ; it++)
            {
                A.push_back(i,*it,0.00);
            }
            row_indices.clear();
        }
#else
        int number_of_threads = omp_get_max_threads();
        vector<unsigned int> matrix_partition;
        CreatePartition(number_of_threads, indices.size(), matrix_partition);
        for( int k=0; k<number_of_threads; ++k )
        {
            #pragma omp parallel
            if( omp_get_thread_num() == k )
            {
                for( std::size_t i = matrix_partition[k]; i < matrix_partition[k+1]; ++i )
                {
                    std::vector<std::size_t>& row_indices = indices[i];
                    std::sort(row_indices.begin(), row_indices.end());

                    for(std::vector<std::size_t>::iterator it = row_indices.begin(); it != row_indices.end() ; ++it)
                    {
                        A.push_back(i, *it, 0.00);
                    }

                    row_indices.clear();
                }
            }
        }
#endif
    }

    //**********AUXILIARY FUNCTION**************************************************************
    //******************************************************************************************
    inline void AddUnique(std::vector<std::size_t>& v, const std::size_t& candidate) const
    {
        std::vector<std::size_t>::iterator i = v.begin();
        std::vector<std::size_t>::iterator endit = v.end();
        while ( i != endit && (*i) != candidate)
        {
            ++i;
        }
        if( i == endit )
        {
            v.push_back(candidate);
        }

    }

    //**********AUXILIARY FUNCTION**************************************************************
    //******************************************************************************************
    inline void CreatePartition(unsigned int number_of_threads,const int number_of_rows, vector<unsigned int>& partitions) const
    {
        partitions.resize(number_of_threads+1);
        int partition_size = number_of_rows / number_of_threads;
        partitions[0] = 0;
        partitions[number_of_threads] = number_of_rows;
        for(unsigned int i = 1; i<number_of_threads; i++)
            partitions[i] = partitions[i-1] + partition_size ;
    }

private:

    struct SpatialKey
    {
        public:
            SpatialKey(int ix, int iy, int iz) : x(ix), y(iy), z(iz) {}
            bool operator<(const SpatialKey& rOther) const
            {
                if(x == rOther.x)
                {
                    if(y == rOther.y)
                    {
                        return z < rOther.z;
                    }
                    else
                        return y < rOther.y;
                }
                else
                    return x < rOther.x;
            }
            int kx() const {return x;}
            int ky() const {return y;}
            int kz() const {return z;}
        private:
            int x, y, z;
    };

    ModelPart::ElementsContainerType mpElements;
    int mEchoLevel;
    double mDx, mDy, mDz;
    std::map<SpatialKey, std::set<std::size_t> > mBinElements;
    std::set<std::size_t> mActiveNodes;
    std::map<std::size_t, std::size_t> mNodeRowId;

};//Class VariableFastTransferUtility

}//namespace Kratos.

#endif /* KRATOS_VARIABLE_FAST_TRANSFER_UTILITY  defined */
