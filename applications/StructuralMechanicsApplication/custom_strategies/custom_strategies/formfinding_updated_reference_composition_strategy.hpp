// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Long Chen
//
//

#if !defined(KRATOS_FORMFINDING_UPDATED_REFERENCE_COMPOSITION_STRATEGY)
#define KRATOS_FORMFINDING_UPDATED_REFERENCE_COMPOSITION_STRATEGY


/* System includes */
//#include<limits>
//#include<iostream>
//#include<iomanip>

/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
//#include "solving_strategies/convergencecriterias/convergence_criteria.h"

//default builder and solver


//default builder and solver
//#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

namespace Kratos
{
template<class TSparseSpace,
    class TDenseSpace, // = DenseSpace<double>,
    class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
    >
class FormfindingUpdatedReferenceCompositionStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions 
    ///@{
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

    // Counted pionter of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(FormfindingUpdatedReferenceCompositionStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TSchemeType TSchemeType;

    //typedef typename BaseType::DofSetType DofSetType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    ///@}
    ///@name Life Cycle

    ///@{

    /**
    * Constructor.
    */

    FormfindingUpdatedReferenceCompositionStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
        )
    {
        mResidualBasedNewtonRaphsonStrategy = ResidualBasedNewtonRaphsonStrategy(
            model_part, 
            pScheme,
            pNewLinearSolver,
            pNewConvergenceCriteria,
            MaxIterations, 
            CalculateReactions,
            ReformDofSetAtEachStep,
            MoveMeshFlag);
    };

    FormfindingUpdatedReferenceCompositionStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
        )
    {
        mResidualBasedNewtonRaphsonStrategy = ResidualBasedNewtonRaphsonStrategy(
            model_part,
            pScheme,
            pNewLinearSolver,
            pNewConvergenceCriteria,
            pNewBuilderAndSolver,
            MaxIterations,
            CalculateReactions,
            ReformDofSetAtEachStep,
            MoveMeshFlag);
    };

    virtual ~FormfindingUpdatedReferenceCompositionStrategy()
    {
        Clear();    // is it needed?
    };

    void SetScheme(typename TSchemeType::Pointer pScheme)
    {
        mResidualBasedNewtonRaphsonStrategy.SetScheme(pScheme);
    };

    typename TSchemeType::Pointer GetScheme()
    {
        return mpScheme;
    };

    // Set and Get the BuilderAndSolver

    void SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer);


private:
    ResidualBasedNewtonRaphsonStrategy mResidualBasedNewtonRaphsonStrategy;

};



}

#endif /* KRATOS_FORMFINDING_UPDATED_REFERENCE_COMPOSITION_STRATEGY defined*/
