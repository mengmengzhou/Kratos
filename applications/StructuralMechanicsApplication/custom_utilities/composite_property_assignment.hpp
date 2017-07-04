//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson
//                   Philipp Bucher
//

#if !defined(KRATOS_COMPOSITE_PROPERTY_ASSIGNMENT_HPP_INCLUDED )
#define  KRATOS_COMPOSITE_PROPERTY_ASSIGNMENT_HPP_INCLUDED

// System includes
#include <iterator>
#include <math.h>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes
#include "structural_mechanics_application_variables.h"

namespace Kratos {

	///@name Kratos Globals
	///@{

	///@}
	///@name Type Definitions
	///@{
	typedef ModelPart::ElementsContainerType::iterator ElementsIteratorType;

	///@}
	///@name  Enum's
	///@{

	///@}
	///@name  Functions
	///@{

	///@}
	///@name Kratos Classes
	///@{

	/// Transfer eigenvectors to solution step variables for GiD output or solution initialization.
	/**
	* Example Python Code:
	* # Eigenvectors are first computed and stored in the nodal variable EIGENVECTOR_MATRIX.
	* for step in range(NumEigenvalues):
	*   main_model_part.ProcessInfo[TIME] = float(step+1)
	*   EigenvectorToSolutionStepVariableTransferUtility().Transfer(main_model_part,step,0)
	*   gid_output.PrintOutput()
	*/
	class CompositePropertyAssignment {
	public:
		///@name Type Definitions
		///@{

		KRATOS_CLASS_POINTER_DEFINITION(CompositePropertyAssignment);

		///@}
		///@name Life Cycle
		///@{

		CompositePropertyAssignment() {}

		virtual ~CompositePropertyAssignment() {}

		///@}
		///@name Operators
		///@{

		///@}
		///@name Operations
		///@{

		void Execute(ModelPart& rSubModelpart, Vector3 GlobalFiberDirection, ProcessInfo& rCurrentProcessInfo)
		{
			// Check to see if the composite orientation assignment has already
			// been performed on the current modelPart
			// Just look at the first element to save time
			ElementsIteratorType& firstElement = rSubModelpart.ElementsBegin();
			Properties elementProperties = (*firstElement).GetProperties();

			if (elementProperties.Has(ORTHOTROPIC_ORIENTATION_ASSIGNMENT))
			{
				// the composite orientation assignment has already been done
			}
			else
			{
				// perform the composite orientation assignment
				compositeOrientationAssignment(rSubModelpart, 
					GlobalFiberDirection, rCurrentProcessInfo);
			}
		}

		///@}
		///@name Access
		///@{

		///@}
		///@name Inquiry
		///@{

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
		void compositeOrientationAssignment(ModelPart& rSubModelpart, Vector3 GlobalFiberDirection, ProcessInfo& rCurrentProcessInfo)
		{
			// Declare working variables
			Matrix LCSOrientation;
			Vector3 localGlobalFiberDirection;
			Vector3 localAxis1 = Vector3(3, 0.0);
			double cosTheta, theta;
			Properties::Pointer pElementProps;

			for (auto& element : rSubModelpart.Elements())
			{
				// get current element properties
				pElementProps = element.pGetProperties();


				// get local orientation of GlobalFiberDirection
				element.Calculate(LOCAL_ELEMENT_ORIENTATION, LCSOrientation, rCurrentProcessInfo);
				localGlobalFiberDirection = prod(LCSOrientation, GlobalFiberDirection);
				localGlobalFiberDirection[2] = 0.0; // remove z component


				// get element local axis 1 vector
				for (size_t i = 0; i < 3; i++)
				{
					localAxis1[i] = LCSOrientation(0, i);
				}


				// compute current angle 'theta' between local axis 1 and 
				// localGlobalFiberDirection
				cosTheta = inner_prod(localAxis1, localGlobalFiberDirection);
				cosTheta /= (std::sqrt(inner_prod(localAxis1, localAxis1)));
				cosTheta /= (std::sqrt(inner_prod(localGlobalFiberDirection, localGlobalFiberDirection)));
				theta = acos(cosTheta); // current angle between ax1 and fiber
				theta *= -1.0; // rotation necessary to align


				// set required rotation in element properties.
				// rotation is performed in the element's InitializeSolutionStep
				// call
				pElementProps = element.pGetProperties();
				pElementProps->SetValue(ORTHOTROPIC_ORIENTATION_ASSIGNMENT, theta);


				// add option to write out angles so they don't have to be computed next time
					// or maybe this should be a separate python call

			}// sub-modelpart element loop
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

		///@}

	}; // class CompositePropertyAssignment

	   ///@}

	   ///@name Type Definitions
	   ///@{

	   ///@}

}
// namespace Kratos
#endif  // KRATOS_COMPOSITE_PROPERTY_ASSIGNMENT_HPP_INCLUDED defined
