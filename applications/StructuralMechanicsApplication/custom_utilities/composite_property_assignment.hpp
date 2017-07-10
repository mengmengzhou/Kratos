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
			const ElementsIteratorType& firstElement = rSubModelpart.ElementsBegin();
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
			Vector localAxis1 = ZeroVector(3);
			double cosTheta, theta;
			Properties::Pointer pElementProps;

			std::cout << "Global fiber direction for model part '" << rSubModelpart.Name() << "' is: " << GlobalFiberDirection << std::endl;

			// Check incoming GlobalFiberDirection is valid
			if (inner_prod(GlobalFiberDirection, GlobalFiberDirection) == 0.0)
			{
				KRATOS_ERROR <<
					"Defined global fiber direction for subModelPart " <<
					rSubModelpart.Name() << " has zero length" << std::endl;
			}

			GlobalFiberDirection /= std::sqrt(inner_prod(GlobalFiberDirection, GlobalFiberDirection));

			for (auto& element : rSubModelpart.Elements())
			{
				// get current element properties
				pElementProps = element.pGetProperties();

				// get local orientation of GlobalFiberDirection
				element.Calculate(LOCAL_ELEMENT_ORIENTATION, LCSOrientation, rCurrentProcessInfo);


				// Project general global orientation onto element plane
				// global cartesian coords.
				Vector3& A = GlobalFiberDirection;
				Vector3 B; // global cartesian element normal
				for (size_t i = 0; i < 3; i++)
				{
					B[i] = LCSOrientation(2, i);
				}
				double magB = std::sqrt(inner_prod(B, B));

				Vector temp1 = MathUtils<double>::CrossProduct(A, B);
				temp1 /= magB;
				Vector projGlobalFiber = MathUtils<double>::CrossProduct(B, temp1);
				projGlobalFiber /= magB;
				//std::cout << "\nprojected global fiber = " << projGlobalFiber << std::endl;
				projGlobalFiber /= std::sqrt(inner_prod(projGlobalFiber, projGlobalFiber));
				//std::cout << "normalised projected global fiber = " << projGlobalFiber << std::endl;

				localGlobalFiberDirection = prod(LCSOrientation, GlobalFiberDirection);
				//std::cout << "pure local fiber = " << localGlobalFiberDirection << std::endl;
				localGlobalFiberDirection[2] = 0.0;
				localGlobalFiberDirection /= std::sqrt(inner_prod(localGlobalFiberDirection, localGlobalFiberDirection));

				//std::cout << "local fiber = " << localGlobalFiberDirection << std::endl;

				Vector3 localAxis2;

				// get element local axis 1 vector
				for (size_t i = 0; i < 3; i++)
				{
					localAxis1[i] = LCSOrientation(0, i);
					localAxis2[i] = LCSOrientation(1, i);
				}

				//std::cout << "global ax1 = " << localAxis1 << std::endl;
				//std::cout << "global ax2 = " << localAxis2 << std::endl;
				//std::cout << "global ax3 = " << B << std::endl;

				localAxis1 = prod(LCSOrientation, localAxis1);
				localAxis1 /= std::sqrt(inner_prod(localAxis1, localAxis1));
				//std::cout << "flattened lc1 = " << localAxis1 << std::endl;

				localAxis2 = prod(LCSOrientation, localAxis2);
				localAxis2 /= std::sqrt(inner_prod(localAxis2, localAxis2));
				//std::cout << "flattened lc2 = " << localAxis2 << std::endl;

				B = prod(LCSOrientation, B);
				B /= std::sqrt(inner_prod(B, B));
				//std::cout << "flattened lc3 = " << B << std::endl;

				// compute angle 'theta' between local axis 1 and
				// localGlobalFiberDirection
				cosTheta = inner_prod(localAxis1, localGlobalFiberDirection);
				//cosTheta /= (std::sqrt(inner_prod(localAxis1, localAxis1)));
				//cosTheta /= (std::sqrt(inner_prod(localGlobalFiberDirection, localGlobalFiberDirection)));
				theta = std::acos(cosTheta);

				
				//std::cout << "angle in degrees: " << theta / 2.0 / 3.14 * 360 << std::endl;

				// dot between lc2 and localFiberDir
				double dotCheck = inner_prod(localAxis2, localGlobalFiberDirection);
				//std::cout << "lc2 dot localGlobalfiber = " << dotCheck << std::endl;
				if (dotCheck < 0.0)
				{
					// theta is currently negative, flip to positive definition
					theta *= -1.0;
				}


				// Rotate by 180 degrees if we need to
				Matrix localToFiberRotation = Matrix(3, 3, 0.0);
				double c = std::cos(theta);
				double s = std::sin(theta);
				localToFiberRotation(0, 0) = c;
				localToFiberRotation(0, 1) = -s;
				localToFiberRotation(1, 0) = s;
				localToFiberRotation(1, 1) = c;
				localToFiberRotation(2, 2) = 1.0;
				Vector3 temp = prod(localToFiberRotation, localGlobalFiberDirection);

				//std::cout << "temp = " << temp << std::endl;
				double dotCheck1 = inner_prod(temp, localAxis1);
				//std::cout << "temp dot localAxis1 = " << dotCheck << std::endl;
				double dotCheck2 = inner_prod(temp, localAxis2);
				//std::cout << "temp dot localAxis2 = " << dotCheck << std::endl;

				/*
				if (dotCheck1 < 0.0)
				{
					theta += KRATOS_M_PI;
				}
				else if (dotCheck2 < 0.0)
				{
					theta += KRATOS_M_PI;
				}
				*/
				if (std::abs(theta) > 2.0*KRATOS_M_PI)
				{
					theta = 0.0;
				}
				
				bool test_lc1_global = false;
				if (test_lc1_global)
				{
					// rotate localAxis1 to local fiber, then transform to global
					Matrix localToFiberRotation = Matrix(3, 3, 0.0);
					double c = std::cos(theta);
					double s = std::sin(theta);
					localToFiberRotation(0, 0) = c;
					localToFiberRotation(0, 1) = -s;
					localToFiberRotation(1, 0) = s;
					localToFiberRotation(1, 1) = c;
					localToFiberRotation(2, 2) = 1.0;

					Vector rotatedLC1 = prod(localToFiberRotation, localAxis1);
					Vector recoveredGlobalFiber = prod(trans(LCSOrientation), rotatedLC1);
					std::cout << "actual local fiber (z=0) = " << localGlobalFiberDirection << std::endl;
					std::cout << "recovered local fiber = " << rotatedLC1 << std::endl;
					std::cout << "actual global fiber = " << GlobalFiberDirection << std::endl;
					std::cout << "recovered global fiber = " << recoveredGlobalFiber << std::endl;
				}

				// set required rotation in element
				pElementProps = element.pGetProperties();
				pElementProps->SetValue(ORTHOTROPIC_ORIENTATION_ASSIGNMENT, theta);
				element.Calculate(ORTHOTROPIC_ORIENTATION_ASSIGNMENT, theta, rCurrentProcessInfo);

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