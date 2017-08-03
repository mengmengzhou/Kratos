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
		void Execute(ModelPart& rSubModelpart, Vector3 GlobalFiberDirection, Vector3 normalVector, ProcessInfo& rCurrentProcessInfo)
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
					GlobalFiberDirection, normalVector, rCurrentProcessInfo);
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
		void compositeOrientationAssignment(ModelPart& rSubModelpart, Vector3 GlobalFiberDirection, Vector3 normalVector, ProcessInfo& rCurrentProcessInfo)
		{
			// Declare working variables
			Matrix LCSOrientation, R;
			Vector localGlobalFiberDirection, rotation_axis;
			Vector localAxis1 = ZeroVector(3);
			Vector localAxis2 = ZeroVector(3);
			Vector localAxis3 = ZeroVector(3);
			double cosTheta, theta, rotation_angle, dotCheck;
			Properties::Pointer pElementProps;

			// Optional printout of details
			bool printDetails = true;
			if (printDetails)
			{
				std::cout << "Composite orientation assignment details for model part '" << rSubModelpart.Name() << "': " << std::endl;
				std::cout << "\tGlobal fiber direction = \t" << GlobalFiberDirection << std::endl;
				std::cout << "\tProjection direction = \t\t" << normalVector << std::endl;
			}

			// Check incoming GlobalFiberDirection and normalVector are valid
			if (inner_prod(GlobalFiberDirection, GlobalFiberDirection) == 0.0)
			{
				KRATOS_ERROR <<
					"Defined global fiber direction for subModelPart " <<
					rSubModelpart.Name() << " has zero length" << std::endl;
			}
			else if (inner_prod(normalVector, normalVector) == 0.0)
			{
				KRATOS_ERROR <<
					"Defined normal vector for subModelPart " <<
					rSubModelpart.Name() << " has zero length" << std::endl;
			}

			// Normalize
			GlobalFiberDirection /= std::sqrt(inner_prod(GlobalFiberDirection, GlobalFiberDirection));
			normalVector /= std::sqrt(inner_prod(normalVector, normalVector));

			for (auto& element : rSubModelpart.Elements())
			{
				// get current element properties
				pElementProps = element.pGetProperties();

				// get local orientation of GlobalFiberDirection
				element.Calculate(LOCAL_ELEMENT_ORIENTATION, LCSOrientation, rCurrentProcessInfo);

				// get element local axis vectors (global cartesian)
				for (size_t i = 0; i < 3; i++)
				{
					localAxis1[i] = LCSOrientation(0, i);
					localAxis2[i] = LCSOrientation(1, i);
					localAxis3[i] = LCSOrientation(2, i);
				}

				// normalise local axis vectors (global cartesian)
				localAxis1 /= std::sqrt(inner_prod(localAxis1, localAxis1));
				localAxis2 /= std::sqrt(inner_prod(localAxis2, localAxis2));
				localAxis3 /= std::sqrt(inner_prod(localAxis3, localAxis3));

				// Make deep copy of local fiber direction (global cartesian)
				localGlobalFiberDirection = Vector3(GlobalFiberDirection);

				// Flip projection vector such that is it in same dir as LC3
				dotCheck = inner_prod(normalVector, localAxis3);
				Vector tempNormalVector = Vector(normalVector);
				if (dotCheck < 0.0)
				{
					tempNormalVector *= -1.0;
				}

				// get rotation matrix to align element normal with projection vec (global cartesian)
				//rotation_axis = MathUtils<double>::CrossProduct(localAxis3, normalVector); // original
				rotation_axis = MathUtils<double>::CrossProduct(tempNormalVector, localAxis3);
				rotation_angle = inner_prod(tempNormalVector, localAxis3);
				if (abs(rotation_angle) < (1.0 - 1E-6)) // skip if already co-linear
				{
					rotation_angle = std::acos(rotation_angle);
					R = setUpRotationMatrix(rotation_angle, rotation_axis);
					localGlobalFiberDirection = prod(R, localGlobalFiberDirection);
				}


				Vector mytemp = prod((R), tempNormalVector);


				bool debugPrint = false;
				if (element.Id() == 137)
				{
					std::cout << element.Info() << std::endl;
					debugPrint = true;
				}
				
				if (debugPrint)
				{
					std::cout << "\nunrotated normalVector dot localAxis3 = " << inner_prod(tempNormalVector, localAxis3) << std::endl;
					std::cout << "rotated normalVector dot localAxis3 = " << inner_prod(mytemp, localAxis3) << std::endl;
					
					std::cout << "rotated localGlobalFiberDirection dot localAxis3 = " << inner_prod(localGlobalFiberDirection, localAxis3) << std::endl;
					std::cout << "rotated normalVector  = " << mytemp << std::endl;
					std::cout << "localAxis3  = " << localAxis3 << std::endl;
					std::cout << "rotated localGlobalFiberDirection  = " << localGlobalFiberDirection << std::endl;
					
				}
				
				// Correct in-plane transformed global cartesian fiber orientation so that
				// it maximises it's direction in accord with the pure
				// global fiber direction.
				bool rotated_fiber_correction = true;
				if (rotated_fiber_correction)
				{
					// Currently this is only setup to align with +Z by rotating 
					// the in-plane transformed global cartesian fiber orientation
					// such that the X component = 0.

					// This needs to be extended in to general coords.

					std::cout << "Element " << element.Id() << std::endl;
					double tolerance = 1E-6;
					double steps = 16.0;
					double step_size = 2.0*KRATOS_M_PI / steps; // initially 45 degrees
					double central_angle = 0.0; // from current alignment
					double max_dot_prod = -10.0;
					double best_angle = 0.0;
					bool converged = false;
					int iteration_limit = 20;
					int iteration = 0;

					while (converged == false)
					{
						for (size_t angle_step = 0; angle_step < steps; angle_step++)
						{
							double current_angle = best_angle + (angle_step - steps / 2.0)*step_size;
							R = setUpRotationMatrix(current_angle, localAxis3);
							Vector tempFiber = prod(R, localGlobalFiberDirection);
							double current_dot_prod = inner_prod(tempFiber, GlobalFiberDirection);
							current_dot_prod = tempFiber[0];
							if (debugPrint)
							{

								std::cout << current_dot_prod << std::endl;
							}
							if (abs(current_dot_prod) < abs(max_dot_prod))
							{
								max_dot_prod = current_dot_prod;
								best_angle = current_angle;
							}
						}
						step_size /= steps;
						iteration++;
						if (abs(max_dot_prod) < tolerance)
						{
							converged = true;
						}
						if (iteration > iteration_limit)
						{
							converged = true;
						}
					}

					// With the best angle determined, rotate the optimised
					// in-plane transformed global cartesian fiber orientation
					R = setUpRotationMatrix(best_angle, localAxis3);
					localGlobalFiberDirection = prod(R, localGlobalFiberDirection);

					// Make sure it points in the right direction.
					if (inner_prod(localGlobalFiberDirection, GlobalFiberDirection) < 0)
					{
						R = setUpRotationMatrix(KRATOS_M_PI, localAxis3);
						localGlobalFiberDirection = prod(R, localGlobalFiberDirection);
					}
				}
				



				bool global_approach = false;
				if (global_approach)
				{
					// get rotation matrix to align element normal with projection vec (global cartesian)
					//rotation_axis = MathUtils<double>::CrossProduct(localAxis3, normalVector); // original
					rotation_axis = MathUtils<double>::CrossProduct(localGlobalFiberDirection, localAxis1);
					rotation_axis /= std::sqrt(inner_prod(rotation_axis, rotation_axis));
					rotation_angle = inner_prod(localGlobalFiberDirection, localAxis1);
					if (abs(rotation_angle) < (1.0 - 1E-6)) // skip if already co-linear
					{
						rotation_angle = std::acos(rotation_angle);
						//R = setUpRotationMatrix(rotation_angle, rotation_axis);
						//localGlobalFiberDirection = prod(R, localGlobalFiberDirection);
						theta = rotation_angle*-1.0;
					}

					if (debugPrint)
					{
						std::cout << "inner_prod(rotation_axis,localAxis3) = " << inner_prod(rotation_axis, localAxis3) << std::endl;
						std::cout << "rotation axis = " << rotation_axis << std::endl;
						std::cout << "localAxis3 = " << localAxis3 << std::endl;
					}
					
					double dotProd = inner_prod(rotation_axis, localAxis3);
					if (dotProd < (1.0 - 1E-6))
					{
						theta *= -1.0;
					}

					
				}
				else
				{
					// Put everything in local space (local cartesian)
					localGlobalFiberDirection = prod(LCSOrientation, localGlobalFiberDirection);
					localAxis1 = prod(LCSOrientation, localAxis1);
					localAxis2 = prod(LCSOrientation, localAxis2);

					// compute angle 'theta' between local axis 1 and localGlobalFiberDirection (local cartesian)
					cosTheta = inner_prod(localAxis1, localGlobalFiberDirection);
					theta = std::acos(cosTheta);

					// dot between lc2 and localFiberDir (local cartesian)
					dotCheck = inner_prod(localAxis2, localGlobalFiberDirection);
					if (dotCheck < 0.0)
					{
						// theta is currently negative, flip to positive definition
						theta *= -1.0;
					}

					if (debugPrint)
					{
						std::cout << "rotation angle = " << theta / KRATOS_M_PI *180.0 << std::endl;
						std::cout << "localAxis1 = " << localAxis1 << std::endl;
						std::cout << "localAxis2 = " << localAxis2 << std::endl;
						std::cout << "localGlobalFiberDirection = " << localGlobalFiberDirection << std::endl;
					}
					
				}
				
				

				// set required rotation in element
				pElementProps = element.pGetProperties();
				pElementProps->SetValue(ORTHOTROPIC_ORIENTATION_ASSIGNMENT, theta);
				element.Calculate(ORTHOTROPIC_ORIENTATION_ASSIGNMENT, theta, rCurrentProcessInfo);

				// add option to write out angles so they don't have to be computed next time
					// or maybe this should be a separate python call
			}// sub-modelpart element loop
		}

		Matrix setUpRotationMatrix(double angle, Vector& rotation_axis)
		{
			Matrix rotationMatrix(3, 3, 0.0);

			double u = rotation_axis[0];
			double v = rotation_axis[1];
			double w = rotation_axis[2];

			double L = (u*u + v * v + w * w);
			double u2 = u * u;
			double v2 = v * v;
			double w2 = w * w;

			rotationMatrix(0, 0) = (u2 + (v2 + w2) * cos(angle)) / L;
			rotationMatrix(0, 1) = (u * v * (1 - cos(angle)) - w * sqrt(L) * sin(angle)) / L;
			rotationMatrix(0, 2) = (u * w * (1 - cos(angle)) + v * sqrt(L) * sin(angle)) / L;
			//rotationMatrix(0,3) = 0.0;

			rotationMatrix(1, 0) = (u * v * (1 - cos(angle)) + w * sqrt(L) * sin(angle)) / L;
			rotationMatrix(1, 1) = (v2 + (u2 + w2) * cos(angle)) / L;
			rotationMatrix(1, 2) = (v * w * (1 - cos(angle)) - u * sqrt(L) * sin(angle)) / L;
			//rotationMatrix(1,3) = 0.0;

			rotationMatrix(2, 0) = (u * w * (1 - cos(angle)) - v * sqrt(L) * sin(angle)) / L;
			rotationMatrix(2, 1) = (v * w * (1 - cos(angle)) + u * sqrt(L) * sin(angle)) / L;
			rotationMatrix(2, 2) = (w2 + (u2 + v2) * cos(angle)) / L;
			//rotationMatrix(2,3) = 0.0;

			//rotationMatrix(3,0) = 0.0;
			//rotationMatrix(3,1) = 0.0;
			//rotationMatrix(3,2) = 0.0;
			//rotationMatrix(3,3) = 1.0;

			return rotationMatrix;
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