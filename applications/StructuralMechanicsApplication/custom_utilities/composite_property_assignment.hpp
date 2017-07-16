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
			Vector localGlobalFiberDirection;
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

			// Normalize
			GlobalFiberDirection /= std::sqrt(inner_prod(GlobalFiberDirection, GlobalFiberDirection));

			for (auto& element : rSubModelpart.Elements())
			{
				// get current element properties
				pElementProps = element.pGetProperties();

				// get local orientation of GlobalFiberDirection
				element.Calculate(LOCAL_ELEMENT_ORIENTATION, LCSOrientation, rCurrentProcessInfo);

				// get element local axis vectors
				Vector3 localAxis2, localAxis3;
				for (size_t i = 0; i < 3; i++)
				{
					localAxis1[i] = LCSOrientation(0, i);
					localAxis2[i] = LCSOrientation(1, i);
					localAxis3[i] = LCSOrientation(2, i);
				}


				std::cout << "localAxis1 still in global cartesian!" << std::endl;
				//localAxis1 = prod(LCSOrientation, localAxis1);
				localAxis1 /= std::sqrt(inner_prod(localAxis1, localAxis1));
				//std::cout << "flattened lc1 = " << localAxis1 << std::endl;

				std::cout << "localAxis2 still in global cartesian!" << std::endl;
				//localAxis2 = prod(LCSOrientation, localAxis2);
				localAxis2 /= std::sqrt(inner_prod(localAxis2, localAxis2));
				//std::cout << "flattened lc2 = " << localAxis2 << std::endl;

				std::cout << "localAxis3 still in global cartesian!" << std::endl;
				//localAxis3 = prod(LCSOrientation, localAxis3);
				localAxis3 /= std::sqrt(inner_prod(localAxis3, localAxis3));


				// Get local fiber direction in element space
				//Vector3 transformedGlobalFiber = ZeroVector(3);
				//computeTransformedFiber(LCSOrientation, GlobalFiberDirection, projectionDir, transformedGlobalFiber);
				//localGlobalFiberDirection = prod(LCSOrientation, transformedGlobalFiber);

				std::cout << "localGlobalFiberDirection still in global cartesian!" << std::endl;
				localGlobalFiberDirection = GlobalFiberDirection;
				//localGlobalFiberDirection = prod(LCSOrientation, GlobalFiberDirection);
				//localGlobalFiberDirection[2] = 0.0;
				localGlobalFiberDirection /= std::sqrt(inner_prod(localGlobalFiberDirection, localGlobalFiberDirection));

				// flatten localGlobalFiberDirection to make it in-plane with the element
				Vector Y_normal = Vector(3, 0.0);
				Y_normal[1] = 1.0;
				Vector rotation_axis = MathUtils<double>::CrossProduct(localAxis3, Y_normal);
				double rotation_angle = inner_prod(Y_normal, localAxis3);
				rotation_angle = std::acos(rotation_angle);
				
				//rotation_angle *= -1.0; // this rotates localGlobalFiberDirection to localAxis 3
				// therefore we need to rotate by 90deg - rotation angle
				
				//rotation_angle = KRATOS_M_PI / 2.0 - rotation_angle;
				
				Matrix R = setUpRotationMatrix(rotation_angle, rotation_axis[0], rotation_axis[1], rotation_axis[2]);


				localGlobalFiberDirection = prod(R, localGlobalFiberDirection);

				
				localAxis1 = prod(LCSOrientation, localAxis1);
				localAxis2 = prod(LCSOrientation, localAxis2);
				localGlobalFiberDirection = prod(LCSOrientation, localGlobalFiberDirection);

				// compute angle 'theta' between local axis 1 and localGlobalFiberDirection
				cosTheta = inner_prod(localAxis1, localGlobalFiberDirection);
				theta = std::acos(cosTheta);
				//std::cout << "angle in degrees: " << theta / 2.0 / 3.14 * 360 << std::endl;



				// dot between lc2 and localFiberDir
				double dotCheck = inner_prod(localAxis2, localGlobalFiberDirection);
				if (dotCheck < 0.0)
				{
					// theta is currently negative, flip to positive definition
					theta *= -1.0;
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

		Matrix setUpRotationMatrix(double angle, double u, double v, double w)
		{
			Matrix rotationMatrix(3, 3, 0.0);

			double L = (u*u + v * v + w * w);
			double u2 = u * u;
			double v2 = v * v;
			double w2 = w * w;

			rotationMatrix(0,0) = (u2 + (v2 + w2) * cos(angle)) / L;
			rotationMatrix(0,1) = (u * v * (1 - cos(angle)) - w * sqrt(L) * sin(angle)) / L;
			rotationMatrix(0,2) = (u * w * (1 - cos(angle)) + v * sqrt(L) * sin(angle)) / L;
			//rotationMatrix(0,3) = 0.0;
							  
			rotationMatrix(1,0) = (u * v * (1 - cos(angle)) + w * sqrt(L) * sin(angle)) / L;
			rotationMatrix(1,1) = (v2 + (u2 + w2) * cos(angle)) / L;
			rotationMatrix(1,2) = (v * w * (1 - cos(angle)) - u * sqrt(L) * sin(angle)) / L;
			//rotationMatrix(1,3) = 0.0;
							  
			rotationMatrix(2,0) = (u * w * (1 - cos(angle)) - v * sqrt(L) * sin(angle)) / L;
			rotationMatrix(2,1) = (v * w * (1 - cos(angle)) + u * sqrt(L) * sin(angle)) / L;
			rotationMatrix(2,2) = (w2 + (u2 + v2) * cos(angle)) / L;
			//rotationMatrix(2,3) = 0.0;
							  
			//rotationMatrix(3,0) = 0.0;
			//rotationMatrix(3,1) = 0.0;
			//rotationMatrix(3,2) = 0.0;
			//rotationMatrix(3,3) = 1.0;

			return rotationMatrix;
		}

		void computeTransformedFiber(const Matrix& elementOrientation, const Vector3& originalGlobalFiber, const Vector3& projectionVector, Vector3& transformedFiber)
		{
			Vector3 localAxis1, localAxis2, elementNormal;
			for (size_t i = 0; i < 3; i++)
			{
				localAxis1[i] = elementOrientation(0, i);
				localAxis2[i] = elementOrientation(1, i);
				elementNormal[i] = elementOrientation(2, i);
			}

			// Y axis projection
			// rotation about x and z to get dot prod = 1.0;

			// Y-X plane, rotation about Z
			// Z coord = 0!
			Vector3 XYplaneFiber = Vector3(originalGlobalFiber);
			XYplaneFiber[2] = 0.0;
			normalizeVector3(XYplaneFiber);
			double cosTheta = inner_prod(elementNormal, XYplaneFiber);
			std::cout << cosTheta << std::endl;
			double thetaZ = std::acos(cosTheta);
			// dot between lc1 and localFiberDir
			double dotCheck = inner_prod(localAxis1, XYplaneFiber);
			if (dotCheck < 0.0)
			{
				thetaZ *= -1.0;
			}


			// Y-Z plane, rotation about X
			// X coord = 0!
			Vector3 YZplaneFiber = Vector3(originalGlobalFiber);
			YZplaneFiber[0] = 0.0;
			normalizeVector3(YZplaneFiber);
			cosTheta = inner_prod(elementNormal, YZplaneFiber);
			double thetaX = std::acos(cosTheta);
			// dot between lc2 and localFiberDir
			dotCheck = inner_prod(localAxis2, YZplaneFiber);
			if (dotCheck < 0.0)
			{
				thetaX *= -1.0;
			}
			
			// Make rotation matrices
			Matrix zRot, xRot;
			makeXrotationMatrix(xRot, thetaX);
			makeZrotationMatrix(zRot, thetaZ);
			Matrix totalRotation = prod(xRot, zRot);

			// Give rotated global fiber orientation
			transformedFiber = prod(totalRotation, originalGlobalFiber);
		}

		void normalizeVector3(Vector3& rVector)
		{
			rVector /= std::sqrt(inner_prod(rVector, rVector));
		}

		void makeXrotationMatrix(Matrix& xRotation, double& angleX)
		{
			xRotation.resize(3, 3, 0.0);
			xRotation(0, 0) = 1.0;

			xRotation(1, 1) = std::cos(angleX);
			xRotation(1, 2) = -std::sin(angleX);
			xRotation(2, 1) = std::sin(angleX);
			xRotation(2, 2) = std::cos(angleX);
		}

		void makeYrotationMatrix(Matrix& yRotation, double& angleY)
		{
			yRotation.resize(3, 3, 0.0);
			
			yRotation(0, 0) = std::cos(angleY);
			yRotation(0, 2) = std::sin(angleY);

			yRotation(1, 1) = 1.0;

			yRotation(2, 0) = -std::sin(angleY);
			yRotation(2, 2) = std::cos(angleY);
		}

		void makeZrotationMatrix(Matrix& zRotation, double& angleZ)
		{
			zRotation.resize(3, 3, 0.0);
			

			zRotation(0, 0) = std::cos(angleZ);
			zRotation(0, 1) = -std::sin(angleZ);
			zRotation(1, 0) = std::sin(angleZ);
			zRotation(1, 1) = std::cos(angleZ);

			zRotation(2, 2) = 1.0;
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