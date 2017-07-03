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

		void Execute(ModelPart& rSubModelpart, Vector3 GlobalFiberDirection)
		{
			unsigned int caseID = 0; // 1=tri, 2=quad, 0=error
			// do element check here


			for (auto& element : rSubModelpart.Elements())
			{
				// check element type (1=tri, 2=quad, 0=error)
				checkElementType(element, caseID); //move this out of loop!
				
				// get global orientation of local axes
				// or get element global to local transformation matrix and apply to the fiber vector

				// compute rotation necessary to make loc2 and fiber vectors orthogonal

				// apply rotation through the shell cross section class.
					// probably need to make a public set method in each element

				// add option to write out angles so they don't have to be computed next time
					// or maybe this should be a separate python call

			}// sub-modelpart element loop
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
		void checkElementType(Element& element, unsigned int& caseID)
		{
			element.Initialize();
			Properties elementProps = element.GetProperties();
			std::string elementType = elementProps.GetValue(ELEMENT_TYPE);

			std::string thickTri = "ShellThickElement3D3N";
			std::string thinQuad = "ShellThinElement3D4N";

			if (elementType.compare(thickTri) == 0)
			{
				std::cout << "thick tri!" << std::endl;
				caseID = 1;
			}
			else if (elementType.compare(thinQuad) == 0)
			{
				std::cout << "thin quad!" << std::endl;
				caseID = 2;
			}
			else
			{
				KRATOS_ERROR <<
					"Error: Element type unsupported for use with CompositePropertyAssignment"
					<< std::endl;
			}
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
