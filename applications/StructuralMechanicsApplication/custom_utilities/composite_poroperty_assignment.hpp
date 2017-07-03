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
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"


namespace Kratos
{
  ///@addtogroup ApplicationNameApplication
  ///@{

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
  /** Detail class definition.
  */
  class CompositePropertyAssignment
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of CompositePropertyAssignment
      KRATOS_CLASS_POINTER_DEFINITION(CompositePropertyAssignment);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      CompositePropertyAssignment();

      /// Destructor.
      virtual ~CompositePropertyAssignment();


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      void Execute(ModelPart& rSubModelpart, Vector3 GlobalFiberDirection) {
          for (auto& element : rSubModelpart.Elements())
          {
              element.DoWhaterver();
          }





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
      virtual std::string Info() const;

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const;

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const;


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
      CompositePropertyAssignment& operator=(CompositePropertyAssignment const& rOther);

      /// Copy constructor.
      CompositePropertyAssignment(CompositePropertyAssignment const& rOther);


      ///@}

    }; // Class CompositePropertyAssignment

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    CompositePropertyAssignment& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const CompositePropertyAssignment& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_COMPOSITE_PROPERTY_ASSIGNMENT_HPP_INCLUDED  defined
