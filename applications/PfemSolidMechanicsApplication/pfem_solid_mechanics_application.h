//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_PFEM_SOLID_MECHANICS_APPLICATION_H_INCLUDED )
#define  KRATOS_PFEM_SOLID_MECHANICS_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/constitutive_law.h"

// Core applications
#include "../SolidMechanicsApplication/solid_mechanics_application.h"

//conditions
#include "custom_conditions/skin_multiple_condition.hpp"
#include "custom_conditions/wall_tip_condition.hpp"

#include "custom_conditions/contact_domain_2D_condition.hpp"
#include "custom_conditions/axisym_contact_domain_2D_condition.hpp"

//elements

//constitutive laws
#include "containers/flags.h"
#include "includes/variables.h"


namespace Kratos
{
  ///@name Type	Definitions
  ///@{

  ///@name Kratos Globals
  ///@{ 


  //Define Variables

  //solution
  KRATOS_DEFINE_VARIABLE(int, NUMBER_OF_ACTIVE_CONTACTS );
  KRATOS_DEFINE_VARIABLE(int, NUMBER_OF_STICK_CONTACTS );
  KRATOS_DEFINE_VARIABLE(int, NUMBER_OF_SLIP_CONTACTS );

  //geometrical

  //constitutive law   
  KRATOS_DEFINE_VARIABLE(double, MEAN_ERROR );

  //element

  //thermal

  //mechanical
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( FORCE_CONTACT_NORMAL );
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( FORCE_CONTACT_TANGENT );

  //nodal dofs
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( OFFSET );
  KRATOS_DEFINE_VARIABLE(Vector, BOUNDARY_NORMAL );
  KRATOS_DEFINE_VARIABLE(double, SHRINK_FACTOR );

  //domain definition
  KRATOS_DEFINE_VARIABLE(unsigned int, DOMAIN_LABEL );
  KRATOS_DEFINE_VARIABLE(bool        , RIGID_WALL );
  KRATOS_DEFINE_VARIABLE(double      , WALL_TIP_RADIUS );
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( WALL_REFERENCE_POINT );
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( WALL_VELOCITY );

  //contact condition
  KRATOS_DEFINE_VARIABLE(Condition::Pointer, MASTER_CONDITION );
  KRATOS_DEFINE_VARIABLE(WeakPointerVector< Element >, MASTER_ELEMENTS );
  KRATOS_DEFINE_VARIABLE(WeakPointerVector< Node<3> >, MASTER_NODES );

  //properties
  KRATOS_DEFINE_VARIABLE(bool, PENALTY_CONTACT );
  KRATOS_DEFINE_VARIABLE(bool, FRICTION_ACTIVE );
  KRATOS_DEFINE_VARIABLE(double, TAU_STAB );
  KRATOS_DEFINE_VARIABLE(double, MU_STATIC );
  KRATOS_DEFINE_VARIABLE(double, MU_DYNAMIC );

  //flags
  /* KRATOS_DEFINE_FLAG( FLUID ); */
  /* KRATOS_DEFINE_FLAG( STRUCTURE ); */
  /* KRATOS_DEFINE_FLAG( SOLID ); */
  /* KRATOS_DEFINE_FLAG( RIGID ); */
  /* KRATOS_DEFINE_FLAG( CONTACT ); */
  
  /* KRATOS_DEFINE_FLAG( BOUNDARY ); */
  /* KRATOS_DEFINE_FLAG( FREE_SURFACE ); */
  
  /* KRATOS_DEFINE_FLAG( INTERFACE ); */
  
  /* KRATOS_DEFINE_FLAG( ENGAGED ); */
  /* KRATOS_DEFINE_FLAG( ISOLATED ); */
  
  /* KRATOS_DEFINE_FLAG( REFINE ); */
  /* KRATOS_DEFINE_FLAG( INSERTED ); */
  /* KRATOS_DEFINE_FLAG( RELEASE ); */

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
  class KratosPfemSolidMechanicsApplication : public KratosApplication
  {
  public:


    ///@name Type Definitions
    ///@{
		

    /// Pointer definition of KratosPfemSolidMechanicsApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosPfemSolidMechanicsApplication);


    ///@}
    ///@name Life Cycle 
    ///@{ 

    /// Default constructor.
    KratosPfemSolidMechanicsApplication();

    /// Destructor.
    virtual ~KratosPfemSolidMechanicsApplication(){}


    ///@}
    ///@name Operators 
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual void Register();



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
    virtual std::string Info() const
      {
	return "KratosPfemSolidMechanicsApplication";
      }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
      rOStream << Info();
      PrintData(rOStream);
    }

    ///// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
      KRATOS_WATCH("in my application");
      KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
      rOStream << "Variables:" << std::endl;
      KratosComponents<VariableData>().PrintData(rOStream);
      rOStream << std::endl;
      rOStream << "Elements:" << std::endl;
      KratosComponents<Element>().PrintData(rOStream);
      rOStream << std::endl;
      rOStream << "Conditions:" << std::endl;
      KratosComponents<Condition>().PrintData(rOStream);
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



    //       static const ApplicationCondition  msApplicationCondition; 

    ///@} 
    ///@name Member Variables 
    ///@{ 
    const Condition mCondition2D;
    const Condition mCondition3D;

    const SkinMultipleCondition mSkinMultipleCondition2D;
    const SkinMultipleCondition mSkinMultipleCondition3D;

    const WallTipCondition mWallTipCondition2D;
    const WallTipCondition mWallTipCondition3D;

    const ContactDomain2DCondition   mContactDomain2DCondition;

    const AxisymContactDomain2DCondition    mAxisymContactDomain2DCondition;
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
    KratosPfemSolidMechanicsApplication& operator=(KratosPfemSolidMechanicsApplication const& rOther);

    /// Copy constructor.
    KratosPfemSolidMechanicsApplication(KratosPfemSolidMechanicsApplication const& rOther);


    ///@}    

  }; // Class KratosPfemSolidMechanicsApplication 

  ///@} 


  ///@name Type Definitions       
  ///@{ 


  ///@} 
  ///@name Input and output 
  ///@{ 

  ///@} 


}  // namespace Kratos.

#endif // KRATOS_PFEM_SOLID_APPLICATION_H_INCLUDED  defined 


