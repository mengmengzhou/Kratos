//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//
//  Reference:       This class is adapted from applications/ChimeraApplication
// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "custom_python/add_custom_utilities_to_python.h"

//Processes
#include "custom_utilities/vtk_output.h"
//#include "custom_utilities/multipoint_constraint_data.hpp"

namespace Kratos
{

namespace Python
{


  void  AddCustomUtilitiesToPython()
  {
    using namespace boost::python;


      class_<VtkOutput, boost::noncopyable>("VtkOutput", init< ModelPart&, std::string, Parameters >())
      .def("PrintOutput", &VtkOutput::PrintOutput)
      .def("PrintOutput", &VtkOutput::PrintOutputSubModelPart);      


  }





}  // namespace Python.

} // Namespace Kratos