#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "spaces/ublas_space.h"

#include "custom_python/add_custom_processes_to_python.h"

#include "processes/find_elements_neighbours_process.h"
#include "includes/model_part.h"

namespace Kratos
{

	namespace Python
	{
		

		void AddCustomProcessesToPython()
		{

			using namespace boost::python;

			class_<FindElementalNeighboursProcess, bases<Process>, boost::noncopyable >
				("FindElementalNeighboursProcess", init<ModelPart&, int, unsigned int>())
				.def("Execute", &FindElementalNeighboursProcess::Execute)
				;



		}

	}  // namespace Python.

} // Namespace Kratos