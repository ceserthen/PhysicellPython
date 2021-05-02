
//#include "../../../pybind11/include/pybind11/pybind11.h"
#include <pybind11/pybind11.h>
#include "../../../src/BioFVM/BioFVM.h"

namespace py = pybind11;

PYBIND11_MODULE(biofvmpy, m) 
{
    py::class_<BioFVM::Microenvironment>(m, "Microenvironment")
    .def(py::init<std::string>())
    .def(py::init<>())
    .def("simulate_diffusion_decay", static_cast<void (BioFVM::Microenvironment::*)(double)>(&BioFVM::Microenvironment::simulate_diffusion_decay), "advance the diffusion-decay solver by dt time");
}
