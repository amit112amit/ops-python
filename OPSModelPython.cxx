#include "pybind11/pybind11.h"
#include "pybind11/eigen.h"
#include "OPSModel.h"
#include "LBFGSBWrapper.h"

namespace py = pybind11;
using namespace OPS;

PYBIND11_MODULE(opsmodel, m){
    m.doc() = "Oriented Particle System approach to modelling a point cloud closed shell.";
    py::class_<OPSModel>(m, "OPSModel")
	.def(py::init<>())
	.def("applyKabschAlgorithm", &OPSModel::applyKabschAlgorithm)
	.def("bothMSD", &OPSModel::getMeanSquaredDisplacement)
	.def("compute", &OPSModel::compute)
	.def("constraintSatisfied", &OPSModel::constraintSatisfied)
	.def("generateParallelKicks", &OPSModel::generateParallelKicks)
	.def("initializeFromVTKFile", &OPSModel::initializeFromVTKFile)
	.def("printVTKFile", &OPSModel::printVTKFile)
	.def("polyDataParts", &OPSModel::polyDataParts, "Return coordinates, normals and triangles.",
            py::return_value_policy::move)
	.def("restoreSavedState", &OPSModel::restoreSavedState)
	.def("saveInitialPosition", &OPSModel::saveInitialPosition)
	.def("updateTriangles", &OPSModel::updateTriangles)
	.def("updatePreviousX", &OPSModel::updatePreviousX)
	.def("uzawaUpdate", &OPSModel::uzawaUpdate)
	.def("writeSimulationState", &OPSModel::writeSimulationState)
	.def_property("timestep", &OPSModel::getTimeStep, &OPSModel::setTimeStep)
	.def_property("vtkfilecount", &OPSModel::getVTKFileSuffix, &OPSModel::setVTKFileSuffix)
	.def_property("morseDistance", &OPSModel::getMorseDistance, &OPSModel::setMorseDistance)
	.def_property("morseWellWidth", &OPSModel::getMorseWellWidth, &OPSModel::setMorseWellWidth)
	.def_property("fvk", &OPSModel::getFVK, &OPSModel::setFVK)
	.def_property("constraint", &OPSModel::getConstraint, &OPSModel::setConstraint)
	.def_property("lagrangeCoeff", &OPSModel::getLagrangeCoeff, &OPSModel::setLagrangeCoeff)
	.def_property("penaltyCoeff", &OPSModel::getPenaltyCoeff, &OPSModel::setPenaltyCoeff)
	.def_property("tolerance", &OPSModel::getTolerance, &OPSModel::setTolerance)
	.def_property("brownCoeff", &OPSModel::getBrownCoeff, &OPSModel::setBrownCoeff)
	.def_property("viscosity", &OPSModel::getViscosity, &OPSModel::setViscosity)
	.def_property_readonly("asphericity", &OPSModel::getAsphericity)
	.def_property_readonly("area", &OPSModel::getArea)
	.def_property_readonly("averageEdgeLength", &OPSModel::getAverageEdgeLength)
	.def_property_readonly("averageRadius", &OPSModel::getAverageRadius)
	.def_property_readonly("brownianEnergy", &OPSModel::getBrownianEnergy)
	.def_property_readonly("circularityEnergy", &OPSModel::getCircularityEnergy)
	.def_property_readonly("msd", &OPSModel::getMSD)
	.def_property_readonly("morseEnergy", &OPSModel::getMorseEnergy)
	.def_property_readonly("normalityEnergy", &OPSModel::getNormalityEnergy)
	.def_property_readonly("rmsAngleDeficit", &OPSModel::getRMSAngleDeficit)
	.def_property_readonly("totalEnergy", &OPSModel::getTotalEnergy)
	.def_property_readonly("viscosityEnergy", &OPSModel::getViscosityEnergy)
	.def_property_readonly("volume", &OPSModel::getVolume);

    py::class_<LBFGSBWrapper>(m, "Solver")
	.def(py::init<OPSModel &, int, int, int, double, double>(), py::arg("model"),
            py::arg("m") = 5, py::arg("iprint") = 1000,
            py::arg("maxiter") = 100000, py::arg("factr") = 10.0,
            py::arg("pgtol") = 1e-8)
	.def("solve", &LBFGSBWrapper::solve)
    .def_property("m", &LBFGSBWrapper::getNumHessianCorrections,
            &LBFGSBWrapper::setNumHessianCorrections)
    .def_property("iprint", &LBFGSBWrapper::getPrintCode, &LBFGSBWrapper::setPrintCode)
    .def_property("maxiter", &LBFGSBWrapper::getMaxIterations, &LBFGSBWrapper::setMaxIterations)
    .def_property("factr", &LBFGSBWrapper::getMachineEPSFactor, &LBFGSBWrapper::setMachineEPSFactor)
    .def_property("pgtol", &LBFGSBWrapper::getProjectedGradientTolerance,
            &LBFGSBWrapper::setProjectedGradientTolerance);
}
