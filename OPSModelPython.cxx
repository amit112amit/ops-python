#include "pybind11/pybind11.h"
#include "OPSModel.h"
#include "LBFGSBWrapper.h"

namespace py = pybind11;
using namespace OPS;

PYBIND11_MODULE(opsmodel, m){
    m.doc() = "Oriented Particle System approach to modelling a point cloud closed shell.";
    py::class_<OPSModel>(m, "OPSModel")
	.def(py::init())
	.def("applyKabschAlgorithm", &OPSModel::applyKabschAlgorithm)
	.def("compute", &OPSModel::compute)
	.def("getAsphericity", &OPSModel::getAsphericity)
	.def("getArea", &OPSModel::getArea)
	.def("getAverageEdgeLength", &OPSModel::getAverageEdgeLength)
	.def("getAverageRadius", &OPSModel::getAverageRadius)
	.def("getCircularityEnergy", &OPSModel::getCircularityEnergy)
	.def("getMeanSquaredDisplacement", &OPSModel::getMeanSquaredDisplacement)
	.def("getMSD", &OPSModel::getMSD)
	.def("getMorseEnergy", &OPSModel::getMorseEnergy)
	.def("getNormalityEnergy", &OPSModel::getNormalityEnergy)
	.def("getTimeStep", &OPSModel::getTimeStep)
	.def("getVTKFileSuffix", &OPSModel::getVTKFileSuffix)
	.def("getRMSAngleDeficit", &OPSModel::getRMSAngleDeficit)
	.def("getTotalEnergy", &OPSModel::getTotalEnergy)
	.def("getVolume", &OPSModel::getVolume)
	.def("initializeFromVTKFile", &OPSModel::initializeFromVTKFile)
	.def("printVTKFile", &OPSModel::printVTKFile)
	.def("restoreSavedState", &OPSModel::restoreSavedState)
	.def("setMorseDistance", &OPSModel::setMorseDistance)
	.def("setMorseWellWidth", &OPSModel::setMorseWellWidth)
	.def("setFVK", &OPSModel::setFVK)
	.def("setTimeStep", &OPSModel::setTimeStep)
	.def("setVTKFileSuffix", &OPSModel::setVTKFileSuffix)
	.def("updateTriangles", &OPSModel::updateTriangles)
	.def("updatePreviousX", &OPSModel::updatePreviousX)
	.def("writeSimulationState", &OPSModel::writeSimulationState)
	.def("constraintSatisfied", &OPSModel::constraintSatisfied)
	.def("setConstraint", &OPSModel::setConstraint)
	.def("setLagrangeCoeff", &OPSModel::setLagrangeCoeff)
	.def("setPenaltyCoeff", &OPSModel::setPenaltyCoeff)
	.def("setTolerance", &OPSModel::setTolerance)
	.def("uzawaUpdate", &OPSModel::uzawaUpdate)
	.def("generateParallelKicks", &OPSModel::generateParallelKicks)
	.def("getBrownianEnergy", &OPSModel::getBrownianEnergy)
	.def("setBrownCoeff", &OPSModel::setBrownCoeff)
	.def("getViscosityEnergy", &OPSModel::getViscosityEnergy)
	.def("setViscosity", &OPSModel::setViscosity);

    py::class_<LBFGSBWrapper>(m, "Solver")
	.def(py::init<OPSModel&>())
	.def("setNumHessianCorrections", &LBFGSBWrapper::setNumHessianCorrections)
	.def("setPrintCode", &LBFGSBWrapper::setPrintCode)
	.def("setMaxIterations", &LBFGSBWrapper::setMaxIterations)
	.def("setMachineEPSFactor", &LBFGSBWrapper::setMachineEPSFactor)
	.def("setProjectedGradientTolerance", &LBFGSBWrapper::setProjectedGradientTolerance)
	.def("solve", &LBFGSBWrapper::solve);
}
