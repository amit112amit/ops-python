#include "pybind11/pybind11.h"
#include "pybind11/eigen.h"
#include "LBFGSBWrapper.h"
#include "OPSModel.h"

namespace py = pybind11;

typedef Eigen::Matrix<double_t, Eigen::Dynamic, 3, Eigen::RowMajor> MatrixX3dR;
typedef Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor> MatrixX3iR;

//! Special module level function to expose to Python.
MatrixX3iR convexhull(Eigen::Ref<MatrixX3dR> positions)
{
    Surface_mesh sm;
    auto N = positions.rows();

    // Copy coordinates
    MatrixX3dR points(N, 3);
    points = positions;

    // Project points to unit sphere
    points.colwise().normalize();

    // Insert the projected points in a CGAL vertex_with_info vector
    std::vector<Point_with_info> verts;
    for (auto j = 0; j < N - 1; ++j)
    {
        verts.push_back(std::make_pair(K::Point_3(points(j, 0), points(j, 1), points(j, 2)), j));
    }

    // Find the spherical convex hull
    CGAL::convex_hull_3(verts.begin(), verts.end(), sm, OPS::CH_traits_for_point_with_info());
    auto origids = sm.points();

    MatrixX3iR triangles(sm.number_of_faces(), 3);

    //Iterate over all triangles related to the center vertex
    auto i = 0;
    for (const auto &f : sm.faces())
    {
        auto h = sm.halfedge(f);
        auto f0 = origids[sm.target(h)].second;
        auto f1 = origids[sm.target(sm.next(h))].second;
        auto f2 = origids[sm.target(sm.next(sm.next(h)))].second;
        triangles(i, 0) = f0;
        triangles(i, 1) = f1;
        triangles(i, 2) = f2;
        ++i;
    }

    return triangles;
}

PYBIND11_MODULE(ops, m)
{
    m.doc() = "Oriented Particle System approach for modelling closed shell of point cloud.";
    m.def("sphericalconvexhull", &convexhull, "Project to sphere and make convex hull.");

    py::class_<OPS::OPSModel>(m, "Model")
        .def(py::init<>())
        .def("applyKabschAlgorithm", &OPS::OPSModel::applyKabschAlgorithm)
        .def("bothMSD", &OPS::OPSModel::getMeanSquaredDisplacement)
        .def("compute", &OPS::OPSModel::compute)
        .def("constraintSatisfied", &OPS::OPSModel::constraintSatisfied)
        .def("generateParallelKicks", &OPS::OPSModel::generateParallelKicks)
        .def("initializeFromVTKFile", &OPS::OPSModel::initializeFromVTKFile)
        .def("printVTKFile", &OPS::OPSModel::printVTKFile)
        .def("polyDataParts", &OPS::OPSModel::polyDataParts, "Return coordinates, normals and triangles.",
             py::return_value_policy::move)
        .def("restoreSavedState", &OPS::OPSModel::restoreSavedState)
        .def("saveInitialPosition", &OPS::OPSModel::saveInitialPosition)
        .def("updateTriangles", &OPS::OPSModel::updateTriangles)
        .def("updatePreviousX", &OPS::OPSModel::updatePreviousX)
        .def("uzawaUpdate", &OPS::OPSModel::uzawaUpdate)
        .def("writeSimulationState", &OPS::OPSModel::writeSimulationState)
        .def_property("timestep", &OPS::OPSModel::getTimeStep, &OPS::OPSModel::setTimeStep)
        .def_property("vtkfilecount", &OPS::OPSModel::getVTKFileSuffix, &OPS::OPSModel::setVTKFileSuffix)
        .def_property("morseDistance", &OPS::OPSModel::getMorseDistance, &OPS::OPSModel::setMorseDistance)
        .def_property("morseWellWidth", &OPS::OPSModel::getMorseWellWidth, &OPS::OPSModel::setMorseWellWidth)
        .def_property("fvk", &OPS::OPSModel::getFVK, &OPS::OPSModel::setFVK)
        .def_property("constraint", &OPS::OPSModel::getConstraint, &OPS::OPSModel::setConstraint)
        .def_property("lagrangeCoeff", &OPS::OPSModel::getLagrangeCoeff, &OPS::OPSModel::setLagrangeCoeff)
        .def_property("penaltyCoeff", &OPS::OPSModel::getPenaltyCoeff, &OPS::OPSModel::setPenaltyCoeff)
        .def_property("tolerance", &OPS::OPSModel::getTolerance, &OPS::OPSModel::setTolerance)
        .def_property("brownCoeff", &OPS::OPSModel::getBrownCoeff, &OPS::OPSModel::setBrownCoeff)
        .def_property("viscosity", &OPS::OPSModel::getViscosity, &OPS::OPSModel::setViscosity)
        .def_property_readonly("asphericity", &OPS::OPSModel::getAsphericity)
        .def_property_readonly("area", &OPS::OPSModel::getArea)
        .def_property_readonly("averageEdgeLength", &OPS::OPSModel::getAverageEdgeLength)
        .def_property_readonly("averageRadius", &OPS::OPSModel::getAverageRadius)
        .def_property_readonly("brownianEnergy", &OPS::OPSModel::getBrownianEnergy)
        .def_property_readonly("circularityEnergy", &OPS::OPSModel::getCircularityEnergy)
        .def_property_readonly("msd", &OPS::OPSModel::getMSD)
        .def_property_readonly("morseEnergy", &OPS::OPSModel::getMorseEnergy)
        .def_property_readonly("normalityEnergy", &OPS::OPSModel::getNormalityEnergy)
        .def_property_readonly("rmsAngleDeficit", &OPS::OPSModel::getRMSAngleDeficit)
        .def_property_readonly("totalEnergy", &OPS::OPSModel::getTotalEnergy)
        .def_property_readonly("viscosityEnergy", &OPS::OPSModel::getViscosityEnergy)
        .def_property_readonly("volume", &OPS::OPSModel::getVolume);

    py::class_<OPS::LBFGSBWrapper>(m, "Solver")
        .def(py::init<OPS::OPSModel &, int, int, int, double, double>(), py::arg("model"),
             py::arg("m") = 5, py::arg("iprint") = 1000,
             py::arg("maxiter") = 100000, py::arg("factr") = 10.0,
             py::arg("pgtol") = 1e-8)
        .def("solve", &OPS::LBFGSBWrapper::solve)
        .def_property("m", &OPS::LBFGSBWrapper::getNumHessianCorrections,
                      &OPS::LBFGSBWrapper::setNumHessianCorrections)
        .def_property("iprint", &OPS::LBFGSBWrapper::getPrintCode, &OPS::LBFGSBWrapper::setPrintCode)
        .def_property("maxiter", &OPS::LBFGSBWrapper::getMaxIterations, &OPS::LBFGSBWrapper::setMaxIterations)
        .def_property("factr", &OPS::LBFGSBWrapper::getMachineEPSFactor, &OPS::LBFGSBWrapper::setMachineEPSFactor)
        .def_property("pgtol", &OPS::LBFGSBWrapper::getProjectedGradientTolerance,
                      &OPS::LBFGSBWrapper::setProjectedGradientTolerance);
}
