#if !defined(__OPSMODEL_H__)
#define __OPSMODEL_H__

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <array>
#include <tuple>
#include <math.h>
#include <random>
#include <set>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <vtkDoubleArray.h>
#include <vtkIdFilter.h>
#include <vtkIdList.h>
#include <vtkOctreePointLocator.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkDelaunay3D.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkCellArray.h>
#include <vtkExtractEdges.h>

namespace OPS {

    class OPSModel{
        public:
            typedef std::normal_distribution<double_t> NormD;
            typedef Eigen::Vector3d Vector3d;
            typedef Eigen::Map<Vector3d> MapV3d;
            typedef Eigen::Ref<const Vector3d> RefCV3d;
            typedef Eigen::VectorXd VectorXd;
            typedef Eigen::Map<VectorXd> MapVXd;
            typedef Eigen::Ref<VectorXd> RefVXd;
            typedef Eigen::Matrix3d Matrix3d;
            typedef Eigen::Ref<Matrix3d> RefM3d;
            typedef Eigen::Matrix3Xd Matrix3Xd;
            typedef Eigen::Ref<Matrix3Xd> RefM3Xd;
            typedef Eigen::Map<Matrix3Xd> MapM3Xd;
            typedef Eigen::Matrix<double_t, 3, 4> Matrix3x4d;
            typedef Eigen::Matrix<double_t, 4, 3> Matrix4x3d;
            typedef Eigen::Quaterniond Quaterniond;
            typedef Eigen::AngleAxisd AngleAxisd;
            typedef Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> MatrixX3fR;
            typedef Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor> MatrixX3iR;
            typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
            typedef K::Point_2 Point;
            typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K> Vb;
            typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
            typedef CGAL::Delaunay_triangulation_2<K, Tds> Delaunay;
            typedef Delaunay::Face_circulator Face_circulator;

            OPSModel(){};
            void applyKabschAlgorithm();
            void compute();
            void computeNormals();
            void computeNormals(const VectorXd &rv);
            void diffNormalRotVec();
            void diffNormalRotVec(const VectorXd &rv);
            Eigen::Affine3d find3DAffineTransform(Eigen::Ref<Eigen::Matrix3Xd>,
                    Eigen::Ref<Eigen::Matrix3Xd>);
            double_t getAsphericity();
            double_t getArea();
            double_t getAverageEdgeLength();
            double_t getAverageRadius();
            inline double_t getCircularityEnergy() { return _circEn; }
            inline void getInitialPositions(Matrix3Xd &v) { v = _initialPositions; }
            std::array<double_t, 2> getMeanSquaredDisplacement();
            double_t getMSD();
            inline std::vector<size_t> getInitialNeighbors() { return _initialNearestNeighbor; }
            inline double_t getMorseEnergy() { return _morseEn; }
            inline double_t getNormalityEnergy() { return _normalEn; }
            inline size_t getNumberOfPoints() { return _N; }
            inline size_t getTimeStep() { return _timestep; }
            inline size_t getVTKFileSuffix() { return _vtkfilesuffix; }
            double_t getPointCloudAvgEdgeLen(std::string file);
            double_t getRMSAngleDeficit();
            inline double_t getTotalEnergy() { return (_morseEn + _normalEn + _circEn); }
            double_t getVolume();
            double_t* getF(){ return &_f; }
            double_t* getG(){ return _g.data(); }
            double_t* getX(){ return _x.data(); }
            void initializeFromVTKFile(std::string vtkfile);
            static void initialRotationVector(RefM3Xd pos, RefM3Xd rotVec);
            void printVTKFile(const std::string name);
            std::tuple<MatrixX3fR, MatrixX3fR, MatrixX3iR> polyDataParts();
            void restoreSavedState(std::string stateFile);
            inline void saveInitialPosition() {
                MapM3Xd positions(_x.data(), 3, _N);
                _initialPositions = positions; }
            inline void setInitialNeighbors(std::vector<size_t> &x) {
                _initialNearestNeighbor = x;
            }
            inline void setInitialPositions(Matrix3Xd p) { _initialPositions = p; }
            inline void setMorseDistance(double_t r) { _re = r; }
            inline double_t getMorseDistance() { return _re; }
            inline void setMorseWellWidth(double_t a) { _a = a; }
            inline double_t getMorseWellWidth() { return _a; }
            inline void setFVK(double_t g) { _gamma_inv = 1 / g; }
            inline double_t getFVK() { return 1.0/_gamma_inv; }
            inline void setTimeStep(size_t t) { _timestep = t; }
            inline void setVTKFileSuffix(size_t s) { _vtkfilesuffix = s; }
            void stereoDelaunay();
            void updateTriangles();
            inline void updatePreviousX(){ _pX = _x.head(3 * _N); }
            void writeSimulationState(std::string file);
            // **************** Area constraint methods ***********************//
            inline bool constraintSatisfied() {
                return (std::abs(_value - _constrainedValue) < _tolerance);
            }
            inline void setConstraint(double_t V) { _constrainedValue = V; }
            inline double_t getConstraint() { return _constrainedValue; }
            inline void setLagrangeCoeff(double_t L) { _Lambda_i = L; }
            inline double_t getLagrangeCoeff() { return _Lambda_i; }
            inline void setPenaltyCoeff(double_t K) { _K_i = K; }
            inline double_t getPenaltyCoeff() { return _K_i; }
            inline void setTolerance(double_t t) { _tolerance = t; }
            inline double_t getTolerance() { return _tolerance; }
            inline void uzawaUpdate() {
                _Lambda_i = _Lambda_i - _K_i * (_value - _constrainedValue);
                _K_i *= 10;
            }
            //**************** Brownian Body methods ******************//
            inline void generateParallelKicks() {
                for (auto i = 0; i < 3 * _N; i++) {
                    _xi(i) = _rng(_e2);
                }
            }
            inline std::mt19937 getRandomEngine() { return _e2; }
            inline NormD getRandomGenerator() { return _rng; }
            inline double_t getBrownianEnergy() { return _brownEn; }
            inline void setBrownCoeff(double_t C) { _brownCoeff = C; }
            inline double_t getBrownCoeff() { return _brownCoeff; }
            inline void setRandomEngine(std::mt19937 e) { _e2 = e; }
            inline void setRandomGenerator(NormD r) { _rng = r; }
            // ******************* Viscosity Body methods *******************//
            inline double_t getViscosityEnergy() { return _viscoEn; }
            inline VectorXd getPreviousX() { return _pX; }
            inline void setViscosity(double_t v) { _viscosity = v; }
            inline double_t getViscosity() { return _viscosity; }

        protected:
            bool _updateArea = true, _updateRadius = true, _updateVolume = true;
            double_t _area = 0, _radius = 0, _volume = 0, _morseEn, _normalEn, _circEn,
                     _re = 1.0, _a = 4.62098120373, _gamma_inv = 10.0, _msd, _msd_tgt;
            double_t _brownCoeff = 0.0, _brownEn = 0.0;
            double_t _constrainedValue = 0.0, _K_i = 1000.0, _Lambda_i = 1.0, _value = 0.0;
            double_t _f;
            double_t _tolerance = 1e-10;
            double_t _viscosity = 0.0, _viscoEn = 0.0;
            Delaunay _dt;
            Matrix3Xd _normals, _initialPositions;
            NormD _rng;
            size_t _N, _timestep = 0, _vtkfilesuffix = 0;
            std::vector<std::array<unsigned, 2>> _edges;
            std::vector<std::array<unsigned, 3>> _triangles;
            std::vector<Matrix3d, Eigen::aligned_allocator<Matrix3d>> _diffNormalRV;
            std::vector<size_t> _initialNearestNeighbor;
            std::mt19937 _e2;
            VectorXd _x, _g, _xi, _pX;

            void updateRotationVectors();
    };

} // namespace OPS

#endif //__OPSMODEL_H__
