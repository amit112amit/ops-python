#include "OPSModel.h"

typedef Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> MatrixX3fR;
typedef Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor> MatrixX3iR;

namespace OPS {

    //! Constructor for BrownOPS
    void OPSModel::initializeFromVTKFile(std::string inFile){

        // Read point positions from VTK File
        vtkNew<vtkPolyDataReader> reader;
        reader->SetFileName(inFile.c_str());
        reader->ReadAllVectorsOn();
        reader->Update();
        auto mesh = reader->GetOutput();
        _N = mesh->GetNumberOfPoints();

        // Initialize matrices and vectors
        _x = VectorXd::Zero(6 * _N);
        _g = VectorXd::Zero(6 * _N);
        _diffNormalRV = std::vector<Matrix3d, Eigen::aligned_allocator<Matrix3d>>(
                _N, Matrix3d::Zero());
        _pX = VectorXd::Zero(3 * _N);
        _normals = Matrix3Xd::Zero(3, _N);

        // Read point coordinates from input mesh
        for (auto i = 0; i < _N; ++i)
            mesh->GetPoint(i, &_x(3 * i));

        // Renormalize by the average edge length
        auto L = getPointCloudAvgEdgeLen(inFile);
        _x /= getPointCloudAvgEdgeLen(inFile);
        updatePreviousX();

        // Generate rotation vectors from input point coordinates
        Eigen::Map<Eigen::Matrix3Xd> xpos(_x.data(), 3, _N), xrot(&_x(3 * _N), 3, _N);
        initialRotationVector(xpos, xrot);
        _initialPositions = xpos;

        // Construct triangles and edges
        updateTriangles();

        // Set the initial nearest neighbor map
        void *coords = (void *)_x.data();
        vtkNew<vtkDoubleArray> pointCoords;
        pointCoords->SetVoidArray(coords, 3 * _N, 1);
        pointCoords->SetNumberOfComponents(3);

        vtkNew<vtkPoints> points;
        points->SetData(pointCoords);

        // Construct vtkPolyData
        vtkNew<vtkPolyData> polyData;
        polyData->SetPoints(points);

        // Construct the kd-tree
        vtkNew<vtkOctreePointLocator> octree;
        octree->SetDataSet(polyData);
        octree->BuildLocator();

        for (auto i = 0; i < _N; i++) {
            vtkNew<vtkIdList> neighbors;
            octree->FindClosestNPoints(2, &_x(3 * i), neighbors);
            neighbors->DeleteId(i);
            _initialNearestNeighbor.push_back(neighbors->GetId(0));
        }

        // Brownian body initialization
        std::random_device rd;
        _e2 = std::mt19937(rd());
        _rng = NormD(0., 1.);
        _xi = VectorXd::Zero(3 * _N);
    }


    //! Initialize from saved state file
    void OPSModel::restoreSavedState(std::string stateFile){
        // Read previously stored state to resume the simulation
        ifstream f(stateFile);
        std::string line;
        double_t ignore;
        size_t numsPerLine = 9;

        // Read number of particles
        std::getline(f, line);
        std::getline(f, line);
        std::getline(f, line);
        _N = std::stod(line);

        // Initialize matrices and vectors
        _x = VectorXd::Zero(6 * _N);
        _g = VectorXd::Zero(6 * _N);
        _pX = VectorXd::Zero(3 * _N);
        _diffNormalRV = std::vector<Matrix3d, Eigen::aligned_allocator<Matrix3d>>(
                _N, Matrix3d::Zero());
        _normals = Matrix3Xd::Zero(3, _N);
        _initialPositions = Matrix3Xd::Zero(3, _N);
        _xi = VectorXd::Zero(3 * _N);

        // Read nameSuffix
        std::getline(f, line);
        std::getline(f, line);
        _vtkfilesuffix = std::stod(line);
        _vtkfilesuffix = _vtkfilesuffix > 0 ? _vtkfilesuffix + 1 : 0;
        // Read step number
        std::getline(f, line);
        std::getline(f, line);
        _timestep = std::stod(line) + 1;
        // Read gamma or the FvK number
        std::getline(f, line);
        std::getline(f, line);
        ignore = std::stod(line);
        // Read beta or 1/temperature
        std::getline(f, line);
        std::getline(f, line);
        ignore = std::stod(line);

        // Read random engine state
        std::getline(f, line);
        std::getline(f, line);
        std::istringstream ss(line);
        ss >> _e2;

        // Read random generator state
        std::getline(f, line);
        _rng = NormD(0., 1.);
        f >> _rng;

        // Read neighbors
        size_t id;
        std::getline(f, line);
        std::getline(f, line);
        size_t counter = 0;
        while (counter < _N && std::getline(f, line)) {
            std::istringstream ss(line);
            for (auto i = 0; i < numsPerLine; ++i) {
                ss >> id;
                _initialNearestNeighbor.push_back(id);
                counter++;
            }
        }

        // Read initial position of particles
        double_t xi;
        std::getline(f, line);
        counter = 0;
        while (counter < 3 * _N && std::getline(f, line)) {
            std::istringstream ss(line);
            for (auto i = 0; i < numsPerLine; ++i) {
                ss >> xi;
                _initialPositions(counter++) = xi;
            }
        }

        // Read prevX vector
        std::getline(f, line);
        counter = 0;
        while (counter < 3 * _N && std::getline(f, line)) {
            std::istringstream ss(line);
            for (auto i = 0; i < numsPerLine; ++i) {
                ss >> xi;
                _pX(counter++) = xi;
            }
        }

        // Read x vector
        std::getline(f, line);
        counter = 0;
        while (counter < 6 * _N && std::getline(f, line)) {
            std::istringstream ss(line);
            for (auto i = 0; i < numsPerLine; ++i) {
                ss >> xi;
                _x(counter++) = xi;
            }
        }
        f.close();

        // Calculate the normals based on rotation vectors
        computeNormals();

        // Update the triangle data structure
        updateTriangles();
    }

    //! Function to convert from rotation vectors to normals
    void OPSModel::computeNormals() {
        MapM3Xd rotationVectors(&_x(3 * _N), 3, _N);
        // Assume z-axis of Global Coord Sys is to be rotated
        Quaterniond zaxis(0.0, 0.0, 0.0, 1.0);
        for (auto i = 0; i < _N; ++i) {
            MapV3d normal(&_normals(0, i), 3, 1);
            MapV3d rotVec(&rotationVectors(0, i), 3, 1);
            Vector3d curr =
                (Quaterniond(AngleAxisd(rotVec.norm(), rotVec.normalized())) * zaxis *
                 (Quaterniond(AngleAxisd(rotVec.norm(), rotVec.normalized()))
                  .conjugate()))
                .vec();
            normal << curr(0), curr(1), curr(2);
        }
    }

    //! Function to convert from rotation vectors to normals
    void OPSModel::computeNormals(const VectorXd &x) {
        // Assume z-axis of Global Coord Sys is to be rotated
        Quaterniond zaxis(0.0, 0.0, 0.0, 1.0);
        for (auto i = 0; i < _N; ++i) {
            MapV3d normal(&_normals(0, i), 3, 1);
            Vector3d curr =
                (Quaterniond(AngleAxisd(x.segment<3>(3 * (_N + i)).norm(),
                                        x.segment<3>(3 * (_N + i)).normalized())) *
                 zaxis *
                 (Quaterniond(AngleAxisd(x.segment<3>(3 * (_N + i)).norm(),
                                         x.segment<3>(3 * (_N + i)).normalized()))
                  .conjugate()))
                .vec();
            normal << curr(0), curr(1), curr(2);
        }
    }

    //! Print a VTK file
    void OPSModel::printVTKFile(const std::string name) {
        vtkNew<vtkCellArray> triangles;
        for (const auto &f : _triangles) {
            triangles->InsertNextCell(3);
            for (auto j = 2; j >= 0; --j)
                triangles->InsertCellPoint(f[j]);
        }
        // Extract point coordinates for _polyData from x
        vtkNew<vtkDoubleArray> pointCoords;
        pointCoords->SetVoidArray((void *)_x.data(), 3 * _N, 1);
        pointCoords->SetNumberOfComponents(3);

        vtkNew<vtkPoints> points;
        points->SetData(pointCoords);

        // Convert rotation vectors to point normals
        vtkNew<vtkDoubleArray> pointNormals;
        pointNormals->SetName("PointNormals");
        pointNormals->SetVoidArray((void *)_normals.data(), 3 * _N, 1);
        pointNormals->SetNumberOfComponents(3);

        // Construct vtkPolyData
        vtkNew<vtkPolyData> polyData;
        polyData->SetPoints(points);
        polyData->GetPointData()->SetNormals(pointNormals);
        polyData->SetPolys(triangles);

        vtkNew<vtkPolyDataWriter> writer;
        writer->SetFileName(name.c_str());
        writer->SetInputData(polyData);
        writer->Write();
    }

    //! Returns coordinates, normals and triangles of a polydata
    std::tuple<MatrixX3fR, MatrixX3fR, MatrixX3iR> OPSModel::polyDataParts(){
        MatrixX3fR coords(_N, 3), normals(_N, 3);
        MatrixX3iR triangles(_triangles.size(), 3);

        for (size_t i = 0; i < _triangles.size(); ++i) {
            triangles(i, 2) = _triangles[i][0];
            triangles(i, 1) = _triangles[i][1];
            triangles(i, 0) = _triangles[i][2];
        }

        Eigen::Map<Eigen::VectorXf> coordsvec(coords.data(), 3*_N);
        Eigen::Map<Eigen::VectorXf> normalvec(normals.data(), 3*_N);
        Eigen::Map<Eigen::VectorXd> _normalvec(_normals.data(), 3*_N);
        coordsvec = _x.cast<float>();
        normalvec = _normalvec.cast<float>();

        return std::make_tuple(coords, normals, triangles);
    }

    //! Calculate average edge length as if the particles were on a mesh
    double_t OPSModel::getAverageEdgeLength() {
        MapM3Xd positions(_x.data(), 3, _N);
        double_t avg = 0;
        for (const auto &e : _edges)
            avg += (positions.col(e[1]) - positions.col(e[0])).norm();
        return avg / _edges.size();
    }

    //! Compute the OPSBody energy
    void OPSModel::compute() {

        MapM3Xd positions(_x.data(), 3, _N);
        MapM3Xd posGradient(_g.data(), 3, _N);
        MapM3Xd rotGradient(&_g(3 * _N), 3, _N);

        // Initialize energies and forces to be zero
        _morseEn = 0.0;
        _normalEn = 0.0;
        _circEn = 0.0;
        _f = 0.0;
        _g.setZero(6 * _N);

        computeNormals();
        diffNormalRotVec();

        for (const auto &e : _edges) {
            // Evaluate morse derivatives
            Vector3d rn = (positions.col(e[1]) - positions.col(e[0])).normalized();
            double_t r = (positions.col(e[1]) - positions.col(e[0])).norm();
            double_t exp_1 = exp(-_a * (r - _re));
            double_t exp_2 = exp_1 * exp_1;
            double_t morseEn = exp_2 - 2 * exp_1;
            Vector3d dMdr = 2 * _a * (exp_1 - exp_2) * rn;

            // Evaluate co-normality derivatives
            Vector3d m = _normals.col(e[0]) - _normals.col(e[1]);
            double_t Phi_n = m.squaredNorm();
            Matrix3d M = _diffNormalRV[e[0]];
            Vector3d dPhi_nVi = 2 * M * m;
            Matrix3d N = _diffNormalRV[e[1]];
            Vector3d dPhi_nVj = -2 * N * m;

            // Evaluate co-circularity derivatives
            Vector3d n = _normals.col(e[0]) + _normals.col(e[1]);
            double_t n_dot_rn = n.dot(rn);
            double_t Phi_c = n_dot_rn * n_dot_rn;
            Vector3d dCdr = (2 * n_dot_rn / r) * (n - n_dot_rn * rn);
            Vector3d dPhi_cVi = (2 * n_dot_rn) * M * rn;
            Vector3d dPhi_cVj = (2 * n_dot_rn) * N * rn;

            // Update the energies
            _morseEn += morseEn;
            _normalEn += Phi_n * _gamma_inv;
            _circEn += Phi_c * _gamma_inv;
            _f += morseEn + (Phi_n + Phi_c) * _gamma_inv;

            // Calculate the total derivatives of energy wrt xi, vi and vj
            posGradient.col(e[0]) -= dMdr + dCdr * _gamma_inv;
            posGradient.col(e[1]) += dMdr + dCdr * _gamma_inv;
            rotGradient.col(e[0]) += (dPhi_nVi + dPhi_cVi) * _gamma_inv;
            rotGradient.col(e[1]) += (dPhi_nVj + dPhi_cVj) * _gamma_inv;
        }

        // Prepare area constraint variables
        Eigen::Matrix3Xd grad(3, _N);
        grad.setZero(3, _N);
        _value = 0.0;
        for (const auto &t : _triangles) {
            Vector3d p = positions.col(t[1]) - positions.col(t[0]);
            Vector3d q = positions.col(t[2]) - positions.col(t[0]);
            double_t S = p.cross(q).norm();
            _value += 0.5 * S;
            Vector3d dAdp = (q.dot(q) * p - p.dot(q) * q) / (2 * S);
            Vector3d dAdq = (p.dot(p) * q - p.dot(q) * p) / (2 * S);
            grad.col(t[0]) += -1.0 * (dAdp + dAdq);
            grad.col(t[1]) += dAdp;
            grad.col(t[2]) += dAdq;
        }
        double_t areaDiff = _value - _constrainedValue;
        _f += 0.5 * _K_i * areaDiff * areaDiff - _Lambda_i * areaDiff;
        posGradient += (_K_i * areaDiff - _Lambda_i) * grad;

        // Brownian body and viscosity body calculations
        MapVXd x(_x.data(), 3*_N);
        MapVXd g(_g.data(), 3*_N);
        _brownEn = -1.0 * _brownCoeff * (_xi.dot(x - _pX));
        _viscoEn = 0.5 * _viscosity * ((x - _pX).dot(x - _pX));
        _f += _brownEn + _viscoEn;
        g += -1.0 * _brownCoeff * _xi + _viscosity * (x - _pX);
    }

    //! Compute derivative of the normal wrt Rotation Vector
    void OPSModel::diffNormalRotVec() {
        MapM3Xd rotationVectors(&_x(3 * _N), 3, _N);
        for (auto i = 0; i < _N; ++i) {
            // Read the rotation vector
            Vector3d vi = rotationVectors.col(i);
            double_t s = sin(0.5 * vi.norm()), c = cos(0.5 * vi.norm());
            Quaterniond q(AngleAxisd(vi.norm(), vi.normalized()));
            double_t q0 = q.w(), q1 = q.x(), q2 = q.y(), q3 = q.z();
            Matrix4x3d dpdq;
            dpdq << q2, -q1, q0, q3, -q0, -q1, q0, q3, -q2, q1, q2, q3;
            dpdq = 2 * dpdq;
            Matrix3x4d dqdv;
            dqdv.leftCols(1) = -0.5 * s * vi.normalized();
            dqdv.rightCols(3) =
                (s * Eigen::Matrix3d::Identity() +
                 (0.5 * c - s / vi.norm()) * vi * vi.normalized().transpose()) /
                vi.norm();
            _diffNormalRV[i] = dqdv * dpdq;
        }
    }

    //! Compute derivative of the normal wrt Rotation Vector
    void OPSModel::diffNormalRotVec(const VectorXd &x) {
        for (auto i = 0; i < _N; ++i) {
            // Read the rotation vector
            Vector3d vi = x.segment<3>(3 * (_N + i));
            double_t s = sin(0.5 * vi.norm()), c = cos(0.5 * vi.norm());
            Quaterniond q(AngleAxisd(vi.norm(), vi.normalized()));
            double_t q0 = q.w(), q1 = q.x(), q2 = q.y(), q3 = q.z();
            Matrix4x3d dpdq;
            dpdq << q2, -q1, q0, q3, -q0, -q1, q0, q3, -q2, q1, q2, q3;
            dpdq = 2 * dpdq;
            Matrix3x4d dqdv;
            dqdv.leftCols(1) = -0.5 * s * vi.normalized();
            dqdv.rightCols(3) =
                (s * Eigen::Matrix3d::Identity() +
                 (0.5 * c - s / vi.norm()) * vi * vi.normalized().transpose()) /
                vi.norm();
            _diffNormalRV[i] = dqdv * dpdq;
        }
    }

    //! Calculate rotation vector with given point coordinates
    void OPSModel::initialRotationVector(RefM3Xd pos, RefM3Xd rotVec) {
        // Find unit normal along each point and calculate rotation vector
        // that would map the global z-axis to this unit normal
        for (auto i = 0; i < pos.cols(); ++i) {
            Vector3d x = pos.col(i);
            x.normalize();
            Vector3d axis, cross_prod;
            double_t angle, cross_prod_norm, p3;

            // Cross-product of x with z-axis.
            cross_prod << -x(1), x(0), 0.0;
            cross_prod_norm = cross_prod.norm();
            // Check if x is parallel or anti-parallel to z-axis
            p3 = x(2);
            if (cross_prod_norm < 1e-10) {
                axis << 1.0, 0.0,
                     0.0; // Arbitrarily choose the x-axis
                angle = (p3 > 0.0) ? 0.0 : M_PI;
            } else {
                angle = asin(cross_prod_norm);
                angle = (p3 < 0.0) ? (M_PI - angle) : angle;
                axis = cross_prod.normalized();
            }
            rotVec.col(i) = angle * axis;
        }
    }

    //! Minimize Rigid Body motions by applying Kabsch Algorithm
    void OPSModel::applyKabschAlgorithm() {
        MapM3Xd positions(_x.data(), 3, _N);
        computeNormals();
        Matrix3Xd pseudoNormal(3, _N);
        pseudoNormal = positions + _normals;
        Eigen::Affine3d A;
        MapM3Xd prevX(_pX.data(), 3, _N);
        A = find3DAffineTransform(positions, prevX);

        // Apply the transformation to each column in positions
        for (auto i = 0; i < _N; ++i) {
            Vector3d tempPos = A.linear() * positions.col(i) + A.translation();
            Vector3d tempN = A.linear() * pseudoNormal.col(i) + A.translation();
            positions.col(i) = tempPos;
            _normals.col(i) = (tempN - tempPos).normalized();
        }
        updateRotationVectors();
    }

    //! Update Rotation Vectors as per the current normals e.g. after Kabsch update
    void OPSModel::updateRotationVectors() {
        MapM3Xd rotationVectors(&_x(3 * _N), 3, _N);
        for (auto i = 0; i < _N; ++i) {
            Vector3d x = _normals.col(i);
            Vector3d axis, cross_prod;
            double_t angle, cross_prod_norm, p3;

            // Cross-product of x with z-axis.
            cross_prod << -x(1), x(0), 0.0;
            cross_prod_norm = cross_prod.norm();
            // Check if x is parallel or anti-parallel to z-axis
            p3 = x(2);
            if (cross_prod_norm < 1e-10) {
                axis << 1.0, 0.0,
                     0.0; // Arbitrarily choose the x-axis
                angle = (p3 > 0.0) ? 0.0 : M_PI;
            } else {
                angle = asin(cross_prod_norm);
                angle = (p3 < 0.0) ? (M_PI - angle) : angle;
                axis = cross_prod.normalized();
            }
            rotationVectors.col(i) = angle * axis;
        }
    }

    //! Calculate only the full mean-squared displacement
    double_t OPSModel::getMSD() {
        MapM3Xd positions(_x.data(), 3, _N);
        _msd = 0;
        vtkIdType nn;
        for (auto i = 0; i < _N; ++i) {
            Vector3d xi, xj, diff, xi_diff, xj_diff;
            Vector3d xi0, xj0, xi1, xj1, ni0, nj0;

            nn = _initialNearestNeighbor[i];

            xi0 = _initialPositions.col(i);
            xi1 = positions.col(i);

            xj0 = _initialPositions.col(nn);
            xj1 = positions.col(nn);

            xi_diff = (xi1 - xi0);
            xj_diff = (xj1 - xj0);

            diff = xi_diff - xj_diff;
            _msd += diff.dot(diff);
        }
        _msd /= _N;
        return _msd;
    }

    //! Calculate tangential and full mean-squared displacement
    std::array<double_t, 2> OPSModel::getMeanSquaredDisplacement() {
        MapM3Xd positions(_x.data(), 3, _N);
        _msd = 0;
        _msd_tgt = 0;
        std::array<double_t, 2> msdAll;
        vtkIdType nn;
        // We will subtract off the radial displacement.
        for (auto i = 0; i < _N; ++i) {
            Vector3d xi, xj, diff, xi_diff, xj_diff;
            Vector3d xi0, xj0, xi1, xj1, ni0, nj0;

            nn = _initialNearestNeighbor[i];

            xi0 = _initialPositions.col(i);
            xi1 = positions.col(i);
            ni0 = xi0.normalized();

            xj0 = _initialPositions.col(nn);
            xj1 = positions.col(nn);
            nj0 = xj0.normalized();

            xi_diff = (xi1 - xi0);
            xj_diff = (xj1 - xj0);

            xi = xi_diff - (ni0.dot(xi_diff) * ni0);
            xj = xj_diff - (nj0.dot(xj_diff) * nj0);

            diff = xi_diff - xj_diff;
            _msd += diff.dot(diff);

            diff = xi - xj;
            _msd_tgt += diff.dot(diff);
        }
        _msd /= _N;
        _msd_tgt /= _N;
        msdAll[0] = _msd;
        msdAll[1] = _msd_tgt;
        return msdAll;
    }

    //! Calculate asphericity
    double OPSModel::getAsphericity() {
        MapM3Xd positions(_x.data(), 3, _N);
        double_t asphericity = 0.0;
        double_t R0 = getAverageRadius();
        Eigen::RowVectorXd R(_N);
        R = positions.colwise().norm();
        asphericity += ((R.array() - R0).square()).sum();
        asphericity /= (_N * R0 * R0);
        return asphericity;
    }

    //! Calculate Volume
    double_t OPSModel::getVolume() {
        MapM3Xd positions(_x.data(), 3, _N);
        if (_updateVolume) {
            _volume = 0.0;
            for (const auto &f : _triangles)
                _volume +=
                    0.166666667 * (positions.col(f[0]).dot(
                                positions.col(f[2]).cross(positions.col(f[1]))));
            _updateVolume = false;
        }
        return _volume;
    }

    //! Calculate Area
    double_t OPSModel::getArea() {
        MapM3Xd positions(_x.data(), 3, _N);
        if (_updateArea) {
            _area = 0.0;
            for (const auto &f : _triangles)
                // Calculate area
                _area += 0.5 * (positions.col(f[1]) - positions.col(f[0]))
                    .cross(positions.col(f[2]) - positions.col(f[0]))
                    .norm();
            _updateArea = false;
        }
        return _area;
    }

    //! Get average radius
    double_t OPSModel::getAverageRadius() {
        MapM3Xd positions(_x.data(), 3, _N);
        if (_updateRadius) {
            _radius = 0.0;
            _radius = positions.colwise().norm().sum() / _N;
            _updateRadius = false;
        }
        return _radius;
    }

    //! Get Root Mean Squared Angle Deficit
    double_t OPSModel::getRMSAngleDeficit() {
        MapM3Xd positions(_x.data(), 3, _N);
        // Prepare a vector to store the angle deficit for each vertex
        Eigen::ArrayXd angleDeficit = VectorXd::Constant(_N, 2 * M_PI);
        for (const auto &f : _triangles) {
            // Get coordinates of all vertices
            Vector3d v0 = positions.col(f[0]);
            Vector3d v1 = positions.col(f[1]);
            Vector3d v2 = positions.col(f[2]);

            // Calculate unit vectors along the edges of the triangles
            Vector3d e0 = (v2 - v1).normalized();
            Vector3d e1 = (v0 - v2).normalized();
            Vector3d e2 = (v1 - v0).normalized();

            // Calculate the angles of the triangles
            double_t a0 = std::acos(e2.dot(-e1));
            double_t a1 = std::acos(e0.dot(-e2));
            double_t a2 = std::acos(e1.dot(-e0));

            // Subtract the angles from the corresponding vertices' curvature
            angleDeficit(f[0]) -= a0;
            angleDeficit(f[1]) -= a1;
            angleDeficit(f[2]) -= a2;
        }
        return std::sqrt(angleDeficit.square().mean());
    }

    void OPSModel::stereoDelaunay() {
        MapM3Xd positions(_x.data(), 3, _N);
        // Copy coordinates
        Matrix3Xd points(3, _N);
        points = positions;

        // Project points to unit sphere
        points.colwise().normalize();

        // Reset the center of the sphere to origin by translating
        Vector3d center = points.rowwise().mean();
        points = points.colwise() - center;

        // Rotate all points so that the point in 0th column is along z-axis
        Vector3d c = points.col(0);
        double_t cos_t = c(2);
        double_t sin_t = std::sqrt(1 - cos_t * cos_t);
        Vector3d axis;
        axis << c(1), -c(0), 0.;
        Matrix3d rotMat, axis_cross, outer;
        axis_cross << 0., -axis(2), axis(1), axis(2), 0., -axis(0), -axis(1), axis(0), 0.;

        outer.noalias() = axis * axis.transpose();

        rotMat = cos_t * Matrix3d::Identity() + sin_t * axis_cross + (1 - cos_t) * outer;
        Matrix3Xd rPts(3, _N);
        rPts = rotMat * points; // The points on a sphere rotated

        // Calculate the stereographic projections
        Vector3d p0;
        MapM3Xd l0(&(rPts(0, 1)), 3, _N - 1);
        Matrix3Xd l(3, _N - 1), proj(3, _N - 1);
        p0 << 0, 0, -1;
        c = rPts.col(0);
        l = (l0.colwise() - c).colwise().normalized();
        for (auto j = 0; j < _N - 1; ++j) {
            proj.col(j) = ((p0(2) - l0(2, j)) / l(2, j)) * l.col(j) + l0.col(j);
        }

        // Insert the projected points in a CGAL vertex_with_info vector
        std::vector<std::pair<Point, unsigned>> verts;
        for (auto j = 0; j < _N - 1; ++j) {
            verts.push_back(std::make_pair(Point(proj(0, j), proj(1, j)), j + 1));
        }

        _dt.clear();
        _dt.insert(verts.begin(), verts.end());
    }

    void OPSModel::updateTriangles() {
        stereoDelaunay();
        _edges.clear();
        _triangles.clear();
        // ***************** Our data structure ********************//
        std::set<std::set<unsigned>> tri, edges;

        // Iterate over all vertices and collect first ring neighbors
        for (auto fvi = _dt.all_vertices_begin(); fvi != _dt.all_vertices_end();
                ++fvi) {

            auto vid = _dt.is_infinite(fvi) ? 0 : fvi->info();
            Delaunay::Edge_circulator ec = _dt.incident_edges(fvi), done(ec);

            // Lambda function to get the vertex id for the edge
            auto getVertexId = [](int a, int b) {
                std::set<int> index{0, 1, 2};
                index.erase(a);
                index.erase(b);
                return *index.begin();
            };

            if (ec != 0) {
                do {
                    auto fh = ec->first;
                    auto edgeIndex = getVertexId(fh->index(fvi), ec->second);
                    auto verH = fh->vertex(edgeIndex);
                    auto edgeId = _dt.is_infinite(verH) ? 0 : verH->info();
                    std::set<unsigned> edge{vid, edgeId};
                    auto tryInsertEdge = edges.insert(edge);
                    if (tryInsertEdge.second) {
                        std::array<unsigned, 2> ed{vid, edgeId};
                        _edges.push_back(ed);
                    }
                } while (++ec != done);
            }

            // Iterate over triangles associated with a vertex
            Delaunay::Face_circulator fc = _dt.incident_faces(fvi), done2(fc);
            if (fc != 0) {
                do {
                    std::set<unsigned> t;
                    std::array<unsigned, 3> face;
                    for (auto j = 0; j < 3; ++j) {
                        auto verId =
                            _dt.is_infinite(fc->vertex(j)) ? 0 : fc->vertex(j)->info();
                        t.insert(verId);
                        face[j] = verId;
                    }
                    auto tryInsertFace = tri.insert(t);
                    if (tryInsertFace.second)
                        _triangles.push_back(face);
                } while (++fc != done2);
            }
        }
        _updateRadius = true;
        _updateVolume = true;
        _updateArea = true;
    }

    // Kabsch algorithm. The input 3D points are stored as columns.
    Eigen::Affine3d OPSModel::find3DAffineTransform(Eigen::Ref<Eigen::Matrix3Xd> in,
            Eigen::Ref<Eigen::Matrix3Xd> out) {
        // Default output
        Eigen::Affine3d A;
        A.linear() = Eigen::Matrix3d::Identity(3, 3);
        A.translation() = Eigen::Vector3d::Zero();

        if (in.cols() != out.cols())
            throw "Find3DAffineTransform(): input data mis-match";

        // Find the centroids then shift to the origin
        Eigen::Vector3d in_ctr = Eigen::Vector3d::Zero();
        Eigen::Vector3d out_ctr = Eigen::Vector3d::Zero();
        for (int col = 0; col < in.cols(); col++) {
            in_ctr += in.col(col);
            out_ctr += out.col(col);
        }
        in_ctr /= in.cols();
        out_ctr /= out.cols();
        for (int col = 0; col < in.cols(); col++) {
            in.col(col) -= in_ctr;
            out.col(col) -= out_ctr;
        }
        Eigen::MatrixXd Cov = in * out.transpose();
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov, Eigen::ComputeThinU |
                Eigen::ComputeThinV);

        // Find the rotation
        double_t d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
        if (d > 0)
            d = 1.0;
        else
            d = -1.0;
        Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
        I(2, 2) = d;
        Eigen::Matrix3d R = svd.matrixV() * I * svd.matrixU().transpose();

        // The final transform
        A.linear() = R;
        A.translation() = out_ctr - R * in_ctr;

        return A;
    }

    //! Find average edge lengthf of a point cloud read from VTK file.
    double_t OPSModel::getPointCloudAvgEdgeLen(std::string vtkfile){
        vtkNew<vtkPolyDataReader> rd;
        vtkNew<vtkPoints> newpts;
        vtkNew<vtkDelaunay3D> d3D;
        vtkNew<vtkDataSetSurfaceFilter> dssf;
        vtkNew<vtkExtractEdges> ee;
        vtkNew<vtkIdList> verts;

        rd->SetFileName(vtkfile.c_str());
        d3D->SetInputConnection(rd->GetOutputPort());
        dssf->SetInputConnection(d3D->GetOutputPort());
        ee->SetInputConnection(dssf->GetOutputPort());
        ee->Update();

        auto edges = ee->GetOutput()->GetLines();
        edges->InitTraversal();
        auto avgEdgeLen = 0.0;
        while( edges->GetNextCell(verts) ){
            Vector3d xi, xj;
            ee->GetOutput()->GetPoint( verts->GetId(0), xi.data() );
            ee->GetOutput()->GetPoint( verts->GetId(1), xj.data() );
            avgEdgeLen += (xj-xi).norm();
        }
        avgEdgeLen /= edges->GetNumberOfCells();
        return avgEdgeLen;
    }

    void OPSModel::writeSimulationState(std::string file){
        // Overwrite any existing file
        ofstream f(file.c_str());
        double_t gamma = 1.0 / _gamma_inv;
        double_t beta = 2.0 * _viscosity/(_brownCoeff*_brownCoeff);
        size_t numsPerLine = 9;

        f << "# !!! Auto-generated file DO NOT EDIT !!! #" << std::endl
            << "# Number of particles" << std::endl
            << _N << std::endl
            << "# Latest output VTK file suffix" << std::endl
            << _vtkfilesuffix << std::endl
            << "# Number of time steps completed" << std::endl
            << _timestep << std::endl
            << "# Gamma or the FvK number" << std::endl
            << gamma << std::endl
            << "# Beta or 1/Temperature" << std::endl
            << beta << std::endl
            << "# State of the random number engine" << std::endl
            << _e2 << std::endl
            << "# State of the random number generator" << std::endl
            << _rng << std::endl
            << "# Initial nearest neighbors" << std::endl;

        size_t numLines = _N / numsPerLine;
        size_t startIdx = 0;

        // Write nearest neighbors
        for (auto i = 0; i < numLines; ++i) {
            startIdx = i * numsPerLine;
            for (auto j = 0; j < numsPerLine - 1; ++j) {
                f << _initialNearestNeighbor[startIdx++] << " ";
            }
            f << _initialNearestNeighbor[startIdx++] << std::endl;
        }
        if (_N % numsPerLine > 0) {
            for (; startIdx < _N - 1;) {
                f << _initialNearestNeighbor[startIdx++] << " ";
            }
            f << _initialNearestNeighbor[startIdx] << std::endl;
        }

        // Write initial position vector
        f << "# Particle positions at step 0" << std::endl;
        numLines = (3 * _N) / numsPerLine;
        for (auto i = 0; i < numLines; ++i) {
            startIdx = i * numsPerLine;
            for (auto j = 0; j < numsPerLine - 1; ++j) {
                f << _initialPositions(startIdx++) << " ";
            }
            f << _initialPositions(startIdx++) << std::endl;
        }
        if (3 * _N % numsPerLine > 0) {
            for (; startIdx < 3 * _N - 1;) {
                f << _initialPositions(startIdx++) << " ";
            }
            f << _initialPositions(startIdx) << std::endl;
        }

        // Write prevX vector
        f << "# Particle positions of previous step" << std::endl;
        numLines = (3 * _N) / numsPerLine;
        for (auto i = 0; i < numLines; ++i) {
            startIdx = i * numsPerLine;
            for (auto j = 0; j < numsPerLine - 1; ++j) {
                f << _pX(startIdx++) << " ";
            }
            f << _pX(startIdx++) << std::endl;
        }
        if (3 * _N % numsPerLine > 0) {
            for (; startIdx < 3 * _N - 1;) {
                f << _pX(startIdx++) << " ";
            }
            f << _pX(startIdx) << std::endl;
        }

        // Write x vector
        f << "# Particle positions and rotation vectors of latest step" << std::endl;
        numLines = (6 * _N) / numsPerLine;
        for (auto i = 0; i < numLines; ++i) {
            startIdx = i * numsPerLine;
            for (auto j = 0; j < numsPerLine - 1; ++j) {
                f << _x(startIdx++) << " ";
            }
            f << _x(startIdx++) << std::endl;
        }
        if (6 * _N % numsPerLine > 0) {
            for (; startIdx < 6 * _N - 1;) {
                f << _x(startIdx++) << " ";
            }
            f << _x(startIdx) << std::endl;
        }
        f.close();
    }

} // namespace OPS
