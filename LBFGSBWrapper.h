// Reference:
//  1. Nonlinear quasi-Newton solver using L-BFGS-B code of Jorge Nocedal.
//  http://www.ece.northwestern.edu/~nocedal/lbfgsb.html
//  2. https://github.com/wsklug/voom -> src/Solvers/Lbfgsb.h
//
/////////////////////////////////////////////////////////////////////////

#if !defined(__LBFGSBWRAPPER_H__)
#define __LBFGSBWRAPPER_H__

#include <Eigen/Dense>
#include <iomanip>
#include <ios>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "OPSModel.h"

#ifdef _NO_PRINTING_
#define PRINT(X)                                                               \
  do {                                                                         \
  } while (0)
#define PRINT_CHAR_ARR(X, N)                                                   \
  do {                                                                         \
  } while (0)
#define UNSET(X)                                                               \
  do {                                                                         \
  } while (0)
#else
#define PRINT(X)                                                               \
  do {                                                                         \
    std::cout << X;                                                            \
  } while (0)
#define PRINT_CHAR_ARR(X, N)                                                   \
  do {                                                                         \
    std::cout.write(X, N) << std::endl;                                        \
  } while (0)
#define UNSET(X)                                                               \
  do {                                                                         \
    std::cout.unsetf(X);                                                       \
  } while (0)
#endif

namespace OPS {

//! The l-BFGS-b Wrapper class
class LBFGSBWrapper {
public:
  typedef Eigen::VectorXd Vector_t;
  typedef Eigen::Ref<Vector_t> RefV;
  typedef Eigen::Ref<const Vector_t> RefCV;
  typedef Eigen::Map<Vector_t> MapV;
  typedef Eigen::VectorXi IntVector_t;
  typedef Eigen::Ref<IntVector_t> RefI;
  typedef Eigen::Ref<const IntVector_t> RefCI;

  LBFGSBWrapper(OPSModel &s);
  void setBounds(const RefCI nbd, const RefCV l, const RefCV u);
  void setNumHessianCorrections(size_t m) { _m = m; }
  void setPrintCode(size_t p) { _iprint = p; }
  void setMaxIterations(size_t i) { _maxIterations = i; }
  void setMachineEPSFactor(double_t f) { _factr = f; }
  void setProjectedGradientTolerance(double_t p) { _pgtol = p; }
  void solve();

private:
  void resize(size_t n);
  OPSModel &_model;
  Vector_t _l;
  Vector_t _u;
  Vector_t _wa;
  IntVector_t _nbd;
  IntVector_t _iwa;

  int _n, _m = 5, _iprint = 1000, _maxIterations = 100000, _iterNo;
  double _factr = 10.0, _pgtol = 1e-8;
  double _projg;
};
} // namespace OPS
#endif // __LBFGSBWRAPPER_H__
