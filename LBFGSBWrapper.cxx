#include "LBFGSBWrapper.h"

#if defined(_WIN32)
#if defined(_MSC_VER) || defined(__INTEL_COMPILER)
#define setulb_ SETULB
#endif
#endif

extern "C" void setulb_(int *n, int *m, double_t *x, double_t *l, double_t *u,
	int *nbd, double_t *f, double_t *g, double_t *factr,
	double_t *pgtol, double_t *wa, int *iwa, char *task,
	int *iprint, char *csave, int *lsave, int *isave,
	double_t *dsave);

using std::endl;
using std::left;
using std::right;
using std::scientific;
using std::setw;

namespace OPS {
    //! Constructor
    LBFGSBWrapper(OPSModel &s) : _model(m){
	_n = 6*_model.getNumberOfPoints();
	resize(_n);
	_projg = 0.0;
    }

    void LBFGSBWrapper::resize(size_t n) {
	_n = n;
	_iwa = IntVector_t::Zero(3 * _n);
	_wa = Vector_t::Zero(2 * _m * _n + 5 * _n + 11 * _m * _m + 8 * _m);
	_nbd.resize(_n);
	_l.resize(_n);
	_u.resize(_n);
	_nbd = IntVector_t::Zero(_n);
	_l = Vector_t::Zero(_n);
	_u = Vector_t::Zero(_n);
    }

    void LBFGSBWrapper::setBounds(const RefCI nbd, const RefCV l, const RefCV u) {
	if (nbd.size() != _n || l.size() != _n || u.size() != _n) {
	    std::cout
		<< "LBFGSBWrapper::setBounds(): input arrays are incorrectly sized."
		<< std::endl;
	    return;
	}
	_nbd = nbd;
	_l = l;
	_u = u;
	return;
    }

    void LBFGSBWrapper::solve() {
	// set up misc. arrays and data
	char task[60], csave[60];
	for (int i = 0; i < 60; i++)
	    task[i] = csave[i] = '\0';

	int lsave[4];
	for (int i = 0; i < 4; i++)
	    lsave[i] = 0;

	double_t dsave[29];
	for (int i = 0; i < 29; i++)
	    dsave[i] = 0.0;

	int isave[44];
	for (int i = 0; i < 44; i++)
	    isave[i] = 0;

	MapV g(_model.getG(), _n);

	PRINT("======================================================================"
		"=========="
		<< std::endl
		<< std::endl
		<< "Starting BFGS iterations." << std::endl
		<< std::endl);

	PRINT(setw(14) << scientific << right << "|proj g|" << setw(14) << scientific
		<< right << "|g|" << setw(14) << scientific << right << "f"
		<< setw(14) << right << "iterations" << setw(14) << right
		<< "evaluations" << std::endl
		<< "----------------------------------------------------------"
		"----------------------"
		<< std::endl);
	// We start the iteration by initializing task.

	sprintf(task, "START");

	int iprint = -1;

	_model.compute();

	// ------- the beginning of the loop ----------
	while (true) {

	    // This is the call to the L-BFGS-B code.
	    setulb_(&_n, &_m, _model.getX(), _l.data(), _u.data(), _nbd.data(), _model.getF(),
		    _model.getG(), &_factr, &_pgtol, _wa.data(), _iwa.data(), &(task[0]),
		    &iprint, &(csave[0]), &(lsave[0]), &(isave[0]), &(dsave[0]));
	    if (strncmp(task, "FG", 2) == 0) {
		// The minimization routine has returned to request the
		// function f and gradient g values at the current x.

		_model.compute();

		// Go back to the minimization routine.
		continue;
	    }

	    else if (strncmp(task, "NEW_X", 5) == 0) {
		// stop if maximum number of iterations has been reached
		if (_maxIterations > 0 && isave[29] > _maxIterations) {
		    break;
		}

		// The minimization routine has returned with a new iterate,
		// and we have opted to continue the iteration.
		if (_iprint > 0 && isave[29] % _iprint == 0) {
		    PRINT(setw(14) << scientific << right << dsave[12] << setw(14)
			    << scientific << right << (g.cwiseAbs()).maxCoeff()
			    << setw(14) << scientific << right << *_model.getF() << setw(14)
			    << right << isave[29] << setw(14) << right << isave[33]
			    << std::endl);
		}
		continue;
	    }

	    else if (strncmp(task, "CONV", 4) == 0) {
		_model.compute();
		PRINT_CHAR_ARR(task, 49);
		break;
	    }

	    else if (strncmp(task, "ABNORM", 6) == 0) {
		_model.compute();
		PRINT_CHAR_ARR(task, 49);
		break;
	    }

	    // If task is neither FG nor NEW_X we terminate execution.
	    else {
		PRINT_CHAR_ARR(task, 49);
		break;
	    }

	    // ---------- the end of the loop -------------
	}
	PRINT(setw(14) << scientific << right << dsave[12] << setw(14) << scientific
		<< right << (g.cwiseAbs()).maxCoeff() << setw(14)
		<< scientific << right << *_model.getF() << setw(14) << right << isave[29]
		<< setw(14) << right << isave[33] << std::endl
		<< std::endl
		<< "=========================================================="
		"======================"
		<< std::endl
		<< std::endl);

	UNSET(std::ios_base::scientific);

	// ---------- the end of solve() -------------
	_projg = dsave[12];
	_iterNo = isave[29];

	_model.compute();
    }
} // namespace OPS
