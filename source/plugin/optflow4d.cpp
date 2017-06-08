/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2016 Nils Thuerey
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL)
 * http://www.gnu.org/licenses
 *
 * FlOF - Optical flow plugins and helpers
 *
 ******************************************************************************/

#include <thread>
#include <cmath>
#include <algorithm>
#include <cstring>

#include "levelset.h"
#include "commonkernels.h"
#include "particle.h"
#include "timing.h"
#include "vector4d.h"
#include "grid4d.h"

// for loading slices
#include "fileio.h"
#include "interpol.h"

using namespace std;

namespace Manta
{

// grid helpers

void gridFactor4d(Vec4 s1, Vec4 s2, Vec4 optSize, Vec4 scale, Vec4 &srcFac,
                  Vec4 &retOff);

template <>
void interpolGridTempl<Grid4d<Real> >(Grid4d<Real> &target,
                                      Grid4d<Real> &source)
{
	Vec4 sourceFactor =
	    calcGridSizeFactor4d(source.getSize(), target.getSize());
	KnInterpolateGrid4dTempl<Real>(target, source, sourceFactor,
	                               sourceFactor * 0.5);
}
template <>
void interpolGridTempl<Grid4d<Vec4> >(Grid4d<Vec4> &target,
                                      Grid4d<Vec4> &source)
{
	Vec4 sourceFactor =
	    calcGridSizeFactor4d(source.getSize(), target.getSize());
	KnInterpolateGrid4dTempl<Vec4>(target, source, sourceFactor,
	                               sourceFactor * 0.5);
}

// get index for serialized vec grid
template <class GRID>
IndexInt getIdx(GRID &vel, int DIM, int ti, int tj, int tk, int tt, int d)
{
	return (vel.index(ti, tj, tk, tt) * DIM + d);
}

// declarations

template <class GRID, class GRIDVEC>
void corrVelsOfTempl(GRIDVEC &dst, GRIDVEC &vel, GRID &phiOrg, GRID &phiCurr,
                     GRID &phiTarget, Real threshPhi = 1e10,
                     Real threshNorm = 1e10, Real postVelBlur = 0.,
                     Real resetBndWidth = -1., int maxIter = 100);

// ================================================================================
// misc helpers

template <class GRIDVEC> void ofSrchMax(GRIDVEC &vel)
{
	typedef typename GRIDVEC::BASETYPE VEC;
	Vec4i p(-1);
	VEC v(0.);
	Real max = -1.;
	FOR_IJKT_BND(vel, 0)
	{
		Real l = norm(vel(i, j, k, t));
		if (l > max) {
			max = l;
			p = Vec4i(i, j, k, t);
			v = vel(i, j, k, t);
		}
	}
	debMsg("MAX found at " << p << " = " << v, 1);
}

//! simple blur
template <class S>
void simpleBlurFunc(Grid<S> &a, int iter = 1, Real strength = 1.)
{
	Grid<S> tmp(a.getParent());
	for (int numIt = 0; numIt < iter; numIt++) {
		FOR_IJK_BND(a, 1)
		{
			tmp(i, j, k) = a(i + 1, j, k) + a(i - 1, j, k) + a(i, j + 1, k) +
			               a(i, j - 1, k); // avg
			if (a.is3D()) {
				tmp(i, j, k) += a(i, j, k + 1) + a(i, j, k - 1);
				tmp(i, j, k) *= 1. / 6.;
			} else {
				tmp(i, j, k) *= 1. / 4.;
			}
			tmp(i, j, k) =
			    tmp(i, j, k) * strength +
			    a(i, j, k) * (1. - strength); // optionally reduce strength
		}
		a.swap(tmp);
	}
}
PYTHON() void simpleBlur(Grid<Real> &a, int iter = 1)
{
	simpleBlurFunc<Real>(a, iter);
}
PYTHON() void simpleBlurVec(Grid<Vec3> &a, int iter = 1)
{
	simpleBlurFunc<Vec3>(a, iter);
}

//! gaussian blur func
static inline Real gaussianWeight(Real sigma, Real dSqr)
{
	return exp(-dSqr / (2. * sigma * sigma));
}

KERNEL(fourd, bnd = 1) template <class GridT, class S>
void knGaussianBlur(GridT &a, GridT &tmp, Real sigma, Vec4i &size)
{
	S val(0.);
	Real weight = 0.;
	for (int vt = t - size[3]; vt <= t + size[3]; vt++)
		for (int zk = k - size[2]; zk <= k + size[2]; zk++)
			for (int yj = j - size[1]; yj <= j + size[1]; yj++)
				for (int xi = i - size[0]; xi <= i + size[0]; xi++) {
					if (!a.isInBounds(xi, yj, zk, vt, 0)) continue;

					Real d = (xi - i) * (xi - i) + (yj - j) * (yj - j);
					if (a.is3D()) d += (zk - k) * (zk - k);
					if (a.is4D()) d += (vt - t) * (vt - t);

					Real wcurr = gaussianWeight(sigma, d);
					weight += wcurr;
					val += wcurr * a(xi, yj, zk, vt);
				}

	if (weight > VECTOR_EPSILON) {
		tmp(i, j, k, t) = val / weight;
	} else {
		tmp(i, j, k, t) = a(i, j, k, t);
	}
}

template <class GRID>
void gaussianBlurGeneric(GRID &a, Real sigma = 1., int iter = 1)
{
	typedef typename GRID::BASETYPE S;
	GRID tmp(a.getParent());
	int s = int(1. * sigma + 0.5); // warning - optimization, kernel cut off
	if (s == 0) s = 1;
	Vec4i size(s, s, s, s);
	if (!a.is4D()) size[3] = 0;
	if (!a.is3D()) size[2] = 0;
	for (int numIt = 0; numIt < 2 * iter; numIt++) {
		knGaussianBlur<GRID, S>(a, tmp, sigma, size);
		a.swap(tmp);
	}
}

// minimal CG solver

//! fixed matrix for optical flow problems
template <class T, int DIM> class FixedMatrixOF
{
  protected:
	IndexInt n;           // dimension
	IndexInt indices[32]; // max 32 per row for now...

	std::vector<T> offd;   // nonzero values , laplacian
	std::vector<T> blockd; // nonzero values , block diags

  public:
	int size() const { return n; }

	explicit FixedMatrixOF(int n_ = 0, int dim_ = 3) : n(n_), offd(0), blockd(0)
	{
		offd.resize(n);
		blockd.resize(DIM * n);
		for (IndexInt i = 0; i < n; ++i)
			offd[i] = 0.;
		for (IndexInt i = 0; i < DIM * n; ++i)
			blockd[i] = 0.;
	}

	void setOffd(IndexInt idx, T val) { offd[idx] = val; }
	void setBlockd(IndexInt idx, int d, T val) { blockd[idx * DIM + d] = val; }
	void setOffdIndex(IndexInt l, IndexInt idx) { indices[l] = idx; }

	T blockdVal(IndexInt idx, int l) const { return blockd[idx * DIM + l]; }
	T offdVal(IndexInt idx) const { return offd[idx]; }
	IndexInt offdIndex(int l) const { return indices[l]; }
};

// multiply x with OF matrix
template <class T, int DIM>
void applyMat(const FixedMatrixOF<T, DIM> &mat, const std::vector<T> &x,
              std::vector<T> &result)
{
	assertMsg(mat.size() == x.size(),
	          "Matrix and vector dimensions dont match!");

	for (IndexInt i = 0; i < mat.size(); ++i) {
		T v = 0;
		const T offd = mat.offdVal(i);
		if (offd != 0.) {
			for (int m = (4 - DIM); m < (4 - DIM) + 2 * DIM; ++m) {
				v += offd * x[i + mat.offdIndex(m)];
			}
		}
		const IndexInt blockIdx = (i / DIM) * DIM;
		for (int m = 0; m < DIM; ++m) {
			v += mat.blockdVal(i, m) * x[blockIdx + m];
		}
		result[i] = v;
	}
}

template <class T>
double dotProd(const std::vector<T> &a, const std::vector<T> &b)
{
	double d = 0.;
	for (IndexInt i = 0; i < a.size(); ++i)
		d += a[i] * b[i];
	return d;
}

template <class T> T getMaxNorm(const std::vector<T> &a)
{
	T maxv = a[0];
	for (IndexInt i = 1; i < a.size(); ++i) {
		if (a[i] > maxv) maxv = a[i];
	}
	return maxv;
}

//! compute a = a + b * c
template <class T>
void addScaled(std::vector<T> &a, double b, const std::vector<T> &c)
{
	for (IndexInt i = 0; i < a.size(); ++i)
		a[i] += b * c[i];
}

template <class T> void setConst(std::vector<T> &v, T value)
{
	for (IndexInt i = 0; i < v.size(); ++i)
		v[i] = value;
}

//! cg solver with fixed (hard coded) OF matrix, i.e. without explicit matrix
// construction
template <class T, int DIM> class GridCGOptflow4d
{
  public:
	GridCGOptflow4d(void) {}

	//! run cg solver
	bool solve(const FixedMatrixOF<T, DIM> &matrix, const std::vector<T> &rhs,
	           std::vector<T> &result, T &ret_residual, int &ret_iterations,
	           Real accuracy, int cgMaxIter)
	{
		const int n = matrix.size();
		std::vector<T> srch, res, precond, tmp;
		srch.resize(n);
		tmp.resize(n);
		res.resize(n);

		setConst<T>(result, 0.);
		res = rhs;
		double residual = getMaxNorm(res);
		if (residual < VECTOR_EPSILON) {
			ret_residual = 0;
			ret_iterations = 0;
			return true;
		}
		const double acc = accuracy * residual;
		const double resIni = residual;

		precondInit(matrix, precond);
		precondApply(res, tmp, precond);
		double sigma = dotProd(tmp, res);
		if (sigma == 0 || sigma != sigma) {
			ret_iterations = 0;
			return false;
		}

		srch = tmp;
		for (int iter = 0; iter < cgMaxIter; ++iter) {
			applyMat(matrix, srch, tmp);

			double alpha1 = dotProd(srch, tmp);
			double alpha = sigma / alpha1;
			addScaled(result, alpha, srch);
			addScaled(res, -alpha, tmp);
			residual = getMaxNorm(res);
			ret_residual = residual / resIni;

			if (residual <= acc) {
				ret_iterations = iter + 1;
				return true;
			}
			precondApply(res, tmp, precond);
			double sigmaNew = dotProd(tmp, res);
			double beta = sigmaNew / sigma;
			addScaled(tmp, beta, srch);
			// srch.swap(tmp);
			srch = tmp;
			sigma = sigmaNew;
			ret_iterations = iter + 1;
		}
		ret_residual = residual / resIni;
		return false;
	}

	void precondInit(const FixedMatrixOF<T, DIM> &matrix,
	                 std::vector<T> &precond)
	{
		// diagonal
		precond.resize(matrix.size());
		for (int i = 0; i < matrix.size(); ++i) {
			T diag = matrix.blockdVal(i, i - (i / DIM) * DIM);
			if (diag != 0.)
				precond[i] = 1. / diag;
			else
				precond[i] = 1. / diag;
		}
	}

	void precondApply(const std::vector<T> &x, std::vector<T> &result,
	                  std::vector<T> &precond)
	{
		// diagonal
		for (IndexInt i = 0; i < (IndexInt)result.size(); ++i) {
			result[i] = x[i] * precond[i];
		}
		// turn off:
		// result = x;
	}
};

// ================================================================================
// base OF solve

//! templated, with fixed dimension
template <class GRID, class GRIDVEC, int DIM>
void opticalFlowDim(GRIDVEC &vel, GRID &i0, GRID &i1, GRID *rhsT = NULL,
                    Real wSmooth = 0., Real wEnergy = 0., Real postVelBlur = 0.,
                    Real cgAccuracy = 1e-04, int blurType = 1,
                    Real resetBndWidth = -1.)
{
	typedef typename GRIDVEC::BASETYPE VEC;
	const Real mDt = vel.getParent()->getDt();
	const Real mDx = 1. / vel.getSizeX();

	// const Real mDx2Inv = 1./(mDx*mDx);
	const Real mDx2Inv = 1.; // discretize in local space

	int dim[4] = { vel.getSizeX(), vel.getSizeY(), vel.getSizeZ(),
		           vel.getSizeT() };
	int N = dim[0] * dim[1] * dim[2] * dim[3] * DIM;
	FixedMatrixOF<Real, DIM> Af(N);

	std::vector<Real> rhs(N);
	for (int j = 0; j < N; ++j) {
		rhs[j] = 0.;
	}

	std::vector<Real> result(N);
	for (int j = 0; j < N; ++j) {
		result[j] = 0.;
	}

	static const int nbx[] = { 0, 0, 0, -1, +1, 0, 0, 0 };
	static const int nby[] = { 0, 0, -1, 0, 0, +1, 0, 0 };
	static const int nbz[] = { 0, -1, 0, 0, 0, 0, +1, 0 };
	static const int nbt[] = { -1, 0, 0, 0, 0, 0, 0, +1 };

	debMsg("ofSolveT setting up matrix", 2);
	int cntSkip = 0;

	// copy into structs
	FOR_IJKT_BND(vel, 0)
	{
		for (int d = 0; d < DIM; ++d) {
			IndexInt vIdx = getIdx(vel, DIM, i, j, k, t, d);
			result[vIdx] = vel(i, j, k, t)[d] * mDx;
		}
	}

	Vec4i p(1, 1, 1, 1);
	if (!vel.is3D()) p[2] = 0;
	if (!vel.is4D()) p[3] = 0;
	IndexInt idx = getIdx(vel, DIM, p[0], p[1], p[2], p[3], 0);
	for (int m = 0 + (4 - DIM); m < (4 - DIM) + 2 * DIM; m++) {
		int ti = p[0] + nbx[m];
		int tj = p[1] + nby[m];
		int tk = p[2] + nbz[m];
		int tt = p[3] + nbt[m];
		if (vel.isInBounds(ti, tj, tk, tt, 0)) {
			IndexInt nbidx = getIdx(vel, DIM, ti, tj, tk, tt, 0);
			Af.setOffdIndex(m, nbidx - idx);
		}
	}

	// setup matrix
	FOR_IJKT_BND(vel, 0)
	{
		if (!vel.isInBounds(i, j, k, t, 1)) {
			// add dummy entries, zero rhs, and 1 on diagonal
			for (int d = 0; d < DIM; ++d) {
				IndexInt vIdx = getIdx(vel, DIM, i, j, k, t, d);
				result[vIdx] = 0.;
				Af.setBlockd(vIdx, d, 1.);
			}
			cntSkip += DIM;
			continue;
		}

		// regular entries

		// calculate temporal & spatial gradients
		Real tderiv = (i1(i, j, k, t) - i0(i, j, k, t)) / mDt;

		// use central difference
		Real dxf = 1. / (2. * mDx);
		VEC grad(0.);
		if (true) {
			grad[0] = (i1(i + 1, j, k, t) - i1(i - 1, j, k, t)) * dxf;
			grad[1] = (i1(i, j + 1, k, t) - i1(i, j - 1, k, t)) * dxf;
		}
		if (DIM > 2) {
			grad[2] = (i1(i, j, k + 1, t) - i1(i, j, k - 1, t)) * dxf;
		}
		if (DIM > 3) {
			grad[3] = (i1(i, j, k, t + 1) - i1(i, j, k, t - 1)) * dxf;
		}

		IndexInt vIdx[4] = { 0, 0, 0, 0 };
		for (int d = 0; d < DIM; ++d)
			vIdx[d] = getIdx(vel, DIM, i, j, k, t, d);

		for (int d = 0; d < DIM; ++d) {

			// rhs is gradient dot time derivative
			rhs[vIdx[d]] = -grad[d] * tderiv;

			for (int m = 0 + (4 - DIM); m < (4 - DIM) + 2 * DIM; m++) {
				int ti = std::max(0, std::min(dim[0] - 1, i + nbx[m]));
				int tj = std::max(0, std::min(dim[1] - 1, j + nby[m]));
				int tk = std::max(0, std::min(dim[2] - 1, k + nbz[m]));
				int tt = std::max(0, std::min(dim[3] - 1, t + nbt[m]));

				rhs[vIdx[d]] -= wSmooth *
				                (vel(i, j, k, t)[d] - vel(ti, tj, tk, tt)[d]) *
				                mDx * mDx2Inv;
			}

			// fixed matrix

			Af.setOffd(vIdx[d], -wSmooth * mDx2Inv);

			Real diag = 0.;
			diag += (Real)(2 * DIM) * wSmooth * mDx2Inv;
			diag += wEnergy;

			// add image gradient squared
			for (int m = 0; m < DIM; ++m) {
				if (d == m)
					Af.setBlockd(vIdx[d], m, grad[d] * grad[m] + diag);
				else
					Af.setBlockd(vIdx[d], m, grad[d] * grad[m]);
			}

			// tikhonov prior
			rhs[vIdx[d]] -= wEnergy * vel(i, j, k, t)[d] * mDx;
		}
	}

	debMsg("ofSolveT solving... ", 2);
	MuTime cgTimer;

	// cg solver
	Real cgRes = 1e10;
	int cgIter = -1;
	GridCGOptflow4d<Real, DIM> gridcg;
	gridcg.solve(Af, rhs, result, cgRes, cgIter, cgAccuracy, 1000);
	debMsg("ofSolve fix iterations:" << cgIter << " ("
	                                 << float(cgTimer.update().time / 1000.)
	                                 << "s) ",
	       1);

	// check for diverging solve / nans
	if (cgRes != cgRes) {
		debMsg("Bad news: CG did not converge.", 1);
		// reset solution to zero...
		for (int i = 0; i < N; ++i)
			result[i] = 0.;
	}
	if (cgAccuracy < cgRes) {
		debMsg("Warning: CG did not fully converge. Using result anyway...", 1);
	}

	// copy back
	FOR_IJKT_BND(vel, 0)
	{
		for (int d = 0; d < DIM; ++d) {
			IndexInt vIdx = getIdx(vel, DIM, i, j, k, t, d);
			vel(i, j, k, t)[d] = result[vIdx] / mDx;

			// for debugging, return part of rhs vec if desired
			if ((d == 0) && rhsT) (*rhsT)(i, j, k, t) = rhs[vIdx];
		}
	}

	// optionally blur result
	if (postVelBlur > 0.) {
		switch (blurType) {
		case 1:
			gaussianBlurGeneric<GRIDVEC>(vel,
			                             0.5 * postVelBlur); // gaussian blur
			break;
		default:
			errMsg("NYI");
		}
	}

	// reset outer border (after smoothing!)
	if (resetBndWidth > 0.) {
		int resetBnd = int(resetBndWidth * vel.getSizeX()) + 1;
		FOR_IJKT_BND(vel, 0)
		{
			if (vel.isInBounds(i, j, k, t, resetBnd)) continue;
			vel(i, j, k, t) = VEC(0.);
		}
	}

} // OF

template <class GRID, class GRIDVEC>
void opticalFlowTemplate(GRIDVEC &vel, GRID &i0, GRID &i1, GRID *rhsT = NULL,
                         Real wSmooth = 0., Real wEnergy = 0.,
                         Real postVelBlur = 0., Real cgAccuracy = 1e-04,
                         int blurType = 1, Real resetBndWidth = -1.)
{
	// explicit dimension check
	if (!vel.is3D())
		opticalFlowDim<GRID, GRIDVEC, 2>(vel, i0, i1, rhsT, wSmooth, wEnergy,
		                                 postVelBlur, cgAccuracy, blurType,
		                                 resetBndWidth);
	else if (vel.is4D())
		opticalFlowDim<GRID, GRIDVEC, 4>(vel, i0, i1, rhsT, wSmooth, wEnergy,
		                                 postVelBlur, cgAccuracy, blurType,
		                                 resetBndWidth);
	else if (vel.is3D())
		opticalFlowDim<GRID, GRIDVEC, 3>(vel, i0, i1, rhsT, wSmooth, wEnergy,
		                                 postVelBlur, cgAccuracy, blurType,
		                                 resetBndWidth);
	else
		errMsg("Weird dimension!?");
}

// ================================================================================
// projection step & helpers

// simple (slow) uniform blur
KERNEL(fourd, bnd = 0) template <class GridT, class GridReal, class S>
void knCvExpolBlurOrg(GridT &a, GridT &tmp, GridReal &mark)
{
	Vec4i s(1, 1, 1, 1);
	if (!a.is3D()) s[2] = 0;
	if (!a.is4D()) s[3] = 0;

	S val(0.);
	Real weight = 0.;

	if (mark(i, j, k, t) != 0.) {
		// skip
	} else {
		for (int vt = t - s[3]; vt <= t + s[3]; vt++)
			for (int zk = k - s[2]; zk <= k + s[2]; zk++)
				for (int yj = j - s[1]; yj <= j + s[1]; yj++)
					for (int xi = i - s[0]; xi <= i + s[0]; xi++) {
						if (!a.isInBounds(xi, yj, zk, vt, 0)) continue;
						Real wcurr = 1.;
						weight += wcurr;
						val += wcurr * a(xi, yj, zk, vt);
					}
	}

	if (weight > VECTOR_EPSILON) {
		tmp(i, j, k, t) = val / weight;
	} else {
		tmp(i, j, k, t) = a(i, j, k, t);
	}
}
// optimized versions
KERNEL(fourd, bnd = 1) template <class GridT, class GridReal, class S>
void knCvExpolBlur4d(GridT &a, GridT &tmp, GridReal &mark)
{
	S val(0.);
	if (mark(i, j, k, t) != 0.) return;

	for (int vt = t - 1; vt <= t + 1; vt++)
		for (int zk = k - 1; zk <= k + 1; zk++)
			for (int yj = j - 1; yj <= j + 1; yj++)
				for (int xi = i - 1; xi <= i + 1; xi++) {
					val += a(xi, yj, zk, vt);
				}
	tmp(i, j, k, t) = val * (1. / 81.0);
}
KERNEL(bnd = 1) template <class GridT, class GridReal, class S>
void knCvExpolBlur3d(GridT &a, GridT &tmp, GridReal &mark)
{
	S val(0.);
	if (mark(i, j, k, 0) != 0.) return;
	for (int zk = k - 1; zk <= k + 1; zk++)
		for (int yj = j - 1; yj <= j + 1; yj++)
			for (int xi = i - 1; xi <= i + 1; xi++) {
				val += a(xi, yj, zk, 0);
			}
	tmp(i, j, k, 0) = val * (1. / 27.0);
}
KERNEL(bnd = 1) template <class GridT, class GridReal, class S>
void knCvExpolBlur2d(GridT &a, GridT &tmp, GridReal &mark)
{
	S val(0.);
	if (mark(i, j, k, 0) != 0.) return;
	for (int yj = j - 1; yj <= j + 1; yj++)
		for (int xi = i - 1; xi <= i + 1; xi++) {
			val += a(xi, yj, 0, 0);
		}
	tmp(i, j, k, 0) = val * (1. / 9.0);
}

template <class GRID, class VEC>
inline VEC getNormalInterpTempl(const GRID &phi, const VEC &pos)
{
	const Real h = 0.5;
	VEC n;
	n[0] = (phi.getInterpolated(VEC(pos[0] + h, pos[1], pos[2], pos[3])) -
	        phi.getInterpolated(VEC(pos[0] - h, pos[1], pos[2], pos[3])));
	n[1] = (phi.getInterpolated(VEC(pos[0], pos[1] + h, pos[2], pos[3])) -
	        phi.getInterpolated(VEC(pos[0], pos[1] - h, pos[2], pos[3])));

	if (phi.is3D()) {
		n[2] = (phi.getInterpolated(VEC(pos[0], pos[1], pos[2] + h, pos[3])) -
		        phi.getInterpolated(VEC(pos[0], pos[1], pos[2] - h, pos[3])));
	}

	if (phi.is4D()) {
		n[3] = (phi.getInterpolated(VEC(pos[0], pos[1], pos[2], pos[3] + h)) -
		        phi.getInterpolated(VEC(pos[0], pos[1], pos[2], pos[3] - h)));
	}

	normalize(n);
	return n;
}

template <class GRID, class GRIDVEC, class VEC>
void projectCell(Real &retD, VEC normal, GRIDVEC &dst, GRID &phiOrg,
                 GRID &phiTarget, VEC pos, VEC ijkt, int maxIter)
{
	const Real dt = dst.getParent()->getDt();
	Real step = 0.25;
	if (maxIter < 50) step = 0.5; // reduce step size...
	int lastDir = 0;
	Real targetVal = phiTarget(ijkt[0], ijkt[1], ijkt[2], ijkt[3]);
	Real vn = phiOrg.getInterpolated(pos);
	Real d = (vn - targetVal);

	for (int s = 0; s < maxIter; ++s) {

		VEC pn = pos + (normal * d) * dt;
		if (!phiOrg.isInBounds(pn)) {
			s = maxIter;
			d = 9999.;
			continue; // cause fill in
		}
		vn = phiOrg.getInterpolated(pn);
		if (vn < targetVal) {
			if (lastDir < 0) step *= 0.5;
			d += step;
			lastDir = 1;
		} else {
			if (lastDir > 0) step *= 0.5;
			d -= step;
			lastDir = -1;
		}

		// converged?
		if (step < 1e-01) s = maxIter + 1;
	}
	retD = d;
}

// project each cell along normal
KERNEL(fourd, bnd = 2) template <class GridT, class GridReal, class VEC>
void knProjectCells(GridT &dst, GridT &vel, GridReal &phiOrg,
                    GridReal &phiTarget, GridReal &marker, Real threshPhi,
                    int maxIter)
{
	const Real dt = vel.getParent()->getDt();
	const VEC pos =
	    VEC(i + 0.5f, j + 0.5f, k + 0.5f, t + 0.5f) - vel(i, j, k, t) * dt;
	VEC n = getNormalInterpTempl<GridReal, VEC>(phiOrg, pos);
	Real d = 1e10;

	projectCell<GridReal, GridT, VEC>(d, n, dst, phiOrg, phiTarget, pos,
	                                  VEC(i, j, k, t), maxIter);

	if (fabs(d) > threshPhi) {
		dst(i, j, k, t) = VEC(0.);
		return;
	} // cutoff defo distance

	dst(i, j, k, t) = n * d;
	marker(i, j, k, t) = 1.;
}

//! note - threshNorm currently unused!
template <class GRID, class GRIDVEC>
void corrVelsOfTempl(GRIDVEC &dst, GRIDVEC &vel, GRID &phiOrg, GRID &phiCurr,
                     GRID &phiTarget, Real threshPhi, Real threshNorm,
                     Real postVelBlur, Real resetBndWidth, int maxIter)
{
	typedef typename GRIDVEC::BASETYPE VEC;

	// note - phiCurr is unused! remove?

	const Real blurThreshold = 0.98;
	int doCvExtrapol = 6; // usually enough
	doCvExtrapol = vel.getMaxAbsValue() + 4;
	debMsg("Extrapol distance increased to " << doCvExtrapol, 1);

	GRID marker(vel.getParent());
	GRIDVEC tmp(vel.getParent());

	int singleStep = 0;
	// do only one step with given blur , if <threshold
	if (postVelBlur < blurThreshold) singleStep = 1;

	// while( (postVelBlur>=1.0) || singleStep)
	while ((postVelBlur >= blurThreshold) || singleStep) {
		debMsg("Performing projection, params dist: " << threshPhi << ", blur:" << postVelBlur, 1);

		dst.setConst(VEC(0.));
		marker.setConst(0);
		knProjectCells<GRIDVEC, GRID, VEC>(dst, vel, phiOrg, phiTarget, marker,
		                                   threshPhi, maxIter);

		// same as for OF:
		debMsg("Performing extrapolation:" << doCvExtrapol, 1);
		if (doCvExtrapol) {
			for (int i = 0; i < doCvExtrapol; ++i) {
				// debMsg("corrVelsOf fill in iter "<<i,1);
				tmp.copyFrom(dst);
				if (dst.is4D())
					knCvExpolBlur4d<GRIDVEC, GRID, VEC>(dst, tmp, marker);
				else if (dst.is3D())
					knCvExpolBlur3d<GRIDVEC, GRID, VEC>(dst, tmp, marker);
				else
					knCvExpolBlur2d<GRIDVEC, GRID, VEC>(dst, tmp, marker);
				dst.swap(tmp);
			}
		}

		// optionally blur result
		if (postVelBlur > VECTOR_EPSILON) {
			// gaussian blur (in line with OF call, apply with half sigma)
			gaussianBlurGeneric<GRIDVEC>(dst, 0.5 * postVelBlur);
		}

		// reset outer border (after smoothing!)
		const int resetBnd = int(resetBndWidth * dst.getSizeX()) + 1;
		if (resetBndWidth > 0.) dst.setBound(VEC(0.), resetBnd);

		// directly add to target
		vel.addScaled(dst, VEC(-1.));

		postVelBlur *= 0.5;

		// note - no reset of dst grid anymore, only done in beginning
		if (singleStep > 0) singleStep--;

	} // blur iter
}
PYTHON() void corrVelsOf3d(Grid<Vec3> &dst, Grid<Vec3> &vel, Grid<Real> &phiOrg,
                           Grid<Real> &phiCurr, Grid<Real> &phiTarget,
                           Real threshPhi = 1e10, Real threshNorm = 1e10,
                           Real postVelBlur = 0., Real resetBndWidth = -1.,
                           int maxIter = 100)
{
	corrVelsOfTempl<Grid<Real>, Grid<Vec3> >(
	    dst, vel, phiOrg, phiCurr, phiTarget, threshPhi, threshNorm,
	    postVelBlur, resetBndWidth, maxIter);
}

// ================================================================================

// some code duplication (unfortunately) to account for centered vel grids (mac
// grids used in advection.cpp)
// helper functions for centered advection, almost identical to 4d version (no
// higher order interpolation)

KERNEL(bnd = 1) template <class T>
void semiLagrangeCent3d(Grid<Vec3> &vel, Grid<T> &dst, Grid<T> &src, Real dt)
{
	Vec3 pos = Vec3(i + 0.5f, j + 0.5f, k + 0.5f) - vel(i, j, k) * dt;
	dst(i, j, k) = src.getInterpolated(pos);
}

template <class GRID> void advectSemiLagrangeCent3d(Grid<Vec3> &vel, GRID &orig)
{
	typedef typename GRID::BASETYPE T;
	GRID fwd(vel.getParent());
	semiLagrangeCent3d<T>(vel, fwd, orig, vel.getParent()->getDt());
	orig.swap(fwd);
}

PYTHON() void advectCent3d(Grid<Vec3> &vel, GridBase *grid)
{
	if (grid->getType() & GridBase::TypeReal) {
		advectSemiLagrangeCent3d<Grid<Real> >(vel, *((Grid<Real> *)grid));
	} else if (grid->getType() & GridBase::TypeVec3) {
		advectSemiLagrangeCent3d<Grid<Vec3> >(vel, *((Grid<Vec3> *)grid));
	} else
		errMsg("AdvectSemiLagrange3d: Grid Type is not supported (only Real, "
		       "Vec3)");
}

// automatically substep for given cfl condition (also converts to mac grid)
// using centered velocity grid
static void advectCflHelperCentered(Real cfl, Grid<Vec3> &vel, GridBase *grid)
{
	FluidSolver *solver = vel.getParent();
	Real orgDt = solver->getDt();
	Real maxVel = vel.getMaxValue() * orgDt;
	int steps = int(maxVel / cfl) + 1;

	solver->mDt = orgDt / Real(steps);
	for (int s = 0; s < steps; ++s) {
		advectCent3d(vel, grid);
	}
	solver->mDt = orgDt;
}

PYTHON() void advectSemiLagrangeCfl(Real cfl, FlagGrid &flags, Grid<Vec3> &vel,
                                    GridBase *grid, int order = 1,
                                    Real velFactor = 1., int orderSpace = 1)
{
	Grid<Vec3> velTmp(vel.getParent());
	velTmp.copyFrom(vel); // copy & scale...
	velTmp.multConst(Vec3(velFactor));
	advectCflHelperCentered(cfl, velTmp, grid);
}

// template glue code - choose 3d/4d advection based on template arguments
// calls python 3d version (advectSemiLagrangeCfl) for 3d grids, or custom 4d
// version
template <class GRIDFLAG, class GRIDVEC, class GRIDBASE>
void advectCflTemplate(Real cfl, GRIDFLAG *flags, GRIDVEC &vel, GRIDBASE *grid,
                       int order = 1, Real velFactor = 1., int orderSpace = 1)
{
	errMsg("advectCflTemplate - Only valid for specific instantiations");
}

template <>
void advectCflTemplate<FlagGrid, Grid<Vec3>, GridBase>(
    Real cfl, FlagGrid *flags, Grid<Vec3> &vel, GridBase *grid, int order,
    Real velFactor, int orderSpace)
{
	advectSemiLagrangeCfl(cfl, *flags, vel, grid, order, velFactor, orderSpace);
}

// ================================================================================
// multi scale version

// error metric - calculate difference between two input grids
template <class GRID>
Real calcLsDiffTempl(GRID &i0, GRID &i1, GRID *out = NULL, Real correction = 1.,
                     int bnd = 0)
{
	double accu = 0.;
	FOR_IJKT_BND(i0, bnd)
	{
		if ((i0(i, j, k, t) < 0. && i1(i, j, k, t) < 0.) ||
		    (i0(i, j, k, t) >= 0. && i1(i, j, k, t) >= 0.)) {
			// zero error , nothing to do
			if (out) (*out)(i, j, k, t) = 0.;
		} else {
			// use correction factor here to account for possible LS rescaling
			Real d = fabs(i0(i, j, k, t) - i1(i, j, k, t)) * correction;
			if (d > 1.) d = 1.; // limit difference
			accu += d;
			if (out) (*out)(i, j, k, t) = d;
		}
	}

	int sx = i0.getSizeX();
	sx -= 2 * bnd;
	int sy = i0.getSizeY();
	sy -= 2 * bnd;
	int sz = i0.getSizeZ();
	if (i0.is3D()) sz -= 2 * bnd;
	int st = i0.getSizeT();
	if (i0.is4D()) st -= 2 * bnd;
	accu *= 1000.; // debug, factor for more readable numbers...
	if (i0.getSizeT() > 1) accu *= 1000.; // even more for 4d
	accu *= 1. / double(sx * sy * sz * st);
	return accu;
}
PYTHON() Real
    calcLsDiff3d(Grid<Real> &i0, Grid<Real> &i1, Grid<Real> *out = NULL,
                 Real correction = 1., int bnd = 0)
{
	return calcLsDiffTempl<Grid<Real> >(i0, i1, out, correction, bnd);
}

// not nice: separate types... pass central 3d/4d type?
template <class GRIDINT, class GRID, class GRIDVEC>
Real opticalFlowMultiscaleTemplate(
    GRIDVEC &vel, GRID &i0, GRID &i1, GRID *rhsT = NULL, // debug info
    Real wSmooth = 0., Real wEnergy = 0., int level = 0, Real postVelBlur = 0.,
    Real cgAccuracy = 1e-04, int blurType = 1, Real cfl = 999,
    int orderTime = 1, int orderSpace = 1, Real resetBndWidth = -1.,
    int multiStep = 1, int projSizeThresh = 9999, int minGridSize = 10,
    bool doFinalProject = false )
{
	int dim[4] = { vel.getSizeX(), vel.getSizeY(), vel.getSizeZ(), 1 };
	typedef typename GRID::BASETYPE_GRID GRIDBASE;
	typedef typename GRIDVEC::BASETYPE VEC;
	FluidSolver *solver = vel.getParent();

	debMsg("Solving FlOF " << vel.getSize(), 1);
	if (solver->getDt() != 1.0) errMsg("Invalid, only dt 1 for now!");
	const int resetBnd =
	    resetBndWidth > 0 ? int(resetBndWidth * i0.getSizeX()) + 1 : 0;

	// parameters for projection step
	const Real projMaxDist = 4.; // tried 10, doesnt matter too much?
	const Real projMaxIter = 40.;
	const Real projNormThresh = -9.; // unused
	// rescale error for downscaled sdf values
	const Real lsDiffFac = 0.1 / 20.;

	// not needed on lowest level, only when advection from previous level is
	// applied
	// thus, on lowest level i0warped == i0
	GRID i0warped(solver);
	i0warped.copyFrom(i0);

	GRIDINT flagsDumm(solver);
	FOR_IJKT_BND(flagsDumm, 0)
	{
		flagsDumm(i, j, k, t) = FlagGrid::TypeFluid;
	} // dummy flags
	// pre-OF error, not overly interesting across levels - all the same...
	// only used for final eval
	Real errPreOf = calcLsDiffTempl<GRID>(i0, i1, NULL, lsDiffFac, resetBnd);

	// in practice - fairly small grids still work
	if (dim[0] > minGridSize) {

		FluidSolver small(
		    Vec3i(dim[0] / 2, dim[1] / 2, solver->is3D() ? dim[2] / 2 : 1),
		    solver->is3D() ? 3 : 2,
		    solver->has4D() ? (solver->getFourthDim() / 2) : -1);

		GRIDVEC velSm(&small);
		GRID i0Sm(&small);
		GRID i1Sm(&small);
		GRID tmpDebug(&small); // debug

		// downsample , seems to work better than with blur
		interpolGridTempl<GRID>(i0Sm, i0);
		interpolGridTempl<GRID>(i1Sm, i1);

		interpolGridTempl<GRIDVEC>(velSm, vel);
		velSm.multConst(VEC(0.5));

		// recurse
		opticalFlowMultiscaleTemplate<GRIDINT, GRID, GRIDVEC>(
		    velSm, i0Sm, i1Sm, &tmpDebug, wSmooth, wEnergy, level + 1,
		    postVelBlur, cgAccuracy, blurType, cfl, orderTime, orderSpace,
		    resetBndWidth, multiStep, projSizeThresh, minGridSize,
		    doFinalProject );

		// upsample
		interpolGridTempl<GRIDVEC>(vel, velSm);
		vel.multConst(VEC(2.));

		if (rhsT) rhsT->copyFrom(i0warped);
	}

	// pre-warp
	advectCflTemplate<GRIDINT, GRIDVEC, GRIDBASE>(
	    cfl, &flagsDumm, vel, &i0warped, orderTime, 1., orderSpace);
	i0warped.setBoundNeumann(0);

	// error before running OF on this level
	Real errCurr =
	    calcLsDiffTempl<GRID>(i0warped, i1, NULL, lsDiffFac, resetBnd);

	// use gradient projection for sizes > projSizeThreshold
	bool doProject = false;
	if (i0.getSizeX() > projSizeThresh) {
		doProject = true;
		multiStep = 1; // automatically does multiple steps...
		if (doFinalProject) doFinalProject = false; // turn off
	}

	// multi step OF
	const int MAX_STEPS = 10;
	assertMsg(multiStep < MAX_STEPS, "Too many of substeps!");
	if (multiStep > 1) {
		GRIDVEC *vs[MAX_STEPS];
		for (int of = 0; of < multiStep; ++of)
			vs[of] = NULL;

		GRIDVEC *vs2[MAX_STEPS];
		for (int of = 0; of < multiStep; ++of)
			vs2[of] = new GRIDVEC(i0.getParent());
		GRIDVEC tmpVel(vel); // velocity acculumator

		Real velBlur = postVelBlur;
		Real errLast = errCurr;
		int ofStepsCurr = multiStep;
		for (int of = 0; of < ofStepsCurr; ++of) {
			// n'th solve
			vs[of] = new GRIDVEC(i0.getParent()); // zero init
			// vs[of] = new GRIDVEC(tmpVel); // pass previous velocity for
			// smoothness, slightly worse in practice

			opticalFlowTemplate<GRID, GRIDVEC>(
			    *vs[of], i0warped, i1, rhsT, wSmooth, wEnergy, velBlur,
			    cgAccuracy, blurType, resetBndWidth);

			velBlur *= (3. / 4.); // reduce while iterating
			if (velBlur < 2.) velBlur = 2.;

			// prepare final advection , note - lots of copies necessary...
			// backup copy
			for (int k = of; k >= 0; --k) {
				vs2[k]->copyFrom(*vs[k]);
			}
			// readvect (align for current iteration of OF)
			for (int k = of - 1; k >= 0; --k) {
				for (int l = 0; l < k; ++l) {
					advectCflTemplate<GRIDINT, GRIDVEC, GRIDBASE>(
					    cfl, &flagsDumm, *vs2[k], vs2[l], orderTime, 1.,
					    orderSpace);
				}
			}
			// sum up
			tmpVel.copyFrom(vel);
			for (int k = of; k >= 0; --k) {
				tmpVel.add(*vs2[k]);
			}
			// advect
			GRID i0warp2(solver); // note, could be removed...
			i0warp2.copyFrom(i0);
			advectCflTemplate<GRIDINT, GRIDVEC, GRIDBASE>(
			    cfl, &flagsDumm, tmpVel, &i0warp2, orderTime, 1., orderSpace);
			i0warp2.setBoundNeumann(0);

			// calc error for curr version
			Real errC =
			    calcLsDiffTempl<GRID>(i0warp2, i1, NULL, lsDiffFac, resetBnd);
			debMsg("Current error, s" << i0warp2.getSizeX() << " " << of
			                          << " = " << errC << " , fac "
			                          << (errC / errLast), 1);
			i0warped.copyFrom(i0warp2);
			// stop?
			if (of > 0 && (errC / errLast) > 0.95) {
				// discard result, and stop
				vs[of]->setConst(VEC(0.));
				ofStepsCurr = of + 1;
			}

			errLast = errC;
		}

		// final alignment of completed vels, advect velocities forward with
		// others
		for (int of = ofStepsCurr - 1; of >= 0; --of) {
			for (int l = 0; l < of; ++l) {
				advectCflTemplate<GRIDINT, GRIDVEC, GRIDBASE>(
				    cfl, &flagsDumm, *vs[of], vs[l], orderTime, 1., orderSpace);
			}
		}

		for (int of = 0; of < ofStepsCurr; ++of) {
			vel.add(*vs[of]);
			delete vs[of];
		}
		for (int of = 0; of < multiStep; ++of) {
			delete vs2[of];
		}

		// end multi step
	} else {
		// regular, single solve
		GRIDVEC velCurr(i0.getParent()); // zero init
		// GRIDVEC velCurr( vel ); // pass previous velocity for smoothness,
		// slightly worse in practice
		if (!doProject) {
			opticalFlowTemplate<GRID, GRIDVEC>(
			    velCurr, i0warped, i1, rhsT, wSmooth, wEnergy, postVelBlur,
			    cgAccuracy, blurType, resetBndWidth);
			vel.add(velCurr);
		} else {
			// we need a second tmp grid to accumulate projection multi step
			// note - no optical flow solve here, only project! testing
			// version... not a good alternative in practice
			GRIDVEC velTmp2(i0.getParent());
			Real blur = postVelBlur;
			debMsg("Velocity correction only enabled , no of solve! " << blur, 1);
			corrVelsOfTempl<GRID, GRIDVEC>(velCurr, velTmp2, i0warped, i0warped,
			                               i1, projMaxDist, projNormThresh,
			                               blur, resetBndWidth, projMaxIter);
			vel.add(velTmp2);
			velCurr.multConst(VEC(-1.)); // "dirty" accumulation for now...
			vel.add(velCurr);
		}
	}

	// final projection step?
	if ((level == 0) && doFinalProject) {
		GRIDVEC velCurr(i0.getParent());

		Real finalProjBlur = postVelBlur * 0.5; // reduce blur for this step
		finalProjBlur = 2.;                     // old
		finalProjBlur = 4.;                     // new default
		corrVelsOfTempl<GRID, GRIDVEC>(velCurr, vel, i0, i0, i1, projMaxDist,
		                               projNormThresh, finalProjBlur,
		                               resetBndWidth, projMaxIter);
	}

	// readvect final , calculate error for highest level
	if (level == 0) {
		i0warped.copyFrom(i0);
		advectCflTemplate<GRIDINT, GRIDVEC, GRIDBASE>(
		    cfl, &flagsDumm, vel, &i0warped, orderTime, 1., orderSpace);
		i0warped.setBoundNeumann(0);
		Real errFinal =
		    calcLsDiffTempl<GRID>(i0warped, i1, NULL, lsDiffFac, resetBnd);

		// calculate and output per level error values
		// return this error to caller
		errCurr = errFinal;
		debMsg("Final error=" << errFinal << ", factor "
		                        << (errFinal / (errPreOf + VECTOR_EPSILON))
		                        << " ", 1);
	}

	return errCurr;
}

PYTHON() void opticalFlowMultiscale3d(
    Grid<Vec3> &vel, Grid<Real> &i0, Grid<Real> &i1,
    Grid<Real> *rhsT = NULL, // debug info
    Real wSmooth = 0., Real wEnergy = 0., int level = 0, Real postVelBlur = 0.,
    Real cgAccuracy = 1e-04, int blurType = 1, Real cfl = 999,
    int orderTime = 1, int orderSpace = 1, Real resetBndWidth = -1.,
    int multiStep = 1, int projSizeThresh = 9999, int minGridSize = 10,
    bool doFinalProject = false) 
{
	opticalFlowMultiscaleTemplate<FlagGrid, Grid<Real>, Grid<Vec3> >(
	    vel, i0, i1, rhsT, wSmooth, wEnergy, level, postVelBlur, cgAccuracy,
	    blurType, cfl, orderTime, orderSpace, resetBndWidth, multiStep,
	    projSizeThresh, minGridSize, doFinalProject );
}

// ================================================================================

// smaller helpers

// copy 2d slices into 3d grid
PYTHON() void placeGrid2d(Grid<Real> &src, Grid<Real> &dst, int dstz)
{
	FOR_IJK(src)
	{
		if (!dst.isInBounds(Vec3i(i, j, dstz))) continue;
		dst(i, j, dstz) = src(i, j, 0);
	}
}
// and back 3d plane into 2d buffer
PYTHON() void getSliceFrom3d(Grid<Real> &src, int srcz, Grid<Real> &dst,
                             int dstz = 0)
{
	const int bnd = 0;
	for (int j = bnd; j < src.getSizeY() - bnd; j++)
		for (int i = bnd; i < src.getSizeX() - bnd; i++) {
			if (!dst.isInBounds(Vec3i(i, j, dstz))) continue;
			dst(i, j, dstz) = src(i, j, srcz);
		}
}
// repeat start frame at srcz, n times
PYTHON() void repeatStartFrame(Grid<Real> &phi, int srcz, int n = 0,
                               int bnd = 0)
{
	errMsg("Use repeatFrame instead!");
}

PYTHON() void repeatFrame(Grid<Real> &phi, Real srcz, Real range = 0,
                          int bnd = 0)
{
	// repeat 1st
	for (int k = int(srcz - range + 0.5); k <= int(srcz); k++) {
		for (int j = bnd; j < phi.getSizeY() - bnd; j++)
			for (int i = bnd; i < phi.getSizeX() - bnd; i++) {
				phi(i, j, k) =
				    phi.getInterpolated(Vec3(i, j, srcz) + Vec3(0.5));
			}
	}
}

// helper to "remove" levelset and init const values
PYTHON() void constOverride(Grid<Real> &phi, Real val = 1.)
{
	FOR_IJK(phi)
	{
		if (phi(i, j, k) > 0.)
			phi(i, j, k) = -val;
		else
			phi(i, j, k) = val;
	}
}

// ================================================================================
// 4d code

// ================================================================================
// 4d helper functions

// corresponding place/slice funcs for 3d/4d

PYTHON() void repeatFrame4d(Grid4d<Real> &phi, Real srct, Real range = 0,
                            int bnd = 0)
{
	// repeat 1st
	for (int t = int(srct - range + 0.5); t <= int(srct); t++) {
		for (int k = bnd; k < phi.getSizeZ() - bnd; ++k) 
			for (int j = bnd; j < phi.getSizeY() - bnd; j++)
				for (int i = bnd; i < phi.getSizeX() - bnd; i++) {
					phi(i, j, k, t) =
					    phi.getInterpolated(Vec4(i, j, k, srct) + Vec4(0.5));
				}
	}
}
PYTHON() void repeatStartFrame4d(Grid4d<Real> &phi, int srct, int n = 0,
                                 int bnd = 0)
{
	errMsg("Use repeatFrame instead!");
}

//! First order advection for 4d grids

KERNEL(fourd, bnd = 1) template <class T>
void semiLagrange4d(Grid4d<Vec4> &vel, Grid4d<T> &dst, Grid4d<T> &src, Real dt)
{
	Vec4 pos =
	    Vec4(i + 0.5f, j + 0.5f, k + 0.5f, t + 0.5f) - vel(i, j, k, t) * dt;
	dst(i, j, k, t) = src.getInterpolated(pos);
}

template <class GRID>
void advectSemiLagrange4d(Grid4d<Vec4> &vel, GRID &orig, Real dtFac)
{
	typedef typename GRID::BASETYPE T;
	GRID fwd(vel.getParent());
	semiLagrange4d<T>(vel, fwd, orig, vel.getParent()->getDt() * dtFac);
	orig.swap(fwd);
}

PYTHON() void advect4d(Grid4d<Vec4> &vel, Grid4dBase *grid, Real dtFac = 1.)
{
	if (grid->getType() & Grid4dBase::TypeReal) {
		advectSemiLagrange4d<Grid4d<Real> >(vel, *((Grid4d<Real> *)grid),
		                                    dtFac);
	} else if (grid->getType() & Grid4dBase::TypeVec3) {
		advectSemiLagrange4d<Grid4d<Vec3> >(vel, *((Grid4d<Vec3> *)grid),
		                                    dtFac);
	} else if (grid->getType() & Grid4dBase::TypeVec4) {
		advectSemiLagrange4d<Grid4d<Vec4> >(vel, *((Grid4d<Vec4> *)grid),
		                                    dtFac);
	} else
		errMsg("AdvectSemiLagrange4d: Grid Type is not supported (only Real, "
		       "Vec3, Vec4)");
}

// almost identical to 3d version
static void advectCflHelper4d(Real cfl, Grid4d<Vec4> &vel, Grid4dBase *grid)
{
	FluidSolver *solver = vel.getParent();
	Real orgDt = solver->getDt();
	Real maxVel = vel.getMaxValue() * orgDt;
	int steps = int(maxVel / cfl) + 1;

	solver->mDt = orgDt / Real(steps);
	for (int s = 0; s < steps; ++s) {
		advect4d(vel, grid);
	}
	solver->mDt = orgDt;
}

static const Vec4i nbs4d[8] = { Vec4i(-1, 0, 0, 0), Vec4i(1, 0, 0, 0),
	                            Vec4i(0, -1, 0, 0), Vec4i(0, 1, 0, 0),
	                            Vec4i(0, 0, -1, 0), Vec4i(0, 0, 1, 0),
	                            Vec4i(0, 0, 0, -1), Vec4i(0, 0, 0, 1) };

//! simple extrapolation step 4d

KERNEL(fourd, bnd = 1) template <class S>
void knExtrap4dLsSimple(Grid4d<S> &val, int distance, Grid4d<int> &tmp,
                        const int d, S direction)
{
	const int dim = 4;
	if (tmp(i, j, k, t) != 0) return;

	// copy from initialized neighbors
	Vec4i p(i, j, k, t);
	int nbs = 0;
	S avg = 0.;
	for (int n = 0; n < 2 * dim; ++n) {
		if (tmp(p + nbs4d[n]) == d) {
			avg += val(p + nbs4d[n]);
			nbs++;
		}
	}

	if (nbs > 0) {
		tmp(p) = d + 1;
		val(p) = avg / nbs + direction;
	}
}

KERNEL(fourd, bnd = 1) template <class S>
void knSetRemaining4d(Grid4d<S> &phi, Grid4d<int> &tmp, S set)
{
	if (tmp(i, j, k, t) != 0) return;
	phi(i, j, k, t) = set;
}

PYTHON() void extrap4dLsSimple(Grid4d<Real> &phi, int distance = 4,
                               bool inside = false)
{
	Grid4d<int> tmp(phi.getParent());
	tmp.clear();
	const int dim = 4;

	// by default, march outside
	Real direction = 1.;
	if (!inside) {
		FOR_IJKT_BND(phi, 1)
		{
			if (phi(i, j, k, t) < 0.) {
				tmp(i, j, k, t) = 1;
			}
		}
	} else {
		direction = -1.;
		FOR_IJKT_BND(phi, 1)
		{
			if (phi(i, j, k, t) > 0.) {
				tmp(i, j, k, t) = 1;
			}
		}
	}
	// + first layer around
	FOR_IJKT_BND(phi, 1)
	{
		Vec4i p(i, j, k, t);
		if (tmp(p)) continue;
		for (int n = 0; n < 2 * dim; ++n) {
			if (tmp(p + nbs4d[n]) == 1) {
				tmp(i, j, k, t) = 2;
				n = 2 * dim;
			}
		}
	}

	// extrapolate for distance
	for (int d = 2; d < 1 + distance; ++d) {
		knExtrap4dLsSimple<Real>(phi, distance, tmp, d, direction);
	}

	// set all remaining cells to max
	knSetRemaining4d<Real>(phi, tmp, Real(direction * (distance + 2)));
}

PYTHON() void extrapolateVec4Simple(Grid4d<Vec4> &vel, Grid4d<Real> &phi,
                                    int distance)
{
	const int dim = 4;
	Grid4d<int> tmp(vel.getParent());
	tmp.clear();

	// mark all inside
	FOR_IJKT_BND(vel, 1)
	{
		if (phi(i, j, k, t) < 0.) {
			tmp(i, j, k, t) = 1;
		}
	}
	// + first layer outside
	FOR_IJKT_BND(vel, 1)
	{
		Vec4i p(i, j, k, t);
		if (tmp(p)) continue; // skip any non-initialized ones
		for (int n = 0; n < 2 * dim; ++n) {
			if (tmp(p + nbs4d[n]) == 1) {
				tmp(i, j, k, t) = 2;
				n = 2 * dim;
			}
		}
	}

	for (int d = 2; d < distance + 1; ++d) {
		knExtrap4dLsSimple<Vec4>(vel, distance, tmp, d, Vec4(0.));
	}
	knSetRemaining4d<Vec4>(vel, tmp, Vec4(0.));
}

// ================================================================================
// load 3d slices & place into 4d grid

// check no fourd kernel?
KERNEL(bnd = 1)
void knLoadPlaceGridIpol(Grid<Real> &dummy, Grid4d<Real> &phi, Grid<Real> &tmp,
                         const Vec3 &sourceFactor3, const Vec3 &offset3,
                         const int slice, const Real wc)
{
	Vec3 p = Vec3(i, j, k) * sourceFactor3 + offset3;
	phi(i, j, k, slice) += wc * tmp.getInterpolated(Vec3(p[0], p[1], p[2]));
}

KERNEL(fourd, bnd = 1)
void knLoadPlaceGridRescale(Grid4d<Real> &phi, const std::vector<Real> &w)
{
	if (w[t] <= 0.) return;
	phi(i, j, k, t) *= w[t];
}

// read & place a time slices from file, interpol to dst size, then do advection lookup in source-LS
// overall, similar to interpolateGrid - place grid with offset & size
// fileIdx range is rescaled to target size.t , if both are -1 assume 0..size.t  range
PYTHON() void loadPlaceGrid4d(std::string fname, Grid4d<Real> &phi, Vec4 offset,
                              Vec4 scale, int fileIdxStart = -1,
                              int fileIdxEnd = -1, int debugSkipLoad = 999999,
                              Real spread = 1., Vec4 overrideSize = Vec4(-1.),
                              Real overrideTimeOff = 0.,
                              int overrideGoodRegion = 0,
                              Real loadTimeScale = 1.,
                              bool rescaleSdfValues = false,
                              Real sdfIsoOff = 0., Real repeatStartFrame = 0.)
{
	// init size , get from first data file (index 0)
	Vec3i dim(0);
	int dimT = 0;
	char fnCur[1024];
	snprintf(fnCur, 1024, fname.c_str(), fileIdxStart);
	getUniFileSize(fnCur, dim[0], dim[1], dim[2], &dimT);
	if (dim[0] < 1 || dim[1] < 1 || dim[2] < 1) {
		errMsg("Invalid src size " << dim << " from " << fname << " ");
	}
	debMsg("Found size " << dim << "," << dimT << " in " << fname, 1);
	FluidSolver solver(Vec3i(dim[0], dim[1], dim[2]), 3, 1); // fourd=1 , one time slice at a time!
	Grid<Real> tmp(&solver);

	// workaround for looping over dst region in 3d, content doesnt matter!
	Grid<Real> dummy(phi.getParent());

	// rename to same as in loadAdvectTimeSlice
	Vec4i v1(dim[0], dim[1], dim[2], 1);
	Vec4 defoSize = toVec4(phi.getSize());
	if (overrideSize.x > 0.) defoSize = overrideSize;
	Vec4 defoOffset = offset;

	// note - size & rescaling identical to loadAdvectTimeSlice
	Vec4 sourceFactor, off2 = defoOffset;
	gridFactor4d(toVec4(v1), defoSize, Vec4(-1.), scale, sourceFactor, off2);
	Vec3 sourceFactor3 = Vec3(sourceFactor.x, sourceFactor.y, sourceFactor.z);
	Vec3 offset3 = Vec3(off2.x, off2.y, off2.z);

	// rescale sdf values if necessary
	Real valueScale = 1.;
	if (rescaleSdfValues) {
		valueScale = 1. / sourceFactor.x; // already determines spatial scaling
	}

	std::vector<Real> w(phi.getSizeT()); // can stay actual size! keep weights
	                                     // in local coordinates
	for (int slice = 0; slice < (int)w.size(); ++slice) {
		w[slice] = 0.;
	}

	// reset for accumulation
	if (overrideGoodRegion <= 0) {
		phi.clear();
	} else {
		memset(&phi(0, 0, 0, overrideGoodRegion), 0,
		       sizeof(Real) * phi.getSizeX() * phi.getSizeY() * phi.getSizeZ() *
		           (phi.getSizeT() - overrideGoodRegion));
	}

	// use 0 -> gridsize.t if unspecified
	if (fileIdxStart < 0) fileIdxStart = 0;
	if (fileIdxEnd < 0) fileIdxEnd = defoSize.t;

	// repeat first/start frame like on 4d load
	int repeatOff = 0;
	if (repeatStartFrame > 0) {
		repeatOff = int((fileIdxEnd - fileIdxStart) * repeatStartFrame + 1.0);
		debMsg("loadPlaceGrid4d: repeating first frame "
		           << fileIdxStart << " for " << (repeatOff), 1);
	}

	// loop over file range
	for (int fileid = fileIdxStart - repeatOff;
	     fileid < std::min(fileIdxEnd, debugSkipLoad); ++fileid) {
		// load only on demand
		int fileidClamp = fileid > fileIdxStart ? fileid : fileIdxStart;
		snprintf(fnCur, 1024, fname.c_str(), fileidClamp);
		bool didLoad = false; // to activate tmp.load( fnCur);

		// note - duration has to be adjusted by time scaling...
		Real duration = 1. * (defoSize.t * scale.t / loadTimeScale) /
		                Real(fileIdxEnd - fileIdxStart);
		Real dstt =
		    offset.t + overrideTimeOff + duration * (fileid - fileIdxStart);
		debMsg("loadPlaceGrid4d "
		           << fileid << " file:" << fnCur << ", to " << dstt << " for "
		           << duration << ", spread" << spread << "   off" << defoOffset
		           << " s" << defoSize << " scale=" << scale << ", file idx "
		           << fileIdxStart << " " << fileIdxEnd, 2);
		for (int slice = int(dstt - spread * duration) + 0;
		     slice < int(dstt + spread * duration) + 2; ++slice) {
			// for first call, init all (overrideGoodRegion=0), then continually
			// update later part of domain
			if (slice < overrideGoodRegion || slice >= phi.getSizeT()) continue;
			Real wc = 1. - fabs(Real(slice) - dstt) / (duration * spread);
			if (wc < VECTOR_EPSILON) continue;
			debMsg("             placing #" << fileid << " -> t:" << slice
			                                << ", w=" << wc, 5);
			w[slice] += wc;

			// note - dummy has to match 3d size of phi 4d grid!
			if (!didLoad) {
				tmp.load(fnCur);

				// fix up level sets...
				if (rescaleSdfValues && sdfIsoOff != 0.) {
					// note - fix outer side, but no 0.1*sizeX offset necessary
					// here! fulll data...
					int bnd = 1;
					FOR_IJK_BND(tmp, bnd)
					{
						tmp(i, j, k) += sdfIsoOff;
					}                      // offset only inner values
					tmp.setBound(1., bnd); // set outside value at border
				}
				didLoad = true;
			}
			knLoadPlaceGridIpol(dummy, phi, tmp, sourceFactor3, offset3, slice,
			                    wc);
		}
	}
	// normalize contributions
	for (int slice = 0; slice < (int)w.size(); ++slice) {
		if (w[slice] < 1e-03) {
			w[slice] = 0.;
		} else {
			w[slice] = 1. / w[slice];
			w[slice] *= valueScale; // scale SDF values if necessary
		}
	}
	knLoadPlaceGridRescale(phi, w);
}

// shift grid forward, to refresh partial load
PYTHON() void shiftForwGrid4d(Grid4d<Real> &phi, int overrideGoodRegion = 0)
{
	Vec4i dim = phi.getSize();
	FluidSolver solver(Vec3i(dim[0], dim[1], dim[2]), 3, 1); // fourd=1!
	Grid<Real> tmp(&solver);
	const int shift = dim.t - overrideGoodRegion;

	// loop in blocks of shift size
	debMsg("Shifting for " << dim.t << " , region " << overrideGoodRegion
	                       << " , by " << shift, 1);
	for (int slice = 0; slice < overrideGoodRegion; slice += shift) {
		int len = shift;
		if (slice + shift >= overrideGoodRegion)
			len = overrideGoodRegion - slice;

		memcpy(&phi(0, 0, 0, slice), &phi(0, 0, 0, slice + shift),
		       sizeof(Real) * dim.x * dim.y * dim.z * len);
	}
}

// ================================================================================
// load vel & advect

// special semi lagrange lookup for sliced velocity load
KERNEL(fourd, bnd = 1)
void knSemiLagrangeLookupSlice4d(Grid4d<Vec4> &velSlice, Grid<Real> &dst,
                                 Grid4d<Real> &src, Real time, Real dt)
{
	Vec4 pos4 = Vec4(i + 0.5f, j + 0.5f, k + 0.5f, time + 0.5f) -
	            velSlice(i, j, k, 0) * dt;
	dst(i, j, k, 0) = src.getInterpolated(pos4);
}

// interpolate velocity (note - not really parallel)
KERNEL(bnd = 1)
void loadAdvReInterpol(Grid4d<Vec4> &vdst, Grid4d<Vec4> &v1,
                       const Vec3 &sourceFactor3, const Vec3 &offset3,
                       const Vec4 &defoFactor, const Vec3i &dim)
{
	// interpolate in 3d space
	Vec3 pos3 = (Vec3(i, j, k) * sourceFactor3 + offset3);
	vdst(i, j, k, 0) = interpol<Vec4>(&v1[0], dim, v1.getStrideZ(), pos3);

	// scale magnitude
	vdst(i, j, k, 0) *= defoFactor;
}

// special semi lagrange lookup for sliced velocity load , interpolate velocity
// on the fly
// note, skips 10 by default! assuming at least 100^3 with 10% skip
KERNEL(bnd = 10)
void knSemiLagrangeLookupSlice4d_WithVel(Grid<Real> &dst, Grid4d<Real> &src,
                                         Real time, Real dt, Grid4d<Vec4> &v1,
                                         const Vec3 &sourceFactor3,
                                         const Vec3 &offset3,
                                         const Vec4 &defoFactor,
                                         const Vec3i &dim, int bordSkip)
{
	if (!dst.isInBounds(Vec3i(i, j, k), bordSkip)) return;

	// interpolate velocity in 3d space
	Vec3 pos3 = (Vec3(i, j, k) * sourceFactor3 + offset3);
	Vec4 vdst = interpol<Vec4>(&v1[0], dim, v1.getStrideZ(), pos3);
	vdst *= defoFactor;

	Vec4 pos4 = Vec4(i + 0.5f, j + 0.5f, k + 0.5f, time + 0.5f) - vdst * dt;
	dst(i, j, k, 0) = src.getInterpolated(pos4);
}

// read & init 2 slices from file, interpol to dst size, then do advection lookup in source-LS
// (slower version without caching and global vars; dummyID unused - just for compatiblity)
// note - some params like thirdAlpha and bordSkip not supported! only for
// compatbility with loadAdvectTimeSlice_opt
PYTHON() void loadAdvectTimeSlice(
    int dummyID, std::string fname, Grid<Real> &dst, Grid4d<Real> &phi,
    Real time, Real blendAlpha, Real loadTimeScale, Vec4 defoOffset,
    Vec4 defoScale, Vec4 defoFactor, Vec4 overrideSize = Vec4(-1.),
    Real overrideTimeOff = 0., // for shifted / partial loading
    Grid<Vec3> *debugVel = NULL, Grid<Real> *debugVelT = NULL,
    bool zeroVel = false, Real thirdAlpha = 0., int bordSkip = 1,
    Real fourthAlpha = 0., float defoAniFac = 1.)
{
	// read from file, then round
	Vec3i dim(0);
	int dimT = 0;
	getUniFileSize(fname, dim[0], dim[1], dim[2], &dimT);
	if (dim[0] < 1 || dim[1] < 1 || dim[2] < 1) {
		errMsg("Invalid src size " << dim << " from " << fname << " ");
	}
	debMsg("Found size " << dim << "," << dimT << " in " << fname, 1);

	// create custom sized 4d data (1 time slice)
	FluidSolver solver(Vec3i(dim[0], dim[1], dim[2]), 3, 1); // one time slice
	Grid4d<Vec4> v1(&solver), v2(&solver);

	Vec4 dim4(dim[0], dim[1], dim[2], dimT);
	Vec4 sourceFactor(0.), off2 = defoOffset;
	Vec4i defoSize = phi.getSize();
	if (overrideSize.x > 0.) {
		// truncate to ints
		defoSize = Vec4i(overrideSize.x, overrideSize.y, overrideSize.z,
		                 overrideSize.t);
	}
	gridFactor4d(toVec4(dim4), toVec4(defoSize), Vec4(-1.), defoScale,
	             sourceFactor, off2);
	Vec3 sourceFactor3 = Vec3(sourceFactor.x, sourceFactor.y, sourceFactor.z);
	Vec3 offset3 = Vec3(off2.x, off2.y, off2.z);

	// load 2 slices (time is given in dst units, convert to src)
	Real srcTime = (time)*sourceFactor.t * loadTimeScale + off2.t - 0.5; 
	int t = (int)srcTime;
	int tp1 = t + 1;
	Real tw = srcTime - (Real)t;
	t = std::min(t, dimT - 1); // make sure we dont go too far...
	tp1 = std::min(tp1, dimT - 1);
	debMsg("Reading velocity at t=" << t << " w=" << tw << "    sf "
	                                << sourceFactor << " " << off2, 1);
	readGrid4dUni<Vec4>(fname, NULL, t, &v1);
	readGrid4dUni<Vec4>(fname, NULL, tp1, &v2);

	// linearly blend in time
	v1.multConst(Vec4(1. - tw));
	v2.multConst(Vec4(tw));
	v1.add(v2);
	if (zeroVel) {
		debMsg("Debug - zeroing deformation!", 1);
		v1.setConst(Vec4(0.));
	}

	// re-interpol to size of dst?  problem - Grid<Vec4> not supported!
	FluidSolver solver2(dst.getSize(), 3, 1); // one time slice
	Grid4d<Vec4> vdst(&solver2);
	vdst.clear();

	// possibly: interpolate on the fly -> knSemiLagrangeLookupSlice4d_WithVel
	Real mint = 1e10, maxt = -1e10;
	FOR_IJK_BND(vdst, 1)
	{
		// interpolate in 3d space
		Vec3 pos3 = (Vec3(i, j, k) * sourceFactor3 + offset3);
		vdst(i, j, k, 0) = interpol<Vec4>(&v1[0], dim, v1.getStrideZ(), pos3);

		// scale magnitude
		vdst(i, j, k, 0) *= defoFactor * defoAniFac;
		Real val = vdst( i, j, k, 0)[3]; 
		mint = std::min(val, mint);
		maxt = std::max(val, maxt);
	}

	// return vel data for debugging
	if (debugVel) {
		FOR_IJK(*debugVel) { (*debugVel)(i, j, k) = toVec3(vdst(i, j, k, 0)); }
	}
	if (debugVelT) {
		FOR_IJK(*debugVelT) { (*debugVelT)(i, j, k) = vdst(i, j, k, 0)[3]; }
	}

	// run special advect for lookup
	debMsg("Performing SL lookup, alpha=" << blendAlpha << ". T-range: " << mint
	                                      << " to " << maxt, 1);
	knSemiLagrangeLookupSlice4d(vdst, dst, phi, time + overrideTimeOff,
	                            blendAlpha);
}

PYTHON() void calcObfDiff(Grid<Real> &phi1, Grid<Real> &phi2,
                          Grid<Real> &phiDiff, Grid<Vec3> &vel1,
                          Grid<Vec3> &vel2, Grid<Real> &velt1,
                          Grid<Real> &velt2, Grid<Real> &velDiff, int bnd)
{
	FOR_IJK_BND(phiDiff, bnd)
	{
		phiDiff(i, j, k) = fabs(phi1(i, j, k) - phi2(i, j, k));
		Vec4 d;
		d[0] = vel1(i, j, k)[0] - vel2(i, j, k)[0];
		d[1] = vel1(i, j, k)[1] - vel2(i, j, k)[1];
		d[2] = vel1(i, j, k)[2] - vel2(i, j, k)[2];
		d[3] = velt1(i, j, k) - velt2(i, j, k);
		velDiff(i, j, k) = norm(d);
	}
}

// current max nr of defo volumes
static const int MAXDV = 10;

// min vel window of 20% by default, max'ed by overrideFactor/partialload from script
static const Real DEFOVOLWIDTH = 0.2;

//! optimized version with global vars
class LoadAdvectData
{
  public:
	Vec3i dim;
	int dimT;
	FluidSolver *solver;
	Grid4d<Vec4> *v1;
	Grid4d<Vec4> *v2;
	Grid4d<Vec4> *tmp;
	int lastT;
	void *fileHandle1;
	std::string fname1;
	// vol defo data
	bool useDefoVols;
	FluidSolver *defovolSolver;
	int defovolOff;
	int numDv;
	Grid4d<Vec4> *defovol[MAXDV];
	void *fileHandle[MAXDV];
	std::string fname[MAXDV];
	Grid4d<Vec4> *defovolTmp;
	bool doAligned;

	LoadAdvectData()
	    : dim(0), dimT(0), solver(NULL), v1(NULL), v2(NULL), tmp(NULL),
	      lastT(-1), fileHandle1(NULL), fname1("unin1"), useDefoVols(true),
	      defovolSolver(NULL), defovolOff(0), defovolTmp(NULL), doAligned(false)
	{
		for (int i = 0; i < MAXDV; ++i) {
			defovol[i] = NULL;
			fileHandle[i] = NULL;
			fname[i] = std::string("unin2");
		}
	};

	void updateDefoVol(int t)
	{
		// slice / volume defo load
		Grid4d<Vec4> &vt = (*this->tmp);
		int currt = this->defovol[0]->getSizeT() / 2;
		this->defovolOff = t - currt;
		for (int dv = 0; dv < this->numDv; ++dv) {
			Grid4d<Vec4> &sl = (*this->defovol[dv]);
			if (this->lastT == t) {
				// all good
			} else if (this->lastT + 1 == t) {
				debMsg("Defovol " << dv << " reading single defovol vel at "
				                  << t, 1);
				// shift all back
				int tl = 0;
				for (; tl < sl.getSizeT() - 1; ++tl) {
					FOR_IJK_BND(vt, 0)
					{
						sl(i, j, k, tl) = sl(i, j, k, tl + 1);
					}
				}
				tl = sl.getSizeT() - 1;
				int it = std::min(std::max(t - currt + tl, 0),
				                  dimT - 1); // make sure we dont go too far...

				readGrid4dUni<Vec4>(this->fname[dv], NULL, it, &vt,
				                    &this->fileHandle[dv]);
				FOR_IJK_BND(vt, 0) { sl(i, j, k, tl) = vt(i, j, k, 0); } // copy
			} else {
				debMsg("Defovol " << dv << " reading all defovol vels at " << t, 1);
				for (int tl = 0; tl < sl.getSizeT(); ++tl) {
					int it =
					    std::min(std::max(t - currt + tl, 0),
					             dimT - 1); // make sure we dont go too far...
					readGrid4dUni<Vec4>(this->fname[dv], NULL, it, &vt); // tmp
					FOR_IJK_BND(vt, 0)
					{
						sl(i, j, k, tl) = vt(i, j, k, 0);
					} // copy
				}
			}
		} // dv
	}
};

// global storage of advection data sets (with integer ID)
static std::map<int, LoadAdvectData *> LATS;

//! initialize data for load&advect, note - lots of unused params, just for
// convenience
PYTHON() void loadAdvectTimeSlice_OptInit(int ID, std::string fname1,
                                          bool useDefoVols, bool doAligned,
                                          Real partialLoadFac = DEFOVOLWIDTH)
{
	// read from file, then round
	Vec3i dim(0);
	int dimT = 0;
	getUniFileSize(fname1, dim[0], dim[1], dim[2], &dimT);
	if (dim[0] < 1 || dim[1] < 1 || dim[2] < 1) {
		errMsg("Invalid src size " << dim << " from " << fname1 << " ");
	}
	debMsg("Found size " << dim << "," << dimT << " in " << fname1, 1);
	if (!LATS[ID]) {
		LATS[ID] = new LoadAdvectData();
	}
	LoadAdvectData &lats = *LATS[ID];

	// create custom sized 4d data (1 time slice)
	lats.solver =
	    new FluidSolver(Vec3i(dim[0], dim[1], dim[2]), 3, 1); // one time slice
	lats.v1 = new Grid4d<Vec4>(lats.solver);
	lats.v2 = new Grid4d<Vec4>(lats.solver);
	lats.tmp = new Grid4d<Vec4>(lats.solver);
	lats.dim = dim;
	lats.dimT = dimT;
	lats.fname1 = fname1;
	lats.useDefoVols = useDefoVols;

	if (lats.useDefoVols) {
		Real defoVolWidth = std::max(DEFOVOLWIDTH, partialLoadFac);
		lats.defovolSolver =
		    new FluidSolver(Vec3i(dim[0], dim[1], dim[2]), 3,
		                    int(dimT * defoVolWidth)); // vel slices
		lats.defovolTmp = new Grid4d<Vec4>(lats.defovolSolver);
		lats.defovol[0] = new Grid4d<Vec4>(lats.defovolSolver);
		lats.fname[0] = fname1;
		debMsg("Created defovol solver , " << int(dimT * 0.2) << " , defo0 "
		                                   << fname1,
		       3);
		lats.numDv = 1;
		lats.doAligned = doAligned;
	}
}
PYTHON() void loadAdvectTimeSlice_OptAdd(int ID, std::string fname)
{
	if (!LATS[ID]) {
		return;
	}
	LoadAdvectData &lats = *LATS[ID];
	// add another deformation data set volume
	if (lats.useDefoVols) {
		lats.defovol[lats.numDv] = new Grid4d<Vec4>(lats.defovolSolver);
		lats.fname[lats.numDv] = fname;
		lats.numDv++;
		assertMsg((lats.numDv < MAXDV), "Too many defovolumes loaded!");
		debMsg("Added defovol " << lats.numDv << " file , " << fname, 3);
	}
}
//! cleanup load advect data
PYTHON() void loadAdvectTimeSlice_Finish(int ID)
{
	if (!LATS[ID]) {
		return;
	}
	LoadAdvectData &lats = *LATS[ID];
	delete lats.v1;
	delete lats.v2;
	delete lats.tmp;
	delete lats.solver;
	if (lats.defovolSolver) {
		for (int i = 0; i < lats.numDv; ++i) {
			delete lats.defovol[i];
		}
		delete lats.defovolTmp;
		delete lats.defovolSolver;
	}
	readGrid4dUniCleanup(&lats.fileHandle1);
	delete LATS[ID];
	LATS[ID] = NULL;
}
PYTHON() void loadAdvectTimeSlice_OptRun(
    int ID, std::string fname, Grid<Real> &dst, Grid4d<Real> &phi, Real time,
    Real blendAlpha, Real loadTimeScale, Vec4 defoOffset, Vec4 defoScale,
    Vec4 defoFactor, Vec4 overrideSize = Vec4(-1.),
    Real overrideTimeOff = 0., // for shifted / partial loading
    Grid<Vec3> *debugVel = NULL, Grid<Real> *debugVelT = NULL,
    bool zeroVel = false, Real thirdAlpha = 0., int bordSkip = 1,
    Real fourthAlpha = 0., float defoAniFac = 1.)
{
	assertMsg(LATS[ID], "Load-advect data id " << ID << " not initialized!");
	LoadAdvectData &lats = *LATS[ID];
	Vec3i dim(lats.dim);
	int dimT = lats.dimT;
	Grid4d<Vec4> &v1 = (*lats.v1), &v2 = (*lats.v2), &vt = (*lats.tmp);

	// same code as for regular loadAdvectTimeSlice! up to lastT check
	Vec4 dim4(dim[0], dim[1], dim[2], dimT);
	Vec4 sourceFactor(0.), off2 = defoOffset;
	Vec4i defoSize = phi.getSize();
	if (overrideSize.x > 0.) {
		defoSize = Vec4i(overrideSize.x, overrideSize.y, overrideSize.z,
		                 overrideSize.t);
	}
	gridFactor4d(toVec4(dim4), toVec4(defoSize), Vec4(-1.), defoScale,
	             sourceFactor, off2);
	Vec3 sourceFactor3 = Vec3(sourceFactor.x, sourceFactor.y, sourceFactor.z);
	Vec3 offset3 = Vec3(off2.x, off2.y, off2.z);

	// load 2 slices (time is given in dst units, convert to src)
	Real srcTime = (time)*sourceFactor.t * loadTimeScale + off2.t - 0.5; 
	int t = (int)srcTime;
	int tp1 = t + 1;
	Real tw = srcTime - (Real)t;
	t = std::min(t, dimT - 1); // make sure we dont go too far...
	tp1 = std::min(tp1, dimT - 1);

	if (!lats.useDefoVols) {
		if (lats.lastT == t) {
			// all good
		} else if (lats.lastT + 1 == t) {
			debMsg("Reading single velocity at "
			           << t << " w " << tw << "    sf " << sourceFactor << " "
			           << off2,
			       1);
			v1.copyFrom(v2);
			readGrid4dUni<Vec4>(lats.fname1, NULL, tp1, lats.v2,
			                    &lats.fileHandle1);
		} else {
			debMsg("Reading double velocity at "
			           << t << " w " << tw << "    sf " << sourceFactor << " "
			           << off2,
			       1);
			readGrid4dUni<Vec4>(lats.fname1, NULL, t, lats.v1,
			                    &lats.fileHandle1);
			readGrid4dUni<Vec4>(lats.fname1, NULL, tp1, lats.v2,
			                    &lats.fileHandle1);
		}

		// linearly blend in time , use 3rd temp array here, dont overwrite v1
		// (possibly remove for on-the-fly interpol)
		vt.setConst(Vec4(0.));
		vt.addScaled(v1, Vec4(1. - tw));
		vt.addScaled(v2, Vec4(tw));

	} else {
		debMsg("Updating defo vol at " << t << " w " << tw << "    sf "
		                               << sourceFactor << " " << off2,
		       1);
		lats.updateDefoVol(t);
		if (lats.numDv == 2) {
			// get velocity from volume, accumulate
			if (!lats.doAligned) {

				// aligned , w/o backmove
				debMsg("Adding aligned defovol 0 alpha " << blendAlpha << " , "
				                                         << thirdAlpha,
				       1);
				Grid4d<Vec4> &dvol1 = (*lats.defovol[0]);
				Grid4d<Vec4> &dvol2 = (*lats.defovol[1]);
				FOR_IJK_BND(vt, 0)
				{
					Vec4 p(i, j, k, srcTime - lats.defovolOff);
					Vec4 v2 = dvol2.getInterpolated(p + Vec4(0.5));
					vt(i, j, k, 0) =
					    blendAlpha * dvol1.getInterpolated(p + Vec4(0.5) - v2) +
					    thirdAlpha * v2;
				}

			} else {

				// align & move "back"
				debMsg("Adding aligned2 defovols, alphas "
				           << blendAlpha << " , " << thirdAlpha,
				       1);
				Grid4d<Vec4> &dvol1 = (*lats.defovol[0]);
				Grid4d<Vec4> &dvol2 = (*lats.defovol[1]);
				Grid4d<Vec4> &dvt = (*lats.defovolTmp);
				FOR_IJKT_BND(dvt, 0)
				{
					Vec4 p(i, j, k, t); 
					Vec4 v2 = dvol2.getInterpolated(p + Vec4(0.5));
					dvt(i, j, k, t) =
					    blendAlpha * dvol1.getInterpolated(p + Vec4(0.5) - v2) +
					    thirdAlpha * v2;
				}
				FOR_IJK_BND(vt, 0)
				{
					Vec4 p(i, j, k, srcTime - lats.defovolOff);
					Vec4 v1 =
					    dvol1.getInterpolated(p + Vec4(0.5)) * -1; // invert!
					vt(i, j, k, 0) = dvt.getInterpolated(
					    p + Vec4(0.5) - (1. - blendAlpha) * v1);
				}
			}
		} else if (lats.numDv == 3) {
			// 3 defos, aligned , w/o backmove
			debMsg("Adding 3 aligned defovols, alphas " << blendAlpha << " , "
			                                            << thirdAlpha << " , "
			                                            << fourthAlpha,
			       1);
			Grid4d<Vec4> &dvol1 = (*lats.defovol[0]);
			Grid4d<Vec4> &dvol2 = (*lats.defovol[1]);
			Grid4d<Vec4> &dvol3 = (*lats.defovol[2]);
			FOR_IJK_BND(vt, 0)
			{
				Vec4 p(i, j, k, srcTime - lats.defovolOff);
				Vec4 v3 = dvol3.getInterpolated(p + Vec4(0.5));
				Vec4 v2 = dvol2.getInterpolated(p + Vec4(0.5) - v3);
				Vec4 v1 = dvol1.getInterpolated(p + Vec4(0.5) - v3 - v2);
				vt(i, j, k, 0) =
				    blendAlpha * v1 + thirdAlpha * v2 + fourthAlpha * v3;
			}
		} else {
			errMsg("Code currently only supports 2 deformation volumes");
		}

		// reset alpha value, already accumulated in modified vt values
		blendAlpha = 1.;
	}

	if (zeroVel) {
		debMsg("Debug - zeroing deformation!", 1);
		vt.setConst(Vec4(0.));
	}

	lats.lastT = t;
	MuTime start;

	debMsg("Performing comb-opt SL lookup, alpha=" << blendAlpha << " ", 1);
	// assertMsg( bordSkip>10, "Warning - dont use for small sizes...");
	if (bordSkip < 10) debMsg("Warning - dont use for small sizes...", 1);
	knSemiLagrangeLookupSlice4d_WithVel(
		dst, phi, time + overrideTimeOff, blendAlpha, vt, sourceFactor3,
		offset3, defoFactor * defoAniFac, dim, bordSkip);
}

// ================================================================================
// instantiate 4d templates

PYTHON() void opticalFlow4d(Grid4d<Vec4> &vel, Grid4d<Real> &i0,
                            Grid4d<Real> &i1, Grid4d<Real> *rhsT = NULL,
                            Real wSmooth = 0., Real wEnergy = 0.,
                            Real postVelBlur = 0., Real cgAccuracy = 1e-04,
                            int blurType = 1, Real resetBndWidth = -1.)
{
	opticalFlowTemplate<Grid4d<Real>, Grid4d<Vec4> >(
	    vel, i0, i1, rhsT, wSmooth, wEnergy, postVelBlur, cgAccuracy, blurType,
	    resetBndWidth);
}

PYTHON() void corrVelsOf4d(Grid4d<Vec4> &dst, Grid4d<Vec4> &vel,
                           Grid4d<Real> &phiOrg, Grid4d<Real> &phiCurr,
                           Grid4d<Real> &phiTarget, Real threshPhi = 1e10,
                           Real threshNorm = 1e10, Real postVelBlur = 0.,
                           Real resetBndWidth = -1., int maxIter = 100)
{
	corrVelsOfTempl<Grid4d<Real>, Grid4d<Vec4> >(
	    dst, vel, phiOrg, phiCurr, phiTarget, threshPhi, threshNorm,
	    postVelBlur, resetBndWidth, maxIter);
}

PYTHON() Real
    calcLsDiff4d(Grid4d<Real> &i0, Grid4d<Real> &i1, Grid4d<Real> *out = NULL,
                 Real correction = 1., int bnd = 0)
{
	return calcLsDiffTempl<Grid4d<Real> >(i0, i1, out, correction, bnd);
}

// ugly , just for testing , copied from calcLsDiffTempl
template <class GRID>
Real calcSmokeDiffTempl(GRID &i0, GRID &i1, GRID *out = NULL,
                        Real correction = 1., int bnd = 0)
{
	double accu = 0.;
	FOR_IJKT_BND(i0, bnd)
	{
		Real d = fabs(i0(i, j, k, t) - i1(i, j, k, t)) * correction;
		accu += d;
	}
	int sx = i0.getSizeX();
	sx -= 2 * bnd;
	int sy = i0.getSizeY();
	sy -= 2 * bnd;
	int sz = i0.getSizeZ();
	if (i0.is3D()) sz -= 2 * bnd;
	int st = i0.getSizeT();
	if (i0.is4D()) st -= 2 * bnd;
	accu *= 1000.; // debug, factor for more readable numbers...
	if (i0.getSizeT() > 1) accu *= 1000.; // even more for 4d
	accu *= 1. / double(sx * sy * sz * st);
	return accu;
}
PYTHON() Real
    calcSmokeDiff4d(Grid4d<Real> &i0, Grid4d<Real> &i1,
                    Grid4d<Real> *out = NULL, Real correction = 1., int bnd = 0)
{
	return calcSmokeDiffTempl<Grid4d<Real> >(i0, i1, out, correction, bnd);
}

template <>
void advectCflTemplate<Grid4d<int>, Grid4d<Vec4>, Grid4dBase>(
    Real cfl, Grid4d<int> *flags, Grid4d<Vec4> &vel, Grid4dBase *grid,
    int order, Real velFactor, int orderSpace)
{
	Grid4d<Vec4> velTmp(vel.getParent());
	velTmp.copyFrom(vel);
	velTmp.multConst(Vec4(velFactor));
	// note - flags and order are not used here!
	advectCflHelper4d(cfl, velTmp, grid);
}

PYTHON() void opticalFlowMultiscale4d(
    Grid4d<Vec4> &vel, Grid4d<Real> &i0, Grid4d<Real> &i1,
    Grid4d<Real> *rhsT = NULL, // debug info
    Real wSmooth = 0., Real wEnergy = 0., int level = 0, Real postVelBlur = 0.,
    Real cgAccuracy = 1e-04, int blurType = 1, Real cfl = 999,
    int orderTime = 1, int orderSpace = 1, Real resetBndWidth = -1.,
    int multiStep = 1, int projSizeThresh = 9999, int minGridSize = 10,
    bool doFinalProject = false )
{
	opticalFlowMultiscaleTemplate<Grid4d<int>, Grid4d<Real>, Grid4d<Vec4> >(
	    vel, i0, i1, rhsT, wSmooth, wEnergy, level, postVelBlur, cgAccuracy,
	    blurType, cfl, orderTime, orderSpace, resetBndWidth, multiStep,
	    projSizeThresh, minGridSize, doFinalProject);
}

} // namespace
