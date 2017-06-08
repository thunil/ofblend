/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL)
 * http://www.gnu.org/licenses
 *
 * Use this file to test new functionality
 *
 ******************************************************************************/

#include "levelset.h"
#include "commonkernels.h"
#include "particle.h"
#include "grid4d.h"
#include <cmath>

using namespace std;

namespace Manta
{

// two simple example kernels

KERNEL(idx, reduce = +)
    returns(double sum = 0) double reductionTest(const Grid<Real> &v)
{
	sum += v[idx];
}

KERNEL(idx, reduce = min)
    returns(double sum = 0) double minReduction(const Grid<Real> &v)
{
	if (sum < v[idx]) sum = v[idx];
}

// ... add own test code ...

//! rescale range of inputs for optical flow
template <class GridT> KERNEL(fourd) void knRemapLsValues(GridT &phi, Real max)
{
	Real &a = phi(i, j, k, t);
	a *= -1;
	if (a > max) a = max;
	if (a < -max) a = -max;
	a /= max;
}

//! remap values of an SDF [0..max] to the [0..1] range, and invert (inside
//positive!)
PYTHON() void remapLsValues(Grid<Real> &grid, Real max)
{
	knRemapLsValues<Grid<Real> >(grid, max);
}

//! remap values of 4d grid
PYTHON() void remapLsValues4d(Grid4d<Real> &grid, Real max)
{
	knRemapLsValues<Grid4d<Real> >(grid, max);
}

//! helper for tests
PYTHON() void debugMaskVel(Grid4d<Real> &phi, Grid4d<Vec4> &vel)
{
	FOR_IJKT_BND(phi, 0)
	{
		if (phi(i, j, k, t) > 0.) continue;
		vel(i, j, k, t) = Vec4(0.);
	}
}
PYTHON() void debugSetVel(Grid4d<Real> &phi, Grid4d<Vec4> &vel)
{
	FOR_IJKT_BND(phi, 0)
	{
		// if(phi(i,j,k,t) > 0.) continue;
		vel(i, j, k, t) = Vec4(phi(i, j, k, t));
	}
}
PYTHON() void debugZeroVel(Grid4d<Vec4> &vel, int comp)
{
	FOR_IJKT_BND(vel, 0) { vel(i, j, k, t)[comp] = 0.; }
}
PYTHON() void scalePartPos(BasicParticleSystem &parts, Real scale)
{
	for (int i = 0; i < parts.size(); ++i) {
		parts[i].pos *= scale;
	}
}

// simple blur , for time interpol

KERNEL(bnd = 1) template <class S>
void knSimpleBlurSpecial(Grid<S> &a, Grid<S> &tmp, Real thresh, int bord)
{
	if (!a.isInBounds(Vec3i(i, j, k), bord)) {
		tmp(i, j, k) = a(i, j, k);
		return;
	}
	if (a(i, j, k) < thresh) {
		tmp(i, j, k) = a(i, j, k);
		return;
	}

	tmp(i, j, k) =
	    a(i + 1, j, k) + a(i - 1, j, k) + a(i, j + 1, k) + a(i, j - 1, k);
	if (a.is3D()) {
		tmp(i, j, k) += a(i, j, k + 1) + a(i, j, k - 1);
		tmp(i, j, k) *= 1. / 6.;
	} else {
		tmp(i, j, k) *= 1. / 4.;
	}
}
// smoothen only positive values, for SDF interpol (2d/3d)
template <class S>
void templ_simpleBlurSpecial(Grid<S> &a, int iter = 1, Real thresh = 0.,
                             int bord = 0)
{
	// Real strength = 1.;
	Grid<S> tmp(a.getParent());
	for (int numIt = 0; numIt < iter; numIt++) {
		knSimpleBlurSpecial<S>(a, tmp, thresh, bord);
		a.swap(tmp);
	}
}
PYTHON() void simpleBlurSpecial(Grid<Real> &a, int iter = 1, Real thresh = 0.,
                                int bord = 0)
{
	templ_simpleBlurSpecial<Real>(a, iter, thresh, bord);
}

PYTHON() void deleteTopParts(BasicParticleSystem &parts, Real maxHeight)
{
	for (IndexInt idx = 0; idx < (IndexInt)parts.size(); idx++) {
		if (!parts.isActive(idx)) continue;

		if (parts.getPos(idx)[1] > maxHeight) {
			parts.kill(idx); // out of domain, remove
		}
	}
}

PYTHON() void limitPvel(BasicParticleSystem &parts,
                        ParticleDataImpl<Vec3> &pvel, Real max = 1.)
{
	for (IndexInt idx = 0; idx < (IndexInt)parts.size(); idx++) {
		if (!parts.isActive(idx)) continue;
		Vec3 v = pvel[idx];
		Real vl = norm(v);
		if (vl > max) {
			v *= max / vl;
			pvel[idx] = v;
		}
	}
}

//! kernel to set velocity components of mac grid to value for a boundary of w
//cells
KERNEL() void knSetBoundaryMAC(Grid<Vec3> &grid, Vec3 value, int w)
{
	if (i <= w || i >= grid.getSizeX() - w || j <= w - 1 ||
	    j >= grid.getSizeY() - 1 - w ||
	    (grid.is3D() && (k <= w - 1 || k >= grid.getSizeZ() - 1 - w)))
		grid(i, j, k).x = value.x;
	if (i <= w - 1 || i >= grid.getSizeX() - 1 - w || j <= w ||
	    j >= grid.getSizeY() - w ||
	    (grid.is3D() && (k <= w - 1 || k >= grid.getSizeZ() - 1 - w)))
		grid(i, j, k).y = value.y;
	if (i <= w - 1 || i >= grid.getSizeX() - 1 - w || j <= w - 1 ||
	    j >= grid.getSizeY() - 1 - w ||
	    (grid.is3D() && (k <= w || k >= grid.getSizeZ() - w)))
		grid(i, j, k).z = value.z;
}

//! only set normal velocity components of mac grid to value for a boundary of w
//cells
KERNEL() void knSetBoundaryMACNorm(Grid<Vec3> &grid, Vec3 value, int w)
{
	if (i <= w || i >= grid.getSizeX() - w) grid(i, j, k).x = value.x;
	if (j <= w || j >= grid.getSizeY() - w) grid(i, j, k).y = value.y;
	if ((grid.is3D() && (k <= w || k >= grid.getSizeZ() - w)))
		grid(i, j, k).z = value.z;
}

//! set velocity components of mac grid to value for a boundary of w cells
//(optionally only normal values)
PYTHON() void setBoundMAC(MACGrid &vel, Vec3 value, int boundaryWidth,
                          bool normalOnly = false)
{
	if (!normalOnly)
		knSetBoundaryMAC(vel, value, boundaryWidth);
	else
		knSetBoundaryMACNorm(vel, value, boundaryWidth);
}

// init

PYTHON() Real debugGridAvg4d(Grid4d<Real> &phi, int brd = 0)
{
	double cnt = 0., accu = 0.;
	FOR_IJKT_BND(phi, brd)
	{
		accu += phi(i, j, k, t);
		cnt += 1.;
	}
	return (Real)(accu * 1000000. / cnt);
}

PYTHON() Real debugVelAvg4d(Grid4d<Vec4> &v, int brd = 0)
{
	double cnt = 0., accu = 0.;
	FOR_IJKT_BND(v, brd)
	{
		accu += norm(v(i, j, k, t));
		cnt += 1.;
	}
	return (Real)(accu * 1. / cnt);
}

PYTHON() void initVecFromScalar(Grid4d<Real> &source,
                                Grid4d<Vec4> &target) // helper
{
	FOR_IJKT_BND(source, 0)
	{
		Real v = source(i, j, k, t);
		target(i, j, k, t) = Vec4(v, v, v, v);
	}
}

PYTHON() void initTestCheckerboard(Grid4d<Real> &val, Grid4d<Vec4> *vec = NULL,
                                   int brd = 0)
{
	// init 4 patches per direction
	Real num = 4.;
	FOR_IJKT_BND(val, brd)
	{
		int ci = int(i / (val.getSizeX() / num));
		int cj = int(j / (val.getSizeY() / num));
		int ck = int(k / (val.getSizeZ() / num));
		int ct = int(t / (val.getSizeT() / num));
		// accu += norm( v(i, j, k, t) ); cnt += 1.;
		Real v = -1.;
		if ((ci + cj + ck + ct) % 2 == 1) v = 1.;

		val(i, j, k, t) = (v);
		if (vec) (*vec)(i, j, k, t) = Vec4(v);
	}
}

} // namespace
