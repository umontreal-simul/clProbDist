/* This file is part of clProbDist.
 *
 * Copyright 2015-2016  Pierre L'Ecuyer, Universite de Montreal and Advanced Micro Devices, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Authors:
 *
 *   Nabil Kemerchou <kemerchn@iro.umontreal.ca>    (2015)
 *   David Munger <mungerd@iro.umontreal.ca>        (2015)
 *   Pierre L'Ecuyer <lecuyer@iro.umontreal.ca>     (2015)
 *
 */

#pragma once
#ifndef CLPROBDIST_PRIVATE_GAMMADIST_CH
#define CLPROBDIST_PRIVATE_GAMMADIST_CH

#ifndef _CLPROBDIST_GAMMA_OBJ_MEM
#define _CLPROBDIST_GAMMA_OBJ_MEM
#endif

#define signum(a)    (((a) < 0) ? -1 : (((a) > 0) ? 1 : 0))

typedef struct _clprobdistGammaParams{
	cl_double alpha;
	cl_double lambda;
	cl_double logFactor;
} clprobdistGammaParams;

struct _clprobdistGamma {
	clprobdistContinuous contDistObj;
	clprobdistGammaParams params;
};

constant const cl_double clprobdistGamma_ALIM = 1.0E5;

constant const cl_double clprobdistGammaXBIG = 100.0;
constant const cl_double clprobdistGammaXBIGM = 1000.0;

// clprobdistGammaEPSARRAY[j]: Epsilon required for j decimal degits of precision
constant const cl_double clprobdistGammaEPSARRAY[] = {
	0.5, 0.5E-1, 0.5E-2, 0.5E-3, 0.5E-4, 0.5E-5, 0.5E-6, 0.5E-7, 0.5E-8,
	0.5E-9, 0.5E-10, 0.5E-11, 0.5E-12, 0.5E-13, 0.5E-14, 0.5E-15, 0.5E-16,
	0.5E-17, 0.5E-18, 0.5E-19, 0.5E-20, 0.5E-21, 0.5E-22, 0.5E-23, 0.5E-24,
	0.5E-25, 0.5E-26, 0.5E-27, 0.5E-28, 0.5E-29, 0.5E-30, 0.5E-31, 0.5E-32,
	0.5E-33, 0.5E-34, 0.5E-35
};

//Private functions
/*
* This is the function  1 + (1 - x*x + 2*x*log(x)) / ((1 - x)*(1 - x))
*/
static cl_double clprobdistGamma_mybelog(cl_double x)
{
	if (x < 1.0e-30)
		return 0.0;
	if (x > 1.0e30)
		return 2.0*(log(x) - 1.0) / x;
	if (x == 1.0)
		return 1.0;

	cl_double t = 1.0 - x;
	if (x < 0.9 || x > 1.1) {
		cl_double w = (t + x*log(x)) / (t*t);
		return 2.0 * w;
	}

	// For x near 1, use a series expansion to avoid loss of precision.
	cl_double term;
	cl_double EPS = 1.0e-12;
	cl_double tpow = 1.0;
	cl_double sum = 0.5;
	int j = 3;
	do {
		tpow *= t;
		term = tpow / (j * (j - 1));
		sum += term;
		j++;
	} while (fabs(term / sum) > EPS);
	return 2.0*sum;
}

//************************************
// Root Finder Function argument type
//************************************
typedef struct _clprobdistRootFinderFunc {
	int d;
	cl_double alp, u;
} clprobdistRootFinderFunc;
cl_double clprobdistGammaCDF_1(cl_double alpha, int d, cl_double x, clprobdistStatus* err);
cl_double clprobdistGammaCDF_2(cl_double alpha, int d, cl_double x, clprobdistStatus* err);
cl_double clprobdistRootFinderFuncEvaluate(cl_double x, clprobdistRootFinderFunc* f, clprobdistStatus * err) {

	return f->u - clprobdistGammaCDF_1(f->alp, f->d, x, err);
}

//*******************************
//Rootfinder
//*******************************
constant const cl_double clprobdistGamma_MINVAL = DBL_MIN; // 5.0e-308;
constant const cl_double clprobdistRootFinderDBL_EPSILON = 2.220446049250313E-16;

/**
* Computes a root <SPAN CLASS="MATH"><I>x</I></SPAN> of the function in <TT>f</TT> using the
*     Brent-Dekker method. The interval <SPAN CLASS="MATH">[<I>a</I>, <I>b</I>]</SPAN> must contain the root <SPAN CLASS="MATH"><I>x</I></SPAN>.
*     The calculations are done with an approximate relative precision
*     <TT>tol</TT>.  Returns <SPAN CLASS="MATH"><I>x</I></SPAN> such that <SPAN CLASS="MATH"><I>f</I> (<I>x</I>) = 0</SPAN>.
*
* @param a left endpoint of initial interval
*
*    @param b right endpoint of initial interval
*
*    @param f the function which is evaluated
*
*    @param tol accuracy goal
*
*    @return the root <SPAN CLASS="MATH"><I>x</I></SPAN>
*
*/
cl_double clprobdistRootFinderBrentDekker(cl_double a, cl_double b, clprobdistRootFinderFunc* f, cl_double tol, clprobdistStatus* err) {
	cl_double EPS = 0.5E-15;
	int MAXITER = 120;    // Maximum number of iterations
	cl_double c, d, e;
	cl_double fa, fb, fc;
	//bool DEBUG = false;
	if (err) *err = CLPROBDIST_SUCCESS;

	// Special case I = [b, a]
	if (b < a) {
		cl_double ctemp = a;
		a = b;
		b = ctemp;
	}

	// Initialization
	fa = clprobdistRootFinderFuncEvaluate(a, f, err);
	fb = clprobdistRootFinderFuncEvaluate(b, f, err);
	c = a;
	fc = fa;
	d = e = b - a;
	tol += EPS + clprobdistRootFinderDBL_EPSILON; // in case tol is too small

	if (fabs(fc) < fabs(fb)) {
		a = b;
		b = c;
		c = a;
		fa = fb;
		fb = fc;
		fc = fa;
	}

	int i;
	for (i = 0; i < MAXITER; i++) {
		cl_double s, p, q, r;
		cl_double tol1 = tol + 4.0 * clprobdistRootFinderDBL_EPSILON * fabs(b);
		cl_double xm = 0.5 * (c - b);
		/*if (DEBUG) {
			cl_double error = fabs(fa - fb);
			*err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): [a, b] = [%g, %g]   fa = %g,   fb = %g   |fa - fb| = %.2g%n", __func__, a, b, fa, fb, error);
		}*/

		if (fabs(fb) <= clprobdistGamma_MINVAL) {
			return b;
		}
		if (fabs(xm) <= tol1) {
			if (fabs(b) > clprobdistGamma_MINVAL)
				return b;
			else
				return 0;
		}

		if ((fabs(e) >= tol1) && (fabs(fa) > fabs(fb))) {
			if (a != c) {
				// Inverse quadratic interpolation
				q = fa / fc;
				r = fb / fc;
				s = fb / fa;
				p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
				q = (q - 1.0) * (r - 1.0) * (s - 1.0);
			}
			else {
				// Linear interpolation
				s = fb / fa;
				p = 2.0 * xm * s;
				q = 1.0 - s;
			}

			// Adjust signs
			if (p > 0.0)
				q = -q;
			p = fabs(p);

			// Is interpolation acceptable ?
			if (((2.0 * p) >= (3.0 * xm * q - fabs(tol1 * q)))
				|| (p >= fabs(0.5 * e * q))) {
				d = xm;
				e = d;
			}
			else {
				e = d;
				d = p / q;
			}
		}
		else {
			// Bisection necessary
			d = xm;
			e = d;
		}

		a = b;
		fa = fb;
		if (fabs(d) > tol1)
			b += d;
		else if (xm < 0.0)
			b -= tol1;
		else
			b += tol1;
		fb = clprobdistRootFinderFuncEvaluate(b,f, err);
		if (fb * (signum(fc)) > 0.0) {
			c = a;
			fc = fa;
			d = e = b - a;
		}
		else {
			a = b;
			b = c;
			c = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}
	}

	if (i >= MAXITER)
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): WARNING: root finding does not converge", __func__);
	return b;
}

/**
* Computes a root <SPAN CLASS="MATH"><I>x</I></SPAN> of the function in <TT>f</TT> using the
*     <SPAN  CLASS="textit">bisection</SPAN> method. The interval <SPAN CLASS="MATH">[<I>a</I>, <I>b</I>]</SPAN> must contain the root <SPAN CLASS="MATH"><I>x</I></SPAN>.
*     The calculations are done with an approximate relative precision
*     <TT>tol</TT>.  Returns <SPAN CLASS="MATH"><I>x</I></SPAN> such that <SPAN CLASS="MATH"><I>f</I> (<I>x</I>) = 0</SPAN>.
*
* @param a left endpoint of initial interval
*
*    @param b right endpoint of initial interval
*
*    @param f the function which is evaluated
*
*    @param tol accuracy goal
*
*    @return the root <SPAN CLASS="MATH"><I>x</I></SPAN>
*
*/
cl_double clprobdistRootFinderBisection(cl_double a, cl_double b, clprobdistRootFinderFunc* f, cl_double tol, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;
	// Case I = [b, a]
	if (b < a) {
		cl_double ctemp = a;
		a = b;
		b = ctemp;
	}
	cl_double xa = a;
	cl_double xb = b;
	cl_double ya = clprobdistRootFinderFuncEvaluate(a,f, err);
	cl_double x = 0, y = 0;
	int MAXITER = 1200;   // Maximum number of iterations
	//bool DEBUG = false;
	tol += clprobdistRootFinderDBL_EPSILON; // in case tol is too small

	//if (DEBUG) printf("\niter              xa                   xb              f(x)\n");

	cl_bool fini = CL_FALSE;
	int i = 0;
	while (!fini) {
		x = (xa + xb) / 2.0;
		y = clprobdistRootFinderFuncEvaluate(x,f, err);
		if ((fabs(y) <= clprobdistGamma_MINVAL) ||
			(fabs(xb - xa) <= tol * fabs(x)) ||
			(fabs(xb - xa) <= clprobdistGamma_MINVAL)) {
			if (fabs(x) > clprobdistGamma_MINVAL)
				return x;
			else
				return 0;
		}
		if (y * ya < 0.0)
			xb = x;
		else
			xa = x;
		++i;
		//if (DEBUG)	printf("%3d    %18.12g     %18.12g    %14.4g%n",i, xa, xb, y);
		if (i > MAXITER) {
			if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): ***** bisection:  SEARCH DOES NOT CONVERGE", __func__);
			fini = CL_TRUE;
		}
	}
	return x;
}


//Static functions

cl_double clprobdistGammaDensity(cl_double alpha, cl_double lambda, int decprec, cl_double x, clprobdistStatus* err)
{
	if (err) *err = CLPROBDIST_SUCCESS;
	if (alpha <= 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): alpha <= 0", __func__);
		return -1;
	}

	if (lambda <= 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): lambda <= 0", __func__);
		return -1;
	}


	if (x <= 0)
		return 0.0;
	cl_double z = alpha * log(lambda*x) - lambda*x - lgamma(alpha);
	if (z > -clprobdistGammaXBIGM)
		return exp(z) / x;
	else
		return 0.0;
}

cl_double clprobdistGammaStdDeviation(cl_double alpha, cl_double lambda, int decprec, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;
	if (alpha <= 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): alpha <= 0", __func__);
		return -1;
	}

	if (lambda <= 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): lambda <= 0", __func__);
		return -1;
	}

	return (sqrt(alpha) / lambda);
}
cl_double clprobdistGammaInverseCDF_1(cl_double alpha, int d, cl_double u, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;
	if (alpha <= 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): alpha <= 0", __func__);
		return -1;
	}

	if (u > 1.0 || u < 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s():  u not in [0,1]", __func__);
		return -1;
	}

	if (u <= 0.0)
		return 0;

	if (u >= 1.0)
		return DBL_MAX;

	if (d <= 0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): d <= 0", __func__);
		return -1;
	}

	if (d > 15)
		d = 15;

	cl_double EPS = pow(10.0, -d);
	cl_double sigma = clprobdistGammaStdDeviation(alpha, 1.0, 0, err);
	cl_double x = clprobdistNormalInverseCDF(alpha, sigma, u, err);

	if (x < 0.)
		x = 0.;
	cl_double cdf = clprobdistGammaCDF_1(alpha, d, x, err);

	cl_double xmax;
	if (alpha < 1.0)
		xmax = 100.0;
	else
		xmax = alpha + 40.0 * sigma;

	clprobdistRootFinderFunc f;
	f.d = d;
	f.alp = alpha;
	f.u = u;
	//clprobdistRootFinderFunc f = new clprobdistRootFinderFunc(alpha, d, u);

	if (u <= 1.0e-8 || alpha <= 1.5) {
		if (cdf < u)
			return clprobdistRootFinderBisection(x, xmax, &f, EPS, err);
		else
			return clprobdistRootFinderBisection(0, x, &f, EPS, err);
	}
	else {
		if (cdf < u)
			return clprobdistRootFinderBrentDekker(x, xmax, &f, EPS, err);
		else
			return clprobdistRootFinderBrentDekker(0, x, &f, EPS, err);
	}
}

cl_double clprobdistGammaMean(cl_double alpha, cl_double lambda, int decprec, clprobdistStatus* err) {

	if (err) *err = CLPROBDIST_SUCCESS;
	if (alpha <= 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): alpha <= 0", __func__);
		return -1;
	}

	if (lambda <= 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): lambda <= 0", __func__);
		return -1;
	}

	return alpha / lambda;

}
cl_double clprobdistGammaVariance(cl_double alpha, cl_double lambda, int decprec, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;
	if (alpha <= 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): alpha <= 0", __func__);
		return -1;
	}

	if (lambda <= 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): lambda <= 0", __func__);
		return -1;
	}

	return alpha / (lambda * lambda);
}

cl_double clprobdistGammaComplCDF_1(cl_double alpha, int d, cl_double x, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;

	if (alpha <= 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): alpha <= 0", __func__);
		return -1;
	}

	if (d <= 0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): d <= 0", __func__);
		return -1;
	}

	if (x <= 0.0)
		return 1.0;
	if (1.0 == alpha)
		return clprobdistExponentialComplCDF(1.0, x, err);

	if (alpha >= 70.0) {
		if (x >= alpha * clprobdistGammaXBIG)
			return 0.0;
	}
	else {
		if (x >= clprobdistGammaXBIGM)
			return 0.0;
	}

	if (alpha >= clprobdistGamma_ALIM) {
		cl_double d2 = x + 1.0 / 3.0 - alpha - 0.02 / alpha;
		cl_double S = alpha - 1.0 / 2.0;
		cl_double w = clprobdistGamma_mybelog(S / x);
		cl_double y = d2 * sqrt(w / x);
		return clprobdistStdNormalComplCDF(y, err);
	}

	if (x <= 1.0 || x < alpha)
		return 1.0 - clprobdistGammaCDF_2(alpha, d, x, err);


	cl_double V[6] = { 0 };
	cl_double EPS = clprobdistGammaEPSARRAY[d];
	cl_double RENORM = 1.0E100;
	cl_double R, dif;
	int i;
	cl_double factor = exp(alpha*log(x) - x - lgamma(alpha));

	cl_double A = 1.0 - alpha;
	cl_double B = A + x + 1.0;
	cl_double term = 0.0;
	V[0] = 1.0;
	V[1] = x;
	V[2] = x + 1.0;
	V[3] = x * B;
	cl_double res = V[2] / V[3];

	do {
		A += 1.0;
		B += 2.0;
		term += 1.0;
		V[4] = B * V[2] - A * term * V[0];
		V[5] = B * V[3] - A * term * V[1];
		if (V[5] != 0.0) {
			R = V[4] / V[5];
			dif = fabs(res - R);
			if (dif <= EPS*R){
				return factor*res;
			}

			res = R;
		}
		for (i = 0; i < 4; i++)
			V[i] = V[i + 2];
		if (fabs(V[4]) >= RENORM) {
			for (i = 0; i < 4; i++)
				V[i] /= RENORM;
		}
	} while (CL_TRUE);
}
cl_double clprobdistGammaComplCDF(cl_double alpha, cl_double lambda, int d, cl_double x, clprobdistStatus* err) {
	return clprobdistGammaComplCDF_1(alpha, d, lambda*x, err);
}

cl_double clprobdistGammaCDF_1(cl_double alpha, int d, cl_double x, clprobdistStatus* err){
	cl_double ret = clprobdistGammaCDF_2(alpha, d, x, err);
	if (ret < 0.0)
		return 1.0 - clprobdistGammaComplCDF_1(alpha, d, x, err);
	return ret;
}

cl_double clprobdistGammaCDF_2(cl_double alpha, int d, cl_double x, clprobdistStatus* err){
	if (err) *err = CLPROBDIST_SUCCESS;

	if (x <= 0.0)
		return 0.0;

	if (1.0 == alpha)
		return clprobdistExponentialCDF(1.0, x, err);

	if (alpha > 10.0) {
		if (x > alpha * 10.0)
			return 1.0;
	}
	else {
		if (x > clprobdistGammaXBIG)
			return 1.0;
	}

	if (alpha >= clprobdistGamma_ALIM) {
		cl_double d2 = x + 1.0 / 3.0 - alpha - 0.02 / alpha;
		cl_double S = alpha - 1.0 / 2.0;
		cl_double w = clprobdistGamma_mybelog(S / x);
		cl_double y = d2 * sqrt(w / x);
		return clprobdistStdNormalCDF(y, err);
	}

	if (x <= 1.0 || x < alpha) {
		cl_double factor, z, rn, term;
		factor = exp(alpha*log(x) - x - lgamma(alpha));
		cl_double EPS = clprobdistGammaEPSARRAY[d];
		z = 1.0;
		term = 1.0;
		rn = alpha;
		do {
			rn += 1.0;
			term *= x / rn;
			z += term;
		} while (term >= EPS * z);
		return z*factor / alpha;
	}

	// this is a workaround to avoid birecusivity on GPU's
	// by returning -1.0, we let the caller know it must call ComplCDF()
	return -1.0;
}

// Dynamic functions
cl_double clprobdistGammaCDF(cl_double alpha, cl_double lambda, int d, cl_double x, clprobdistStatus* err){
	return clprobdistGammaCDF_1(alpha, d, lambda*x, err);
}

cl_double clprobdistGammaCDFWithObject(_CLPROBDIST_GAMMA_OBJ_MEM const clprobdistGamma* distObj, cl_double x, clprobdistStatus* err)
{
	return clprobdistGammaCDF(distObj->params.alpha,distObj->params.lambda,distObj->contDistObj.decPrec,x,err);

}

cl_double clprobdistGammaInverseCDF(cl_double alpha, cl_double lambda, int d, cl_double u, clprobdistStatus* err)
{

	return clprobdistGammaInverseCDF_1(alpha, d, u, err) / lambda;
}
cl_double clprobdistGammaInverseCDFWithObject(_CLPROBDIST_GAMMA_OBJ_MEM const clprobdistGamma* distObj, cl_double u, clprobdistStatus* err)
{
	return clprobdistGammaInverseCDF_1(distObj->params.alpha, distObj->contDistObj.decPrec, u, err) / distObj->params.lambda;
}

cl_double clprobdistGammaComplCDFWithObject(_CLPROBDIST_GAMMA_OBJ_MEM const clprobdistGamma* distObj, cl_double x, clprobdistStatus* err) {
	return clprobdistGammaComplCDF(distObj->params.alpha, distObj->params.lambda, distObj->contDistObj.decPrec, x, err);
}
cl_double clprobdistGammaDensityWithObject(_CLPROBDIST_GAMMA_OBJ_MEM const clprobdistGamma* distObj, cl_double x, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;

	if (x <= 0)
		return 0.0;

	cl_double z = distObj->params.logFactor + (distObj->params.alpha - 1.0) * log(x) - distObj->params.lambda * x;

	if (z > -clprobdistGammaXBIGM)
		return exp(z);
	else
		return 0.0;
}
cl_double clprobdistGammaMeanWithObject(_CLPROBDIST_GAMMA_OBJ_MEM const clprobdistGamma* distObj, clprobdistStatus* err) {
	return clprobdistGammaMean(distObj->params.alpha, distObj->params.lambda, distObj->contDistObj.decPrec, err);
}
cl_double clprobdistGammaVarianceWithObject(_CLPROBDIST_GAMMA_OBJ_MEM const clprobdistGamma* distObj, clprobdistStatus* err) {

	return clprobdistGammaVariance(distObj->params.alpha, distObj->params.lambda, distObj->contDistObj.decPrec, err);
}
cl_double clprobdistGammaStdDeviationWithObject(_CLPROBDIST_GAMMA_OBJ_MEM const clprobdistGamma* distObj, clprobdistStatus* err) {
	return clprobdistGammaStdDeviation(distObj->params.alpha, distObj->params.lambda, distObj->contDistObj.decPrec, err);

}

cl_double clprobdistGammaGetAlpha(_CLPROBDIST_GAMMA_OBJ_MEM const clprobdistGamma* distObj, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;
	return distObj->params.alpha;
}
cl_double clprobdistGammaGetLambda(_CLPROBDIST_GAMMA_OBJ_MEM const clprobdistGamma* distObj, clprobdistStatus* err) {
	
	if (err) *err = CLPROBDIST_SUCCESS;
	return distObj->params.lambda;
}

#endif
