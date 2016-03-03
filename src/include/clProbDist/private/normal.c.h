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
#ifndef CLPROBDIST_PRIVATE_NORMALDIST_CH
#define CLPROBDIST_PRIVATE_NORMALDIST_CH

#ifndef _CLPROBDIST_NORMAL_OBJ_MEM
#define _CLPROBDIST_NORMAL_OBJ_MEM
#endif

struct _clprobdistNormal {
	clprobdistContinuous contDistObj;
	cl_double mu;
	cl_double sigma;
};

//Private variables and Functions
constant cl_double clprobdistNormal_SQRT2PI = 2.50662827463100050; // Sqrt(2*Pi)
constant cl_double clprobdistNormal_SQRT2 = 1.4142135623730951;

constant const cl_double clprobdistNormalXBIG = 100.0;
constant const cl_double clprobdistNormalXBIGM = 1000.0;

/*
* The precision of cl_double is 16 decimals; we shall thus use coeffmax = 24
* coefficients. But the approximation is good to 30 decimals of precision
* with 44 coefficients.
*/
constant int clprobdistNormal_COEFFMAX = 24;

constant cl_double clprobdistNormal_AbarF[] = {
	6.10143081923200418E-1,
	-4.34841272712577472E-1,
	1.76351193643605501E-1,
	-6.07107956092494149E-2,
	1.77120689956941145E-2,
	-4.32111938556729382E-3,
	8.54216676887098679E-4,
	-1.27155090609162743E-4,
	1.12481672436711895E-5,
	3.13063885421820973E-7,
	-2.70988068537762022E-7,
	3.07376227014076884E-8,
	2.51562038481762294E-9,
	-1.02892992132031913E-9,
	2.99440521199499394E-11,
	2.60517896872669363E-11,
	-2.63483992417196939E-12,
	-6.43404509890636443E-13,
	1.12457401801663447E-13,
	1.7281533389986098E-14,
	-4.2641016949424E-15,
	-5.4537197788E-16,
	1.5869760776E-16,
	2.08998378E-17,
	-0.5900E-17
};

static cl_double clprobdistNormalEvalCheby(constant cl_double* a, int n, cl_double x, clprobdistStatus* err)  {
	if (fabs(x) > 1.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): Chebychev polynomial evaluated at x outside [-1, 1]", __func__);
		return 0;
	}

	cl_double xx = 2.0*x;
	cl_double b0 = 0.0;
	cl_double b1 = 0.0;
	cl_double b2 = 0.0;
	for (int j = n; j >= 0; j--) {
		b2 = b1;
		b1 = b0;
		b0 = (xx*b1 - b2) + a[j];
	}

	if (err) *err = CLPROBDIST_SUCCESS;
	return (b0 - b2) / 2.0;
}

//***********************************************************************
// Implementation of static functions
//***********************************************************************
cl_double clprobdistStdNormalDensity(cl_double x, clprobdistStatus* err)  {
	if (err) *err = CLPROBDIST_SUCCESS;
	return exp(-0.5*x*x) / clprobdistNormal_SQRT2PI;
}

cl_double clprobdistNormalDensity(cl_double mu, cl_double sigma, cl_double x, clprobdistStatus* err)  {
	if (sigma <= 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): sigma <= 0", __func__);
		return 0;
	}
	else{
		cl_double z = (x - mu) / sigma;
		return exp(-0.5*z*z) / (clprobdistNormal_SQRT2PI*sigma);
	}
}


constant cl_double clprobdistNormal_NORMAL2_A[] = {
	6.10143081923200417926465815756e-1,
	-4.34841272712577471828182820888e-1,
	1.76351193643605501125840298123e-1,
	-6.0710795609249414860051215825e-2,
	1.7712068995694114486147141191e-2,
	-4.321119385567293818599864968e-3,
	8.54216676887098678819832055e-4,
	-1.27155090609162742628893940e-4,
	1.1248167243671189468847072e-5,
	3.13063885421820972630152e-7,
	-2.70988068537762022009086e-7,
	3.0737622701407688440959e-8,
	2.515620384817622937314e-9,
	-1.028929921320319127590e-9,
	2.9944052119949939363e-11,
	2.6051789687266936290e-11,
	-2.634839924171969386e-12,
	-6.43404509890636443e-13,
	1.12457401801663447e-13,
	1.7281533389986098e-14,
	-4.264101694942375e-15,
	-5.45371977880191e-16,
	1.58697607761671e-16,
	2.0899837844334e-17,
	-5.900526869409e-18,
	-9.41893387554e-19
	/*,     2.14977356470e-19,
	4.6660985008e-20,
	-7.243011862e-21,
	-2.387966824e-21,
	1.91177535e-22,
	1.20482568e-22,
	-6.72377e-25,
	-5.747997e-24,
	-4.28493e-25,
	2.44856e-25,
	4.3793e-26,
	-8.151e-27,
	-3.089e-27,
	9.3e-29,
	1.74e-28,
	1.6e-29,
	-8.0e-30,
	-2.0e-30
	*/
};

cl_double clprobdistStdNormalCDF(cl_double x, clprobdistStatus* err) {
	/*
	* Returns P[X < x] for the normal distribution.
	* As in J. L. Schonfelder, Math. of Computation, Vol. 32,
	* pp 1232--1240, (1978).
	*/

	cl_double t, r;

	if (x <= -clprobdistNormalXBIG)
		return 0.0;
	if (x >= clprobdistNormalXBIG)
		return 1.0;

	x = -x / clprobdistNormal_SQRT2;
	if (x < 0) {
		x = -x;
		t = (x - 3.75) / (x + 3.75);
		r = 1.0 - 0.5 * exp(-x * x) * clprobdistNormalEvalCheby(clprobdistNormal_NORMAL2_A, clprobdistNormal_COEFFMAX, t, err);
	}
	else {
		t = (x - 3.75) / (x + 3.75);
		r = 0.5 * exp(-x * x) * clprobdistNormalEvalCheby(clprobdistNormal_NORMAL2_A, clprobdistNormal_COEFFMAX, t, err);
	}
	return r;
}

cl_double clprobdistNormalCDF(cl_double mu, cl_double sigma, cl_double x, clprobdistStatus* err) {
	if (sigma <= 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): sigma <= 0", __func__);
		return 0;
	}
	else
		return clprobdistStdNormalCDF((x - mu) / sigma, err);
}


cl_double clprobdistStdNormalComplCDF(cl_double x, clprobdistStatus* err) {
	/*
	* Returns P[X >= x] = 1 - F (x) where F is the normal distribution by
	* computing the complementary distribution directly; it is thus more
	* precise in the tail.
	*/

	cl_double KK = 5.30330085889910643300;      // 3.75 Sqrt (2)
	cl_double y, t;
	int neg;

	if (x >= clprobdistNormalXBIG)
		return 0.0;
	if (x <= -clprobdistNormalXBIG)
		return 1.0;

	if (x >= 0.0)
		neg = 0;
	else {
		neg = 1;
		x = -x;
	}

	t = (x - KK) / (x + KK);
	y = clprobdistNormalEvalCheby(clprobdistNormal_AbarF, 24, t, err);
	y = y * exp(-x * x / 2.0) / 2.0;

	if (err) *err = CLPROBDIST_SUCCESS;
	if (neg == 1)
		return 1.0 - y;
	else
		return y;
}

cl_double clprobdistNormalComplCDF(cl_double mu, cl_double sigma, cl_double x, clprobdistStatus* err) {
	if (sigma <= 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): sigma <= 0", __func__);
		return 0;
	}
	else
		return clprobdistStdNormalComplCDF((x - mu) / sigma, err);
}

constant cl_double clprobdistNormal_InvP1[] = {
	0.160304955844066229311E2,
	-0.90784959262960326650E2,
	0.18644914861620987391E3,
	-0.16900142734642382420E3,
	0.6545466284794487048E2,
	-0.864213011587247794E1,
	0.1760587821390590
};

constant cl_double clprobdistNormal_InvQ1[] = {
	0.147806470715138316110E2,
	-0.91374167024260313396E2,
	0.21015790486205317714E3,
	-0.22210254121855132366E3,
	0.10760453916055123830E3,
	-0.206010730328265443E2,
	0.1E1
};

constant cl_double clprobdistNormal_InvP2[] = {
	-0.152389263440726128E-1,
	0.3444556924136125216,
	-0.29344398672542478687E1,
	0.11763505705217827302E2,
	-0.22655292823101104193E2,
	0.19121334396580330163E2,
	-0.5478927619598318769E1,
	0.237516689024448000
};

constant cl_double clprobdistNormal_InvQ2[] = {
	-0.108465169602059954E-1,
	0.2610628885843078511,
	-0.24068318104393757995E1,
	0.10695129973387014469E2,
	-0.23716715521596581025E2,
	0.24640158943917284883E2,
	-0.10014376349783070835E2,
	0.1E1
};

constant cl_double clprobdistNormal_InvP3[] = {
	0.56451977709864482298E-4,
	0.53504147487893013765E-2,
	0.12969550099727352403,
	0.10426158549298266122E1,
	0.28302677901754489974E1,
	0.26255672879448072726E1,
	0.20789742630174917228E1,
	0.72718806231556811306,
	0.66816807711804989575E-1,
	-0.17791004575111759979E-1,
	0.22419563223346345828E-2
};

constant cl_double clprobdistNormal_InvQ3[] = {
	0.56451699862760651514E-4,
	0.53505587067930653953E-2,
	0.12986615416911646934,
	0.10542932232626491195E1,
	0.30379331173522206237E1,
	0.37631168536405028901E1,
	0.38782858277042011263E1,
	0.20372431817412177929E1,
	0.1E1
};

cl_double clprobdistStdNormalInverseCDF(cl_double u, clprobdistStatus* err) {
	/*
	* Returns the inverse of the cdf of the normal distribution.
	* Rational approximations giving 16 decimals of precision.
	* J.M. Blair, C.A. Edwards, J.H. Johnson, "Rational Chebyshev
	* approximations for the Inverse of the Error Function", in
	* Mathematics of Computation, Vol. 30, 136, pp 827, (1976)
	*/

	int i;
	cl_bool negatif;
	cl_double y, z, v, w;
	cl_double x = u;

	if (u < 0.0 || u > 1.0)
	{
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): u is not in [0, 1]", __func__);
		return -1;
	}
	if (u <= 0.0)
		return DBL_MIN;// Double.NEGATIVE_INFINITY;
	if (u >= 1.0)
		return DBL_MAX; // Double.POSITIVE_INFINITY;

	// Transform x as argument of InvErf
	x = 2.0 * x - 1.0;
	if (x < 0.0) {
		x = -x;
		negatif = CL_TRUE;
	}
	else {
		negatif = CL_FALSE;
	}

	if (x <= 0.75) {
		y = x * x - 0.5625;
		v = w = 0.0;
		for (i = 6; i >= 0; i--) {
			v = v * y + clprobdistNormal_InvP1[i];
			w = w * y + clprobdistNormal_InvQ1[i];
		}
		z = (v / w) * x;

	}
	else if (x <= 0.9375) {
		y = x * x - 0.87890625;
		v = w = 0.0;
		for (i = 7; i >= 0; i--) {
			v = v * y + clprobdistNormal_InvP2[i];
			w = w * y + clprobdistNormal_InvQ2[i];
		}
		z = (v / w) * x;

	}
	else {
		if (u > 0.5)
			y = 1.0 / sqrt(-log(1.0 - x));
		else
			y = 1.0 / sqrt(-log(2.0 * u));
		v = 0.0;
		for (i = 10; i >= 0; i--)
			v = v * y + clprobdistNormal_InvP3[i];
		w = 0.0;
		for (i = 8; i >= 0; i--)
			w = w * y + clprobdistNormal_InvQ3[i];
		z = (v / w) / y;
	}

	if (negatif) {
		if (u < 1.0e-105) {
			cl_double RACPI = 1.77245385090551602729;
			w = exp(-z * z) / RACPI;  // pdf
			y = 2.0 * z * z;
			v = 1.0;
			cl_double term = 1.0;
			// Asymptotic series for erfc(z) (apart from exp factor)
			for (i = 0; i < 6; ++i) {
				term *= -(2 * i + 1) / y;
				v += term;
			}
			// Apply 1 iteration of Newton solver to get last few decimals
			z -= u / w - 0.5 * v / z;

		}
		return -(z * clprobdistNormal_SQRT2);

	}
	else
		return z * clprobdistNormal_SQRT2;
}

cl_double clprobdistNormalInverseCDF(cl_double mu, cl_double sigma, cl_double u, clprobdistStatus* err)  {
	if (sigma <= 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): sigma <= 0", __func__);
		return 0;
	}
	else
		return mu + sigma*clprobdistStdNormalInverseCDF(u, err);
}



cl_double clprobdistNormalMean(cl_double mu, cl_double sigma, clprobdistStatus* err) {
	if (sigma <= 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): sigma <= 0", __func__);
		return 0;
	}
	else{
		if (err) *err = CLPROBDIST_SUCCESS;
		return mu;
	}
}

cl_double clprobdistNormalVariance(cl_double mu, cl_double sigma, clprobdistStatus* err) {
	if (sigma <= 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): sigma <= 0", __func__);
		return 0;
	}
	else{
		if (err) *err = CLPROBDIST_SUCCESS;
		return (sigma * sigma);
	}
}

cl_double clprobdistNormalStdDeviation(cl_double mu, cl_double sigma, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;
	return sigma;
}

//***********************************************************************
// Implementation with Dist object
//***********************************************************************

cl_double clprobdistNormalDensityWithObject(_CLPROBDIST_NORMAL_OBJ_MEM const clprobdistNormal* dist, cl_double x, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;
	cl_double z = (x - dist->mu) / dist->sigma;
	return exp(-0.5*z*z) / (clprobdistNormal_SQRT2PI*dist->sigma);
}

cl_double clprobdistNormalCDFWithObject(_CLPROBDIST_NORMAL_OBJ_MEM const clprobdistNormal* dist, cl_double x, clprobdistStatus* err) {
	return clprobdistStdNormalCDF((x - dist->mu) / dist->sigma, err);
}

cl_double clprobdistNormalComplCDFWithObject(_CLPROBDIST_NORMAL_OBJ_MEM const clprobdistNormal* dist, cl_double x, clprobdistStatus* err) {
	return clprobdistStdNormalComplCDF((x - dist->mu) / dist->sigma, err);
}

cl_double clprobdistNormalInverseCDFWithObject(_CLPROBDIST_NORMAL_OBJ_MEM const clprobdistNormal* dist, cl_double u, clprobdistStatus* err) {
	return dist->mu + dist->sigma * clprobdistStdNormalInverseCDF(u, err);
}

cl_double clprobdistNormalMeanWithObject(_CLPROBDIST_NORMAL_OBJ_MEM const clprobdistNormal* dist, clprobdistStatus* err) {
	return clprobdistNormalMean(dist->mu, dist->sigma, err);
}

cl_double clprobdistNormalVarianceWithObject(_CLPROBDIST_NORMAL_OBJ_MEM const clprobdistNormal* dist, clprobdistStatus* err) {
	return clprobdistNormalVariance(dist->mu, dist->sigma, err);
}

cl_double clprobdistNormalStdDeviationWithObject(_CLPROBDIST_NORMAL_OBJ_MEM const clprobdistNormal* dist, clprobdistStatus* err) {
	return clprobdistNormalStdDeviation(dist->mu, dist->sigma, err);
}

cl_double clprobdistNormalGetMu(_CLPROBDIST_NORMAL_OBJ_MEM const clprobdistNormal* dist, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;
	return dist->mu;
}

cl_double clprobdistNormalGetSigma(_CLPROBDIST_NORMAL_OBJ_MEM const clprobdistNormal* dist, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;
	return dist->sigma;
}
#endif
