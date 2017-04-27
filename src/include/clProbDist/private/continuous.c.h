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
#ifndef CLPROBDIST_PRIVATE_CONTINUOUSDIST_CH
#define CLPROBDIST_PRIVATE_CONTINUOUSDIST_CH


#if 0
/***********************************
* Getters and Setters
***********************************/
cl_int clprobdistContinuousGetDecPrec(clprobdistContinuous* distObj, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;
	return distObj->decPrec;
}

cl_double clprobdistContinuousGetXinf(clprobdistContinuous* distObj, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;
	return distObj->supportA;
}

cl_double clprobdistContinuousGetXsup(clprobdistContinuous* distObj, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;
	return distObj->supportB;
}

clprobdistStatus clprobdistContinuousSetXinf(cl_double xa, clprobdistContinuous* distObj) {
	distObj->supportA = xa;
	return CLPROBDIST_SUCCESS;
}

clprobdistStatus clprobdistContinuousSetXsup(cl_double xb, clprobdistContinuous* distObj) {
	distObj->supportB = xb;
	return CLPROBDIST_SUCCESS;
}



/***********************************
* Abstract functions
***********************************/
cl_double clprobdistContinuousCDF(cl_double x, clprobdistStatus* err) {
	if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): Abstract function", __func__);
	return -1;
}

clprobdistStatus clprobdistContinuousFindInterval(cl_double u, cl_double* iv, clprobdistContinuous* distObj, clprobdistStatus* err) {
	// Finds an interval [a, b] that certainly contains x defined as
	// u = cdf(x). The result is written in iv[0] = a and iv[1] = b.

	if (u > 1.0 || u < 0.0)
		return clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): u not in [0, 1]", __func__);

	const cl_double XLIM = DBL_MAX / 2.0;
	const cl_double B0 = 8.0;

	cl_double b = B0;
	while (b < XLIM && u > clprobdistContinuousCDF(b, err))
		b *= 2.0;
	if (b > B0) {
		iv[0] = b / 2.0;
		iv[1] = fmin(b, distObj->supportB);
		return CLPROBDIST_SUCCESS;
	}

	cl_double a = -B0;
	while (a > -XLIM && u < clprobdistContinuousCDF(a, err))
		a *= 2.0;
	if (a < -B0) {
		iv[1] = a / 2.0;
		iv[0] = fmax(a, distObj->supportA);
		return CLPROBDIST_SUCCESS;
	}
	iv[0] = fmax(a, distObj->supportA);
	iv[1] = fmin(b, distObj->supportB);

	return CLPROBDIST_SUCCESS;
}

cl_double clprobdistContinuousDensity(cl_double x, clprobdistStatus* err) {
	if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): Abstract function", __func__);
	return -1;
}

cl_double clprobdistContinuousGetMean(clprobdistStatus* err) {
	if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): Abstract function", __func__);
	return -1;
}

cl_double clprobdistContinuousGetVariance(clprobdistStatus* err) {
	if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): Abstract function", __func__);
	return -1;
}

cl_double clprobdistContinuousGetStdDeviation(clprobdistStatus* err) {
	if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): Abstract function", __func__);
	return -1;
}

/***********************************
* Default implementations of functions
***********************************/
cl_double clprobdistContinuousComplCDF(cl_double x, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;
	return 1.0 - clprobdistContinuousCDF(x, err);

}

cl_double clprobdistContinuousInverseBrent(cl_double a, cl_double b, cl_double u, cl_double tol, clprobdistContinuous* distObj, clprobdistStatus* err)  {
	if (u > 1.0 || u < 0.0) {
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): u not in [0, 1]", __func__);
		return -1;
	}
	if (b < a) {
		cl_double ctemp = a;   a = b;   b = ctemp;
	}
	if (u <= 0.0) {
		//println("********** WARNING,  inverseBrent:   u = 0");
		return distObj->supportA;
	}
	if (u >= 1.0) {
		//println("********** WARNING,  inverseBrent:   u = 1");
		return distObj->supportB;
	}
	int MAXITER = 50;      // Maximum number of iterations
	tol += clprobdistEPSARRAY[distObj->decPrec] + DBL_EPSILON;    // in case tol is too small

	cl_double ua = clprobdistContinuousCDF(a, err) - u;
	if (ua > 0.0) {
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): u < cdf(a)", __func__);
		return -1;
	}

	cl_double ub = clprobdistContinuousCDF(b, err) - u;
	if (ub < 0.0) {
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): u > cdf(b)", __func__);
		return -1;
	}

	/*bool DEBUG = false;
	if (DEBUG) {
	String ls = System.getProperty("line.separator");
	System.out.println(
	"-------------------------------------------------------------"
	+ ls + "u = " + PrintfFormat.g(20, 15, u));
	System.out.println
	(ls + "iter           b                  c               F(x) - u" + ls);
	}*/
	// Initialize
	cl_double c = a;
	cl_double uc = ua;
	cl_double len = b - a;
	cl_double t = len;
	if (fabs(uc) < fabs(ub)) {
		a = b; b = c; c = a;
		ua = ub; ub = uc; uc = ua;
	}
	int i;
	for (i = 0; i < MAXITER; ++i) {
		cl_double tol1 = tol + 4.0 * DBL_EPSILON * fabs(b);
		cl_double xm = 0.5 * (c - b);
		/*if (DEBUG) {
		System.out.println(PrintfFormat.d(3, i) + "  " +
		PrintfFormat.g(18, decPrec, b) + "  " +
		PrintfFormat.g(18, decPrec, c) + "  " +
		PrintfFormat.g(14, 4, ub));
		}*/
		if (fabs(ub) == 0.0 || (fabs(xm) <= tol1)) {
			if (b <= distObj->supportA) return distObj->supportA;
			if (b >= distObj->supportB) return distObj->supportB;
			return b;
		}

		cl_double s, p, q, r;
		if ((fabs(t) >= tol1) && (fabs(ua) > fabs(ub))) {
			if (a == c) {
				// linear interpolation
				s = ub / ua;
				q = 1.0 - s;
				p = 2.0 * xm * s;
			} else {
				// quadratic interpolation
				q = ua / uc;
				r = ub / uc;
				s = ub / ua;
				p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
				q = (q - 1.0) * (r - 1.0) * (s - 1.0);
			}
			if (p > 0.0)
				q = -q;
			p = fabs(p);

			// Accept interpolation?
			if ((2.0 * p >= (3.0 * xm * q - fabs(q * tol1))) ||
			        (p >= fabs(0.5 * t * q))) {
				len = xm;
				t = len;
			} else {
				t = len;
				len = p / q;
			}

		} else {
			len = xm;
			t = len;
		}

		a = b;
		ua = ub;
		if (fabs(len) > tol1)
			b += len;
		else if (xm < 0.0)
			b -= tol1;
		else
			b += tol1;

		ub = clprobdistContinuousCDF(b, err) - u;

		if (ub * (uc / fabs(uc)) > 0.0) {
			c = a;
			uc = ua;
			len = b - a;
			t = len;
		} else if (fabs(uc) < fabs(ub)) {
			a = b; b = c; c = a;
			ua = ub; ub = uc; uc = ua;
		}
	}
	if (i >= MAXITER) {
		/*	String lineSep = System.getProperty("line.separator");
		System.out.println(lineSep +"*********** inverseBrent:   no convergence after " + MAXITER +" iterations");*/
	}
	return b;
}

cl_double clprobdistContinuousInverseBisection(cl_double u, clprobdistContinuous* distObj, clprobdistStatus* err) {
	int MAXITER = 100;              // Maximum number of iterations
	cl_double EPSILON = clprobdistEPSARRAY[distObj->decPrec];  // Absolute precision

	if (u > 1.0 || u < 0.0) {
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): u not in [0, 1]", __func__);
		return -1;
	}

	if (distObj->decPrec > DBL_DIG)	{
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): decPrec too large", __func__);
		return -1;
	}

	if (distObj->decPrec <= 0)	{
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): decPrec <= 0", __func__);
		return -1;
	}


	cl_double x = 0.0;
	if (u <= 0.0) {
		x = distObj->supportA;
		return x;
	}
	if (u >= 1.0) {
		x = distObj->supportB;
		return x;
	}

	cl_double iv[2] = { 0, 0 };
	clprobdistContinuousFindInterval(u, iv, distObj, err);
	cl_double xa = iv[0];
	cl_double xb = iv[1];

	cl_double ya = clprobdistContinuousCDF(xa, err) - u;
	cl_double y;

	cl_bool fini = CL_FALSE;
	int i = 0;
	while (!fini) {
		x = (xa + xb) / 2.0;

		y = clprobdistContinuousCDF(x, err) - u;
		if ((y == 0.0) ||
		        (fabs((xb - xa) / (x + DBL_EPSILON)) <= EPSILON)) {
			fini = CL_TRUE;
		} else if (y * ya < 0.0) {
			xb = x;
		} else
			xa = x;
		++i;

		if (i > MAXITER) {
			fini = CL_TRUE;
		}
	}
	return x;
}

cl_double clprobdistContinuousInverseCDF(cl_double u, clprobdistContinuous* distObj, clprobdistStatus* err) {

	cl_double iv[2] = { 0, 0 };

	clprobdistContinuousFindInterval(u, iv, distObj, err);

	return clprobdistContinuousInverseBrent(iv[0], iv[1], u, clprobdistEPSARRAY[distObj->decPrec], distObj, err);
}
#endif

#endif