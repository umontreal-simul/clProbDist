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

/* @file clProbDist.c
* @brief Implementation of functions defined in clProbDist.h and private.h
*/

#include <clProbDist/clProbDist.h>
#include "private.h"
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#ifdef _MSC_VER
#include <io.h>
#else
#include <unistd.h>
#endif

#define CASE_ERR_(code,msg) case code: base = msg; break
#define CASE_ERR(code)      CASE_ERR_(CLPROBDIST_ ## code, MSG_ ## code)

extern char errorString[1024];

const char* clprobdistGetErrorString()
{
	return errorString;
}


static char lib_path_default1[] = "/usr";
static char lib_path_default1_check[] = "/usr/include/clProbDist/clProbDist.h";
static char lib_path_default2[] = ".";
static char lib_path_default2_check[] = "./include/clProbDist/clProbDist.h";

const char* clprobdistGetLibraryRoot()
{
	const char* lib_path = getenv("CLPROBDIST_ROOT");

	if (lib_path != NULL && lib_path[0] != 0)
			return lib_path;

	// check if lib_path_default1_check exists
	if (
#ifdef _MSC_VER
	_access(lib_path_default1_check, 0) != -1
#else
	access(lib_path_default1_check, F_OK) != -1
#endif
	)
		return lib_path_default1;

	// last resort
	if (
#ifdef _MSC_VER
	_access(lib_path_default2_check, 0) != -1
#else
	access(lib_path_default2_check, F_OK) != -1
#endif
	)
		return lib_path_default2;

	return NULL;
}


static char lib_includes[1024];

const char* clprobdistGetLibraryDeviceIncludes(cl_int* err)
{
	int nbytes;
	const char* root = clprobdistGetLibraryRoot();

	if (err)
			*err = CLPROBDIST_SUCCESS;

	if (root == NULL) {
		if (err)
			*err = clprobdistSetErrorString(CLPROBDIST_INVALID_ENVIRONMENT, "environment variable CLPROBDIST_ROOT not set");
		return NULL;
	}
#ifdef _MSC_VER
	nbytes = sprintf_s(
#else
	nbytes = snprintf(
#endif
		lib_includes,
		sizeof(lib_includes),
		"-I\"%s/include\"",
		root);

#ifdef _MSC_VER
	if (nbytes < 0) {
#else
	if (nbytes >= sizeof(lib_includes)) {
#endif
		if (err)
			*err = clprobdistSetErrorString(CLPROBDIST_OUT_OF_RESOURCES, "value of CLPROBDIST_ROOT too long (max = %u)", sizeof(lib_includes) - 16);
		return NULL;
	}
	return lib_includes;
	}
