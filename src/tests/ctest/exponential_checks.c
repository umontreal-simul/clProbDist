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
 *   David Munger <mungerd@iro.umontreal.ca>        (2015)
 *   Nabil Kemerchou <kemerchn@iro.umontreal.ca>    (2015)
 *   Pierre L'Ecuyer <lecuyer@iro.umontreal.ca>     (2015)
 *
 *
 ***********************************************************************
 Copyright (c) 2015 Advanced Micro Devices, Inc. 
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without 
 modification, are permitted provided that the following conditions 
 are met:
 
 1. Redistributions of source code must retain the above copyright 
 notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright 
 notice, this list of conditions and the following disclaimer in the 
 documentation and/or other materials provided with the distribution.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
 HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
 ***********************************************************************
 */

#define CTEST_PROBDIST_TYPE    Exponential
#define CTEST_PROBDIST_TYPE_C  EXPONENTIAL
#define CTEST_PROBDIST_HEADER exponential
#define CTEST_PROBDIST_PARAM_COUNT 1
#define CTEST_PROBDIST_PARAM_DEF cl_double p1;
#include "mangle.h"

typedef struct { CTEST_PROBDIST_PARAM_DEF } DistParams;

DistParams CTEST_MANGLE(inputDistParams)[] = {
  { 1e-9  },
  { 1.0   },
  { 1e9 }
};

cl_double CTEST_MANGLE(inputCDFValues)[] = {
  1e-9,
  1e-6,
  1e-3,
  0.01,
  0.1,
  0.5,
  0.9,
  0.99,
  1 - 1e-3,
  1 - 1e-6,
  1 - 1e-9
};

MemType CTEST_MANGLE(objMemTypes)[] = {
  MEM_PARAMS,
  MEM_GLOBAL,
  MEM_CONSTANT,
  MEM_LOCAL,
  MEM_PRIVATE
};

// The following values were obtained from SSJ:
cl_double CTEST_MANGLE(expectedInverseCDFValues)
  [CTEST_ARRAY_SIZE(CTEST_MANGLE(inputDistParams))]
  [CTEST_ARRAY_SIZE(CTEST_MANGLE(inputCDFValues))] = {
  { 1.00000000050000000, 1000.00050000033330, 1000500.33358353340, 10050335.8535014410, 105360515.657826300, 693147180.559945200, 2302585092.99404570, 4605170185.98809000, 6907755278.98213600, 13815510557.9355160, 20723265865.2283400 },
  { 1.00000000050000010e-09, 1.00000050000033340e-06, 0.00100050033358353350, 0.0100503358535014420, 0.105360515657826310, 0.693147180559945300, 2.30258509299404600, 4.60517018598809000, 6.90775527898213600, 13.8155105579355180, 20.7232658652283420 },
  { 1.00000000050000020e-18, 1.00000050000033350e-15, 1.00050033358353350e-12, 1.00503358535014410e-11, 1.05360515657826310e-10, 6.93147180559945300e-10, 2.30258509299404600e-09, 4.60517018598809000e-09, 6.90775527898213600e-09, 1.38155105579355170e-08, 2.07232658652283400e-08 }
};
cl_double CTEST_MANGLE(expectedDensityValues)
  [CTEST_ARRAY_SIZE(CTEST_MANGLE(inputDistParams))]
  [CTEST_ARRAY_SIZE(CTEST_MANGLE(inputCDFValues))] = {
  { 9.99999999000000100e-10, 9.99999000000000000e-10, 9.99000000000000000e-10, 9.90000000000000000e-10, 9.00000000000000100e-10, 5.00000000000000000e-10, 9.99999999999999800e-11, 1.00000000000000140e-11, 1.00000000000000120e-12, 1.00000000002875620e-15, 9.99999971718069600e-19 },
  { 0.999999999000000000, 0.999999000000000000, 0.999000000000000000, 0.990000000000000000, 0.900000000000000000, 0.500000000000000000, 0.0999999999999999800, 0.0100000000000000140, 0.00100000000000000100, 1.00000000002875600e-06, 9.99999971718069600e-10 },
  { 999999999.000000000, 999999000.000000000, 999000000.000000000, 990000000.000000000, 900000000.000000000, 500000000.000000000, 99999999.9999999900, 10000000.0000000150, 1000000.00000000100, 1000.00000002875610, 0.999999971718069500 }
};

#include "checks.c.h"
