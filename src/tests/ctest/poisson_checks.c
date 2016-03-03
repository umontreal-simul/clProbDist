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

#define CTEST_PROBDIST_TYPE   Poisson
#define CTEST_PROBDIST_TYPE_C POISSON
#define CTEST_PROBDIST_HEADER poisson
#define CTEST_PROBDIST_PARAM_COUNT 1
#define CTEST_PROBDIST_DISCRETE 1
#define CTEST_PROBDIST_PARAM_DEF cl_double p1;
// The Poisson dist sometimes rely on approximations
#define CTEST_REL_TOL 1e-12
#include "mangle.h"

typedef struct { CTEST_PROBDIST_PARAM_DEF } DistParams;

DistParams CTEST_MANGLE(inputDistParams)[] = {
  { 1.0    },
  { 10.0   },
  { 50.0   },
  { 150.0  },
  { 300.0  },
  { 1000.0 }
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
  MEM_LOCAL
};

// The following values were obtained from SSJ:
cl_int CTEST_MANGLE(expectedInverseCDFValues)
  [CTEST_ARRAY_SIZE(CTEST_MANGLE(inputDistParams))]
  [CTEST_ARRAY_SIZE(CTEST_MANGLE(inputCDFValues))] = {
  {                  0,                  0,                  0,                  0,                  0,                  1,                  2,                  4,                  5,                  9,                 11 },
  {                  0,                  0,                  2,                  3,                  6,                 10,                 14,                 18,                 21,                 28,                 34 },
  {                 14,                 20,                 30,                 34,                 41,                 50,                 59,                 67,                 73,                 87,                 98 },
  {                 83,                 96,                114,                122,                134,                150,                166,                179,                189,                212,                229 },
  {                202,                221,                248,                260,                278,                300,                322,                341,                355,                386,                410 },
  {                816,                853,                904,                927,                960,               1000,               1041,               1074,               1099,               1154,               1195 }
};
cl_double CTEST_MANGLE(expectedDensityValues)
  [CTEST_ARRAY_SIZE(CTEST_MANGLE(inputDistParams))]
  [CTEST_ARRAY_SIZE(CTEST_MANGLE(inputCDFValues))] = {
  { 0.367879441171442300, 0.367879441171442300, 0.367879441171442300, 0.367879441171442300, 0.367879441171442300, 0.367879441171442300, 0.183939720585721140, 0.0153283100488100940, 0.00306566200976201900, 1.01377711963029740e-06, 9.21615563300270300e-09 },
  { 4.53999297624848400e-05, 4.53999297624848400e-05, 0.00226999648812424220, 0.00756665496041414050, 0.0630554580034511800, 0.125110035721133270, 0.0520771044460261660, 0.00709110899319528600, 0.000888610149523218700, 1.48906740991700640e-06, 1.53776714207131960e-09 },
  { 1.35035393235192670e-09, 7.56051491837190800e-07, 0.000677198457150210000, 0.00380269460244236700, 0.0262190624351560700, 0.0563250063251908160, 0.0241258465748114760, 0.00358358081694166950, 0.000456843531868277570, 5.91351469405618400e-07, 6.45606463450500900e-10 },
  { 7.50408128976223800e-10, 5.81055924561978900e-07, 0.000334802972926895550, 0.00221015265348031900, 0.0142089375548418500, 0.0325554094568365550, 0.0135692580550197070, 0.00213041477695351780, 0.000269115952914787400, 3.25067613709759250e-07, 4.49439483559981000e-10 },
  { 3.84337548345669500e-10, 2.83394166462109060e-07, 0.000210029142768276800, 0.00151342419780077860, 0.0104600660933882300, 0.0230265461491873280, 0.0101111760118445210, 0.00147776414404461190, 0.000181497454257799480, 2.52902321578535940e-07, 2.78712248886699100e-10 },
  { 1.97545172513642100e-10, 1.56523912987383600e-07, 0.000113306390183136550, 0.000852943819568859200, 0.00572234128634273200, 0.0126146113487215060, 0.00539528353494033150, 0.000840558034942925800, 0.000104488386665294580, 1.46427495762500900e-07, 1.97230004313202000e-10 }
};
cl_double CTEST_MANGLE(expectedCDFValues)
  [CTEST_ARRAY_SIZE(CTEST_MANGLE(inputDistParams))]
  [CTEST_ARRAY_SIZE(CTEST_MANGLE(inputCDFValues))] = {
  { 0.367879441171442300, 0.367879441171442300, 0.367879441171442300, 0.367879441171442300, 0.367879441171442300, 0.735758882342884600, 0.919698602928605800, 0.996340153172656300, 0.999405815182418300, 0.999999888574521700, 0.999999999168389200 },
  { 4.53999297624848400e-05, 4.53999297624848400e-05, 0.00276939571551157540, 0.0103360506759257160, 0.130141420882482930, 0.583039750192985400, 0.916541527065337200, 0.992813495396145700, 0.999300349487665200, 0.999999235533359200, 0.999999999393967700 },
  { 1.85680233651023700e-09, 1.23518722187100030e-06, 0.00159402731860628900, 0.0107814591643343260, 0.112289062553116980, 0.537516690853147500, 0.907734948041066900, 0.991121004183620300, 0.999113567411716900, 0.999999245925228900, 0.999999999353929600 },
  { 1.65129210220285830e-09, 1.56687095880790100e-06, 0.00129899516913011300, 0.0105625116801581750, 0.101243974586825740, 0.521697179707476900, 0.909383377505140500, 0.990582052783946100, 0.999065454779913900, 0.999999252131021500, 0.999999999175605200 },
  { 1.15362263576127200e-09, 1.04223410646353050e-06, 0.00112442471756642160, 0.0100710239375584420, 0.106227073983067510, 0.515348757262923400, 0.901959352659438200, 0.990692156389725600, 0.999101143081445800, 0.999999158740829900, 0.999999999263733200 },
  { 1.04963971946566310e-09, 1.02735972464810400e-06, 0.00108942507680500140, 0.0102821571484102950, 0.105256720664464310, 0.508409367168506400, 0.904674696358401400, 0.990168693070978100, 0.999037369594133500, 0.999999090341006800, 0.999999999017994000 }
};

#include "checks.c.h"
