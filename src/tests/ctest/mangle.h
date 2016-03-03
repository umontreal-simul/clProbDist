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

#pragma once
#ifndef CTEST_COMMON_H
#define CTEST_COMMON_H

#ifdef __APPLE__
#include <Opencl/cl.h>
#else
#include <CL/cl.h>
#endif


#ifndef CTEST_PROBDIST_TYPE
#error "CTEST_PROBDIST_TYPE undefined"
#endif

#ifndef CTEST_PROBDIST_HEADER
#error "CTEST_PROBDIST_HEADER undefined"
#endif

#ifndef CTEST_PROBDIST_PARAM_COUNT
#error "CTEST_PROBDIST_PARAM_COUNT undefined"
#endif

#define PROBDIST_HOST_HEADER   <clProbDist/CTEST_PROBDIST_HEADER.h>
#define PROBDIST_DEVICE_HEADER <clProbDist/CTEST_PROBDIST_HEADER.clh>

#define _PROBDIST_MANGLE(ident)               _PROBDIST_MANGLE_(ident,CTEST_PROBDIST_TYPE)
#define _PROBDIST_MANGLE_(ident,PROBDIST)     _PROBDIST_MANGLE__(ident,PROBDIST)
#define _PROBDIST_MANGLE__(ident,probdist)    clprobdist ## probdist ## ident

#define CTEST_MANGLE(ident)        CTEST_MANGLE_(ident,CTEST_PROBDIST_TYPE)
#define CTEST_MANGLE_(ident,probdist)   CTEST_MANGLE__(ident,probdist)
#define CTEST_MANGLE__(ident,probdist)  ctest ## probdist ## _ ## ident

// utility macros
#define CTEST_ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0]))

#if CTEST_PROBDIST_PARAM_COUNT == 1
#define CTEST_DIST_PARAMS(params) (params).p1
#elif CTEST_PROBDIST_PARAM_COUNT == 2
#define CTEST_DIST_PARAMS(params) (params).p1, (params).p2
#elif CTEST_PROBDIST_PARAM_COUNT == 3
#define CTEST_DIST_PARAMS(params) (params).p1, (params).p2, (params).p3
#else
#error "unsupported value of CTEST_PROBDIST_PARAM_COUNT"
#endif

// for size types
#if defined(_MSC_VER) || defined(__MINGW32__)
    #define SIZE_T_FORMAT   "Iu"
#elif defined(__GNUC__)
    #define SIZE_T_FORMAT   "zu"
#else
    #define SIZE_T_FORMAT   "lu"
#endif

#define _CTEST_STR_(x)  _CTEST_STR__(x)
#define _CTEST_STR__(x) #x

#define PROBDIST_TYPE_S          _CTEST_STR_(CTEST_PROBDIST_TYPE)
#define PROBDIST_TYPE_SC         _CTEST_STR_(CTEST_PROBDIST_TYPE_C)
#define PROBDIST_DEVICE_HEADER_S _CTEST_STR_(PROBDIST_DEVICE_HEADER)

#define clprobdistObject                    _PROBDIST_MANGLE()
#define clprobdistCreate                    _PROBDIST_MANGLE(Create)
#define clprobdistDestroy                   _PROBDIST_MANGLE(Destroy)
#define clprobdistDensityWithObject         _PROBDIST_MANGLE(DensityWithObject)
#define clprobdistProbWithObject            _PROBDIST_MANGLE(ProbWithObject)
#define clprobdistCDFWithObject             _PROBDIST_MANGLE(CDFWithObject)
#define clprobdistComplCDFWithObject        _PROBDIST_MANGLE(ComplCDFWithObject)
#define clprobdistInverseCDFWithObject      _PROBDIST_MANGLE(InverseCDFWithObject)
#define clprobdistMeanWithObject            _PROBDIST_MANGLE(MeanWithObject)
#define clprobdistVarianceWithObject        _PROBDIST_MANGLE(VarianceWithObject)
#define clprobdistStdDeviationWithObject    _PROBDIST_MANGLE(StdDeviationWithObject)
#define clprobdistDensity                   _PROBDIST_MANGLE(Density)
#define clprobdistProb                      _PROBDIST_MANGLE(Prob)
#define clprobdistCDF                       _PROBDIST_MANGLE(CDF)
#define clprobdistComplCDF                  _PROBDIST_MANGLE(ComplCDF)
#define clprobdistInverseCDF                _PROBDIST_MANGLE(InverseCDF)
#define clprobdistMean                      _PROBDIST_MANGLE(Mean)
#define clprobdistVariance                  _PROBDIST_MANGLE(Variance)
#define clprobdistStdDeviation              _PROBDIST_MANGLE(StdDeviation)

#include PROBDIST_HOST_HEADER

// defined in util.c
extern int ctestVerbose;

/*! @brief Enumeration for memory types
 */
typedef enum MemType_ { MEM_PARAMS, MEM_PRIVATE, MEM_LOCAL, MEM_CONSTANT, MEM_GLOBAL } MemType;

#endif
