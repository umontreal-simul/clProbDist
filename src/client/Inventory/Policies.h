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
 *
 */

#ifndef POLICIES_SS_H
#define POLICIES_SS_H

#include "Types.h"
#include "../common.h"

typedef struct OnePolicyData_ {
	int n;
	int n1;
	int m;
	int s;
	int S;
	int selectedDist;
	double lambda;
	double mu;
	double sigma;
	ExecType execType;
	simResult* SimResults;
} OnePolicyData;

typedef struct SeveralPoliciesData_ {
	int n;
	int n1;
	int m;
	int* s;
	int* S;
	int P;
	int selectedDist;
	double lambda;
	double mu;
	double sigma;
	ExecOption optionType;
	simResult * SimResults;
	bool useCreator;
} SeveralPoliciesData;

int one_Policy(cl_context context, cl_device_id device, cl_command_queue queue, void* data_);
int several_Policies(cl_context context, cl_device_id device, cl_command_queue queue, void* data_);

#endif
