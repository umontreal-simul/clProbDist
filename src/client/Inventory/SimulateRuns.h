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

#ifndef SIMULATERUNS_SS_H
#define SIMULATERUNS_SS_H

#include <clRNG/mrg32k3a.h>
#include <clProbDist/poisson.h>
#include <clProbDist/normal.h>

#include "../common.h"

//Device Simulation
clrngStatus inventorySimulateRunsGPU(cl_context context, cl_device_id device, cl_command_queue queue, 
	                                 int m, int* s, int* S, int P, int n, int n1, int n2, ExecOption Option, 
									 const char * kernelName, 
									 size_t streamsBufSize, clrngMrg32k3aStream * streams_demand, clrngMrg32k3aStream * streams_order, 
									 int selectedDist, double* lambdaArr, 
									 size_t poissonBufSize, clprobdistPoisson* poissonDist, size_t normalBufSize, clprobdistNormal* normalDist,
									 double *stat, simResult * results);

//CPU Simulation
clrngStatus inventorySimulateRunsCPU(int m, int s, int S, int n, 
	                                double* lambdaArr,clrngMrg32k3aStream* streams_demand, clrngMrg32k3aStream* streams_order, 
									int selectedDist, clprobdistPoisson* poissonDist, clprobdistNormal* normalDist,
									double *stat_profit, ExecType execType, simResult * results);



#endif