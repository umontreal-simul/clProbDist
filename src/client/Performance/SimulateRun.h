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

#ifndef SIMULATERUN_SS_H
#define SIMULATERUN_SS_H

#include <clRNG/mrg32k3a.h>
#include <clProbDist/poisson.h>
#include <clProbDist/normal.h>

#include "../common.h"

//Device Simulation
void performanceOnGPU_Poisson(cl_context context, cl_device_id device, cl_command_queue queue, const char * kernelName,
	int n, int n1, double lambda, size_t streamsBufSize, clrngMrg32k3aStream* streams, ExecType execType);

void performanceOnGPU_NoramlDist(cl_context context, cl_device_id device, cl_command_queue queue, const char * kernelName,
	int n, int n1, size_t streamsBufSize, clrngMrg32k3aStream* streams, size_t normalBufSize, clprobdistNormal* normalDist);

void performanceOnGPU_PoissonGlobalLocal(cl_context context, cl_device_id device, cl_command_queue queue, const char * kernelName,
	int n, int n1, double lambda, size_t streamsBufSize, clrngMrg32k3aStream* streams, ExecType execType);

//CPU Simulation
void performanceOnCPU_Poisson(int n, double lambda, clrngMrg32k3aStream* stream, ExecType execType);

#endif