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

#ifndef INVENTORY_TYPES_SS_H
#define INVENTORY_TYPES_SS_H

typedef enum ExecType_ {
	basic = 1,
	Case_a = 2, //Simulate n runs on n workitems using two streams and their substreams
	Case_b = 3, //Simulate n runs on n workitems using n distinct streams
	Case_c = 4, //Simulate n2 runs on n1 workitmes using 2n1 streams and their substreams
	Case_d = 5  //Simulate n2 runs on n1 workitems using 2n2 streams
} ExecType;

typedef enum ExecOption_ {
	OnePolicy = 0,
	Option1 = 1, // Simulate n = n1.n2 runs for p policies, with n1 work items and n2 runs per work item
	Option2 = 2  // Simulate n = n1.n2 runs for p policies, with n1p work items and n2 runs per work item,
} ExecOption;

#endif
