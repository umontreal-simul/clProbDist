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
#ifndef CLPROBDIST_PRIVATE_H
#define CLPROBDIST_PRIVATE_H

/*! @brief Set the current error string
*
*  The error string will be constructed based on the error code \c err and on
*  the optional message \c msg.
*
*  @param[in]  err     Error code.
*  @param[in]  msg     Additional error message (format string).  Can be `NULL`.
*  @param[in]  ...     Additional arguments for the format string.
*  @return     The value of err (for convenience).
*/
clprobdistStatus clprobdistSetErrorString(cl_int err, const char* msg, ...);

void arraycopy(cl_double* src, int index_s, cl_double *dest, int index_d, int size);

#endif
