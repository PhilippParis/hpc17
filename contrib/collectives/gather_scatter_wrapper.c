/*  ReproMPI Benchmark
 *
 *  Copyright 2015 Alexandra Carpen-Amarie, Sascha Hunold
 Research Group for Parallel Computing
 Faculty of Informatics
 Vienna University of Technology, Austria

 <license>
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 </license>
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include "mpi.h"
#include "collective_ops/collectives.h"
#include "gather_scatter_implem.h"

/***************************************/
// MPI_Gather
inline void execute_MY_Gather(collective_params_t* params) {
  MY_Gather(params->sbuf, params->count, params->datatype, params->rbuf, params->count, params->datatype,
      params->root, MPI_COMM_WORLD);
}


void initialize_data_MY_Gather(const basic_collective_params_t info, const long count, collective_params_t* params) {
  // allocate the send and receive buffers for Gather
  initialize_data_Gather(info, count, params);

  // additional initializations
  init_MY_Gather(params->count, params->datatype, params->count, params->datatype, params->root, MPI_COMM_WORLD);
}

void cleanup_data_MY_Gather(collective_params_t* params) {

  cleanup_MY_Gather(params->count, params->datatype, params->count, params->datatype, params->root, MPI_COMM_WORLD);
  // free send and receive buffers
  cleanup_data_Gather(params);
}

/***************************************/

/***************************************/
// Scatter
inline void execute_MY_Scatter(collective_params_t* params) {
  MY_Scatter(params->sbuf, params->count, params->datatype, params->rbuf, params->count, params->datatype,
      params->root, MPI_COMM_WORLD);
}


void initialize_data_MY_Scatter(const basic_collective_params_t info, const long count, collective_params_t* params) {
  // allocate the send and receive buffers for Scatter
  initialize_data_Scatter(info, count, params);

  // additional initializations
  init_MY_Scatter(params->count, params->datatype, params->count, params->datatype, params->root, MPI_COMM_WORLD);
}


void cleanup_data_MY_Scatter(collective_params_t* params) {

  cleanup_MY_Scatter(params->count, params->datatype, params->count, params->datatype, params->root, MPI_COMM_WORLD);
  cleanup_data_Scatter(params);
}

/***************************************/


