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
#include <time.h>
#include "mpi.h"

#include "../../collectives/gather_scatter_implem.h"
#include "buffer_handling.h"
#include "mpi_colls_check.h"


const collective_ops_t collective_calls[] = {
    { "MY_Gather", &MY_Gather, &init_MY_Gather, &cleanup_MY_Gather, CHECK_RES_ROOT_ONLY },
    { "MY_Scatter", &MY_Scatter, &init_MY_Scatter, &cleanup_MY_Scatter, CHECK_RES_ALL_PROCS }
};
static const int N_CALLS = 2;

int find_coll_functions(const char* name, collective_ops_t* coll_ops) {
  int index = -1;
  int i;
  int found = 0;

  for (i = 0; i < N_CALLS; i++) {
    if (strcmp(name, collective_calls[i].name) == 0) {
      index = i;
      break;
    }
  }
  if (index >= 0) {
    *coll_ops = collective_calls[index];
    found = 1;
  }
  return found;
}



void test_Scatter(const mpicall_implem_t mpicall_implem, basic_collective_params_t* params, void **resbuf_ref,
    void **resbuf_test, size_t *resbuf_size) {
  void *sendbuf, *recvbuf1, *recvbuf2;
  MPI_Aint sendtype_extent, recvtype_extent, lb;
  int rank, size;

  MPI_Comm_rank(params->comm, &rank);
  MPI_Comm_size(params->comm, &size);
  MPI_Type_get_extent(params->sendtype, &lb, &sendtype_extent);
  MPI_Type_get_extent(params->recvtype, &lb, &recvtype_extent);

  // count is the number of elements on each process
  sendbuf = NULL;
  if (rank == params->root) { // sendbuf only required at the root
    sendbuf = (char*) calloc(size * params->sendcount * sendtype_extent, sizeof(char));
    set_buffer_random(sendbuf, size * params->sendcount, params->sendtype);
  }

  recvbuf1 = (char*) calloc(params->recvcount * recvtype_extent, sizeof(char));
  recvbuf2 = (char*) calloc(params->recvcount * recvtype_extent, sizeof(char));

  MPI_Scatter(sendbuf, params->sendcount, params->sendtype, recvbuf1, params->recvcount, params->recvtype, params->root,
      params->comm);
  mpicall_implem(sendbuf, params->sendcount, params->sendtype, recvbuf2, params->recvcount, params->recvtype,
      params->root, params->comm);

  if (rank == params->root) {
    free(sendbuf);
  }

  *resbuf_ref = recvbuf1;
  *resbuf_test = recvbuf2;
  *resbuf_size = params->recvcount * recvtype_extent;
}

void test_Gather(const mpicall_implem_t mpicall_implem, basic_collective_params_t* params, void **resbuf_ref,
    void **resbuf_test, size_t *resbuf_size) {
  void *sendbuf, *recvbuf1, *recvbuf2;
  MPI_Aint sendtype_extent, recvtype_extent, lb;
  int rank, size;

  MPI_Comm_rank(params->comm, &rank);
  MPI_Comm_size(params->comm, &size);
  MPI_Type_get_extent(params->sendtype, &lb, &sendtype_extent);
  MPI_Type_get_extent(params->recvtype, &lb, &recvtype_extent);

  // count is the number of elements on each process
  sendbuf = (char*) calloc(params->sendcount * sendtype_extent, sizeof(char));

  recvbuf1 = NULL;
  recvbuf2 = NULL;
  if (rank == params->root) { // recv buffers only needed at the root
    recvbuf1 = (char*) calloc(size * params->recvcount * recvtype_extent, sizeof(char));
    recvbuf2 = (char*) calloc(size * params->recvcount * recvtype_extent, sizeof(char));
  }

  set_buffer_random(sendbuf, params->sendcount, params->sendtype);

  MPI_Gather(sendbuf, params->sendcount, params->sendtype, recvbuf1, params->recvcount, params->recvtype, params->root,
      params->comm);
  mpicall_implem(sendbuf, params->sendcount, params->sendtype, recvbuf2, params->recvcount, params->recvtype,
      params->root, params->comm);

  free(sendbuf);

  *resbuf_ref = recvbuf1;
  *resbuf_test = recvbuf2;
  if (rank == params->root) { // recv buffers only needed at the root
    *resbuf_size = size * params->recvcount * recvtype_extent;
  } else {
    *resbuf_size = 0;
  }
}
