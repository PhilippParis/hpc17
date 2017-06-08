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
#include <mpi.h>

#include "test_helpers/test_utils.h"

static const int OUTPUT_RANK = 0;

int main(int argc, char* argv[]) {
  int count = 0;
  int error;
  char *mpicall_name = NULL;
  int rank;
  basic_collective_params_t params;

  srand(1000);

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc < 3) {
    if (rank == OUTPUT_RANK) {
      printf("\nUSAGE: mpirun -np 2 %s mpi_collective count\n\n", argv[0]);
      printf("%-20s - %s\n", "mpi_collective", "MY_Gather or MY_Scatter");
      printf("%-20s - %s\n", "count", "positive integer specifying the send count of the MPI collective");
    }
    MPI_Finalize();
    exit(0);
  }

  mpicall_name = strdup(argv[1]);
  if (mpicall_name == NULL) {
    if (rank == OUTPUT_RANK) {
      fprintf(stderr, "MPI function is null\n");
    }
    MPI_Finalize();
    exit(0);
  }

  error = str_to_positive_int(argv[2], &count);
  if (error) {
    if (rank == OUTPUT_RANK) {
      printf("Invalid count specified: %s \n", argv[2]);
    }
    MPI_Finalize();
    exit(0);
  }

  params.root = 0;
  params.sendtype = MPI_INT;
  params.recvtype = MPI_INT;
  params.comm = MPI_COMM_WORLD;
  params.sendcount = count;
  params.recvcount = count;

  test_collective(mpicall_name, &params);

  free(mpicall_name);

  /* shut down MPI */
  MPI_Finalize();

  return 0;
}
