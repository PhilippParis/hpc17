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

static void test_collective_all_roots(char* name, basic_collective_params_t* params, int size) {
  for (int root = 0; root < size; root++) {
    params->root = root;
    test_collective(name, params);
  }
}



int main(int argc, char* argv[]) {
  int count = 0;
  int error;
  char *mpicall_name = NULL;
  int rank;
  int size;

  srand(1000);

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

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

  // simple type
  {
    if (rank == OUTPUT_RANK) {
      printf(">>>> Simple Type:\n");
    }

    basic_collective_params_t params = {
      .sendtype = MPI_INT,
      .recvtype = MPI_INT,
      .comm = MPI_COMM_WORLD,
      .sendcount = count,
      .recvcount = count
    };

    test_collective_all_roots(mpicall_name, &params, size);
  }

  // contiguous type
  {
    if (rank == OUTPUT_RANK) {
      printf(">>>> Contiguous Type:\n");
    }

    MPI_Datatype type;
    MPI_Type_contiguous(3, MPI_INT, &type);
    MPI_Type_commit(&type);

    basic_collective_params_t params = {
      .sendtype = type,
      .recvtype = type,
      .comm = MPI_COMM_WORLD,
      .sendcount = count,
      .recvcount = count
    };

    test_collective_all_roots(mpicall_name, &params, size);

    MPI_Type_free(&type);
  }

  // contiguous type
  {
      if (rank == OUTPUT_RANK) {
          printf(">>>> Contiguous Type (different count):\n");
      }

      MPI_Datatype send_type;
      MPI_Type_contiguous(4, MPI_INT, &send_type);
      MPI_Type_commit(&send_type);

      MPI_Datatype recv_type;
      MPI_Type_contiguous(2, MPI_INT, &recv_type);
      MPI_Type_commit(&recv_type);

      basic_collective_params_t params = {
          .sendtype = send_type,
          .recvtype = recv_type,
          .comm = MPI_COMM_WORLD,
          .sendcount = count,
          .recvcount = 2 * count
      };

      test_collective_all_roots(mpicall_name, &params, size);

      MPI_Type_free(&send_type);
      MPI_Type_free(&recv_type);
  }

  // vector type
  {
    if (rank == OUTPUT_RANK) {
      printf(">>>> Vector Type:\n");
    }

    MPI_Datatype type;
    MPI_Type_vector(2, 5, 2, MPI_INT, &type);
    MPI_Type_commit(&type);

    basic_collective_params_t params = {
      .sendtype = type,
      .recvtype = type,
      .comm = MPI_COMM_WORLD,
      .sendcount = count,
      .recvcount = count
    };

    test_collective_all_roots(mpicall_name, &params, size);

    MPI_Type_free(&type);
  }

  // vector type
  {
      if (rank == OUTPUT_RANK) {
          printf(">>>> Vector Type (different send/recv type):\n");
      }

      MPI_Datatype send_type;
      MPI_Type_vector(2, 5, 2, MPI_INT, &send_type);
      MPI_Type_commit(&send_type);

      MPI_Datatype recv_type;
      MPI_Type_vector(5, 2, 3, MPI_INT, &recv_type);
      MPI_Type_commit(&recv_type);

      basic_collective_params_t params = {
          .sendtype = send_type,
          .recvtype = recv_type,
          .comm = MPI_COMM_WORLD,
          .sendcount = count,
          .recvcount = count
      };

      test_collective_all_roots(mpicall_name, &params, size);

      MPI_Type_free(&send_type);
      MPI_Type_free(&recv_type);
  }

  // hvector type
  {
    if (rank == OUTPUT_RANK) {
      printf(">>>> HVector Type:\n");
    }

    MPI_Datatype type;
    MPI_Type_create_hvector(2, 5, 2, MPI_INT, &type);
    MPI_Type_commit(&type);

    basic_collective_params_t params = {
      .sendtype = type,
      .recvtype = type,
      .comm = MPI_COMM_WORLD,
      .sendcount = count,
      .recvcount = count
    };

    test_collective_all_roots(mpicall_name, &params, size);

    MPI_Type_free(&type);
  }

  // indexed type
  {
    if (rank == OUTPUT_RANK) {
      printf(">>>> Indexed Type:\n");
    }

    const int blocklens[3] = {2, 3, 1};
    const int displacements[3] = {0, 3, 8};

    MPI_Datatype type;
    MPI_Type_indexed(3, blocklens, displacements, MPI_INT, &type);
    MPI_Type_commit(&type);

    basic_collective_params_t params = {
      .sendtype = type,
      .recvtype = type,
      .comm = MPI_COMM_WORLD,
      .sendcount = count,
      .recvcount = count
    };

    test_collective_all_roots(mpicall_name, &params, size);

    MPI_Type_free(&type);
  }

  // hindexed type
  {
    if (rank == OUTPUT_RANK) {
      printf(">>>> HIndexed Type:\n");
    }

    const int blocklens[3] = {2, 3, 1};
    const MPI_Aint displacements[3] = {0, 3, 8};

    MPI_Datatype type;
    MPI_Type_create_hindexed(3, blocklens, displacements, MPI_INT, &type);
    MPI_Type_commit(&type);

    basic_collective_params_t params = {
      .sendtype = type,
      .recvtype = type,
      .comm = MPI_COMM_WORLD,
      .sendcount = count,
      .recvcount = count
    };

    test_collective_all_roots(mpicall_name, &params, size);

    MPI_Type_free(&type);
  }

  // struct type
  {
    if (rank == OUTPUT_RANK) {
      printf(">>>> Struct Type:\n");
    }

    const int blocklens[2] = {2, 3};
    const MPI_Aint displacements[2] = {0, 64};
    const MPI_Datatype oldtypes[2] = {MPI_INT, MPI_CHAR};

    MPI_Datatype type;
    MPI_Type_create_struct(2, blocklens, displacements, oldtypes, &type);
    MPI_Type_commit(&type);

    basic_collective_params_t params = {
      .sendtype = type,
      .recvtype = type,
      .comm = MPI_COMM_WORLD,
      .sendcount = count,
      .recvcount = count
    };

    test_collective_all_roots(mpicall_name, &params, size);

    MPI_Type_free(&type);
  }

  // subarray type
  {
    if (rank == OUTPUT_RANK) {
      printf(">>>> Subarray Type:\n");
    }

    const int sizes[] = {10, 10};
    const int subsizes[] = {6, 6};
    const int starts[] = {2, 2};

    MPI_Datatype type;
    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &type);
    MPI_Type_commit(&type);

    basic_collective_params_t params = {
      .sendtype = type,
      .recvtype = type,
      .comm = MPI_COMM_WORLD,
      .sendcount = count,
      .recvcount = count
    };

    test_collective_all_roots(mpicall_name, &params, size);

    MPI_Type_free(&type);
  }

  free(mpicall_name);

  /* shut down MPI */
  MPI_Finalize();

  return 0;
}
