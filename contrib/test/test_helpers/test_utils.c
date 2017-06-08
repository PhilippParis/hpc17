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
#include <limits.h>
#include <errno.h>
#include <mpi.h>

#include "mpi_colls_check.h"
#include "buffer_handling.h"
#include "test_utils.h"


static char* test_res_str[] = {
    [TEST_PASSED] = "Passed",
    [TEST_FAILED] = "FAILED"
};
static const int OUTPUT_RANK = 0;


void test_collective(const char* mpicall_name, basic_collective_params_t* params) {
  void *resbuf_ref, *resbuf_test;
  size_t resbuf_size;
  test_res_codes_t test_result;
  int rank;
  collective_ops_t mpicall_funcs;
  int found = 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (mpicall_name == NULL) {
    if (rank == OUTPUT_RANK) {
      fprintf(stderr, "MPI function is null\n");
    }
    MPI_Finalize();
    exit(0);
  }

  found = find_coll_functions(mpicall_name, &mpicall_funcs);
  if (found == 0) {
    if (rank == OUTPUT_RANK) {
      fprintf(stderr, "MPI function \"%s\" is invalid\n", mpicall_name);
    }
    MPI_Finalize();
    exit(0);
  }
  mpicall_funcs.mpicall_init(params->sendcount, params->sendtype, params->recvcount, params->recvtype, params->root,
      params->comm);

  if (strcmp(mpicall_name, "MY_Gather") == 0) {
    // test an implementation of MPI_Gather
    test_Gather(mpicall_funcs.mpicall_func, params, &resbuf_ref, &resbuf_test, &resbuf_size);

  } else if (strcmp(mpicall_name, "MY_Scatter") == 0) {
    // test an implementation of MPI_Scatter
    test_Scatter(mpicall_funcs.mpicall_func, params, &resbuf_ref, &resbuf_test, &resbuf_size);
  } else {
    if (rank == OUTPUT_RANK) {
      fprintf(stderr, "MPI function \"%s\" is unknown\n", mpicall_name);
    }
    MPI_Finalize();
    exit(0);
  }

  test_result = check_results(params, resbuf_ref, resbuf_test, resbuf_size, mpicall_funcs.check_type);
  print_test_result(mpicall_name, params, test_result);
  mpicall_funcs.mpicall_cleanup(params->sendcount, params->sendtype, params->recvcount, params->recvtype, params->root,
      params->comm);

  // cleanup result buffers
  free(resbuf_ref);
  free(resbuf_test);

}



test_res_codes_t check_results(basic_collective_params_t* params, void *resbuf_ref, void *resbuf_test,
    size_t resbuf_size, test_res_check_t check_type) {

  int rank, size;
  test_res_codes_t local_test_result = TEST_FAILED;
  test_res_codes_t test_result;
  int *test_res_all_procs;
  int i;

  MPI_Comm_rank(params->comm, &rank);
  MPI_Comm_size(params->comm, &size);

  test_res_all_procs = (int*) calloc(size, sizeof(int));

  local_test_result = TEST_PASSED;
  if (check_type == CHECK_RES_ROOT_ONLY) {
    if (rank == params->root) { // only check receive buffers on the test's root process
      if (!identical(resbuf_ref, resbuf_test, resbuf_size)) {
        local_test_result = TEST_FAILED;
      }
    }
  } else {  // each process has to check its local receive buffers
    if (!identical(resbuf_ref, resbuf_test, resbuf_size)) {
      local_test_result = TEST_FAILED;
    }
  }

  // make sure all processes return the same test result
  MPI_Allgather(&local_test_result, 1, MPI_INT, test_res_all_procs, 1, MPI_INT, params->comm);

  test_result = TEST_PASSED;
  for (i = 0; i < size; i++) {
    if (test_res_all_procs[i] == TEST_FAILED) {
      test_result = TEST_FAILED;
      break;
    }
  }

  free(test_res_all_procs);
  return test_result;
}

void print_test_result(const char* mpicall_name, basic_collective_params_t* params, test_res_codes_t test_res) {
  int rank;
  char tmps[50];
  MPI_Comm_rank(params->comm, &rank);

  if (rank == OUTPUT_RANK) {
    sprintf(tmps, "%s (count=%d)", mpicall_name, params->sendcount);
    printf("%-50s ... %-10s\n", tmps, test_res_str[test_res]);
  }
}


int str_to_positive_int(const char *str, int* result) {
  char *endptr;
  int error = 0;
  long res;

  errno = 0;
  res = strtol(str, &endptr, 10);

  /* Check for various possible errors */
  if ((errno == ERANGE && (res == LONG_MAX || res == LONG_MIN)) || (errno != 0 && res == 0)) {
    error = 1;
  }
  if (endptr == str) {  // no digits parsed
    error = 1;
  }
  if (endptr != str + strlen(str)) {  // not all characters were parsed
    error = 1;
  }
  if (res < 0 || res > INT_MAX) {  // obtained number is not a positive int
    error = 1;
  }

  if (!error) {
    *result = (int)res;
  }
  else {
    *result = 0;
  }

  return error;
}



