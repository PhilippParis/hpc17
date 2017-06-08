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

#ifndef REPROMPI_MPI_COLLS_VERIFICATION_H_
#define REPROMPI_MPI_COLLS_VERIFICATION_H_

typedef enum test_res_check {
  CHECK_RES_ALL_PROCS,
  CHECK_RES_ROOT_ONLY
} test_res_check_t;

typedef struct basic_collparams {
  int root;
  MPI_Datatype sendtype, recvtype;
  int sendcount, recvcount;
  MPI_Comm comm;
} basic_collective_params_t;


typedef int (*mpicall_implem_t)(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount,
    MPI_Datatype recvtype, int root, MPI_Comm comm);
typedef void (*mpicall_aux_func_t)(int sendcount, MPI_Datatype sendtype, int recvcount, MPI_Datatype recvtype, int root,
    MPI_Comm comm);


typedef struct collective_ops {
  char* name;
  mpicall_implem_t mpicall_func;
  mpicall_aux_func_t mpicall_init;
  mpicall_aux_func_t mpicall_cleanup;
  test_res_check_t check_type;
} collective_ops_t;

int find_coll_functions(const char* name, collective_ops_t* coll_ops);

void test_Scatter(const mpicall_implem_t mpicall_implem, basic_collective_params_t* params, void **resbuf_ref,
    void **resbuf_test, size_t *resbuf_size);
void test_Gather(const mpicall_implem_t mpicall_implem, basic_collective_params_t* params, void **resbuf_ref,
    void **resbuf_test, size_t *resbuf_size);

#endif /* REPROMPI_MPI_COLLS_VERIFICATION_H_ */



