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

#ifndef REPROMPI_TEST_UTILS_H_
#define REPROMPI_TEST_UTILS_H_

#include "mpi_colls_check.h"

typedef enum test_res_codes {
  TEST_PASSED = 0,
  TEST_FAILED = 1
} test_res_codes_t;


void test_collective(const char* mpicall_name, basic_collective_params_t* params);
test_res_codes_t check_results(basic_collective_params_t* params, void *resbuf_ref, void *resbuf_test,
    size_t resbuf_size, test_res_check_t check_type);
void print_test_result(const char* mpicall_name, basic_collective_params_t* params, test_res_codes_t test_res);

int str_to_positive_int(const char *str, int* result);

#endif /* REPROMPI_TEST_UTILS_H_ */



