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

#include "buffer_handling.h"


void set_buffer_random(void *buff, int count, MPI_Datatype datatype) {
  int i;
  MPI_Aint type_extent, lb;
  MPI_Type_get_extent(datatype, &lb, &type_extent);

  for (i = 0; i < count * type_extent; i++) {
    ((char*) buff)[i] = rand() % 100;
  }
}


int identical(void* buffer1, void* buffer2, size_t size) {
  size_t i;
  int identical_bufs = 1;

  for (i = 0; i < size; i++) {
    if (((char*) buffer1)[i] != ((char*) buffer2)[i]) {
      identical_bufs = 0;
      break;
    }
  }
  return identical_bufs;
}


