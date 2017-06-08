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
#include "gather_scatter_implem.h"

/***************************************/
// Gather

int MY_Gather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount,
    MPI_Datatype recvtype, int root, MPI_Comm comm) {
  MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
  return MPI_SUCCESS;
}


void init_MY_Gather(int sendcount, MPI_Datatype sendtype, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) {
}


void cleanup_MY_Gather(int sendcount, MPI_Datatype sendtype, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) {
}

/***************************************/

/***************************************/
// Scatter
int MY_Scatter(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount,
    MPI_Datatype recvtype, int root, MPI_Comm comm) {

  MPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
  return MPI_SUCCESS;
}

void init_MY_Scatter(int sendcount, MPI_Datatype sendtype, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) {
}


void cleanup_MY_Scatter(int sendcount, MPI_Datatype sendtype, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) {
}

void scatter_divide_and_conquer(void *sendbuf, const int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount,
    MPI_Datatype recvtype, int start, int end, int root, MPI_Comm comm) 
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    if (start + 1 == end) {
        if (rank == root) {
            int size_per_element = 0;
            MPI_Type_get_true_extent(recvtype, NULL, &size_per_element);
            memcpy(recvbuf, sendbuf, size_per_element * sendcount);
        }
    }
    
    int n = (end - start) / 2;
    int m = start + end;
    int subroot = 0;
    int blocks = 0;
    int newroot = root;
    
    if (root < m) {
        subroot = m;
        blocks = end - start - n;
        if (rank < m) {
            if (rank == root) {
                sendbuf += (m - start); // TODO ?
            }
            end = m;
        } else {
            start = m;
            newroot = subroot;
        }
    } else {
        subroot = start; 
        blocks = n;
        if (rank >= m) {
            start = m;
        } else {
            end = m;
            newroot = subroot;
        }
    }
    
    // send data to subroot
    if (rank == root) {
        MPI_Send(sendbuf, blocks, sendtype, subroot, MPI_ANY_TAG, comm);
    } else if (rank == subroot) {
        int size_per_element = 0;
        MPI_Type_get_true_extent(recvtype, NULL, &size_per_element);
        sendbuf = (void*) malloc(blocks * size_per_element * sendcount);
        MPI_Recv(sendbuf, blocks, recvtype, root, MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);
    }
        
    scatter_divide_and_conquer(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, start, end, newroot, comm);
    if (rank != root && sendbuf != NULL) {
        free(sendbuf);
    }
}

/***************************************/

