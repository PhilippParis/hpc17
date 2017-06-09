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

void scatter_divide_and_conquer(void *sendbuf, const int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount,
    MPI_Datatype recvtype, int start, int end, int root, MPI_Comm comm) 
{
    
    int rank;
    MPI_Aint size_per_element = 0;
    MPI_Aint lb = 0;
    MPI_Comm_rank(comm, &rank);
    MPI_Type_get_true_extent(recvtype, &lb, &size_per_element);
    
    if (start + 1 == end) {
        memcpy(recvbuf, sendbuf, size_per_element * sendcount);
        return;
    }
    
    int n = (end - start) / 2;
    int m = start + n;
    int subroot = 0;
    int blocks = 0;
    int newroot = root;
    char* shifted_sendbuf = (char*) sendbuf;
    
    if (root < m) {
        subroot = m;
        blocks = end - start - n;
        if (rank < m) {
            if (rank == root) {
                shifted_sendbuf += (m - start) * size_per_element * sendcount;
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
        MPI_Send(shifted_sendbuf, blocks * sendcount, sendtype, subroot, 0, comm);
    } else if (rank == subroot) {
        sendbuf = (void*) malloc(blocks * size_per_element * recvcount);
        MPI_Recv(sendbuf, blocks * recvcount, recvtype, root, 0, comm, MPI_STATUS_IGNORE);
    }
        
    scatter_divide_and_conquer(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, start, end, newroot, comm);
    if (rank != root && sendbuf != NULL) {
        free(sendbuf);
    }
}


// Scatter
int MY_Scatter(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount,
    MPI_Datatype recvtype, int root, MPI_Comm comm) {
    int size = 0;
    MPI_Comm_size(comm, &size);

    scatter_divide_and_conquer(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, 0, size, root, comm);
    return MPI_SUCCESS;
}

void init_MY_Scatter(int sendcount, MPI_Datatype sendtype, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) {
}


void cleanup_MY_Scatter(int sendcount, MPI_Datatype sendtype, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) {
}


/***************************************/

