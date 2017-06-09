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
#include <stdbool.h>
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

void scatter_divide_and_conquer(char *buffer,
                                const int sendcount, const int recvcount,
                                const MPI_Datatype sendtype, const MPI_Datatype recvtype,
                                int start, int end, int root, MPI_Comm comm,
                                const int rank, const MPI_Aint size_per_element)
{
    const int n = (end - start) / 2;
    const int m = start + n;

    if (n == 0) {
        return;
    }

    int subroot = 0;
    int blocks = 0;
    int newroot = root;
    char* sendbuf = buffer;
    char* recvbuf = buffer;
    
    if (root < m) {
        subroot = m;
        blocks = end - start - n;
        if (rank < m) {
            if (rank == root) {
                sendbuf += m * size_per_element * sendcount;
            }
            end = m;
        } else {
            if (rank == subroot) {
                recvbuf += m * size_per_element * sendcount;
            }
            start = m;
            newroot = subroot;
        }
    } else {
        subroot = start; 
        blocks = n;
        if (rank >= m) {
            if (rank == root) {
                sendbuf += start * size_per_element * sendcount;
            }
            start = m;
        } else {
            if (rank == subroot) {
                recvbuf += start * size_per_element * sendcount;
            }
            end = m;
            newroot = subroot;
        }
    }

    // send data to subroot
    if (rank == root) {
        MPI_Send(sendbuf, blocks * sendcount, sendtype, subroot, 0, comm);
    } else if (rank == subroot) {
        MPI_Recv(recvbuf, blocks * recvcount, recvtype, root, 0, comm, MPI_STATUS_IGNORE);
    }

    scatter_divide_and_conquer(buffer, sendcount, recvcount, sendtype, recvtype,
                               start, end, newroot, comm, rank, size_per_element);
}

//#define SCATTER_DEBUG

// Scatter
int MY_Scatter(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
               void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);

    int rank;
    MPI_Comm_rank(comm, &rank);

    MPI_Aint lb;
    MPI_Aint size_per_element;
    MPI_Type_get_true_extent(recvtype, &lb, &size_per_element);

    char* buffer;
    if (rank != root) {
        buffer = (char*)malloc(sendcount * size_per_element * size);
    } else {
        buffer = (char*)sendbuf;
    }

    // =======================================
#ifdef SCATTER_DEBUG
    int* tmp = (int*) sendbuf;

    if (rank == root) {
        for (int p = 0; p < size; ++p) {
            printf("%i:", p);
            for (int i = 0; i < sendcount; ++i) {
                printf("%i, ", *tmp);
                tmp++;
            }
            printf("\n");
        }
    }
#endif
    // ==============================

    scatter_divide_and_conquer(buffer, sendcount, recvcount, sendtype, recvtype,
                               0, size, root, comm, rank, size_per_element);
    memcpy(recvbuf, buffer + sendcount * size_per_element * rank, sendcount * size_per_element);

#ifdef SCATTER_DEBUG
    printf("recv: %i: %i\n", rank, *((int*) recvbuf));
    fflush(stdout);
#endif

    if (rank != root) {
        free(buffer);
    }

    return MPI_SUCCESS;
}

void init_MY_Scatter(int sendcount, MPI_Datatype sendtype, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) {
}


void cleanup_MY_Scatter(int sendcount, MPI_Datatype sendtype, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) {
}


/***************************************/

