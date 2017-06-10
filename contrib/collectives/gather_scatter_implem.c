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

void scatter_divide_and_conquer(char **buffer, unsigned long buffer_offset,
                                const int sendcount, const MPI_Datatype sendtype,
                                int start, int end, const int root, const MPI_Comm comm,
                                const int rank, const MPI_Aint size_per_element)
{
    const int n = (end - start) / 2;
    const int m = start + n;

    if (n == 0) {
        return;
    }

    int subroot;
    int blocks;
    int newroot;
    unsigned long send_offset = 0;

    if (root < m) {
        subroot = m;
        blocks = end - start - n;
        if (rank < m) {
            if (rank == root) {
                send_offset = m * size_per_element * sendcount;
            }
            end = m;
            newroot = root;
        } else {
            if (rank == subroot) {
                buffer_offset = m * size_per_element * sendcount;
            }
            start = m;
            newroot = subroot;
        }
    } else {
        subroot = start;
        blocks = n;
        if (rank >= m) {
            if (rank == root) {
                send_offset = start * size_per_element * sendcount;
            }
            start = m;
            newroot = root;
        } else {
            if (rank == subroot) {
                buffer_offset = start * size_per_element * sendcount;
            }
            end = m;
            newroot = subroot;
        }
    }

    // send data to subroot
    if (rank == root) {
        assert(*buffer != NULL);
        MPI_Send(*buffer + send_offset - buffer_offset, blocks * sendcount, sendtype, subroot, 0, comm);
    } else if (rank == subroot) {
        assert(*buffer == NULL);
        *buffer = (char*)malloc(sendcount * size_per_element * blocks);
        MPI_Recv(*buffer, blocks * sendcount, sendtype, root, 0, comm, MPI_STATUS_IGNORE);
    }

    scatter_divide_and_conquer(buffer, buffer_offset, sendcount, sendtype, start, end, newroot,
                               comm, rank, size_per_element);
}

// Scatter
int MY_Scatter(const void* sendbuf, const int sendcount, const MPI_Datatype sendtype,
               void* recvbuf, const int recvcount, const MPI_Datatype recvtype,
               const int root, const MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);

    int rank;
    MPI_Comm_rank(comm, &rank);

    if (((rank == root) && (sendcount == 0)) || ((rank != root) && (recvcount == 0))) {
        // nothing to do
        return MPI_SUCCESS;
    }

    if ((rank != root) && (recvbuf == MPI_IN_PLACE)) {
        // only root can use MPI_IN_PLACE
        return MPI_ERR_BUFFER;
    }

    MPI_Aint send_lb;
    MPI_Aint send_size_per_element;
    MPI_Type_get_extent(sendtype, &send_lb, &send_size_per_element);

    MPI_Aint recv_lb;
    MPI_Aint recv_size_per_element;
    MPI_Type_get_extent(recvtype, &recv_lb, &recv_size_per_element);

    char* buffer = NULL;
    if (rank == root) {
        buffer = (char*)sendbuf;
    }

    scatter_divide_and_conquer(&buffer, 0, sendcount, sendtype, 0, size, root,
                               comm, rank, send_size_per_element);

    if (rank == root) {
        if (recvbuf != MPI_IN_PLACE) {
            memset(recvbuf, 0, recvcount * recv_size_per_element);
            MPI_Sendrecv(buffer + sendcount * send_size_per_element * rank,
                         sendcount, sendtype, rank, 0, recvbuf, recvcount, recvtype,
                         rank, 0, comm, MPI_STATUS_IGNORE);
        }
    } else if (buffer != NULL) {
        memset(recvbuf, 0, recvcount * recv_size_per_element);
        MPI_Sendrecv(buffer, sendcount, sendtype,
                     rank, 0, recvbuf, recvcount, recvtype, rank, 0, comm,
                     MPI_STATUS_IGNORE);
        free(buffer);
    }

    return MPI_SUCCESS;
}

void init_MY_Scatter(int sendcount, MPI_Datatype sendtype, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) {
}


void cleanup_MY_Scatter(int sendcount, MPI_Datatype sendtype, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) {
}


/***************************************/

