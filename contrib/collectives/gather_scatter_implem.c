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

static int min(int a, int b)
{
    return (a < b) ? a : b;
}

static int to_virtual_rank(int rank, int root, int size)
{
    return (rank - root + size) % size;
}

static int to_real_rank(int vrank, int root, int size)
{
    return (vrank + root) % size;
}

static int get_parent_vrank(int vrank, int size)
{
    int d = 1;
    while((vrank & d) != d && (d < size)) {
        d <<= 1;
    }
    return vrank & ~d;
}

/***************************************/
// Gather

static int binominal_tree_gather(const char* sendbuf, const int sendcount, const MPI_Datatype sendtype,
                                 char* recvbuf, const int recvcount, const MPI_Datatype recvtype,
                                 const int root, const MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);

    int rank;
    MPI_Comm_rank(comm, &rank);
    const int vrank = to_virtual_rank(rank, root, size);

    MPI_Aint send_lb;
    MPI_Aint send_size_per_element;
    MPI_Type_get_extent(sendtype, &send_lb, &send_size_per_element);

    MPI_Aint recv_lb;
    MPI_Aint recv_size_per_element;
    MPI_Type_get_extent(recvtype, &recv_lb, &recv_size_per_element);

    if (((rank == root) && (recvcount == 0)) || ((rank != root) && (sendcount == 0))) {
        // nothing to do
        return MPI_SUCCESS;
    }

    if (size == 1) {
        // root local gather
        return MPI_Sendrecv(sendbuf, sendcount, sendtype, rank, 0,
                            recvbuf, recvcount, recvtype, rank, 0,
                            comm, MPI_STATUS_IGNORE);
    }

    char* tmpbuffer = NULL;

    int d = 1;
    while (((vrank & d)) != d && (d < size)) {
        if ((vrank + d) < size) {
            if (tmpbuffer == NULL) {
                tmpbuffer = (char*)malloc(sendcount * send_size_per_element * size);

                // copy local data into temp buffer because send will use this temp buffer
                memcpy(tmpbuffer + sendcount * send_size_per_element * vrank,
                       sendbuf, sendcount * send_size_per_element);
            }

            const int src_vrank = vrank + d;
            const int blocks = min(d, size - src_vrank); // each node local data counts 1 block
            MPI_Recv(tmpbuffer + sendcount * send_size_per_element * src_vrank,
                     blocks * sendcount, sendtype, to_real_rank(src_vrank, root, size),
                     0, comm, MPI_STATUS_IGNORE);
        }
        d <<= 1;
    }

    if (rank != root) {
        const int dst_vrank = vrank - d;
        const int blocks = min(d, size - vrank); // each node local data counts 1 block

        if (blocks == 1) {
            // leaf node
            assert(tmpbuffer == NULL);
            MPI_Send(sendbuf, sendcount, sendtype, to_real_rank(dst_vrank, root, size),
                     0, comm);
        } else {
            assert(tmpbuffer != NULL);
            MPI_Send(tmpbuffer + sendcount * send_size_per_element * vrank,
                     blocks * sendcount, sendtype, to_real_rank(dst_vrank, root, size),
                     0, comm);
            free(tmpbuffer);
        }
    } else {
        // root node
        assert(tmpbuffer != NULL);
        MPI_Sendrecv(tmpbuffer, sendcount * (size - rank), sendtype, rank, 0,
                    recvbuf + recv_size_per_element * recvcount * rank,
                    recvcount * (size - rank), recvtype, rank, 0,
                    comm, MPI_STATUS_IGNORE);
        MPI_Sendrecv(tmpbuffer + send_size_per_element * sendcount * (size - rank),
                    sendcount * rank, sendtype, rank, 0,
                    recvbuf, recvcount * rank, recvtype, rank, 0,
                    comm, MPI_STATUS_IGNORE);
        free(tmpbuffer);
    }

    return MPI_SUCCESS;
}

int MY_Gather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount,
    MPI_Datatype recvtype, int root, MPI_Comm comm) {
  return binominal_tree_gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
}


void init_MY_Gather(int sendcount, MPI_Datatype sendtype, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) {
}


void cleanup_MY_Gather(int sendcount, MPI_Datatype sendtype, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) {
}

/***************************************/

/***************************************/

void scatter_divide_and_conquer(char **buffer, unsigned long buffer_offset,
                                const int sendcount, const MPI_Datatype sendtype,
                                void* recvbuf, const int recvcount, const MPI_Datatype recvtype,
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
        if (blocks == 1) {
            // no forwarding required, receive directly into recvbuf -> done
            MPI_Recv(recvbuf, recvcount, recvtype, root, 0, comm, MPI_STATUS_IGNORE);
            return;
        } else {
            // forwarding required, receive into temp buffer
            *buffer = (char*)malloc(sendcount * size_per_element * blocks);
            MPI_Recv(*buffer, blocks * sendcount, sendtype, root, 0, comm, MPI_STATUS_IGNORE);
        }
    }

    scatter_divide_and_conquer(buffer, buffer_offset, sendcount, sendtype,
                               recvbuf, recvcount, recvtype, start, end, newroot,
                               comm, rank, size_per_element);
}


static void create_wrapping_datatype(int* blocklens, int* displacements, int sendcount, int vrank, int size, int root)
{
    if (to_real_rank(vrank, root, size) >= size) {
        return;
    }
    
    *blocklens = sendcount;
    *displacements = to_real_rank(vrank, root, size) * sendcount;    
    
    int d = 1;
    while((vrank & d) != d && ((vrank | d) < size)) {
        create_wrapping_datatype(++blocklens, ++displacements, sendcount, vrank | d, size, root);
        d <<= 1;
    }
}



static int binominal_tree_scatter(const char* sendbuf, const int sendcount, const MPI_Datatype sendtype,
                                 char* recvbuf, const int recvcount, const MPI_Datatype recvtype,
                                 const int root, const MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);

    int rank;
    MPI_Comm_rank(comm, &rank);
    const int vrank = to_virtual_rank(rank, root, size);

    MPI_Aint send_lb;
    MPI_Aint send_size_per_element;
    MPI_Type_get_extent(sendtype, &send_lb, &send_size_per_element);

    MPI_Aint recv_lb;
    MPI_Aint recv_size_per_element;
    MPI_Type_get_extent(recvtype, &recv_lb, &recv_size_per_element);

    if (((rank == root) && (recvcount == 0)) || ((rank != root) && (sendcount == 0))) {
        // nothing to do
        return MPI_SUCCESS;
    }

    if (size == 1) {
        // root local scatter
        return MPI_Sendrecv(sendbuf, sendcount, sendtype, rank, 0,
                            recvbuf, recvcount, recvtype, rank, 0,
                            comm, MPI_STATUS_IGNORE);
    }
    
    if (rank == root) {
        
        for (int i = 0; i < size; ++i) {
            //printf("expected at %i: %i\n", i, ((int*)sendbuf)[i]);
        }
        
        // root sends to children
        int d = 1;
        while (d < size) {
            const int blocks = min(d, size - d);
            const int real_recv = to_real_rank(d, root, size);
            
            MPI_Datatype type;
            int blocklens[blocks];
            int displacements[blocks];
            create_wrapping_datatype(blocklens, displacements, sendcount, d, size, root);
            MPI_Type_indexed(blocks, blocklens, displacements, sendtype, &type);
            MPI_Type_commit(&type);
            
            /*
            printf("blocklens: ");
            for (int i=0; i < blocks; ++i) {
                printf("%i, ", blocklens[i]);
            }
            printf("\n");
            printf("displacements: ");
            for (int i=0; i < blocks; ++i) {
                printf("%i, ", displacements[i]);
            }
            printf("\n");
            */
            
            //printf("send %i blocks from %i to %i: %i\n", blocks, rank, real_recv, ((int*)(sendbuf + real_recv * send_size_per_element * sendcount))[0]);
            MPI_Send(sendbuf, 1, type, real_recv ,0, comm);
            MPI_Type_free(&type);
            
            d <<= 1;
        }
        
        // send block to itself
        *recvbuf = (char*)malloc(recv_size_per_element * recvcount);
        MPI_Sendrecv(sendbuf + sendcount * send_size_per_element * rank,
                     sendcount, sendtype, rank, 0, recvbuf, recvcount, recvtype,
                     rank, 0, comm, MPI_STATUS_IGNORE);
        
    } else {
        // receive from parent
        const int parent_vrank = get_parent_vrank(vrank, size);
        const int blocks = min(vrank - parent_vrank, size - vrank);
        const int real_sender = to_real_rank(parent_vrank, root, size);
        //printf("recv %i blocks at %i from %i \n", blocks, rank, real_sender);
        *recvbuf = (char*)malloc(blocks * recv_size_per_element * recvcount);
        MPI_Recv(recvbuf, blocks * recvcount, recvtype, real_sender, 0, comm, MPI_STATUS_IGNORE);

        // forward to children
        if (blocks > 1) {
            int d = 1;
            while((vrank & d) != d && (d < size)) {
                const int vrecv = vrank | d;
                const int shifted_vrecv = vrecv - vrank;
                const int real_recv = to_real_rank(vrecv, root, size);
                const int blocks = min(shifted_vrecv, size - vrecv);
                
                //printf("forward %i blocks at %i to %i: %i\n", blocks, rank, real_recv, ((int*)(recvbuf + shifted_vrecv * send_size_per_element * sendcount))[0]);
                
                MPI_Send(recvbuf + shifted_vrecv * send_size_per_element * sendcount, blocks * sendcount, sendtype, real_recv ,0, comm);
                d <<= 1;
            }
        }
    }
    printf("%i: %i\n", rank, *((int*)recvbuf));
    return MPI_SUCCESS;
}


// Scatter
int MY_Scatter(const void* sendbuf, const int sendcount, const MPI_Datatype sendtype,
               void* recvbuf, const int recvcount, const MPI_Datatype recvtype,
               const int root, const MPI_Comm comm)
{
    return binominal_tree_scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
    
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

    if (recvbuf != MPI_IN_PLACE) {
        memset(recvbuf, 0, recvcount * recv_size_per_element);
    }

    scatter_divide_and_conquer(&buffer, 0, sendcount, sendtype,
                               recvbuf, recvcount, recvtype, 0, size, root,
                               comm, rank, send_size_per_element);

    if (rank == root) {
        if (recvbuf != MPI_IN_PLACE) {
            MPI_Sendrecv(buffer + sendcount * send_size_per_element * rank,
                         sendcount, sendtype, rank, 0, recvbuf, recvcount, recvtype,
                         rank, 0, comm, MPI_STATUS_IGNORE);
        }
    } else if (buffer != NULL) {
        // used a temp buffer for forwarding -> get the relevant data from the
        // temp buffer and free it
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

