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

/***************************************/
// binomial tree helper functions

static int binominal_tree_virtual_rank(int rank, int root, int size)
{
    return (rank - root + size) % size;
}

static int binominal_tree_real_rank(int vrank, int root, int size)
{
    return (vrank + root) % size;
}

static int binominal_tree_parent_rank(int vrank, int size)
{
    int d = 1;
    while((vrank & d) != d && (d < size)) {
        d <<= 1;
    }
    return vrank & ~d;
}

/***************************************/
// Gather (Divide-And-Conquer)

static void gather_divide_and_conquer_inner(const char* sendbuf, const int sendcount,
                                            const MPI_Datatype sendtype,
                                            const MPI_Aint send_size_per_element,
                                            char* recvbuf, const int recvcount,
                                            const MPI_Datatype recvtype,
                                            const MPI_Aint recv_size_per_element,
                                            char *tmpbuf, int tmpbuf_offset,
                                            const MPI_Datatype tmptype,
                                            const int tmp_block_size,
                                            int start, int end, const int root,
                                            const MPI_Comm comm, const int rank,
                                            const bool is_gather_root)
{
    const int n = (end - start) / 2;
    const int m = start + n;

    if (n == 0) {
        return;
    }

    int subroot;
    int blocks;
    int newroot;
    int recv_offset = 0;

    if (root < m) {
        subroot = m;
        blocks = end - start - n;
        if (rank < m) {
            if (rank == root) {
                recv_offset = m;
            }
            end = m;
            newroot = root;
        } else {
            if (rank == subroot) {
                tmpbuf_offset = m;
            }
            start = m;
            newroot = subroot;
        }
    } else {
        subroot = start;
        blocks = n;
        if (rank >= m) {
            if (rank == root) {
                recv_offset = start;
            }
            start = m;
            newroot = root;
        } else {
            if (rank == subroot) {
                tmpbuf_offset = start;
            }
            end = m;
            newroot = subroot;
        }
    }

    // subroot nodes which are responsible for more than 1 block are "forwarding
    // nodes" and thus require a temporary buffer.
    bool free_tmpbuf_on_return = false;
    if ((rank == subroot) && (blocks > 1) && (tmpbuf == NULL)) {
        assert(!is_gather_root);
        tmpbuf = (char*)malloc(blocks * tmp_block_size);
        // copy sendbuf data to the beginning of the tmpbuf (because of the
        // tmpbuf_offset the subroot local data is always at the beginning of
        // the tmpbuf)
        int position = 0;
        MPI_Pack(sendbuf, sendcount, sendtype, tmpbuf, tmp_block_size, &position, comm);
        free_tmpbuf_on_return = true;
    }

    gather_divide_and_conquer_inner(sendbuf, sendcount, sendtype, send_size_per_element,
                                    recvbuf, recvcount, recvtype, recv_size_per_element,
                                    tmpbuf, tmpbuf_offset, tmptype, tmp_block_size,
                                    start, end, newroot, comm, rank, is_gather_root);

    // get data from subroot
    if (rank == root) {
        if (is_gather_root) {
            // no forwarding required -> receive directly into recvbuf
            assert(recvbuf != NULL);
            assert(tmpbuf == NULL);
            MPI_Recv(recvbuf + recv_offset * recv_size_per_element * recvcount,
                     blocks * recvcount, recvtype, subroot, 0, comm, MPI_STATUS_IGNORE);
        } else {
            // forwarding required -> receive blocks into tmpbuf
            assert(tmpbuf != NULL);
            MPI_Recv(tmpbuf + (recv_offset - tmpbuf_offset) * tmp_block_size,
                     blocks * tmp_block_size, tmptype, subroot, 0, comm, MPI_STATUS_IGNORE);
        }
    } else if (rank == subroot) {
        if (blocks == 1) {
            // no forwarding required -> send sendbuf
            assert(sendbuf != NULL);
            assert(tmpbuf == NULL);
            MPI_Send(sendbuf, sendcount, sendtype, root, 0, comm);
        } else {
            // forwarding required -> send blocks from tmpbuf
            assert(tmpbuf != NULL);
            MPI_Send(tmpbuf, blocks * tmp_block_size, tmptype, root, 0, comm);
        }
    }

    if (free_tmpbuf_on_return) {
        free(tmpbuf);
    }
}

static int gather_divide_and_conquer(const void* sendbuf,
                                     const int sendcount, const MPI_Datatype sendtype,
                                     void* recvbuf,
                                     const int recvcount, const MPI_Datatype recvtype,
                                     const int root, const MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);

    int rank;
    MPI_Comm_rank(comm, &rank);

    if (((rank == root) && (recvcount == 0)) || ((rank != root) && (sendcount == 0))) {
        // nothing to do
        return MPI_SUCCESS;
    }

    if ((rank != root) && (sendbuf == MPI_IN_PLACE)) {
        // only root can use MPI_IN_PLACE
        return MPI_ERR_BUFFER;
    }

    MPI_Aint send_lb;
    MPI_Aint send_size_per_element;
    MPI_Type_get_extent(sendtype, &send_lb, &send_size_per_element);

    // one block contains all data of one node
    int tmp_block_size;
    MPI_Pack_size(sendcount, sendtype, comm, &tmp_block_size);

    MPI_Aint recv_lb;
    MPI_Aint recv_size_per_element;
    MPI_Type_get_extent(recvtype, &recv_lb, &recv_size_per_element);

    if ((rank == root) && (sendbuf != MPI_IN_PLACE)) {
        // copy root sendbuf into recvbuf
        memset(recvbuf, 0, recvcount * recv_size_per_element * size);
        MPI_Sendrecv(sendbuf, sendcount, sendtype, rank, 0,
                     recvbuf + recvcount * recv_size_per_element * rank,
                     recvcount, recvtype, rank, 0, comm, MPI_STATUS_IGNORE);
    }

    gather_divide_and_conquer_inner(sendbuf, sendcount, sendtype, send_size_per_element,
                                    recvbuf, recvcount, recvtype, recv_size_per_element,
                                    NULL, 0, MPI_BYTE, tmp_block_size,
                                    0, size, root, comm, rank, rank == root);

    return MPI_SUCCESS;
}

/***************************************/
// Gather (Binominal Tree)

static int gather_binominal_tree(const char* sendbuf, const int sendcount,
                                 const MPI_Datatype sendtype,
                                 char* recvbuf, const int recvcount,
                                 const MPI_Datatype recvtype,
                                 const int root, const MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);

    int rank;
    MPI_Comm_rank(comm, &rank);
    const int vrank = binominal_tree_virtual_rank(rank, root, size);

    if (((vrank == 0) && (recvcount == 0)) || ((vrank != 0) && (sendcount == 0))) {
        // nothing to do
        return MPI_SUCCESS;
    }

    if ((vrank != 0) && (sendbuf == MPI_IN_PLACE)) {
        // only root can use MPI_IN_PLACE
        return MPI_ERR_BUFFER;
    }

    MPI_Aint send_lb;
    MPI_Aint send_size_per_element;
    MPI_Type_get_extent(sendtype, &send_lb, &send_size_per_element);

    MPI_Aint recv_lb;
    MPI_Aint recv_size_per_element;
    MPI_Type_get_extent(recvtype, &recv_lb, &recv_size_per_element);

    if ((vrank == 0) && (sendbuf != MPI_IN_PLACE)) {
        // copy root sendbuf into recvbuf
        memset(recvbuf, 0, recvcount * recv_size_per_element * size);
        MPI_Sendrecv(sendbuf, sendcount, sendtype, rank, 0,
                     recvbuf + recvcount * recv_size_per_element * rank,
                     recvcount, recvtype, rank, 0, comm, MPI_STATUS_IGNORE);
    }

    if (size == 1) {
        // nothing more to do
        return MPI_SUCCESS;
    }

    // Use direct receive when root is 0, because in this case no reordering will
    // be required.
    const bool use_direct_recv = (root == 0);

    char* tmpbuf = NULL;

    int d = 1;
    while (((vrank & d)) != d && (d < size)) {
        if ((vrank + d) < size) {
            const int src_vrank = vrank + d;
            const int blocks = min(d, size - src_vrank); // each node local data counts 1 block

            if ((vrank == 0) && use_direct_recv) {
                // receive directly into recvbuf, no reordering required
                const int src_rank = binominal_tree_real_rank(src_vrank, root, size);
                MPI_Recv(recvbuf + recvcount * recv_size_per_element * src_rank,
                         blocks * recvcount, recvtype, src_rank, 0, comm, MPI_STATUS_IGNORE);
            } else {
                if (tmpbuf == NULL) {
                    // determine the number of expected blocks
                    int max_blocks;
                    if (vrank == 0) {
                        // root node
                        max_blocks = size;
                    } else {
                        // forwarding nodes, use max. number of childrens to
                        // estimate the max. number of blocks
                        max_blocks = min(1 << __builtin_ctz(vrank), size - vrank);
                    }
                    tmpbuf = (char*)malloc(sendcount * send_size_per_element * max_blocks);

                    if (vrank != 0) {
                        // copy sendbuf data to the beginning of the tmpbuf (because
                        // of the tmpbuf offset of vrank the subroot local data is
                        // always at the beginning of the tmpbuf)
                        memcpy(tmpbuf, sendbuf, sendcount * send_size_per_element);
                    }
                }

                // tmpbuf has an offset of vrank thus src_rank-vrank
                MPI_Recv(tmpbuf + sendcount * send_size_per_element * (src_vrank - vrank),
                         blocks * sendcount, sendtype,
                         binominal_tree_real_rank(src_vrank, root, size),
                         0, comm, MPI_STATUS_IGNORE);
            }
        }
        d <<= 1;
    }

    if (vrank != 0) {
        const int dst_vrank = vrank - d;
        const int blocks = min(d, size - vrank); // each node local data counts 1 block

        if (blocks == 1) {
            // leaf node
            assert(tmpbuf == NULL);
            MPI_Send(sendbuf, sendcount, sendtype,
                     binominal_tree_real_rank(dst_vrank, root, size), 0, comm);
        } else {
            // forwarding node
            assert(tmpbuf != NULL);
            MPI_Send(tmpbuf, blocks * sendcount, sendtype,
                     binominal_tree_real_rank(dst_vrank, root, size), 0, comm);
            free(tmpbuf);
        }
    } else {
        // root node
        if (!use_direct_recv) {
            // copy from tmpbuf to recvbuf and reorder the data, note that the
            // root local data is already in recvbuf
            assert(tmpbuf != NULL);
            MPI_Sendrecv(tmpbuf + sendcount * send_size_per_element,
                         sendcount * (size - rank - 1), sendtype, rank, 0,
                         recvbuf + recv_size_per_element * recvcount * (rank + 1),
                         recvcount * (size - rank - 1), recvtype, rank, 0,
                         comm, MPI_STATUS_IGNORE);
            MPI_Sendrecv(tmpbuf + send_size_per_element * sendcount * (size - rank),
                         sendcount * rank, sendtype, rank, 0,
                         recvbuf, recvcount * rank, recvtype, rank, 0,
                         comm, MPI_STATUS_IGNORE);
            free(tmpbuf);
        }
    }

    return MPI_SUCCESS;
}

int MY_Gather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount,
    MPI_Datatype recvtype, int root, MPI_Comm comm) {
    /*return gather_binominal_tree(sendbuf, sendcount, sendtype, recvbuf,
                                 recvcount, recvtype, root, comm);*/
    return gather_divide_and_conquer(sendbuf, sendcount, sendtype,
                                     recvbuf, recvcount, recvtype, root, comm);
}


void init_MY_Gather(int sendcount, MPI_Datatype sendtype, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) {
}


void cleanup_MY_Gather(int sendcount, MPI_Datatype sendtype, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) {
}

/***************************************/





/***************************************/
// Scatter (Divide-And-Conquer)

static void scatter_divide_and_conquer_inner(char **buffer, unsigned long buffer_offset,
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

    scatter_divide_and_conquer_inner(buffer, buffer_offset, sendcount, sendtype,
                                     recvbuf, recvcount, recvtype, start, end, newroot,
                                     comm, rank, size_per_element);
}

static int scatter_divide_and_conquer(const void* sendbuf, const int sendcount, const MPI_Datatype sendtype,
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

    if (recvbuf != MPI_IN_PLACE) {
        memset(recvbuf, 0, recvcount * recv_size_per_element);
    }

    scatter_divide_and_conquer_inner(&buffer, 0, sendcount, sendtype,
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


/***************************************/
// Scatter (binomial-tree)

static int binominal_tree_scatter(const char* sendbuf, const int sendcount, const MPI_Datatype sendtype,
                                 char* recvbuf, const int recvcount, const MPI_Datatype recvtype,
                                 const int root, const MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);

    int rank;
    MPI_Comm_rank(comm, &rank);
    const int vrank = binominal_tree_virtual_rank(rank, root, size);

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
    
    char* tmpbuffer = NULL;
    
    if (rank == root) {
        // root sends to children
        int d = 1;
        while (d < size) {
            const int blocks = min(d, size - d);
            const int real_recv = binominal_tree_real_rank(d, root, size);
            
            const int message_size = blocks * send_size_per_element * sendcount;
            const int sendbuf_offset = real_recv * send_size_per_element * sendcount;
            const int sendbuf_size = size * send_size_per_element * sendcount;
            const int sendbuf_remaining_bytes = sendbuf_size - sendbuf_offset;
            
            if (message_size > sendbuf_remaining_bytes) {
                // data wraps around sendbuf -> copy data to tmpbuffer before sending
                tmpbuffer = (char*)malloc(message_size);
                memcpy(tmpbuffer, sendbuf + sendbuf_offset, sendbuf_remaining_bytes);
                memcpy(tmpbuffer + sendbuf_remaining_bytes, sendbuf, message_size - sendbuf_remaining_bytes);
                
                MPI_Send(tmpbuffer, blocks * sendcount, sendtype, real_recv ,0, comm);
                free(tmpbuffer);
            } else {
                MPI_Send(sendbuf + sendbuf_offset, blocks * sendcount, sendtype, real_recv ,0, comm);
            }
            
            d <<= 1;
        }
        
        // copy node relevant data to receive buffer
        MPI_Sendrecv(sendbuf + sendcount * send_size_per_element * rank,
                     sendcount, sendtype, rank, 0, recvbuf, recvcount, recvtype,
                     rank, 0, comm, MPI_STATUS_IGNORE);
        
    } else {
        const int parent_vrank = binominal_tree_parent_rank(vrank, size);
        const int blocks = min(vrank - parent_vrank, size - vrank);
        const int real_sender = binominal_tree_real_rank(parent_vrank, root, size);
        
        if (blocks > 1) {
            // receive more than one block -> use tmpbuffer
            tmpbuffer = (char*)malloc(blocks * send_size_per_element * sendcount);
            MPI_Recv(tmpbuffer, blocks * sendcount, sendtype, real_sender, 0, comm, MPI_STATUS_IGNORE);
            
            // forward to children
            int d = 1;
            while((vrank & d) != d && ((vrank | d) < size)) {
                const int vrecv = vrank | d;
                const int shifted_vrecv = vrecv - vrank;
                const int real_recv = binominal_tree_real_rank(vrecv, root, size);
                const int blocks = min(shifted_vrecv, size - vrecv);
                
                MPI_Send(tmpbuffer + shifted_vrecv * send_size_per_element * sendcount, blocks * sendcount, sendtype, real_recv ,0, comm);
                d <<= 1;
            }
            
            // copy node relevant data to receive buffer
            MPI_Sendrecv(tmpbuffer, sendcount, sendtype, rank, 0, recvbuf, recvcount, recvtype,
                     rank, 0, comm, MPI_STATUS_IGNORE);
            free(tmpbuffer);
        } else {
            // direct receive
            MPI_Recv(recvbuf, blocks * recvcount, recvtype, real_sender, 0, comm, MPI_STATUS_IGNORE);
        }
    }
    return MPI_SUCCESS;
}

// Scatter
int MY_Scatter(const void* sendbuf, const int sendcount, const MPI_Datatype sendtype,
               void* recvbuf, const int recvcount, const MPI_Datatype recvtype,
               const int root, const MPI_Comm comm)
{
    //return scatter_divide_and_conquer(sendbuf, sendcount, sendtype, recvbuf,
    //                                  recvcount, recvtype, root, comm);
    
    return binominal_tree_scatter(sendbuf, sendcount, sendtype, recvbuf, 
                                  recvcount, recvtype, root, comm);
}

void init_MY_Scatter(int sendcount, MPI_Datatype sendtype, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) {
}


void cleanup_MY_Scatter(int sendcount, MPI_Datatype sendtype, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) {
}


/***************************************/

