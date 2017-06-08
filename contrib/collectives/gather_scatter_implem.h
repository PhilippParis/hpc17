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

#ifndef GATHER_SCATTER_IMPLEM_H_
#define GATHER_SCATTER_IMPLEM_H_

/***************************************/
// Gather
int MY_Gather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount,
    MPI_Datatype recvtype, int root, MPI_Comm comm);

void init_MY_Gather(int sendcount, MPI_Datatype sendtype, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);

void cleanup_MY_Gather(int sendcount, MPI_Datatype sendtype, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
/***************************************/

/***************************************/
// Scatter
int MY_Scatter(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount,
    MPI_Datatype recvtype, int root, MPI_Comm comm);

void init_MY_Scatter(int sendcount, MPI_Datatype sendtype, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);

void cleanup_MY_Scatter(int sendcount, MPI_Datatype sendtype, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);

/***************************************/

#endif /* GATHER_SCATTER_IMPLEM_H_ */
