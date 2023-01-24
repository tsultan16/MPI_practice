program test_MPI_RMA                                                                                                                                                      

use mpi                                                                                                                                                     


integer, parameter :: m=10

 

INTEGER i,map(m),comm,ierr,nproc,rank

REAL*8 A(m),B(m)

                                                                                                                                                           

call MPI_INIT(ierr)

call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)

call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)


PRINT*,"nproc=",nproc,"rank=",rank

 

do i=1,m

   A(i)=rank*0.1d0

   B(i)=0.d0

   if (rank==0) then

      map(i)=1

   else

      map(i)=0

   endif

   print*, 'node', rank, A(i), map(i)

enddo

 

call MPI_BARRIER(MPI_COMM_WORLD, ierr)

 

call SUM(A, B, map, m, MPI_COMM_WORLD, nproc)

 

call MPI_FINALIZE(ierr)

 

end program test_MPI_RMA


SUBROUTINE SUM(A, B, map, m, comm, p)

USE MPI

 

INTEGER m, map(m), comm, p, sizeofreal, win, ierr

REAL*8 A(m), B(m)

 

CALL MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, sizeofreal, ierr)

 

CALL MPI_WIN_CREATE(B, m*sizeofreal, sizeofreal, MPI_INFO_NULL,  &
                    comm, win, ierr)

 

CALL MPI_WIN_FENCE(0, win, ierr)

DO i=1,m

  j = map(i)/m

  k = MOD(map(i),m)

  CALL MPI_ACCUMULATE(A(i), 1, MPI_DOUBLE_PRECISION, j, k, 1, MPI_DOUBLE_PRECISION,   &
                      MPI_SUM, win, ierr)

END DO

CALL MPI_WIN_FENCE(0, win, ierr)

 

CALL MPI_WIN_FREE(win, ierr)

RETURN

END
