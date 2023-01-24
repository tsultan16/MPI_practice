PROGRAM read_file

IMPLICIT NONE

INTEGER, PARAMETER :: nranks = 10
REAL(8), ALLOCATABLE :: array(:)
INTEGER :: i, arr_size

! set array size in bytes
arr_size = 0
DO i = 1, nranks
    arr_size = arr_size + i 
END DO


ALLOCATE(array(arr_size))


OPEN(UNIT=10, FILE = 'output_parallel.dat', FORM = 'UNFORMATTED', ACCESS = 'STREAM', STATUS = 'OLD')

DO i = 1, arr_size
    READ(10) array(i)
END DO

CLOSE(UNIT=10)

PRINT*,''
PRINT*,'Data read from parallel file = ',array
PRINT*,''

DEALLOCATE(array)

END PROGRAM read_file