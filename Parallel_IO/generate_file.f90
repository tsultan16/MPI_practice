PROGRAM generate_file

IMPLICIT NONE

REAL*8, ALLOCATABLE :: array(:)
INTEGER :: i, arr_size

! set array size in bytes
arr_size = 20

ALLOCATE(array(arr_size))

DO i = 1, arr_size
    array(i) = i
END DO

OPEN(UNIT=10, FILE = 'output.dat', FORM = 'UNFORMATTED', ACCESS = 'STREAM')

DO i = 1, arr_size
    WRITE(10) array(i)
END DO

CLOSE(UNIT=10)

PRINT*,''
PRINT*,'Sample unformatted binary file generated. File size (bytes) =',8*arr_size
PRINT*,''

DEALLOCATE(array)

END PROGRAM generate_file