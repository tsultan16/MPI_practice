! Overview:  This program demonstrates pack and unpack subroutines that can be used to
!            communicate derived data types across MPI ranks. For this example, we use  
!            a "particle" data type and communicate data for lots of particles across neighboring ranks.
!


MODULE PARTICLE_TYPE

    TYPE particle
        INTEGER :: id = 0, x = 0, y = 0, z = 0
        TYPE(particle), POINTER :: next => null()
    END TYPE particle

    TYPE cell_head
         TYPE(particle), POINTER :: p => null()
    END TYPE cell_head

END MODULE PARTICLE_TYPE


PROGRAM datapack

USE MPI
USE PARTICLE_TYPE
IMPLICIT NONE


INTEGER :: comm2d, ierr, ndim, myrank, numprocs(1), dims(2), coord(2), req(4), &
           coltype, status(MPI_STATUS_SIZE), source, tag, destination, buffer_size
LOGICAL :: isperiodic(2), reorder 
INTEGER, PARAMETER :: N = 6, nb = 1, maxSteps = 2, maxparticles = 20, nvars = 4, max_buffer_particles = 5
REAL(8) :: ti, tf
INTEGER, ALLOCATABLE ::  N_particles(:,:)
INTEGER :: ftot(N+2,N+2)
INTEGER :: ilow, ihi, jlow, jhi, i, j, iter, np, p(2), my_particle_num, bndry
TYPE(cell_head), ALLOCATABLE :: cell_list(:,:)
TYPE(particle), POINTER :: current
REAL(8), ALLOCATABLE :: MPI_buffer_out(:), MPI_buffer_in(:) 


! MPI initialization
CALL MPI_INIT(ierr)

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs(1), ierr)  

np = numprocs(1)

! create a new cartesian communicator for our 2d domain decomposition

isperiodic(1) = .FALSE. ! periodic boundaries
isperiodic(2) = .FALSE. ! periodic boundaries
reorder = .TRUE. ! allow MPI rank reordering 
ndim = 2 ! 2d decomposition
dims(1) = 2 !SQRT(REAL(np))    ! number of ranks in x direction
dims(2) = 2 !SQRT(REAL(np)) ! number of ranks in y direction


p(1) = dims(1)
p(2) = dims(2)

CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndim, dims, isperiodic, reorder, comm2d, ierr)

CALL MPI_COMM_RANK(comm2d, myrank, ierr)

! get rank coordinate
CALL MPI_CART_COORDS(comm2d, myrank, ndim, coord, ierr) 


! compute boundary indices for this rank
CALL compute_bound(ilow, ihi, jlow, jhi)


!PRINT*,''
!PRINT*,'My rank, coordinate, ilow, ihi, jlow, jhi =',myrank, coord, ilow, ihi, jlow, jhi
!PRINT*,''

! initialize particle distribution
ALLOCATE(cell_list(ilow-nb:ihi+nb, jlow-nb:jhi+nb))
ALLOCATE(N_particles(ilow-nb:ihi+nb, jlow-nb:jhi+nb))
N_Particles = 0
CALL init_1d(N_particles, cell_list, ilow, ihi, jlow, jhi)

!allocate memory for MPI buffers
buffer_size = 1+nvars*max_buffer_particles
ALLOCATE(MPI_buffer_in(1:buffer_size),MPI_buffer_out(1:buffer_size))


! create MPI derived data type for column buffer
CALL MPI_TYPE_VECTOR(1+jhi-jlow,            &  ! count 
                     1,                     &  ! block length
                     3+ihi-ilow,            &  ! stride
                     MPI_DOUBLE_PRECISION,  &  ! oldtype
                     coltype,               &  ! newtype
                     ierr)
CALL MPI_TYPE_COMMIT(coltype, ierr)
 

CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
ti = MPI_WTIME()

! Jacobi iteration loop
iter = 0

DO i = 1,maxSteps

    !PRINT*,'Rank # = ',myrank,' Iteration #',iter
    ! exchange ghost point data with neighbor ranks 
    !PRINT*,'Rank # = ',myrank,' Starting Exchange.'

    CALL advect2d(N_particles, cell_list, ilow, ihi, jlow, jhi)

    DO bndry = 1, 4
    
        ! clear MPI buffers
        MPI_buffer_in = 0.0
        MPI_buffer_out = 0.0

        ! pack up outgoing particle buffer
        CALL pack_MPI_buffer_out(MPI_buffer_out, buffer_size, N_particles, cell_list, bndry, ilow, ihi, jlow, jhi)

        ! MPI communications
        CALL exchange_bndry_particles(MPI_buffer_out, MPI_buffer_in, buffer_size, bndry, ilow, ihi, jlow, jhi)

        ! unpack incoming particle buffer
        CALL unpack_MPI_buffer_in(MPI_buffer_in, N_particles, cell_list, bndry, ilow, ihi, jlow, jhi)
        !PRINT*,'Rank # = ',myrank,' Exchange Complete.'
    END DO
   
END DO


CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
tf = MPI_WTIME()

!PRINT*,''
!PRINT*,'Rank #',myrank, ' Done.'
!PRINT*,'Num processes = ', numprocs
!PRINT*,''

IF(myrank .EQ. 0) THEN
    PRINT*,''
    PRINT*,'Time elapsed(s) = ',tf-ti
    PRINT*,''
END IF


!PRINT*,'Rank# = ',myrank,' f = '
!DO j = jhi+1,jlow-1,-1
!    DO i = ilow-1, ihi+1
!        WRITE(*,FMT='(1i4)',ADVANCE='no') N_particles(i,j)
!        !WRITE(*,FMT='(1L4)',ADVANCE='no') ASSOCIATED(cell_list(i,j)%p)
!    END DO
!    PRINT*,''
!END DO


CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! I/O: Send flattened array back to rank 0
CALL prepare_u(N_particles)

!CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

!PRINT*,''
!PRINT*,'Rank #',myrank, ' Done.'
!PRINT*,''

DEALLOCATE(N_particles)



! free MPI derived type
CALL MPI_TYPE_FREE(coltype, ierr)
! MPI termination
CALL MPI_FINALIZE(ierr)


CONTAINS


SUBROUTINE compute_bound(ilow, ihi, jlow, jhi)

    INTEGER, INTENT(INOUT) :: ilow, ihi, jlow, jhi

    ilow = 0+coord(1)*(N+2)/dims(1)
    ihi = ilow-1+(N+2)/dims(1)
    jlow = 0+coord(2)*(N+2)/dims(2)
    jhi = jlow-1+(N+2)/dims(2)

END SUBROUTINE compute_bound


SUBROUTINE init_1d(f, cell_list, ilow, ihi, jlow, jhi)
    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi
    INTEGER, INTENT(INOUT) :: f(ilow-nb:ihi+nb, jlow-nb:jhi+nb)
    TYPE(cell_head), INTENT(INOUT) :: cell_list(ilow-nb:ihi+nb, jlow-nb:jhi+nb)
    TYPE(particle), POINTER :: current
    INTEGER :: i, j, k, my_start_id, my_end_id, counter

    my_particle_num = 0
    f(:,:) = 0
    DO i = ilow, ihi
        DO j = jlow, jhi
            f(i,j) = fn(i,j)
            IF(f(i,j) .GT. 0) my_particle_num = my_particle_num + f(i,j) 
        END DO
    END DO

    ! compute start ID for my particles
    IF(myrank .EQ. 0) THEN

       my_start_id = 1
       my_end_id = my_start_id + my_particle_num - 1
       CALL MPI_SEND(my_end_id, 1, MPI_INTEGER, myrank+1, myrank, MPI_COMM_WORLD, ierr)

    ELSE IF(myrank .EQ. numprocs(1)-1) THEN

       CALL MPI_RECV(my_start_id, 1, MPI_INTEGER, myrank-1, myrank-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
       my_start_id = my_start_id + 1
       my_end_id = my_start_id + my_particle_num - 1
       
    ELSE

       CALL MPI_RECV(my_start_id, 1, MPI_INTEGER, myrank-1, myrank-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
       my_start_id = my_start_id + 1
       my_end_id = my_start_id + my_particle_num - 1
       CALL MPI_SEND(my_end_id, 1, MPI_INTEGER, myrank+1, myrank, MPI_COMM_WORLD, ierr) 
   
    END IF

    !PRINT*,'my_start_id, my_end_id =',my_start_id,my_end_id

    counter = 0
    ! allocate memory for particles
    DO i = ilow, ihi
        DO j = jlow, jhi

            IF(f(i,j) .GT. 0) THEN
       
               
                ALLOCATE(cell_list(i,j)%p)
                current => cell_list(i,j)%p
                ! define particle state
                current%id = my_start_id+counter
                current%x = i
                current%y = j 
                counter = counter + 1

                DO k = 1, f(i,j)-1
       
                    ALLOCATE(current%next)
                    ! define particle state
                    current=>current%next
                    current%id = my_start_id+counter
                    current%x = i
                    current%y = j 
                    counter = counter + 1

                END DO        

            END IF

        END DO
    END DO

    !PRINT*,'Particle memory allocation complete.'

    DO i = ilow, ihi
        DO j = jlow, jhi
    
        IF(f(i,j) .GT. 0) THEN
            !PRINT*,'CELL# ',i,j,', N=',f(i,j)
            !PRINT*,''
            current => cell_list(i,j)%p
            
            DO WHILE(ASSOCIATED(current))
                !PRINT*,'Tracer ID, x y pos ',current%id,current%x,current%y
                current => current%next
            END DO
            !PRINT*,''
        END IF
        
        END DO
        !PRINT*,''
    END DO

END SUBROUTINE init_1d


! initial particle distribution
FUNCTION fn(i,j) RESULT(fx)
    INTEGER, INTENT(IN) :: i,j
    INTEGER :: fx
    
    IF(i .GE. N/2 .AND. i .LE. 1+N/2 .AND. j .GE. N/2 .AND. j .LE. 1+N/2) THEN
        fx = 5   
    ELSE
        fx = 0
    END IF
   
   

END FUNCTION fn


SUBROUTINE advect2d(N_particles, cell_list, ilow, ihi, jlow, jhi)
    INTEGER, INTENT(INOUT) :: N_particles(ilow-nb:ihi+nb, jlow-nb:jhi+nb)
    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi
    TYPE(cell_head), INTENT(INOUT) :: cell_list(ilow-nb:ihi+nb, jlow-nb:jhi+nb)

    INTEGER :: i,j,k,nout, N_temp(ilow-nb:ihi+nb, jlow-nb:jhi+nb), dest(2)
    TYPE(particle), POINTER :: prev_src, current_src, current_des

    
    N_temp = N_Particles

    ! move half of the particles inside each cell above
    DO i = ilow, ihi
        DO j = jlow, jhi

            dest(:) = (/ i+1, j/)

            nout = N_Particles(i,j)/2
            IF(nout .GT. 0) THEN
                !PRINT*,'CELL #',i,j
               
                DO k = 1, nout
                    current_src => cell_list(i,j)%p
                    prev_src => cell_list(i,j)%p
                    ! remove last particle from source list
                    DO WHILE(ASSOCIATED(current_src%next))
                        prev_src => current_src 
                        current_src => current_src%next
                    END DO

                    !PRINT*,'Moving Particle#,',current_src%id,' into cell#',dest(1),dest(2)

                    ! link particle to destination list 
                    current_des => cell_list(dest(1),dest(2))%p
                    IF(ASSOCIATED(current_des)) THEN ! traverse to last particle in destination list
                        DO WHILE(ASSOCIATED(current_des%next)) 
                            current_des => current_des%next
                        END DO
                        current_des%next => current_src

                    ELSE 
                        cell_list(dest(1),dest(2))%p => current_src
                    END IF

                    ! unlink particle from source cell
                    prev_src%next => null()              
                    current_des => current_des%next

                    ! update particle state
                    current_src%x = dest(1)
                    current_src%y = dest(2) 

                END DO

                ! update particle numbers
                N_temp(i,j) = N_temp(i,j) - nout
                N_temp(dest(1),dest(2)) = N_temp(dest(1),dest(2)) + nout 

            END IF
        END DO
    END DO

    N_Particles = N_temp

END SUBROUTINE advect2d


SUBROUTINE pack_MPI_buffer_out(buffer, buffer_size, N_particles, cell_list, bndry, ilow, ihi, jlow, jhi)

    REAL(8), INTENT(INOUT) :: buffer(:)
    INTEGER, INTENT(IN) :: buffer_size, ilow, ihi, jlow, jhi, bndry
    INTEGER, INTENT(INOUT) :: N_particles(ilow-nb:ihi+nb, jlow-nb:jhi+nb)
    TYPE(cell_head), INTENT(INOUT) :: cell_list(ilow-nb:ihi+nb, jlow-nb:jhi+nb)

    TYPE(particle), POINTER :: current, prev
    INTEGER :: i, j, counter, nout, i_b, i_v, &
               imin, imax, jmin, jmax
    REAL(8) :: particle_data(4)
    LOGICAL :: buffer_full, all_particles_buffered

    !##########################
    !copy particle data into buffer
    !##########################

    buffer_full = .FALSE.

    SELECT CASE(bndry)

    ! left x boundary
    CASE(1) 

        imin = ilow - nb
        imax = imin
        jmin = jlow
        jmax = jhi    

        !PRINT*,'Left x boundary.'

    ! right x boundary
    CASE(2)

        imin = ihi + nb
        imax = imin
        jmin = jlow
        jmax = jhi    

        !PRINT*,'Right x boundary.'

    ! bottom y boundary
    CASE(3)

        imin = ilow
        imax = ihi
        jmin = jlow-nb
        jmax = jmin 

        !PRINT*,'Bottom y boundary.'

    ! top y boundary
    CASE(4)

        imin = ilow 
        imax = ihi
        jmin = jhi+nb
        jmax = jmin    

        !PRINT*,'Top y boundary.'

    END SELECT


    ! find out how many particles in boundary
    nout = 0

    DO i = imin, imax
        DO j = jmin, jmax
            nout = nout + N_particles(i,j) 
        END DO
    END DO

    buffer(1) = nout  ! first entry contains number of particles in buffer 

    IF(nout .GT. max_buffer_particles) THEN
        PRINT*,'WARNING! NOT ENOUGH SPACE IN BUFFER TO SEND ALL PARTICLES! Nout = ',nout
        
        buffer_full = .TRUE.
        all_particles_buffered = .FALSE.     

        STOP

    END IF

    !PRINT*,'NOUT = ',nout
   
    i_b = 2

    IF(nout .GT. 0) THEN

        DO i = imin, imax
            DO j = jmin, jmax

                prev => cell_list(i,j)%p
                current => cell_list(i,j)%p 
      
                DO WHILE(ASSOCIATED(current))

                    !PRINT*,'CELL #',i,j,'. Packing particle #',current%id,' into buffer.'

                    particle_data(1) = current%id   
                    particle_data(2) = current%x
                    particle_data(3) = current%y
                    particle_data(4) = current%z

                    DO i_v = 1,nvars
                        buffer(i_b) = particle_data(i_v)
                        i_b = i_b + 1    
                    END DO

                    current => current%next

                               
                    ! deallocate particle memory
                    DEALLOCATE(prev)  

                    prev => current
    
                END DO

                ! nullify the pointer
                cell_list(i,j)%p => null()
                
            END DO
        END DO    

        ! update particle number in boundary cells
        DO i = imin, imax
            DO j = jmin, jmax
                N_particles(i,j) = 0 
            END DO
        END DO

    END IF

    


END SUBROUTINE pack_MPI_buffer_out



SUBROUTINE unpack_MPI_buffer_in(buffer, N_particles, cell_list, bndry, ilow, ihi, jlow, jhi)
  
    REAL(8), INTENT(INOUT) :: buffer(:)
    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi, bndry
    INTEGER, INTENT(INOUT) :: N_particles(ilow-nb:ihi+nb, jlow-nb:jhi+nb)
    TYPE(cell_head), INTENT(INOUT) :: cell_list(ilow-nb:ihi+nb, jlow-nb:jhi+nb)

    TYPE(particle), POINTER :: current
    INTEGER :: i, j, k, counter, nout, i_b, i_v, &
               imin, imax, jmin, jmax, dest(2)
    REAL(8) :: particle_data(4)
    

    !#################################
    !retrive particle data from buffer
    !#################################

    SELECT CASE(bndry)

    ! left x boundary
    CASE(1) 

        imin = ilow - nb
        imax = imin
        jmin = jlow
        jmax = jhi    

        !PRINT*,'Left x boundary.'

    ! right x boundary
    CASE(2)

        imin = ihi + nb
        imax = imin
        jmin = jlow
        jmax = jhi    

        !PRINT*,'Right x boundary.'

    ! bottom y boundary
    CASE(3)

        imin = ilow
        imax = ihi
        jmin = jlow-nb
        jmax = jmin 

        !PRINT*,'Bottom y boundary.'

    ! top y boundary
    CASE(4)

        imin = ilow 
        imax = ihi
        jmin = jhi+nb
        jmax = jmin    

        !PRINT*,'Top y boundary.'

    END SELECT


    ! find out how many particles in buffer
    nout = buffer(1)

    ! allocate memory for these particles and add insert them into the appropriate cells
    i_b = 2
    DO k = 1, nout 
        particle_data(:) = buffer(i_b:i_b+nvars)  
        i_b = i_b + nvars

        dest(:) = particle_data(2:3)
        !PRINT*,'Incoming Particle#',particle_data(1),' into cell #',dest(1),dest(2)
         

        ! link particle to destination list 
        current => cell_list(dest(1),dest(2))%p
        IF(ASSOCIATED(current)) THEN ! traverse to last particle in destination list
            DO WHILE(ASSOCIATED(current%next)) 
                current => current%next
            END DO
            ALLOCATE(current%next)
            current => current%next             
        ELSE 
            ALLOCATE(cell_list(dest(1),dest(2))%p) 
            current => cell_list(dest(1),dest(2))%p
        END IF

        ! update particle state
        current%id = particle_data(1)    
        current%x = particle_data(2) 
        current%y = particle_data(3) 
        current%z = particle_data(4) 

        ! update particle number array
        N_particles(dest(1),dest(2)) = N_particles(dest(1),dest(2))  + 1 

    END DO

END SUBROUTINE unpack_MPI_buffer_in


SUBROUTINE exchange_bndry_particles(buffer_out, buffer_in, buffer_size, bndry, ilow, ihi, jlow, jhi)

    REAL(8), INTENT(INOUT) :: buffer_out(:), buffer_in(:)
    INTEGER, INTENT(IN) :: buffer_size, bndry, ilow, ihi, jlow, jhi

    INTEGER :: bottom, top, left, right, neighbor_rank

 
    INTEGER :: i,tag,nx
    ! get neighbor ranks
    CALL MPI_CART_SHIFT(comm2d, 0, 1, left, right, ierr) 
    CALL MPI_CART_SHIFT(comm2d, 1, 1, bottom, top, ierr)
     
    SELECT CASE(bndry)

    CASE(1) ! x-
        neighbor_rank = left
        !PRINT*,'Exchanging with left neighbor.' 
    CASE(2) ! x+
        neighbor_rank = right
        !PRINT*,'Exchanging with right neighbor.'
    CASE(3) ! y-
        neighbor_rank = bottom
        !PRINT*,'Exchanging with bottom neighbor.'
    CASE(4) ! y+
        neighbor_rank = top
        !PRINT*,'Exchanging with top neighbor.'
    END SELECT
    
   
    tag = 0
    ! send and receive particle buffers from neighbor

    CALL MPI_SENDRECV(buffer_out,            &  ! send buffer
                      buffer_size,           &  ! send count
                      MPI_DOUBLE_PRECISION,  &  ! data type
                      neighbor_rank,         &  ! dest
                      tag,                   &  ! send tag
                      buffer_in,             &  ! recv buffer
                      buffer_size,           &  ! recv count
                      MPI_DOUBLE_PRECISION,  &  ! data type
                      neighbor_rank,         &  ! source
                      tag,                   &  ! recv tag
                      comm2d, MPI_STATUS_IGNORE, ierr)


END SUBROUTINE exchange_bndry_particles



SUBROUTINE prepare_u(u_new)

    INTEGER, INTENT(IN) :: u_new(:,:)
    INTEGER :: i, j, s_coord(2)
    INTEGER :: u_flat(((N+2)/dims(1))*((N+2)/dims(2)))


    u_flat = 0
    IF(np .GT. 1) THEN 
        IF(myrank .NE. 0) THEN
            tag = myrank
            CALL flatten_array(u_new, u_flat)
            CALL MPI_SEND(u_flat, ((N+2)/dims(1))*((N+2)/dims(2)), MPI_INTEGER, 0, tag, comm2d, ierr)
           
            !PRINT*,'Rank #',myrank,'u_my ='           
            !DO j = 2+(N+2)/np,1,-1
            !    DO i = 1, N+2
            !        WRITE(*,FMT='(1f8.3)',ADVANCE='no') u_new(i,j)
            !    END DO
            !    PRINT*,''
            !END DO
            !PRINT*,''

        ELSE

            !PRINT*,'Rank #',myrank,'u_my ='
            !DO j = 2+(N+2)/np,1,-1
            !    DO i = 1, N+2
            !        WRITE(*,FMT='(1f8.3)',ADVANCE='no') u_new(i,j)
            !    END DO
            !    PRINT*,''
            !END DO
            !PRINT*,''

            CALL flatten_array(u_new, u_flat)
            CALL unflatten_save(u_flat, (/ 0, 0 /)) 
                
            DO i= 1, np-1
                tag = i
                source = i
                CALL MPI_RECV(u_flat, ((N+2)/dims(1))*((N+2)/dims(2)), MPI_INTEGER, source, tag, comm2d, status, ierr)   

                ! get source coordinate
                CALL MPI_CART_COORDS(comm2d, tag, ndim, s_coord, ierr) 

             
                CALL unflatten_save(u_flat, s_coord)
            END DO

            PRINT*,'f ='
            DO j = N+2,1,-1
                DO i = 1, N+2
                    WRITE(*,FMT='(1i4)',ADVANCE='no') ftot(i,j)
                END DO
                PRINT*,''
            END DO

        END IF
    ELSE 
        tag = 0
        
        !PRINT*,'u='
        !DO j = N+2,1,-1
        !    DO i = 1, N+2
        !        u(i,j) = u_new(i,j+1)
        !        !WRITE(*,FMT='(1f8.3)',ADVANCE='no') u(i,j)
        !    END DO
        !    PRINT*,''
        !END DO

    END IF


END SUBROUTINE prepare_u



SUBROUTINE flatten_array(u_this, u_flat)
    INTEGER, INTENT(IN) :: u_this(:,:)
    INTEGER, INTENT(INOUT) :: u_flat(:)
    INTEGER :: i, j, ix
  
    DO i = 2, 2-1+(N+2)/p(1)
        DO j = 2, 2-1+(N+2)/p(2)
            ix = 1+(i-2)+(j-2)*((N+2)/p(1))
            u_flat(ix)  = u_this(i,j)   
        END DO
    END DO


END SUBROUTINE flatten_array



SUBROUTINE unflatten_save(u_flat, s_coord)
    INTEGER, INTENT(IN) :: u_flat(:)
    INTEGER, INTENT(IN) :: s_coord(2)

    INTEGER :: ilow, ihi, jlow, jhi
    
    INTEGER :: i, j, ix

    ilow = 1+s_coord(1)*(N+2)/p(1)
    ihi = ilow-1+(N+2)/p(1)
    jlow = 1+s_coord(2)*(N+2)/p(2)
    jhi = jlow-1+(N+2)/p(2)

    
    DO i = ilow, ihi
        DO j = jlow, jhi
            ix = 1+(i-ilow)+(j-jlow)*(N+2)/p(1)
            ftot(i,j) = u_flat(ix)  
        END DO
    END DO

    
END SUBROUTINE unflatten_save



END PROGRAM datapack
