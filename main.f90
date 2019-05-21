program ns_fd

  use contxt_class
  use   blas_class
  use domain_class
  use matrix_class

  implicit none

  ! Local Declarations
  ! ------------------

  type(geom_struc) :: geom
  type(matr_stor)  :: matr, lap_matr
  type(gauss_stor) :: g_matr, g_lap_matr

  integer  :: nbseg1, nbseg2,  Numline1, Numline2
  real(rp) :: x1, y1
  real(rp) :: x2, y2
  real(rp) :: x3, y3
  real(rp) :: x4, y4

  real(rp), dimension(:), allocatable :: rhs, U, V, Uh, Vh, P
  real(rp) 	:: h,k,nu,mu,Re,T
  integer	:: i,j,nod, it_n
  

  ! Beginning of the main routine
  ! -----------------------------

  write(*,*) "SOLVE THE NAVIER STOKES PROBLEM WITH FINITE DIFFERENCES ..."

  ! Get the mesh
  ! ------------

  read(*,*) x1, y1
  read(*,*) x2, y2
  read(*,*) x3, y3
  read(*,*) x4, y4
  
  
  

  read(*,*) nbseg1, nbseg2
  
  read(*,*) k
  
  read (*,*) T
  
  read(*,*) Re
  
  h = min(1.0_rp/nbseg1,1.0_rp/nbseg2)
  
  nu =  k/(h**2)
  
  mu = k/h
  
  write(*,*) " h = " ,h, "k= ", k, "mu =",  mu, "nu = ",nu
  
  call gm_rect_2d(geom, x1, y1, x2, y2, &
                        x3, y3, x4, y4, nbseg1, nbseg2)

  call gm_print(" The Mesh", geom)

  
  allocate(U(1:geom%nb_vers))
  allocate(V(1:geom%nb_vers))
!-------------------------------  
!  Initial data
!-----------------------------
  U = 0_rp
  V = 0_rp
  
!-------------------------
! First iteration  
!--------------------------------- 
  write(*,*)" ********************************"
  write(*,*)"                ITER k = : ", it_n	

	 
  it_n = 1	 
  
  ! Get the matrix
 
  
  
  call init_matr(matr, geom, U, V, mu, nu, Re, nbseg1+1, nbseg2+1)
  
  call mt_print(" The Matrix", matr)
  
  call csr_to_gauss_matr(matr,g_matr)
  

 ! call g_mt_print(" The Matrix", g_matr)
  


  ! Solve the problem for Uh, Vh, P
  ! -----------------

  allocate(rhs(g_matr%nb_row))
  allocate(Uh(g_matr%nb_row))
  allocate(Vh(g_matr%nb_row)) 
  allocate(P(g_matr%nb_row)) 
  
  !solve for Uh
  call init_rhs_u(rhs, geom, U)

  call gaussian(g_matr,rhs,Uh,50,'Row_op.bin',NumLine1)
  
  !solve for Vh
  
  call init_rhs_v(rhs, geom, V)
  
  call gaussian2(g_matr,rhs,Vh,50,'Row_op.bin',NumLine1)

  !solve for P	
  call init_lap_matr(lap_matr, geom,  nu, nbseg1+1, nbseg2+1)
	
  call init_rhs_p(rhs, geom,  Uh, Vh, nbseg1+1, nbseg2+1)
  
  call csr_to_gauss_matr(lap_matr,g_lap_matr)
  
  call gaussian(g_lap_matr,rhs,P,60,'Row_op_lap.bin',NumLine2)

  ! finding U,V from Uh, Vh, P

  call com_vel(geom, Uh, Vh, P, U, V, k, nbseg1+1, nbseg2+1)
  
  ! Writing data first iteration
 
  open(10, file = "solution.dat")

  
	
  nod = 0
  do j = 1, nbseg2+1
	do i = 1, nbseg1+1
	   nod = nod + 1
	   write(10,"(1x,5e17.9)") geom%ver_coo(:,nod), U(nod), V(nod), P(nod)
   end do
   write(10,*)
  end do
  write(10,*)
  
  
  !~  !solve for testing solution LAPLACE WITH NEUMANN BC
  !~
  !~call init_rhs(rhs, geom, nbseg1+1, nbseg2+1)
  !~
  !~call gaussian2(g_lap_matr,rhs,U,60,'Row_op_lap.txt',NumLine2)
  !~
  !~call graf_gnu(12,"model_sol.dat", geom, U)

  do it_n = 2,floor(T/k)
  
	 
	  !solve for Uh
	  call init_rhs_u(rhs, geom, U)

	  call gaussian2(g_matr,rhs,Uh,50,'Row_op.bin',NumLine1)
	  
	  !solve for Vh
	  
	  call init_rhs_v(rhs, geom, V)
	  
	  call gaussian2(g_matr,rhs,Vh,50,'Row_op.bin',NumLine1)

	  !solve for P	

		
	  call init_rhs_p(rhs, geom,  Uh, Vh, nbseg1+1, nbseg2+1)
	  
	  	  
	  call gaussian2(g_lap_matr,rhs,P,60,'Row_op_lap.bin',NumLine2)

	  ! finding U,V from Uh, Vh, P

	  call com_vel(geom, Uh, Vh, P, U, V, k, nbseg1+1, nbseg2+1)
	  
	  ! Writing data nth iteration	 
 
	if ( mod(it_n,5) .eq. 0 ) then 
	write(*,*)" ********************************"
	write(*,*)"                ITER k = : ", it_n	
    write(*,*)
		write(10,*)" ITER k = : ", it_n		
		  nod = 0
		  do j = 1, nbseg2+1
			do i = 1, nbseg1+1
			   nod = nod + 1
			   write(10,"(1x,5e17.9)") geom%ver_coo(:,nod), U(nod), V(nod), P(nod)
		   end do
		   write(10,*)
		  end do
		  write(10,*)
	end if
	
	
	
  end do
  
  
  
  
  write (*,*) 'Okay . RU gnuplot "plot-sol.gnu" to see the result'
	
   
  
  
  close(10)
  
  deallocate(rhs)
  deallocate(Uh)
  deallocate(Vh)
  deallocate(P)

  deallocate(U)
  deallocate(V)

!~  ! End of the main routine
  ! -----------------------
  

contains
!~
!~  subroutine set_sol(geom,sol)
!~
!~    implicit none
!~
!~    type(geom_struc)      , intent(in) :: geom
!~    real(rp), dimension(:), intent(out) :: sol
!~	
!~	integer i;
!~	
!~	do i=1,geom%nb_vers
!~	 
!~		sol(i) = g_D(geom%ver_coo(1,i),geom%ver_coo(2,i))
!~
!~	end do
!~  end subroutine set_sol
  
  subroutine graf_gnu(fileUit,fname, geom, U)

    implicit none

    character(*)          , intent(in) :: fname
    type(geom_struc)      , intent(in) :: geom
	integer, intent(in) 				:: fileUit
  !  type(matr_stor)       , intent(in) :: matr
    real(rp), dimension(:), intent(in) :: U

    ! Local declarations
    ! ------------------

    integer :: i, j, k

    ! Beginning of internal routine
    ! -----------------------------

    open(fileUit, file = fname)

    k = 0
    do j = 1, nbseg2+1
       do i = 1, nbseg1+1
          k = k + 1
          write(fileUit,"(1x,3e17.9)") geom%ver_coo(:,k), U(k)
       end do
       write(fileUit,*)
    end do

    close(fileUit)

    ! End of internal routine
    ! -----------------------

  end subroutine graf_gnu
  
  !-------------------
 ! Find U,V at next iteration from Uhalf, Vhalf,P
 
  subroutine com_vel(geom, Uh, Vh, P, U, V, k, nb_nod_h, nb_nod_v)

    implicit none
	
    real(rp), dimension(:), intent(out) :: U,V
	real(rp), dimension(:), intent(in) 	:: Uh,Vh,P
    type(geom_struc), intent(in) 		:: geom
	integer         , intent(in) 		:: nb_nod_v, nb_nod_h
	real(rp), intent(in)				:: k
 
	! local decrations
	integer ::  nod
	real(rp):: hx, hy

		
	U(:) = 0.0_rp
	V(:) = 0.0_rp
	
	hx = 1.0_rp / (nb_nod_h-1)
	hy = 1.0_rp / (nb_nod_v-1)
	
	do nod = 1,geom%nb_vers	
		
		if (geom%ver_zon(nod) .eq. 0) then
	
			U(nod) = Uh(nod) - k* (P(nod+1)-P(nod-1))/(2*hx) 
					
			V(nod) = Vh(nod) - k* (P(nod+nb_nod_h) - P(nod-nb_nod_h))/(2*hy)
				
		else	
			select case (geom%ver_zon(nod))
			
				case (1)
					U(nod) = Uh(nod) - k* (P(nod+1)-P(nod-1))/(2*hx) 
					V(nod) = Vh(nod) - k*(-4*P(nod+nb_nod_h) +7* P(nod+2*nb_nod_h)- 3* P(nod+3*nb_nod_h))/(hy)
				case (2)
					U(nod) = Uh(nod) - k* (4*P(nod-1)-7*P(nod-2)+3*P(nod-3))/(hx) 
					V(nod) = Vh(nod) - k*(P(nod+nb_nod_h) -  P(nod-nb_nod_h))/(2*hy)
				case (3)
					U(nod) = Uh(nod) - k* (P(nod+1)-P(nod-1))/(2*hx) 
					V(nod) = Vh(nod) - k*(4*P(nod-nb_nod_h) -7* P(nod-2*nb_nod_h) + 3* P(nod-3*nb_nod_h))/(hy)
				case (4)
					U(nod) = Uh(nod) - k* (-4*P(nod+1)+7*P(nod+2)-3*P(nod+3))/(hx)
					V(nod) = Vh(nod) - k*(P(nod+nb_nod_h) -  P(nod-nb_nod_h))/(2*hy)
				case (5)
					U(nod) = Uh(nod) - k* (-4*P(nod+1)+7*P(nod+2)-3*P(nod+3))/(hx) 
					V(nod) = Vh(nod) - k*(-4*P(nod+nb_nod_h) +7* P(nod+2*nb_nod_h)- 3* P(nod+3*nb_nod_h))/(hy)
				case (6)
					U(nod) = Uh(nod) - k*(4*P(nod-1)-7*P(nod-2)+3*P(nod-3))/(hx) 
					V(nod) = Vh(nod) - k*(-4*P(nod+nb_nod_h) +7* P(nod+2*nb_nod_h)- 3* P(nod+3*nb_nod_h))/(hy)
				case (7)
					U(nod) = Uh(nod) - k*  (4*P(nod-1)-7*P(nod-2)+3*P(nod-3))/(hx) 
					V(nod) = Vh(nod) - k* (4*P(nod-nb_nod_h) -7* P(nod-2*nb_nod_h) + 3* P(nod-3*nb_nod_h))/(hy)
				case (8)
					U(nod) = Uh(nod) - k*  (-4*P(nod+1)+7*P(nod+2)-3*P(nod+3))/(hx) 
					V(nod) = Vh(nod) - k* (4*P(nod-nb_nod_h) -7* P(nod-2*nb_nod_h) + 3* P(nod-3*nb_nod_h))/(hy)
				case default
			  
			end select
		end if	
				
	
	end do
	

  end subroutine com_vel		
    
end program ns_fd
