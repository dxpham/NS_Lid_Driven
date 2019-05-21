module matrix_class

  use contxt_class
  use domain_class
  use function_class

!  use ell_stor_class ! or csr_stor_class
  use csr_stor_class ! or ell_stor_class

contains

!-------------------
! rhs for the V component of the velocity
!----------------------------
  subroutine init_rhs_u(rhs, geom, U)

    implicit none
	
	real(rp), dimension(:), allocatable :: rhs
	real(rp), dimension(:), intent(in) :: U
    type(geom_struc), intent(in) :: geom
 
	! local decrations
	integer ::  nod

	
	rhs(:) = 0.0_rp
	
	do nod = 1,geom%nb_vers	
		
		if (geom%ver_zon(nod) .eq. 0) then
	
			rhs(nod) = U(nod)!f(geom%ver_coo(1,nod),geom%ver_coo(2,nod))       f = 0
				
		else	
			rhs(nod) = gu_D(geom%ver_coo(1,nod),geom%ver_coo(2,nod)) !~ Dirichlet values at the all four boundaries 
		end if	
				
	
	end do
	
	!~! Dirichlet points upper and lower boundaries
	!~do el = 1, nb_nod_h
	!~
	!~	rhs(el) = g_D(geom%ver_coo(1,el),geom%ver_coo(2,el))
	!~	i=geom%nb_vers-(el-1)
	!~	rhs(i) = g_D(geom%ver_coo(1,i),geom%ver_coo(2,i))
	!~		end do
	!~
	!~! Dirichlet points left and right boundaries
	!~
	!~
	!~do i = 1, nb_nod_v-2
	!~	
	!~	el= 1 + i*nb_nod_h
	!~	rhs(el) = g_D(geom%ver_coo(1,el),geom%ver_coo(2,el))
	!~	el= nb_nod_h + i*nb_nod_h
	!~	rhs(el) = g_D(geom%ver_coo(1,el),geom%ver_coo(2,el))
	!~	
	!~end do
	!~
	
!~	! Neuman left and right boundaries
!~	
!~	i=1
!~	do el = 1, nb_nod_v-1
!~
!~		!	rhs(el) = 
!~		
!~	end do
!~	
  end subroutine init_rhs_u

!-------------------
! rhs for the V component of the velocity
!----------------------------
  subroutine init_rhs_v(rhs, geom, V)

    implicit none
	
    real(rp), dimension(:), allocatable :: rhs
	real(rp), dimension(:), intent(in) :: V
    type(geom_struc), intent(in) :: geom
 
	! local decrations
	integer ::  nod

	
	rhs(:) = 0.0_rp
	
	do nod = 1,geom%nb_vers	
		
		if (geom%ver_zon(nod) .eq. 0) then
	
			rhs(nod) = V(nod)!f(geom%ver_coo(1,nod),geom%ver_coo(2,nod))       f = 0
				
		else	
			rhs(nod) = 0_rp  !gu_D(geom%ver_coo(1,nod),geom%ver_coo(2,nod))      V=0 on boundary!~ Dirichlet values at the all four boundaries 
		end if	
				
	
	end do
	

  end subroutine init_rhs_v		

!-------------------
! rhs for the pressure 
!----------------------------
  subroutine init_rhs_p(rhs, geom, U, V, nb_nod_h, nb_nod_v)

    implicit none
	
    real(rp), dimension(:), intent(out) :: rhs
	real(rp), dimension(:), intent(in) 	:: U,V
    type(geom_struc), intent(in) 		:: geom
	integer         , intent(in) 		:: nb_nod_v, nb_nod_h
 
	! local decrations
	integer ::  nod
	real(rp):: hx, hy

	hx = 1.0_rp / (nb_nod_h-1)
	hy = 1.0_rp / (nb_nod_v-1)
		
	rhs(:) = 0.0_rp
	
	do nod = 1,geom%nb_vers	
		
		if (geom%ver_zon(nod) .eq. 0) then
	
			rhs(nod) = (U(nod+1)-U(nod-1))/(2*hx) +&
					(V(nod+nb_nod_h) - V(nod-nb_nod_h))/(2*hy)
				
		else	
			select case (geom%ver_zon(nod))
			
				case (1)
					rhs(nod) = (U(nod+1)-U(nod-1))/(2*hx) + &
					(-4*V(nod+nb_nod_h) +7* V(nod+2*nb_nod_h)- 3* V(nod+3*nb_nod_h))/(hy)
				case (2)
					rhs(nod) = (4*U(nod-1)-7*U(nod-2)+3*U(nod-3))/(hx) + &
						(V(nod+nb_nod_h) -  V(nod-nb_nod_h))/(2*hy)
				case (3)
					rhs(nod) = (U(nod+1)-U(nod-1))/(2*hx) + &
						(4*V(nod-nb_nod_h) -7* V(nod-2*nb_nod_h) + 3* V(nod-3*nb_nod_h))/(hy)
				case (4)
					rhs(nod) = (-4*U(nod+1)+7*U(nod+2)-3*U(nod+3))/(hx) + &
						(V(nod+nb_nod_h) -  V(nod-nb_nod_h))/(2*hy)
				case (5)
					rhs(nod) = (-4*U(nod+1)+7*U(nod+2)-3*U(nod+3))/(hx) &
					+ (-4*V(nod+nb_nod_h) +7* V(nod+2*nb_nod_h)- 3* V(nod+3*nb_nod_h))/(hy)
				case (6)
					rhs(nod) = (4*U(nod-1)-7*U(nod-2)+3*U(nod-3))/(hx) &
					+(-4*V(nod+nb_nod_h) +7* V(nod+2*nb_nod_h)- 3* V(nod+3*nb_nod_h))/(hy)
				case (7)
					rhs(nod) = (4*U(nod-1)-7*U(nod-2)+3*U(nod-3))/(hx) &
					+ (4*V(nod-nb_nod_h) -7* V(nod-2*nb_nod_h) + 3* V(nod-3*nb_nod_h))/(hy)
				case (8)
					rhs(nod) = (-4*U(nod+1)+7*U(nod+2)-3*U(nod+3))/(hx) &
					+ (4*V(nod-nb_nod_h) -7* V(nod-2*nb_nod_h) + 3* V(nod-3*nb_nod_h))/(hy)
				case default
			  
			end select
		end if	
				
	
	end do
	

  end subroutine init_rhs_p		
  
  
  !-------------------
! rhs for the Neumann Laplace problem 
!----------------------------
  subroutine init_rhs(rhs, geom, nb_nod_h, nb_nod_v)

    implicit none
	
    real(rp), dimension(:), intent(out) :: rhs
    type(geom_struc), intent(in) 		:: geom
	integer         , intent(in) 		:: nb_nod_v, nb_nod_h
 
	! local decrations
	integer ::  nod
	real(rp):: hx, hy

	hx = 1.0_rp / (nb_nod_h-1)
	hy = 1.0_rp / (nb_nod_v-1)
		
	rhs(:) = 0.0_rp
	
	do nod = 1,geom%nb_vers	
		
		if (geom%ver_zon(nod) .eq. 0) then
	
			rhs(nod) = f(geom%ver_coo(1,nod),geom%ver_coo(2,nod)) 
				
		else	
			select case (geom%ver_zon(nod))
			
				case (1)
					rhs(nod) =  f(geom%ver_coo(1,nod),geom%ver_coo(2,nod))  + &
					2.0_rp/hy* gy_N(geom%ver_coo(1,nod),geom%ver_coo(2,nod))
				case (2)
					rhs(nod) =  f(geom%ver_coo(1,nod),geom%ver_coo(2,nod))  - &
					2.0_rp/hy* gx_N(geom%ver_coo(1,nod),geom%ver_coo(2,nod))
				case (3)
					rhs(nod) =  f(geom%ver_coo(1,nod),geom%ver_coo(2,nod))  - &
					2.0_rp/hy* gy_N(geom%ver_coo(1,nod),geom%ver_coo(2,nod))
				case (4)
					rhs(nod) =  f(geom%ver_coo(1,nod),geom%ver_coo(2,nod))  + &
					2.0_rp/hx* gx_N(geom%ver_coo(1,nod),geom%ver_coo(2,nod))
				case (5)
					rhs(nod) =  f(geom%ver_coo(1,nod),geom%ver_coo(2,nod))  + &
					2.0_rp/hx* gx_N(geom%ver_coo(1,nod),geom%ver_coo(2,nod)) + &
					2.0_rp/hx* gy_N(geom%ver_coo(1,nod),geom%ver_coo(2,nod))
				case (6)
					rhs(nod) = f(geom%ver_coo(1,nod),geom%ver_coo(2,nod))  - &
					2.0_rp/hx* gx_N(geom%ver_coo(1,nod),geom%ver_coo(2,nod)) + &
					2.0_rp/hx* gy_N(geom%ver_coo(1,nod),geom%ver_coo(2,nod))
				case (7)
					rhs(nod) = f(geom%ver_coo(1,nod),geom%ver_coo(2,nod))  - &
					2.0_rp/hx* gx_N(geom%ver_coo(1,nod),geom%ver_coo(2,nod)) - &
					2.0_rp/hx* gy_N(geom%ver_coo(1,nod),geom%ver_coo(2,nod))
				case (8)
					rhs(nod) = f(geom%ver_coo(1,nod),geom%ver_coo(2,nod))  + &
					2.0_rp/hx* gx_N(geom%ver_coo(1,nod),geom%ver_coo(2,nod)) - &
					2.0_rp/hx* gy_N(geom%ver_coo(1,nod),geom%ver_coo(2,nod))
				case default
			  
			end select
		end if	
				
	
	end do
	

  end subroutine init_rhs	


!-------------------------------
! Set matrix for the velocity navier stokes
!-----------------------------------
	
  subroutine init_matr(matr, geom, Un, Vn, mu, nu, Re, nb_nod_h, nb_nod_v)

    implicit none

    type(matr_stor), intent(out) 		:: matr
    type(geom_struc), intent(in) 		:: geom
	real(rp), dimension(:), intent(in)	:: Un, Vn
	real(rp), intent(in)				:: mu, nu, Re
    integer         , intent(in) 		:: nb_nod_v, nb_nod_h

    ! Local Declarations
    ! ------------------

    integer :: tt_nzc, nb_nzc, nb_lin   !tt_nzc total nonzero const nb_nzc: neighbor nonzero const
    integer :: nod
	real(rp) :: c, e, w, n, s

    ! Beginning of the routine
    ! ------------------------

    ! Matrix allocation
    ! -----------------

    tt_nzc = geom%nb_elms &
           + 2 * (geom%nb_elms -1) &
           + 2 *  (geom%nb_elms-nb_nod_v) 

    nb_nzc = 4

    nb_lin = nb_nod_v * nb_nod_h

    call alloc_matr(matr, tt_nzc, nb_nzc, nb_lin, nb_lin)

    matr%nod_zon => geom%ver_zon
	


    ! Matrix construction
    ! -------------------

	do nod = 1,geom%nb_vers	
		
		if (geom%ver_zon(nod) .eq. 0) then
	
			c = 1.0_rp + 4.0_rp*nu/Re
			
			e = Un(nod)/2*mu - nu/Re
			
			w = -Un(nod)/2*mu - nu/Re
	
			n = Vn(nod)/2*mu - nu/Re
			
			s = -Vn(nod)/2*mu - nu/Re
		
			call fem_add_matr_coef(matr, nod, nod, c) 
			call fem_add_matr_coef(matr, nod, nod + 1, e) 
			call fem_add_matr_coef(matr, nod, nod - 1, w) 
			call fem_add_matr_coef(matr, nod, nod + nb_nod_h, n) 
			call fem_add_matr_coef(matr, nod, nod - nb_nod_h, s) 
		
		else	
			call set_id_row_matr(matr, nod) !~ Dirichlet values at the all four boundaries 
		end if	
				
	
	end do
	
	!~!~ Dirichlet values at the upper and lower boundaries 
	!~
	!~do nod = 1, nb_nod_h
	!~
	!~	call set_id_row_matr(matr, nod)
	!~	call set_id_row_matr(matr, geom%nb_vers-(nod-1))
	!~	
	!~end do
	!~
	!~! Dirichlet points left and right boundaries
	!~
	!~
	!~do i = 1, nb_nod_v-2
	!~	
	!~	nod= 1 + i*nb_nod_h
	!~	call set_id_row_matr(matr, nod)
	!~	nod= nb_nod_h + i*nb_nod_h
	!~	call set_id_row_matr(matr, nod)
	!~	
	!~end do
	
  end subroutine init_matr
  
  !-------------------------------
! Set matrix for the pressure navier stokes - Lap with Neumann bc
!-----------------------------------
  subroutine init_lap_matr(matr, geom, nu, nb_nod_h, nb_nod_v)

    implicit none

    type(matr_stor), intent(out) 		:: matr
    type(geom_struc), intent(in) 		:: geom
	real(rp), intent(in)				::  nu
    integer         , intent(in) 		:: nb_nod_v, nb_nod_h

    ! Local Declarations
    ! ------------------

    integer :: tt_nzc, nb_nzc, nb_lin   !tt_nzc total nonzero const nb_nzc: neighbor nonzero const
    integer ::  nod
	real(rp) :: c, e, w, n, s


    ! Beginning of the routine
    ! ------------------------

    ! Matrix allocation
    ! -----------------

    tt_nzc = geom%nb_elms &
           + 2 * (geom%nb_elms -1) &
           + 2 *  (geom%nb_elms-nb_nod_v)&
		   + 3* (2*nb_nod_v + 2*nb_nod_h)

    nb_nzc = 4

    nb_lin = nb_nod_v * nb_nod_h

    call alloc_matr(matr, tt_nzc, nb_nzc, nb_lin, nb_lin)

    matr%nod_zon => geom%ver_zon
	


    ! Matrix construction
    ! -------------------

	do nod = 1,geom%nb_vers	
			c = -4.0_rp*nu			
			e =  nu				
			w =  nu		
			n =  nu				
			s =  nu
		
		if (geom%ver_zon(nod) .eq. 0) then
	
			
			
			call fem_add_matr_coef(matr, nod, nod, c) 
			call fem_add_matr_coef(matr, nod, nod + 1, e) 
			call fem_add_matr_coef(matr, nod, nod - 1, w) 
			call fem_add_matr_coef(matr, nod, nod + nb_nod_h, n) 
			call fem_add_matr_coef(matr, nod, nod - nb_nod_h, s) 
		
		else	
			! call set_id_row_matr(matr, nod) !~ Dirichlet values at the all four boundaries 
			
			! For Neumann
			
			select case (geom%ver_zon(nod))
			
				case (1)
					call fem_add_matr_coef(matr, nod, nod, c) 
					call fem_add_matr_coef(matr, nod, nod + 1, e) 
					call fem_add_matr_coef(matr, nod, nod - 1, w) 
					call fem_add_matr_coef(matr, nod, nod + nb_nod_h, 2*n) 
					!call fem_add_matr_coef(matr, nod, nod - nb_nod_h, s) 
				case (2)
					call fem_add_matr_coef(matr, nod, nod, c) 
					!call fem_add_matr_coef(matr, nod, nod + 1, e) 
					call fem_add_matr_coef(matr, nod, nod - 1, 2*w) 
					call fem_add_matr_coef(matr, nod, nod + nb_nod_h, n) 
					call fem_add_matr_coef(matr, nod, nod - nb_nod_h, s) 
				case (3)
					call fem_add_matr_coef(matr, nod, nod, c) 
					call fem_add_matr_coef(matr, nod, nod + 1, e) 
					call fem_add_matr_coef(matr, nod, nod - 1, w) 
					!call fem_add_matr_coef(matr, nod, nod + nb_nod_h, n) 
					call fem_add_matr_coef(matr, nod, nod - nb_nod_h, 2*s) 
				case (4)
					call fem_add_matr_coef(matr, nod, nod, c) 
					call fem_add_matr_coef(matr, nod, nod + 1, 2*e) 
					!call fem_add_matr_coef(matr, nod, nod - 1, w) 
					call fem_add_matr_coef(matr, nod, nod + nb_nod_h, n) 
					call fem_add_matr_coef(matr, nod, nod - nb_nod_h, s) 
				case (5)
					call fem_add_matr_coef(matr, nod, nod, c) 
					call fem_add_matr_coef(matr, nod, nod + 1, 2*e) 
					!call fem_add_matr_coef(matr, nod, nod - 1, w) 
					call fem_add_matr_coef(matr, nod, nod + nb_nod_h, 2*n) 
					!call fem_add_matr_coef(matr, nod, nod - nb_nod_h, s) 
				case (6)
					call fem_add_matr_coef(matr, nod, nod, c) 
					!call fem_add_matr_coef(matr, nod, nod + 1, e) 
					call fem_add_matr_coef(matr, nod, nod - 1, 2*w) 
					call fem_add_matr_coef(matr, nod, nod + nb_nod_h, 2*n) 
					!call fem_add_matr_coef(matr, nod, nod - nb_nod_h, s) 
				case (7)
					call fem_add_matr_coef(matr, nod, nod, c) 
					!call fem_add_matr_coef(matr, nod, nod + 1, e) 
					call fem_add_matr_coef(matr, nod, nod - 1, 2*w) 
					!call fem_add_matr_coef(matr, nod, nod + nb_nod_h, n) 
					call fem_add_matr_coef(matr, nod, nod - nb_nod_h, 2*s) 
				case (8)
					call fem_add_matr_coef(matr, nod, nod, c) 
					call fem_add_matr_coef(matr, nod, nod + 1, 2*e) 
					!call fem_add_matr_coef(matr, nod, nod - 1, w) 
					!call fem_add_matr_coef(matr, nod, nod + nb_nod_h, n) 
					call fem_add_matr_coef(matr, nod, nod - nb_nod_h, s) 
				case default
			  
			end select
		end if	
		
	end do
	
	!~!~ Dirichlet values at the upper and lower boundaries 
	!~
	!~do nod = 1, nb_nod_h
	!~
	!~	call set_id_row_matr(matr, nod)
	!~	call set_id_row_matr(matr, geom%nb_vers-(nod-1))
	!~	
	!~end do
	!~
	!~! Dirichlet points left and right boundaries
	!~
	!~
	!~do i = 1, nb_nod_v-1
	!~	
	!~	nod= 1 + i*nb_nod_h
	!~	
	!~	nod= nb_nod_h + i*nb_nod_h
	!~	call set_id_row_matr(matr, nod)
	!~	
	!~end do
	
  end subroutine init_lap_matr

end module matrix_class
