module blas_class

  use contxt_class

!  use ell_stor_class ! OR csr_stor_class
  use csr_stor_class ! OR ell_stor_class

contains

  subroutine solv_jacobi(matr, rhs, sol, eps, itmax)

    implicit none

    type(matr_stor)       , intent(in   ) :: matr
    real(rp), dimension(:), intent(in   ) :: rhs
    real(rp), dimension(:), intent(inout) :: sol
    real(rp)              , intent(in   ) :: eps
    integer               , intent(in   ) :: itmax

    ! Local Declarations

    integer  :: it
    real(rp) :: sol_nrm, des_nrm, res_nr0, res_nrm, errel, erres

    real(rp), dimension(matr%nb_row) :: res, des

    ! Beginning of the routine
    ! ------------------------

    call prod_matr_vect(matr, sol, res)
    res = rhs - res
    call set_bound(matr, res)

    res_nr0 = dot_product(res, res)
    res_nrm = res_nr0

    do it = 1, itmax

       call solv_diag(matr, res, des)

       sol = sol + des

       des_nrm = dot_product(des, des)
       sol_nrm = dot_product(sol, sol)

       errel = des_nrm / sol_nrm
       erres = res_nrm / res_nr0

       write(*,"(' SOLV-JACOBI ', i5, 36e15.7)") &
            it, errel, erres, sol_nrm, des_nrm, res_nrm

       if (errel < eps .and. erres < eps) then
          exit
       end if

       call prod_matr_vect(matr, sol, res)
       res = rhs - res
    !~   call set_bound(matr, res)

       res_nrm = dot_product(res, res)

    end do

  end subroutine solv_jacobi

  subroutine set_bound(matr, vect)

    implicit none
    
    type(matr_stor), intent(in) :: matr
    real(rp), dimension(:)      :: vect

    where (matr%nod_zon > 0)
       vect = 0.0_rp
    end where

  end subroutine set_bound
  
  
    subroutine Cg_Sp_Mat(A, B, Solu, epsCG, i, rho)
    
    implicit none
    
    type(matr_stor),  			intent(in)     			:: A
    real(rp)                , dimension(:), intent(in ) :: B
  !  real(rp)                , dimension(:), intent(in ) :: X0
    real(rp)                , dimension(:), intent(out) :: Solu
	real(rp)                ,               intent(in)	:: epsCG
	!integer                 ,               intent(in)	:: tt_nzc
    real(rp)                ,               intent(out) :: rho
    integer                 ,               intent(out) :: i
     
    ! Declarations locales

    real(rp), dimension(A%nb_row) :: R, P, W
    real(rp)                   :: alpha, beta, norm2B, eps, rho2,dw
    integer                    :: NbIter
    
 !   call push_contxt("GC") 
    
    NbIter = A%tt_nzc

    
    ! write(uprint,*) prefix, NbIter

    !Solu = X0

    call prod_matr_vect(A, Solu, W)
    R = B - W
    P = R
  !  call pr_var("******** R_old:" , R_old)
  
    norm2B =dot_product(B, B)

    eps = epsCG**2

  !  call pr_var("******** norm2B:" , normB)
  !  call pr_var("******** NbPts:" , NbIter)

   ! rho = dot_product(R, R) 
       
    do i = 0, NbIter
	
		call prod_matr_vect(A, P, W)
		rho = dot_product(R, R) 
		dw = dot_product(P,W)
		
		alpha = rho/dw
		Solu = Solu  + alpha * P
		R = R - alpha * W
        rho2 = dot_product(R, R) 
		if (rho2 <= eps*norm2B) then
		rho = sqrt(rho)
	   EXIT
       end if
		beta  = rho2 / rho
		P= R + beta * P

        rho = rho2

        write(*,"(1x,i9,10e13.5)") i, rho, dw, alpha, beta

             
    end do

    if (i == NbIter) then
		write(*,*) " CG does not converge yet, res norm is :", &
            sqrt(rho)
		write(*,*) "         Nb iterations                 = ", i
		
    !~   call pr_var(" CG does not converge yet, res norm is :", &
    !~        sqrt(rho))
    !~   call pr_var("         Nb iterations                 = ", i)
    !~   call pop_contxt() 
      stop "CG des not converge yet..."
    end if

    
    !call pop_contxt()

  end subroutine Cg_Sp_Mat

end module blas_class
