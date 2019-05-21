module csr_stor_class

  use contxt_class

  type matr_stor

     integer :: nb_row			!number of row
     integer :: nb_col 			! number of columns

     integer :: tt_nzc

     integer , dimension(:), pointer :: nrow, ninf    	!k-th entry = index in Col and Val where the k-th row begins
     integer , dimension(:), pointer :: ncol			!column index
     real(rp), dimension(:), pointer :: coef, diag  	!coef off diagonal, diag on diagonal

     ! From the mesh

     integer, dimension(:), pointer :: nod_zon

  end type matr_stor
  
  type row
    
	integer :: l ! length of the row
	integer, dimension(:), pointer :: col
	real(rp), dimension(:), pointer :: coef
		
  end type row 
  
  type gauss_stor  ! for to do gaussian elimination
     
	integer :: nb_row			!number of row
    integer :: nb_col 			! number of columns
	real(rp), dimension(:), pointer:: diag
	type(row), dimension(:), allocatable :: row
  
  end type gauss_stor

contains


  subroutine csr_to_gauss_matr(matr,g_matr)
    
	implicit none
	
	type(matr_stor), intent(in), target :: matr
	type(gauss_stor), intent(out) :: g_matr
	
	! local var
	integer :: i, j
	
	g_matr%nb_row = matr%nb_row 
	g_matr%nb_col = matr%nb_col

	allocate(g_matr%diag(1:g_matr%nb_row))
	g_matr%diag = matr%diag
	
	
	allocate( g_matr%row(1:g_matr%nb_row) )

	do i= 1,g_matr%nb_row
	
		g_matr%row(i)%l = matr%nrow(i+1) - matr%nrow(i)  !plus 1 on the diagonal
		
		allocate(g_matr%row(i)%col(1:matr%nrow(i+1) - matr%nrow(i)))
		g_matr%row(i)%col = matr%ncol(matr%nrow(i) : matr%nrow(i+1)-1)
		
		allocate(g_matr%row(i)%coef(1:matr%nrow(i+1) - matr%nrow(i)))
		g_matr%row(i)%coef = matr%coef(matr%nrow(i) : matr%nrow(i+1)-1)

	end do
  
  end subroutine csr_to_gauss_matr
  
  subroutine gaussian(g_matr,rhs,sol,UnitNum,fname,NumLine)
    
	implicit none
	
	type(gauss_stor), intent(inout) 	:: g_matr
	real(rp), dimension(:), intent(inout)	:: rhs
	real(rp), dimension(:), intent(out) :: sol	
	integer, intent (in) :: UnitNum
	character (*), intent (in) :: fname
	integer, intent (out) :: NumLine
	
	! local var
	integer 							:: i, j, k, n, p
	integer, dimension(:), allocatable 	:: col
	real(rp), dimension(:), allocatable :: coef
	real(rp) 							:: c
	type(row)							:: rowi
	
	
	
	open(UnitNum, file = fname, status='replace',form='unformatted')
	NumLine = 0
	
	do i= 1,g_matr%nb_row-1
	
		do j=i+1,g_matr%nb_row
			
			if ((g_matr%row(j)%l .gt. 0)  .and. (g_matr%row(j)%col(1) .eq. i)) then
			
				!~write(*,*) 'Row i=', i
				!~write(*,*)
				!~call row_print(g_matr%row(i),g_matr%diag(i),i)
				!~write(*,*) 'Row j=', j
				!~write(*,*)
				!~call row_print(g_matr%row(j),g_matr%diag(j),j)	
				!~
			
				n = g_matr%row(i)%l + g_matr%row(j)%l -1 ! the first position of row j is zero 
				allocate(col(1:n)) 
				allocate(coef(1:n)) 
				
				c = g_matr%row(j)%coef(1)  / g_matr%diag(i)
				
				!write(*,*) "CONST c= ", c 
				
				! diagonal value of row j is traeated first
				k = 1  									! Find if row i has a column j.... 
				do while (k .le. g_matr%row(i)%l .and. g_matr%row(i)%col(k) .ne. j) 
					k = k+1	
				end do
				
				rowi%l = g_matr%row(i)%l                   !we'll make change on rowi to make sure g_matr%row(i) is unchanged
				allocate(rowi%col(1:g_matr%row(i)%l))
				allocate(rowi%coef(1:g_matr%row(i)%l))
				rowi%col = g_matr%row(i)%col
				rowi%coef = g_matr%row(i)%coef
			
				if (k .le. rowi%l) then 
					g_matr%diag(j) = g_matr%diag(j) - c*rowi%coef(k)
					
					do p = k, rowi%l -1   ! delete column j out of row i since the main diagonal is done
						rowi%col(p) = rowi%col(p+1)
						rowi%coef(p) = rowi%coef(p+1)
					end do	
					rowi%l = rowi%l -1
					n = n-1
				end if
				
				! treat values on on diagonal of row j
				call gauss_add_rows(rowi,g_matr%row(j),c,col,coef,n) 
				deallocate(rowi%col)
				deallocate(rowi%coef)
				
				g_matr%row(j)%l=n
							
				deallocate(g_matr%row(j)%col)       ! get rid of the old values of row j and replace by new values
				deallocate (g_matr%row(j)%coef)
				
				allocate(g_matr%row(j)%col(1:n))
				allocate(g_matr%row(j)%coef(1:n))
				
				g_matr%row(j)%col(1:n) = col(1:n)
				g_matr%row(j)%coef(1:n) = coef(1:n)
				
				deallocate(col)
				deallocate(coef)
				!~write(*,*)
				!~write(*,*) 'Row j new =', i, 'sum  row j=', j, ' Const c= ', c
				!~write(*,*)
				!~call row_print(g_matr%row(j),g_matr%diag(j),j)
				!~write(*,*) "NEW MATRX"
				!~call g_mt_print(" The Updated Matrix", g_matr)
				! rhs
				
				rhs(j)  = rhs(j) - c*rhs(i)
				
				write(UnitNum) i,j, c
				!write (*,*) " First time i,j,c  =    ", i,j,c
				NumLine = NumLine + 1
			end if
			
			!write(*,*) '(i, j)=', i,j
		end do
		
	end do
	
	close(UnitNum)
	! back sub

	call solv_gauss_Utrian(g_matr, rhs, sol)
  
  end subroutine gaussian
!-----------------------------------
! gaussian 2 takes Utriag resulted from gaussian and the row operation files generated from gausian, do back sub
!----------------------------------------  
  subroutine gaussian2(g_matr,rhs,sol,UnitNum,fname,NumLine)
    
	implicit none
	
	type(gauss_stor), intent(in) 	:: g_matr
	real(rp), dimension(:), intent(inout)	:: rhs
	real(rp), dimension(:), intent(out) :: sol	
	integer, intent (in) :: UnitNum
	integer, intent (in) :: NumLine
	character (*), intent (in) :: fname
	
	! local var
	integer 							:: i, j, k, n, p
	integer, dimension(:), allocatable 	:: col
	real(rp), dimension(:), allocatable :: coef
	real(rp) 							:: c

	
	
	
	open(UnitNum, file = fname,status='old',form='unformatted')
	
	do k= 1,NumLine
	
		read(UnitNum) i,j, c		
		!write (*,*) "i,j,c  =    ", i,j,c
	
		rhs(j)  = rhs(j) - c*rhs(i)
				
	end do
	
	close(UnitNum)
	! back sub

	call solv_gauss_Utrian(g_matr, rhs, sol)
  
  end subroutine gaussian2
  
  
   !------------------------------------------
  ! Solve Upper gauss matrix

  subroutine solv_gauss_Utrian(matr, rhs, sol)


    implicit none

    type(gauss_stor)       , intent(in ) :: matr
    real(rp), dimension(:), intent(in ) :: rhs
    real(rp), dimension(:), intent(out) :: sol

    ! Local Declarations
    ! ------------------

    integer :: i, k
	real(rp):: sum

    ! Beginning of the routine
    ! ------------------------
    sol(matr%nb_row) = rhs(matr%nb_row) / matr%diag(matr%nb_row)

	do i = matr%nb_row-1, 1, -1
		
		
		if (matr%row(i)%l .gt. 0 ) then
			!write (*,"(5a,4i5)") "COL= ",matr%row(i)%col
			sum = 0.0_rp
			do k = 1,matr%row(i)%l	
			  !write(*,*) ' bad i,k,row(i)%l =',i, k, matr%row(i)%l
				sum = sum + matr%row(i)%coef(k) * sol(matr%row(i)%col(k))
				
			end do
			sol(i) = (rhs(i) - sum) / matr%diag(i)
		else 
			sol(i) = rhs(i) / matr%diag(i)	
		end if	
	
    end do

  end subroutine solv_gauss_Utrian
  
  
  subroutine gauss_add_rows(row1,row2,c,col,coef,n) 
	implicit none
	
	type(row), intent(in)					:: row1
	type(row), intent(in)					:: row2
	real(rp), intent(in)					:: c
	integer, dimension(1:n), intent(out)	:: col
	real(rp), dimension(1:n), intent(out)	:: coef
	integer, intent(inout)					:: n
	
	! local var
	integer 							:: k, k1, k2

	if (row1%l .eq. 0) then
		
		col(1:n) = row2%col(2:row2%l)
		coef(1:n) = row2%coef(2:row2%l)
		
	else	
		if (row2%l .ge. 2) then
			
			k1 = 1
			k2 = 2
			k = 1
			do while ( (k1 .le. row1%l) .and. (k2 .le. row2%l) )
				if (row1%col(k1) == row2%col(k2)) then
					n = n -1
					col(k) = row1%col(k1)
					coef(k) = row2%coef(k2) - c*row1%coef(k1)
					k = k+1
					k1 = k1+1
					k2 = k2+1
				else
					if ( row1%col(k1) .lt. row2%col(k2) ) then
						col(k) = row1%col(k1)
						coef(k) = -c*row1%coef(k1)
						k = k+1
						k1 = k1 +1
					else
						col(k) = row2%col(k2)
						coef(k) = row2%coef(k2)
						k = k+1
						k2 = k2 +1
					end if
				end if
			end do
			if (k1 .gt. row1%l) then 
				do while (k2 .le. row2%l) 
					col(k) = row2%col(k2)
					coef(k) = row2%coef(k2)
					k = k+1
					k2 = k2 +1
				end do
			end if
			if (k2 .gt. row2%l) then 
				do while (k1 .le. row1%l) 
					col(k) = row1%col(k1)
					coef(k) = row1%coef(k1)
					k = k+1
					k1 = k1 +1
				end do
			end if
					
					
		else
		    
			col(1:n) = row1%col(1:row1%l)
			coef(1:n) = -c*row1%coef(1:row1%l)
		end if	
	end if	
			
  end subroutine gauss_add_rows
  
  
  
  subroutine dealloc_g_matr(matr)

    implicit none

    type(gauss_stor)     :: matr
    
	! Local Declarations
    ! ------------------
	integer :: i
    ! Beginning of the routine
    ! ------------------------


	do i =1, matr%nb_row
	
		deallocate(matr%row(i)%col)
		deallocate(matr%row(i)%coef)
		
    end do
	
    ! End of the routine
    ! ------------------

  end subroutine dealloc_g_matr
!---------------------------------
! product csr by vector  

  subroutine prod_matr_vect(matr, fact, resu)

    implicit none

    type(matr_stor)       , intent(in ) :: matr
    real(rp), dimension(:), intent(in ) :: fact
    real(rp), dimension(:), intent(out) :: resu

    ! Local Declarations
    ! ------------------

    integer :: i, j, k

    ! Beginning of the routine
    ! ------------------------

    resu = 0.0_rp

    do i = 1, matr%nb_row

       resu(i) = matr%diag(i) * fact(i)

       do k = matr%nrow(i), matr%nrow(i+1)-1
          j       =           matr%ncol(k)
          resu(i) = resu(i) + matr%coef(k) * fact(j)
       end do

    end do

    ! End of the routine
    ! ------------------

  end subroutine prod_matr_vect

!---------------------------------
! Solve lower triang csr matrix

  subroutine solv_Ltrian(matr, rhs, sol)


    implicit none

    type(matr_stor)       , intent(in ) :: matr
    real(rp), dimension(:), intent(in ) :: rhs
    real(rp), dimension(:), intent(out) :: sol

    ! Local Declarations
    ! ------------------

    integer :: i, j, k

    ! Beginning of the routine
    ! ------------------------

    sol = 0.0_rp

    do i = 1, matr%nb_row

       do k = 1, matr%nrow(i), matr%ninf(i)
          j      = matr%ncol(k)
          sol(i) = sol(i) + matr%coef(k) * sol(j)
       end do

       sol(i) = (rhs(i) - sol(i) ) / matr%diag(i) ! Diagonal

    end do

    ! End of the routine
    ! ------------------

  end subroutine solv_Ltrian
  
  !------------------------------------------
  ! Solve Upper csr matrix

  subroutine solv_Utrian(matr, rhs, sol)


    implicit none

    type(matr_stor)       , intent(in ) :: matr
    real(rp), dimension(:), intent(in ) :: rhs
    real(rp), dimension(:), intent(out) :: sol

    ! Local Declarations
    ! ------------------

    integer :: i, j, k

    ! Beginning of the routine
    ! ------------------------

    sol = 0.0_rp

    ! exercice ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ

    ! End of the routine
    ! ------------------

  end subroutine solv_Utrian
  
  !-------------------------------------
  ! Solve diag csr matrix

  subroutine solv_diag(matr, rhs, sol)


    implicit none

    type(matr_stor)       , intent(in ) :: matr
    real(rp), dimension(:), intent(in ) :: rhs
    real(rp), dimension(:), intent(out) :: sol

    ! Local Declarations
    ! ------------------

    ! Beginning of the routine
    ! ------------------------

    sol = rhs / matr%diag ! Diagonal

    ! End of the routine
    ! ------------------

  end subroutine solv_diag
  
  !-------------------------------
  ! allocate csr matrix
  
  subroutine alloc_matr(matr, tt_nzc, nb_nzc, nb_row, nb_col)

    implicit none

    type(matr_stor)     :: matr
    integer, intent(in) :: tt_nzc
    integer, intent(in) :: nb_nzc
    integer, intent(in) :: nb_row
    integer, intent(in) :: nb_col

    ! Local Declarations
    ! ------------------

    ! Beginning of the routine
    ! ------------------------

    matr%tt_nzc = tt_nzc
    matr%nb_row = nb_row
    matr%nb_col = nb_col

    allocate(matr%nrow(nb_row+1))
    allocate(matr%ninf(nb_row  ))
    allocate(matr%ncol(tt_nzc  ))
    allocate(matr%coef(tt_nzc  ))
    allocate(matr%diag(nb_row  ))

    matr%nrow =  1
    matr%ninf =  0
    matr%ncol =  0
    matr%coef = 0.0_rp
    matr%diag = 0.0_rp

    ! End of the routine
    ! ------------------

  end subroutine alloc_matr
  
  subroutine dealloc_matr(matr)

    implicit none

    type(matr_stor)     :: matr
    
	! Local Declarations
    ! ------------------

    ! Beginning of the routine
    ! ------------------------

    matr%tt_nzc = 0
    matr%nb_row = 0
    matr%nb_col = 0

    deallocate(matr%nrow)
    deallocate(matr%ninf)
    deallocate(matr%ncol)
    deallocate(matr%coef)
    deallocate(matr%diag)

    matr%nrow =  0
    matr%ninf =  0
    matr%ncol =  0
    matr%coef = 0.0_rp
    matr%diag = 0.0_rp

    ! End of the routine
    ! ------------------

  end subroutine dealloc_matr

!----------------------------------------
! add value to current value CSR matrix

  subroutine fem_add_matr_coef(matr, i, j, c)

!~  If the coefficient exists, add c to its value. Otherwise, make a new value c at row i column j

    implicit none

    type(matr_stor) :: matr
    integer , intent(in) :: i, j
    real(rp), intent(in) :: c

    ! Local Declarations
    ! ------------------

    integer :: k, l, mx
    logical :: okay

    ! Beginning of the routine
    ! ------------------------

    if ( i == j) then

       matr%diag(i) = matr%diag(i) + c

    else

       mx = matr%nrow(matr%nb_row+1) ! mx = value of current number of coefs of matr

       okay = .false.

     !   write(*,"(1x,36i7)") i, mx, matr%nrow(i), matr%nrow(i+1)

       do k = matr%nrow(i), matr%nrow(i+1)-1

          if (matr%ncol(k) == j) then

             write(*,"(1x,a,2i7)") &
                  "ADD-MATR-COEF: coefficient  exists (updated) : ", i, j

             matr%coef(k) = matr%coef(k) + c

             okay = .true.

             exit

          else if (matr%ncol(k) > j) then

             if (mx < matr%tt_nzc) then
!~				write(*,*) 'Case 1 i,j=', i,j
!~                write(*,"(1x,a4,3i5,/,(20x,6i15))") "AVO ", &
!~                     i, j, k, matr%ncol
				mx = mx + 1
                do l = mx, k+1, -1
                   matr%ncol(l) = matr%ncol(l-1)
                   matr%coef(l) = matr%coef(l-1)
                end do

                matr%ncol(k) = j
                matr%coef(k) = c

                matr%nrow(i+1:matr%nb_row+1) = matr%nrow(i+1:matr%nb_row+1) + 1
                

    !~            write(*,"(1x,a4,3i5,/,(20x,6i15))") "AVO ", &
    !~                 i, j, k, matr%ncol

                okay = .true.

                exit

             else
                write(*,*) "AA-MATR-COEF: Error too many coefficient ..."
                STOP       "AA-MATR-COEF: Error too many coefficient ..."
             end if

          end if

       end do

       if (.not. okay) then

          if (mx < matr%tt_nzc) then

             k = matr%nrow(i+1)
			!~ write(*,*) 
			!~write(*,*) 'Case 2 i,j=', i,j
   !~          write(*,"(1x,a4,3i5,/,(20x,6i15))") "AVN ", &
   !~               i, j, k, matr%ncol
	
             do l = mx, k+1, -1
                matr%ncol(l) = matr%ncol(l-1)
                matr%coef(l) = matr%coef(l-1)
             end do

             matr%ncol(k) = j
             matr%coef(k) = c

             matr%nrow(i+1:matr%nb_row+1) = matr%nrow(i+1:matr%nb_row+1) + 1

    !~         write(*,"(1x,a4,3i5,/,(20x,6i15))") "APN ", &
    !~              i, j, k, matr%ncol

          else
             write(*,*) "AA-MATR-COEF: Error too many coefficient ..."
             STOP       "AA-MATR-COEF: Error too many coefficient ..."
          end if

          !call mt_print(" XXXX", matr)

       end if

    end if

    ! End of the routine
    ! ------------------

  end subroutine fem_add_matr_coef

!---------------------------------------------
! set id row csr

  subroutine set_id_row_matr(matr, i)
	!this sets row i to be 0 ... 0 1 0 ...0 

    implicit none

    type(matr_stor) :: matr
    integer , intent(in) :: i
		
	!~ local declarations
	integer :: k
	real(rp), dimension(matr%tt_nzc-(matr%nrow(i+1)-matr%nrow(i))) :: coef
	
	
	matr%diag(i) = 1.0_rp
	
	if (matr%nrow(i+1)-matr%nrow(i) .gt. 0)	then
		do k = matr%nrow(i), matr%tt_nzc - (matr%nrow(i+1)-matr%nrow(i))
			matr%ncol(k) = matr%ncol(k + (matr%nrow(i+1)-matr%nrow(i)))  
			matr%coef(k) = matr%coef(k + (matr%nrow(i+1)-matr%nrow(i)))  
		end do

		matr%nrow(i+1:matr%nb_row+1) = matr%nrow(i+1:matr%nb_row+1) - (matr%nrow(i+1)-matr%nrow(i))
	end if	
	!~ matr%tt_nzc = matr%tt_nzc - (matr%nrow(i+1)-matr%nrow(i))
	
  end subroutine set_id_row_matr
  
!---------------------------------------------------
! add coefficient csr matrix

  subroutine add_matr_coef(matr, i, j, c)  
  
  !this only set new coefficent in the matrix, it does not add c to the current value if the entry exists.

    implicit none

    type(matr_stor) :: matr
    integer , intent(in) :: i, j
    real(rp), intent(in) :: c

    ! Local Declarations
    ! ------------------

    integer :: k, l, mx
    logical :: okay

    ! Beginning of the routine
    ! ------------------------

    if ( i == j) then

       matr%diag(i) = c

    else

       mx = matr%nrow(matr%nb_row+1)

       okay = .false.

       ! write(*,"(1x,36i7)") i, mx, matr%nrow(i), matr%nrow(i+1)

       do k = matr%nrow(i), matr%nrow(i+1)-1

          if (matr%ncol(k) == j) then

             write(*,"(1x,a,2i7)") &
                  "ADD-MATR-COEF: coefficient  exists (updated) : ", i, j

             matr%coef(k) = c

             okay = .true.

             exit

          else if (matr%ncol(k) > j) then

             if (mx < matr%tt_nzc) then

                !write(*,"(1x,a4,3i5,/,(20x,6i15))") "AVO ", &
                !     i, j, k, matr%ncol

                do l = mx, k+1, -1
                   matr%ncol(l) = matr%ncol(l-1)
                   matr%coef(l) = matr%coef(l-1)
                end do

                matr%ncol(k) = j
                matr%coef(k) = c

                matr%nrow(i+1:matr%nb_row+1) = matr%nrow(i+1:matr%nb_row+1) + 1
                mx = mx + 1

                !write(*,"(1x,a4,3i5,/,(20x,6i15))") "AVO ", &
                !     i, j, k, matr%ncol

                okay = .true.

                exit

             else
                write(*,*) "AA-MATR-COEF: Error too many coefficient ..."
                STOP       "AA-MATR-COEF: Error too many coefficient ..."
             end if

          end if

       end do

       if (.not. okay) then

          if (mx < matr%tt_nzc) then

             k = matr%nrow(i+1)

             !write(*,"(1x,a4,3i5,/,(20x,6i15))") "AVN ", &
             !     i, j, k, matr%ncol

             do l = mx, k+1, -1
                matr%ncol(l) = matr%ncol(l-1)
                matr%coef(l) = matr%coef(l-1)
             end do

             matr%ncol(k) = j
             matr%coef(k) = c

             matr%nrow(i+1:matr%nb_row+1) = matr%nrow(i+1:matr%nb_row+1) + 1

             !write(*,"(1x,a4,3i5,/,(20x,6i15))") "APN ", &
             !     i, j, k, matr%ncol

          else
             write(*,*) "AA-MATR-COEF: Error too many coefficient ..."
             STOP       "AA-MATR-COEF: Error too many coefficient ..."
          end if

          !call mt_print(" XXXX", matr)

       end if

    end if

    ! End of the routine
    ! ------------------

  end subroutine add_matr_coef
!-------------------------------------------------
!  Print csr matrix

  subroutine mt_print(text, matr)

    implicit none

    character(*)   , intent(in) :: text
    type(matr_stor), intent(in) :: matr

    integer :: i, k, n

    write(*,*) text
    
    write(*,*) " Number of rows                  : ", matr%nb_row
    write(*,*) " Number of columns               : ", matr%nb_col
    write(*,*) " Number of non zero coefficients : ", matr%tt_nzc

    do i = 1, matr%nb_row

       n = matr%nrow(i+1) - matr%nrow(i)

       write(*,"(1x, i7,2i4,36i15)") i, matr%ninf(i), n, i, &
            (matr%ncol(k), k = matr%nrow(i), matr%nrow(i+1)-1)

       write(*,"(1x, 7x,8x, 36es15.7)") matr%diag(i), &
            (matr%coef(k), k = matr%nrow(i), matr%nrow(i+1)-1)

    end do

  end subroutine mt_print

!---------------------------------------------------
! Print gauss_stor matrix

  subroutine g_mt_print(text, matr)

    implicit none

    character(*)   , intent(in) :: text
    type(gauss_stor), intent(in) :: matr

    integer :: i, k

    write(*,*) text
    
    write(*,*) " Number of rows                  : ", matr%nb_row
    write(*,*) " Number of columns               : ", matr%nb_col


    do i = 1, matr%nb_row
	
		call row_print(matr%row(i),matr%diag(i),i)
      
    end do

  end subroutine g_mt_print
  
  subroutine row_print(r,diag,i)
	implicit none

    type(row), intent(in) 	:: r
	integer, intent(in)		:: i
	real(rp), intent(in)	:: diag

    integer ::  k
  
	write(*,"(1x, 2i4,36i15)") i,  r%l, i, &
            (r%col(k), k = 1, r%l)

    write(*,"(1x, 8x, 36es15.7)") diag, &
			(r%coef(k), k = 1, r%l)
			
	end subroutine row_print

end module csr_stor_class
