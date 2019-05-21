module function_class

  use contxt_class

  

contains

  real(rp) function f(x,y)
	
	implicit none
	
	real(rp), intent(in) :: x,y
  
	f =  ( 400*(x**2+y**2) - 40 ) * exp(-10*(x**2+y**2)) 
	
  end function f
  
 !------------------------------------------ 
  real(rp) function gu_D(x,y)
	
	implicit none
	
	real(rp), intent(in) :: x,y
	
	if (y .eq. 1.0_rp) then
		gu_D = 1.0_rp
	else 
		gu_D = 0.0_rp
	end if	
  
	!gu_D =   exp(-10*(x**2+y**2))

  end function gu_D
 !------------------------------------------ 
  real(rp) function gx_N(x,y)
	
	implicit none
	
	real(rp), intent(in) :: x,y
  
	gx_N =  -20*x*exp(-10*(x**2+y**2))

  end function gx_N
  !------------------------------------------
  real(rp) function gy_N(x,y)
	
	implicit none
	
	real(rp), intent(in) :: x,y
  
	gy_N =  -20*y*exp(-10*(x**2+y**2))

  end function gy_N
  
 
  
end module function_class
