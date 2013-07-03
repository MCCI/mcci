
module precision
          
  ! Specify the precision we want for real variables, which is the old double
  ! precision on most machines. The KIND value of this precision for this
  ! computer we save in dp. We get the needed number of decimal digits and
  ! range of exponents for double precision from the web page 
  !           http://www.nsc.liu.se/~boein/f77to90/c13.html

  implicit none            ! just say no
	
  integer, parameter :: dp = SELECTED_REAL_KIND(15,307) 

  ! This should be double precision. Selected_Real_Kind(Precision, Exponent)
  ! is a F90 intrinsic, that returns the kind of the real type that has
  ! at least this many decimal digits and at least this exponent range.
  
  ! dp is for the Double Precision. So, we should have (at least) 15 decimal
  ! digits of accuracy, and exponents from (at least) -307 to +307.

  ! Now choose which precision you want

  integer, parameter :: pr = dp
  
end module precision
	
