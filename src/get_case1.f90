subroutine get_case1(idiff1,i_am_mu,i_am_nu,my_case)
  use commonarrays, only: nbft, ntotal, list, my_pair
  use dyn_par
  implicit real*8   (a-h,o-z)
  !common  /aodat/   nsym, nbft, nbpsy(irmax)
  !common  /occupy/  ntotal, n_alpha, n_beta, i_sx2
  !common  /mxcoinc/ list(2,maxocc)
  !common  /pairs/   my_pair(2,maxocc)

!     my_case bit
!     -----------
!        0                  idiff1 doubly occupied in mu
!        1                  idiff1 occurs in nu, too
!        2                  idiff1 doubly occupied in nu 
!        3                  idiff1 occurs in mu, too

  if(idiff1.lt.0.or.idiff1.gt.ntotal)then
   call ldump(i_am_nu,i_am_mu,nu_doubly)
   stop 'get_case1: idiff1 out of range '
  endif
 
  my_case = 0

  if(my_pair(i_am_mu,idiff1).ne.0) my_case=ibset(my_case,0)
  if(my_pair(i_am_nu,idiff1).ne.0) my_case=ibset(my_case,2)

  mu = list(i_am_mu,idiff1)
  nu = list(i_am_nu,idiff1)

!     look for same spatial orbital in other configuration
  if(mu.gt.nbft) mu = mu-nbft
  if(nu.gt.nbft) nu = nu-nbft

  do n=1, ntotal
   if(list(i_am_nu,n).eq.mu .or. list(i_am_nu,n).eq.mu+nbft) my_case=ibset(my_case,1)
   if(list(i_am_mu,n).eq.nu .or. list(i_am_mu,n).eq.nu+nbft) my_case=ibset(my_case,3)
  enddo


!     following isn't necessary but makes code a little easier to read
  if(my_case.eq.0) then
   my_case = 0
  elseif(my_case.eq.12) then
   my_case = 1
  elseif(my_case.eq.15) then
   my_case = 2
  else
   write(*,*)'my_case=',my_case
   call ldump(i_am_nu,i_am_mu,nu_doubly)
   STOP 'my_case1: strange case indeed'
  endif 

  return 
end
