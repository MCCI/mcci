subroutine swap(m,n,ep,i_am_mu,swapped)
  use commonarrays, only: list, nbft, my_pair
  use dyn_par
  implicit real*8   (a-h,o-z)
  logical           swapped
  !common  /aodat/   nsym,nbft,nbpsy(irmax)
  !common  /mxcoinc/ list(2,maxocc)
  !common  /pairs/   my_pair(2,maxocc)

  if(m.eq.n) STOP 'swap: can''t swap yourself!'

  mtemp              = list(i_am_mu,m)
  list(i_am_mu,m)    = list(i_am_mu,n)
  list(i_am_mu,n)    = mtemp

!     if switched a pair, pair index information is ok
  lm = list(i_am_mu,m)
  ln = list(i_am_mu,n)
  if(lm.gt.nbft) lm = lm - nbft
  if(ln.gt.nbft) ln = ln - nbft
  if(lm.ne.ln) then
   mtemp              = my_pair(i_am_mu,m)
   my_pair(i_am_mu,m) = my_pair(i_am_mu,n)
   my_pair(i_am_mu,n) = mtemp
!        pair's index should point to the new location
   if(my_pair(i_am_mu,m).ne.0) my_pair(i_am_mu,my_pair(i_am_mu,m)) = m
   if(my_pair(i_am_mu,n).ne.0) my_pair(i_am_mu,my_pair(i_am_mu,n)) = n
  endif

  ep = -ep

  swapped = .TRUE.

return 
end
