subroutine get_case2(idiff1,idiff2,nu_doubly,i_am_mu,i_am_nu,ep,kck,my_case)
  use commonarrays, only: nbft, ntotal, list, my_pair 
  use dyn_par
  implicit real*8   (a-h,o-z)
  logical           swapped
  !common  /aodat/   nsym, nbft, nbpsy(irmax)
  !common  /occupy/  ntotal, n_alpha, n_beta, i_sx2
  !common  /mxcoinc/ list(2,maxocc)
  !common  /schcase/ icase(16) !doesn't even seem to be used in this subroutine
  !common  /pairs/   my_pair(2,maxocc)

  !     my_case bit
  !     -----------
  !         0          idiff1 orbital doubly occupied in mu
  !         1          idiff1 orbital occurs in nu, too
  !         2          idiff2 orbital doubly occupied in mu 
  !         3          idiff2 orbital occurs in nu, too
  !         4          idiff1 orbital doubly occupied in nu
  !         5          idiff1 orbital occurs in mu, too   
  !         6          idiff2 orbital doubly occupied in nu 
  !         7          idiff2 orbital occurs in mu, too 

  if(idiff1.lt.0.or.idiff1.gt.ntotal .or. idiff2.lt.0.or.idiff2.gt.ntotal) then
     call ldump(i_am_nu,i_am_mu,nu_doubly)
     STOP 'get_case2: idiff1 or idiff2 = 0 '
  endif

  my_case = 0

  if(my_pair(i_am_mu,idiff1).ne.0) my_case=ibset(my_case,0)
  if(my_pair(i_am_mu,idiff2).ne.0) my_case=ibset(my_case,2)
  if(my_pair(i_am_nu,idiff1).ne.0) my_case=ibset(my_case,4)
  if(my_pair(i_am_nu,idiff2).ne.0) my_case=ibset(my_case,6)

  mu_1 = list(i_am_mu,idiff1)
  mu_2 = list(i_am_mu,idiff2)
  nu_1 = list(i_am_nu,idiff1)
  nu_2 = list(i_am_nu,idiff2)

  if(mu_1.gt.nbft) mu_1 = mu_1-nbft
  if(mu_2.gt.nbft) mu_2 = mu_2-nbft
  if(nu_1.gt.nbft) nu_1 = nu_1-nbft
  if(nu_2.gt.nbft) nu_2 = nu_2-nbft

  do n=1, ntotal
     if(list(i_am_nu,n).eq.mu_1 .or.&
          list(i_am_nu,n).eq.mu_1+nbft) my_case=ibset(my_case,1)
     if(list(i_am_nu,n).eq.mu_2 .or.&
          list(i_am_nu,n).eq.mu_2+nbft) my_case=ibset(my_case,3)
     if(list(i_am_mu,n).eq.nu_1 .or.&
          list(i_am_mu,n).eq.nu_1+nbft) my_case=ibset(my_case,5)
     if(list(i_am_mu,n).eq.nu_2 .or.&
          list(i_am_mu,n).eq.nu_2+nbft) my_case=ibset(my_case,7)
  enddo

  if(my_case.eq.0) then
     my_case = 2   ! 2   
  elseif(my_case.eq.80) then
     my_case = 0   ! 0   
  elseif(my_case.eq.85) then
     my_case = 0   ! 0  
  elseif(my_case.eq.178) then
     my_case = 3   ! 3  
  elseif(my_case.eq.48) then
     ispin = 1
     jspin = 1
     if(list(i_am_mu,idiff1).gt.nbft) ispin= 0
     if(list(i_am_mu,idiff2).gt.nbft) jspin= 0
     if(ispin.eq.jspin) then
        my_case = 1  ! 1 
     else
        my_case = 0  ! 0  
     endif
  elseif(my_case.eq.53) then
     my_case = 0     ! 0 
  elseif(my_case.eq.51) then
     ispin = 1
     jspin = 1
     if(list(i_am_mu,idiff1).gt.nbft) ispin = 0
     if(list(i_am_mu,idiff2).gt.nbft) jspin = 0 
     if(ispin.eq.jspin) then
        my_case = 4  ! 4 
     else
        my_case = 0  ! 0   
     endif
  elseif(my_case.eq.60) then
     ispin = 1
     jspin = 1
     if(list(i_am_mu,idiff1).gt.nbft) ispin = 0
     if(list(i_am_mu,idiff2).gt.nbft) jspin = 0 
     if(ispin.eq.jspin) then
        call swap(idiff1,idiff2,ep,i_am_mu,swapped)
        my_case = 4  ! 4
     else
        ispin = 1
        if(list(i_am_mu,my_pair(i_am_mu,idiff2)).gt.nbft) ispin = 0
        jspin = 1
        if(list(i_am_mu,idiff2).gt.nbft) jspin = 0
        kspin = 1
        if(list(i_am_nu,my_pair(i_am_mu,idiff2)).gt.nbft) kspin = 0
        lspin = 1
        if(list(i_am_nu,idiff2).gt.nbft) lspin = 0
        if( kspin.eq.lspin) then
           kck= kck
        elseif(ispin.eq.kspin.and.jspin.eq.lspin) then
           kck= kck+1
        elseif(ispin.eq.lspin.and.jspin.eq.kspin) then
           kck= kck-1
        endif
        call swap(idiff2,my_pair(i_am_mu,idiff2),ep,i_am_mu,swapped)
        call swap(idiff1,idiff2,ep,i_am_mu,swapped)
        my_case = 4  ! 4
     endif
  elseif(my_case.eq.240) then
     ispin = 1
     jspin = 1
     if(list(i_am_mu,idiff1).gt.nbft) ispin = 0
     if(list(i_am_mu,idiff2).gt.nbft) jspin = 0
     if(ispin.eq.jspin) then
        my_case = 1   ! 1
     else
        my_case = 0   ! 0
     endif
  elseif(my_case.eq.245) then
     my_case = 0      ! 0 
  elseif(my_case.eq.243.or.my_case.eq.252) then
     ispin = 1
     jspin = 1
     if(list(i_am_mu,idiff1).gt.nbft) ispin = 0
     if(list(i_am_mu,idiff2).gt.nbft) jspin = 0
     if(ispin.eq.jspin) then
        my_case = 1   ! 1 
     else
        my_case = 0   ! 0
     endif
  elseif(my_case.eq.255) then
     ispin = 1
     jspin = 1
     if(list(i_am_mu,idiff1).gt.nbft) ispin = 0
     if(list(i_am_mu,idiff2).gt.nbft) jspin = 0
     if(ispin.eq.jspin) then
        my_case = 1   ! 1 
     else
        my_case = 5   ! 5
     endif
  else
     write(*,*)'my_case=',my_case
     call ldump(i_am_nu,i_am_mu,nu_doubly)
     STOP 'get_case2: strange case indeed'
  endif

  return 
end subroutine get_case2
