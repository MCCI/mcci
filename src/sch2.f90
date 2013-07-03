subroutine sch2(e,i_am_mu,i_am_nu,nu_doubly,idiff1,idiff2,kck,n_2p,ep)
  use dyn_par
  implicit real*8   (a-h,o-z)
  !common  /mxcoinc/ list(2,maxocc) !none of these blocks seem to be used!
  !common  /aodat/   nsym, nbft, nbpsy(irmax)
  !common  /occupy/  ntotal, n_alpha, n_beta, i_sx2

  !     determine what doubles Slater-Condon-Harris rules apply
  call get_case2(idiff1,idiff2,nu_doubly,i_am_mu,i_am_nu,&
       ep,kck,my_case)

  if(my_case.eq.0) then
     call sch2_0(e,i_am_mu,i_am_nu,idiff1,idiff2,kck,n_2p,ep)
  elseif(my_case.eq.1) then
     call sch2_1(e,i_am_mu,i_am_nu,idiff1,idiff2,kck,n_2p,ep)
  elseif(my_case.eq.2) then
     call sch2_2(e,i_am_mu,i_am_nu,idiff1,idiff2,kck,n_2p,ep)
  elseif(my_case.eq.3) then
     call sch2_3(e,i_am_mu,i_am_nu,&
          idiff1,idiff2,nu_doubly,kck,n_2p,ep)
  elseif(my_case.eq.4) then
     call sch2_4(e,i_am_mu,i_am_nu,idiff1,idiff2,kck,n_2p,ep)
  elseif(my_case.eq.5) then
     call sch2_5(e,i_am_mu,i_am_nu,idiff1,idiff2,kck,n_2p,ep)
  else
     write(*,*)'my_case',my_case
     call ldump(i_am_nu,i_am_mu,nu_doubly)
     STOP 'sch2: no conditions were met'
  endif

  return 
end subroutine sch2
