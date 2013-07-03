subroutine sch2_3(e,i_am_mu,i_am_nu,idiff1,idiff2,nu_doubly,kck,n_2p,ep)
  use commonarrays, only: list, nbft, e2ints, ntotal, i_sx2, ipoint
  use dyn_par
  implicit real*8    (a-h,o-z)
  !common  /aodat/   nsym, nbft, nbpsy(irmax)
  !common  /mxcoinc/ list(2,maxocc)
  !common  /ints/    e1ints(max1), e2ints(max2)
  !common  /occupy/  ntotal, n_alpha, n_beta, i_sx2
  !common  /pointer/ ipoint(max2)

  !     two orbital difference Slater Condon Harris rules
  !     (one orbital difference forced into a two orbital difference
  !     by spin ordering constraints)

  iorb    = list(i_am_mu,idiff1)
  jorb    = list(i_am_mu,idiff2)
  korb    = list(i_am_nu,idiff1)
  lorb    = list(i_am_nu,idiff2)

  if(iorb.gt.nbft) iorb = iorb - nbft
  if(jorb.gt.nbft) jorb = jorb - nbft
  if(korb.gt.nbft) korb = korb - nbft
  if(lorb.gt.nbft) lorb = lorb - nbft

  if(iorb .ge. korb) then
     ik = ipoint(iorb) + korb
  else
     ik = ipoint(korb) + iorb
  endif

  if(jorb .ge. lorb) then
     jl = ipoint(jorb) + lorb
  else	
     jl = ipoint(lorb) + jorb
  endif

  if(ik.ge.jl) then
     ikjl = ipoint(ik) + jl
  else 
     ikjl = ipoint(jl) + ik
  endif

  sc= ck(i_sx2,n_2p,kck)
  e = sc*e2ints(ikjl)

  ispin = 1
  if(list(i_am_mu,idiff2).gt.nbft) ispin = 0

  kspin = 1
  if(list(i_am_nu,idiff2).gt.nbft) kspin = 0

  idiff1_spin = 1
  if(list(i_am_mu,idiff1).gt.nbft) idiff1_spin = 0

  do n=nu_doubly+1, ntotal

     iorb = list(i_am_nu,n)
     if(iorb.gt.nbft) iorb = iorb-nbft

     jspin =1
     if(list(i_am_mu,n).gt.nbft) jspin = 0

     lspin =1
     if(list(i_am_nu,n).gt.nbft) lspin = 0

     !        sum here is spin restricted based upon spin in mu
     if(jspin.eq.idiff1_spin) then

        if( iorb .ge. korb) then
           ik = ipoint(iorb) + korb
        else
           ik = ipoint(korb) + iorb
        endif

        if( jorb .ge. iorb) then
           ji = ipoint(jorb) + iorb
        else	
           ji = ipoint(iorb) + jorb
        endif

        if(ik.ge.ji) then
           ikji = ipoint(ik) + ji
        else 
           ikji = ipoint(ji) + ik
        endif

        if( kspin.eq.lspin ) then
           sc= ck(i_sx2,n_2p,kck)
        elseif(ispin.eq.kspin.and.jspin.eq.lspin) then
           sc= ck(i_sx2,n_2p,kck+1)
        elseif(ispin.eq.lspin.and.jspin.eq.kspin) then
           sc= ck(i_sx2,n_2p,kck-1)
        endif

        e = e + sc*e2ints(ikji)

     endif
  enddo

  e =  ep*e

  return
end subroutine sch2_3
