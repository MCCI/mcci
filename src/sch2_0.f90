subroutine sch2_0(e,i_am_mu,i_am_nu,idiff1,idiff2,kck,n_2p,ep)
  use commonarrays, only: list, nbft, e2ints, ipoint, i_sx2
  use dyn_par
  implicit real*8    (a-h,o-z)
  !common  /mxcoinc/ list(2,maxocc)
  !common  /aodat/   nsym, nbft, nbpsy(irmax)
  !common  /ints/    e1ints(max1), e2ints(max2)
  !common  /pointer/ ipoint(max2)
  !common  /occupy/  ntotal, n_alpha, n_beta, i_sx2

  iorb    = list(i_am_mu,idiff1)
  jorb    = list(i_am_mu,idiff2)
  korb    = list(i_am_nu,idiff1)
  lorb    = list(i_am_nu,idiff2)

  if(iorb.gt.nbft) iorb = iorb - nbft
  if(jorb.gt.nbft) jorb = jorb - nbft
  if(korb.gt.nbft) korb = korb - nbft
  if(lorb.gt.nbft) lorb = lorb - nbft

  if(iorb.ge.korb) then
     ik = ipoint(iorb) + korb
  else
     ik = ipoint(korb) + iorb
  endif

  if(jorb.ge.lorb) then
     jl = ipoint(jorb) + lorb
  else
     jl = ipoint(lorb) + jorb
  endif

  if(ik.ge.jl) then
     ikjl = ipoint(ik) + jl
  else
     ikjl = ipoint(jl) + ik
  endif

  sc = ck(i_sx2,n_2p,kck)

  e = ep*sc*e2ints(ikjl)

  return
end subroutine sch2_0
