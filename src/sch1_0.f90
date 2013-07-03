subroutine sch1_0(e,i_am_mu,i_am_nu,nu_doubly,idiff1,kck,n_2p,ep)
  use commonarrays, only: nbft, e1ints, e2ints, ntotal, i_sx2, ipoint, list
  use dyn_par
  implicit real*8    (a-h,o-z)
  !common  /aodat/   nsym, nbft, nbpsy(irmax)
  !common  /ints/    e1ints(max1), e2ints(max2)
  !common  /occupy/  ntotal, n_alpha, n_beta, i_sx2
  !common  /pointer/ ipoint(max2)
  !common  /mxcoinc/ list(2,maxocc)

  !     one orbital difference Slater Condon Harris rules

  jorb = list(i_am_mu,idiff1)
  lorb = list(i_am_nu,idiff1)

  jspin = 1
  lspin = 1

  if(jorb.gt.nbft) then
     jorb  = jorb - nbft
     jspin = 0
  endif

  if(lorb.gt.nbft) then
     lorb  = lorb - nbft
     lspin = 0
  endif

  if(jorb .ge. lorb) then
     jl = ipoint(jorb) + lorb
  else
     jl = ipoint(lorb) + jorb
  endif

  e = e1ints(jl)

  do n = 2, nu_doubly, 2

     iorb = list(i_am_nu,n)
     if(iorb.gt.nbft) iorb = iorb-nbft

     ii = ipoint(iorb) + iorb

     if(ii.ge.jl) then
        iijl = ipoint(ii) + jl
     else 
        iijl = ipoint(jl) + ii
     endif

     e = e + 2.0*e2ints(iijl)          ! "Coulomb terms"

     if(iorb .ge. lorb) then
        il = ipoint(iorb) + lorb
     else
        il = ipoint(lorb) + iorb
     endif

     if(jorb .ge. iorb) then
        ji = ipoint(jorb) + iorb
     else	
        ji = ipoint(iorb) + jorb
     endif

     if( il .ge. ji ) then
        ilji = ipoint(il) + ji
     else
        ilji = ipoint(ji) + il
     endif

     e = e - e2ints(ilji)        ! "Exchange terms"

  enddo

  sc1 = ck(i_sx2,n_2p,kck)

  e = sc1*e

  do n=nu_doubly+1, ntotal

     if(n.ne.idiff1) then

        ispin = 1
        iorb = list(i_am_mu,n)
        if(iorb.gt.nbft) then
           ispin = 0
           iorb  = iorb - nbft
        endif

        kspin = 1
        if(list(i_am_nu,n).gt.nbft) kspin = 0

        ii = ipoint(iorb) + iorb

        if(ii.ge.jl) then
           iijl = ipoint(ii) + jl
        else 
           iijl = ipoint(jl) + ii
        endif

        e = e + sc1*e2ints(iijl)          ! "Coulomb terms"

        if( iorb .ge. lorb) then
           il = ipoint(iorb) + lorb
        else	
           il = ipoint(lorb) + iorb
        endif

        if( jorb .ge. iorb) then
           ji = ipoint(jorb) + iorb
        else
           ji = ipoint(iorb) + jorb
        endif

        if(il.ge.ji) then
           ilji = ipoint(il) + ji
        else 
           ilji = ipoint(ji) + il
        endif

        !           need a new Sanibel coefficient
        if( ispin.eq.jspin .or. kspin.eq.lspin ) then
           sc2 = sc1
        elseif(ispin.eq.kspin .and. jspin.eq.lspin) then		
           sc2 =  ck(i_sx2,n_2p,kck+1)
        elseif(ispin.eq.lspin .and. jspin.eq.kspin) then		
           sc2 =  ck(i_sx2,n_2p,kck-1)
        endif

        e = e - sc2*e2ints(ilji) 

     endif

  enddo

  e = ep*e

  return
end subroutine sch1_0
