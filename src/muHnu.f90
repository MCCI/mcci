subroutine muHnu(ici,jci,e,n_2p,s_overlap)
  use commonarrays, only: nbft, icij, nword, i_sx2
  use dyn_par
  implicit real*8    (a-h,o-z)
  integer           ijsingles(iword), ijdoubles(iword),ij_sd(iword)
  !common  /aodat/   nsym, nbft, nbpsy(irmax)
  !common  /config/  icij(2,iword,maxc), nword
  !common  /occupy/  ntotal, n_alpha, n_beta, i_sx2

  e = 0.0
  s_overlap = 0.0

  !     find number of differences and if matrix element = 0, return
  do n=1,nword 

     isingles = ieor(icij(1,n,ici),icij(2,n,ici))
     jsingles = ieor(icij(1,n,jci),icij(2,n,jci))

     idoubles = iand(icij(1,n,ici),icij(2,n,ici))
     jdoubles = iand(icij(1,n,jci),icij(2,n,jci))

     !        all singles
     ijsingles(n) = ieor(isingles,jsingles)
     !        doubles only in one configuration
     ijdoubles(n) = ieor(idoubles,jdoubles)
     ij_sd(n)     = iand(ijsingles(n),ijdoubles(n))

  enddo

  !     sum number of difference orbitals
  ndiff = 0
  do j = 0, nbft-1
     n = j/int_bits+1
     jshift = j - (n-1)*int_bits
     itest = 0
     if(btest(ijdoubles(n),jshift)) itest = 1
     ndiff = ndiff+2*itest
     itest = 0
     if(btest(ijsingles(n),jshift)) itest = 1
     ndiff = ndiff+  itest
     itest = 0
     if(btest(ij_sd(n),jshift)) itest = 1
     ndiff = ndiff-2*itest
  enddo
  if(ndiff.gt.4) return                     ! matrix element = 0

  !     note ndiff can change due to restrictions on the orderings,
  !     in any event, 4 --> 2, 2 --> 1 , 0 --> 0
  call reorder(ici,jci,i_am_mu,i_am_nu,nu_doubly,ndiff,idiff1,idiff2,kck,n_2p,ep)

  if(ndiff.gt.2) return

  !     Slater Condon Harris rules
  if( ndiff .eq. 0) then
     s_overlap = ck(i_sx2,n_2p,kck)
     call sch0(e,i_am_mu,i_am_nu,nu_doubly,kck,n_2p,ep)       ! no orbital differences
  elseif( ndiff .eq. 1) then
     call sch1(e,i_am_mu,i_am_nu,nu_doubly,idiff1,kck,n_2p,ep) 		! one orbital difference
  elseif( ndiff .eq. 2) then
     call sch2(e,i_am_mu,i_am_nu,nu_doubly,idiff1,idiff2,kck,n_2p,ep) 		! two orbital differences
  else
     write(*,*) 'ndiff', ndiff
     STOP 'muHnu: no conditions were met'
  endif

      return
      end
