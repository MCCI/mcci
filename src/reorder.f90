subroutine reorder(ici,jci,i_am_mu,i_am_nu,nu_doubly,ndiff,idiff1,idiff2,kck,n_2p,ep)
  use commonarrays, only: icij, nword, nbft, ntotal, n_alpha, n_beta, list, my_pair
  use dyn_par
  implicit real*8   (a-h,o-z)
  integer           idoubles(iword), jdoubles(iword)
  logical           spin1, spin2, swapped
  !common  /config/  icij(2,iword,maxc), nword
  !common  /aodat/   nsym, nbft, nbpsy(irmax)
  !common  /occupy/  ntotal, n_alpha, n_beta, i_sx2
  !common  /mxcoinc/ list(2,maxocc)
  !common  /pairs/   my_pair(2,maxocc)

  !     convert from "dictionary ordering" to "maximal coincidence ordering" 

  !     the references to configurations as mu or nu refers to
  !     Harris, J Chem Phys, v. 46, p. 2769 (1967)

  !     initialize list and pair lists to 0
  do n=1,ntotal
     list(1,n) = 0
     list(2,n) = 0
     my_pair(1,n) = 0
     my_pair(2,n) = 0
  enddo

  !     obtain lists of doubly occupied orbitals
  k = 0
  l = 0

  do n=1,nword
     idoubles(n) = iand(icij(1,n,ici),icij(2,n,ici))
     jdoubles(n) = iand(icij(1,n,jci),icij(2,n,jci))
  enddo

  do j=0,nbft-1

     n = j/int_bits +1
     jshift = j - (n-1)*int_bits

     if( btest( idoubles(n),jshift ) ) then
        k = k + 1
        list(1,k) = j+1
        my_pair(1,k) = k+1
        k = k+1
        list(1,k) = j+1 + nbft
        my_pair(1,k) = k-1
     endif

     if( btest( jdoubles(n),jshift ) ) then 
        l = l + 1
        list(2,l) = j+1
        my_pair(2,l) = l+1
        l = l+1
        list(2,l) = j+1 + nbft
        my_pair(2,l) = l-1
     endif

  enddo

  kdoubly = k
  ldoubly = l
  !     obtain lists of singly occupied orbitals
  do j=0, nbft-1

     n = j/int_bits + 1
     jshift = j - (n-1)*int_bits
     spin1   = btest(icij(1,n,ici), jshift)
     spin2   = btest(icij(2,n,ici), jshift)
     if( ( (.not.spin1) .and. spin2 ) .or. ( spin1 .and. (.not.spin2)) ) then 
        k = k+1
        if(spin1) list(1,k) = j+1
        if(spin2) list(1,k) = j+1 + nbft
     endif

     spin1   = btest(icij(1,n,jci), jshift)
     spin2   = btest(icij(2,n,jci), jshift)
     if( spin1 .NEQV. spin2 ) then 
        l = l+1
        if(spin1) list(2,l) = j+1
        if(spin2) list(2,l) = j+1 + nbft
     endif

  enddo

  !     these should simply be the number of occupied orbitals
  if(k.ne.ntotal .or. l.ne.ntotal) then
     write(*,*)'nttl nalpha nbeta',ntotal,n_alpha,n_beta
     write(*,*)'kdoubly ldoubly',kdoubly,ldoubly
     write(*,*)'k l',k,l
     call ldump(1,2,0)
     STOP 'error reorder: wrong particle number'
  endif

  if(kdoubly.gt.ldoubly) then
     i_am_nu= 1
     i_am_mu= 2
  elseif(ldoubly.gt.kdoubly) then
     i_am_nu= 2
     i_am_mu= 1
  elseif( kdoubly .eq. ldoubly) then

     !        in case of the tie, choose i_am_nu as follows

     !        singly occupied orbs in "1" occuring as doubly occupied orbs in "2"
     ks_ld = 0
     do k =kdoubly+1, ntotal
        do l= 1, ldoubly
           if(list(1,k) .eq. list(2,l) ) ks_ld = ks_ld + 1
        enddo
     enddo

     !        singly occupied orbs in "2" occuring as doubly occupied orbs in "1"
     ls_kd = 0
     do l =ldoubly+1, ntotal
        do k= 1, kdoubly
           if(list(2,l) .eq. list(1,k)) ls_kd = ls_kd + 1
        enddo
     enddo

     !        the number of doubly occupied in nu occuring as singles in mu .ge.
     !                                         mu                        nu
     if(ls_kd .ge. ks_ld) then
        i_am_nu = 1
        i_am_mu = 2
     elseif(ks_ld .gt. ls_kd) then
        i_am_nu = 2
        i_am_mu = 1
     endif

  endif

  if( i_am_nu .eq. 1) then
     nu_doubly = kdoubly
  elseif( i_am_nu .eq. 2) then
     nu_doubly = ldoubly
  endif

  !     at this point, orbitals in list are in their starting permutation-
  !     all exchanges must be counted after this point
  ep = 1.0

  !     first try to get best ordering while treating spin orbitals
  !     as distinct
  n = 0
  do while (n.lt.ntotal)
     n = n+1
     if(list(i_am_mu,n).ne.list(i_am_nu,n))then
        do m=1, ntotal
           if(n.ne.m) then
              if(list(i_am_mu,m).eq.list(i_am_nu,n))then
                 call swap(m,n,ep,i_am_mu,swapped)
                 !		     n = n - 1
              endif
           endif
        enddo
     endif
  enddo

  !     now try to get best ordering while considering space orbitals
  n = 0
  do while (n.lt.ntotal)
     n = n + 1
     if(list(i_am_mu,n) .ne. list(i_am_nu,n)) then
        do m=1,ntotal
           if(n.ne.m) then
              if(list(i_am_mu,m).ne.list(i_am_nu,m))then
                 l1 = list(i_am_mu,m)
                 !                    change to opposite spin orbital
                 if(l1.gt.nbft) then
                    l1 = l1-nbft
                 else
                    l1 = l1+nbft
                 endif
                 if(l1.eq.list(i_am_nu,n))then
                    call swap(m,n,ep,i_am_mu,swapped)
                    !			n = n-1
                 endif
              endif
           endif
        enddo
     endif
  enddo

  idiff1 = 0
  idiff2 = 0
  ndiff = 0
  do n=1,ntotal
     l1=list(i_am_nu,n)
     l2=list(i_am_mu,n)
     if(l1.gt.nbft) l1=l1-nbft
     if(l2.gt.nbft) l2=l2-nbft
     if(l1.ne.l2)then
        ndiff = ndiff + 1
        if(idiff1.eq.0) then
           idiff1 = n
        elseif(idiff1.ne.0) then
           idiff2 = n
        endif
     endif
  enddo

  !     if there is still a problem with a one orbital difference, fix it
  if(ndiff.eq.1 .and. idiff1.le.nu_doubly) then

     ispin = 1
     if(list(i_am_mu,idiff1).gt.nbft) ispin = 0

     kspin = 1
     if(list(i_am_nu,idiff1).gt.nbft) kspin = 0

     if(ispin.ne.kspin) then

        swapped = .FALSE.

        n = my_pair(i_am_mu,idiff1)
        if(n.ne.0 .and. n.gt.nu_doubly) then
           m = idiff1
           call swap(m,n,ep,i_am_mu,swapped)
        endif

        n=nu_doubly
        do while(n.lt.ntotal .and. .not. swapped)

           n = n + 1

           jspin = 1
           if(list(i_am_mu,n).gt.nbft) jspin = 0

           if(ispin.ne.jspin) then
              m = idiff1
              call swap(m,n,ep,i_am_mu,swapped)
           endif

        enddo
        if(.not.swapped) then
           call ldump(i_am_nu,i_am_mu,nu_doubly)
           STOP 'reorder: error with single difference'
        endif
     endif
  endif

  !     now let's look at ndiff = 2
  if(ndiff.eq.2 .and. idiff1.le.nu_doubly.and. idiff2.gt.nu_doubly) then

     !        if I am doubly occupied in mu, I should be idiff1
     if(my_pair(i_am_mu,idiff1).eq.0 .and.my_pair(i_am_mu,idiff2).ne.0) then
        n = idiff1
        m = idiff2
        call swap(m,n,ep,i_am_mu,swapped)
     endif

     ispin = 1
     if(list(i_am_mu,idiff1).gt.nbft) ispin = 0

     kspin = 1
     if(list(i_am_nu,idiff1).gt.nbft) kspin = 0

     if(ispin.ne.kspin) then

        swapped = .FALSE.

        n = my_pair(i_am_mu,idiff1)
        if(n.ne.0.and.n.gt.nu_doubly) then
           m = idiff1
           call swap(m,n,ep,i_am_mu,swapped)
        endif

        kspin = 1
        if(list(i_am_mu,idiff2).gt.nbft) kspin = 0

        if(ispin.ne.kspin .and. .not.swapped) then

           n= idiff1
           m= idiff2

           call swap(m,n,ep,i_am_mu,swapped)

        endif

        if(.not.swapped) then
           ndiff = 3
           return
        endif

     endif
  endif

  if(ndiff.eq.2 .and. idiff1.le.nu_doubly.and. idiff2.le.nu_doubly) then

     ispin = 1
     if(list(i_am_mu,idiff1).gt.nbft) ispin = 0

     jspin = 1
     if(list(i_am_mu,idiff2).gt.nbft) jspin = 0

     kspin = 1
     if(list(i_am_nu,idiff1).gt.nbft) kspin = 0

     lspin = 1
     if(list(i_am_nu,idiff2).gt.nbft) lspin = 0

     if(ispin.ne.kspin .and. jspin.eq.lspin) then
        if(my_pair(i_am_mu,idiff1).eq.0 .and.my_pair(i_am_mu,idiff2).ne.0) then
           n=idiff1
           m=idiff2
           call swap(m,n,ep,i_am_mu,swapped)
           itmp = ispin
           ispin = jspin
           jspin = itmp
        endif
     endif

     if(ispin.eq.kspin .and. jspin.ne.lspin) then
        if(my_pair(i_am_mu,idiff1).ne.0 .and.my_pair(i_am_mu,idiff2).eq.0) then
           n=idiff1
           m=idiff2
           call swap(m,n,ep,i_am_mu,swapped)
           itmp = ispin
           ispin = jspin
           jspin = itmp
        endif
     endif

     swapped = .FALSE.

     !        A everything spin matched
     if(ispin.eq.kspin .and. jspin.eq.lspin) swapped = .true. ! already okay

     !        B both differences aren't spin matched, and are different to each other
     if(ispin.ne.jspin .and.ispin.ne.kspin .and. jspin.ne.lspin .and. .not.swapped) then
        n= idiff1
        m= idiff2
        call swap(m,n,ep,i_am_mu,swapped)
     endif

     !        C both difference orbitals aren't spin matched, but can be swapped with their pairs 
     if(ispin.eq.jspin .and.ispin.ne.kspin .and. jspin.ne.lspin .and. .not.swapped) then

        n = my_pair(i_am_mu,idiff1)
        m = my_pair(i_am_mu,idiff2)

        if(n.ne.0 .and. n.gt.nu_doubly .and.m.ne.0 .and. m.gt.nu_doubly) then

           n= my_pair(i_am_mu,idiff1)
           m= idiff1
           call swap(m,n,ep,i_am_mu,swapped)

           n= my_pair(i_am_mu,idiff2)
           m= idiff2
           call swap(m,n,ep,i_am_mu,swapped)

        endif
     endif

     !        D one difference orbital isn't spin matched, but can be swapped with its pair
     if(ispin.ne.kspin .and. jspin.eq.lspin .and. .not.swapped) then ! ispin,jspin equal or not
        n = my_pair(i_am_mu,idiff1) 
        if(n.ne.0 .and. n.gt.nu_doubly) then 
           m = idiff1
           call swap(m,n,ep,i_am_mu,swapped)
        endif
     endif

     !        E same as D, but other difference orbital
     if(ispin.eq.kspin .and. jspin.ne.lspin .and. .not.swapped) then ! ispin,jspin equal or not
        n = my_pair(i_am_mu,idiff2)
        if(n.ne.0.and.n.gt.nu_doubly) then 
           m = idiff2
           call swap(m,n,ep,i_am_mu,swapped)
        endif
     endif

     if(.not.swapped) then
        ndiff = 3
        return
     endif

  endif


  idiff1 = 0
  idiff2 = 0
  ndiff = 0
  do n=1,ntotal
     l1=list(i_am_nu,n)
     l2=list(i_am_mu,n)
     if(l1.gt.nbft) l1=l1-nbft
     if(l2.gt.nbft) l2=l2-nbft
     if(l1.ne.l2)then
        ndiff = ndiff + 1
        if(idiff1.eq.0) then
           idiff1 = n
        elseif(idiff1.ne.0) then
           idiff2 = n
        endif
     endif
  enddo

  if(ndiff.ge.3) then
     write(*,*)'after third pass'
     write(*,*)'idiff1 idiff2 total#',idiff1,idiff2,ndiff
     call ldump(i_am_nu,i_am_mu,nu_doubly)
     STOP 'reorder: found a difference after ordering?'
  endif

  !     determine k value for spin mismatch
  kck = 0
  do n=nu_doubly+1,ntotal

     spin1   = .FALSE.
     spin2   = .FALSE.

     if( list(1,n) .le. nbft ) spin1 = .TRUE.
     if( list(2,n) .le. nbft ) spin2 = .TRUE.

     if( spin1 .NEQV. spin2 ) kck = kck + 1

  enddo

  if( mod(kck,2) .ne. 0 .and. test) then
     call ldump(i_am_nu,i_am_mu,nu_doubly)
     STOP 'reorder: error in spin exchanges'
  endif

  kck = kck/2       ! counting number of exchanges to pair like spins

  n_2p = ntotal - nu_doubly   ! number of singly occupied orbs. in nu

  return 
end subroutine reorder
