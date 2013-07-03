subroutine doubles(iprev,ilength,ispin,seed)
  use commonarrays, only: icij,nword,ntotal,n_alpha,n_beta,i_sx2,&
                          nsym,nbft,nbpsy,irrep,nfreeze,ifreeze,nactive,iactive
  use dyn_par
  implicit real*8   (a-h,o-z)
  !common  /config/  icij(2,iword,maxc), nword
  !common  /occupy/  ntotal, n_alpha, n_beta, i_sx2
  !common  /aodat/   nsym, nbft, nbpsy(irmax)
  !common  /sym/     irrep(0:irmax-1)
  !common  /froze/   nfreeze, ifreeze(maxocc)
  !common  /active/  nactive, iactive(maxocc)

  nhalf= 0
  if(ispin.eq.1) then
     nhalf= n_alpha
  elseif(ispin.eq.2) then
     nhalf= n_beta
  elseif(ispin.eq.3) then
     nhalf1= n_alpha
     nhalf2= n_beta
  endif

  if(ispin.eq.1.or.ispin.eq.2) then

     !        randomly pick two orbitals from the current configuration
     it = 0
10   call random(seed,rand)
     iorb1 = int(rand*(nhalf-nfreeze))+1
20   call random(seed,rand)
     iorb2 = int(rand*(nhalf-nfreeze))+1
     if(iorb1.le.0.or.iorb1.gt.(nhalf-nfreeze)) then
        write(*,*) 'iorb1',iorb1
        STOP 'Error doubles: random 1'
     endif
     if(iorb2.le.0.or.iorb2.gt.(nhalf-nfreeze)) then
        write(*,*)'iorb2',iorb2
        STOP 'Error doubles: random 1'
     endif

     it = it + 1
     if(it.gt.500) return ! got stuck here

     if(iorb1.eq.iorb2) then
        goto 20 
     endif

     icount = 0
     do j=1, nactive
        ja=iactive(j)-1
        n = ja/int_bits + 1
        jshift = ja - (n-1)*int_bits
        if( btest(icij(ispin,n,ilength),jshift) ) icount=icount+1	! BW
        if (icount.eq.iorb1) then
           j1shift = jshift
           j1 = ja
           n1 = n
           iorb1 =-1 ! reset iorb1 to avoid changing j1 on next loop
        elseif (icount.eq.iorb2) then
           j2shift = jshift
           j2 = ja
           n2 = n
           iorb2 =-1 ! reset iorb2 to avoid changing j2 on next loop
        endif
     enddo

     isym1 =0
     isym2=0
     !        orbitals j1 and j2 belong to which symmetry groups?
     ilo  = 0
     do k=1, nsym
        ihi = ilo + nbpsy(k)-1
        if(j1.ge.ilo .and. j1.le.ihi) isym1=k-1
        if(j2.ge.ilo .and. j2.le.ihi) isym2=k-1
        ilo = ilo + nbpsy(k)
     enddo

     irprod = ieor( irrep(isym1), irrep(isym2) )		! BW

     !        randomly generate two orbitals with the same symmetry product 
     it = 0
30   call random(seed,rand)

     j3 = int(rand*nactive)+1
     if(j3.le.0 .or. j3.gt. nactive) STOP 'Error in doubles selection'
     j3 = iactive(j3)-1

     it = it+1
     if(it.gt.500) return ! got stuck here

     !        is/was this orbital occupied in the previous configuration?
     n3 = j3/int_bits+1
     j3shift = j3 -(n3-1)*int_bits
     if(btest(icij(ispin,n3,iprev),j3shift)) goto 30      	! BW

40   call random(seed,rand)
     j4 = int(rand*nactive)+1
     if(j4.le.0 .or. j4.gt. nactive) STOP 'Error in doubles selection'
     j4 = iactive(j4)-1

     it = it+1
     if(it.gt.500) return ! got stuck here

     if(j3.eq.j4) goto 40   

     !        is/was this orbital occupied in the previous configuration?
     n4 = j4/int_bits+1
     j4shift = j4 -(n4-1)*int_bits
     if(btest(icij(ispin,n4,iprev),j4shift)) goto 40	! BW

     isym3=0
     isym4=0
     !        orbitals j3 and j4 belong to which symmetry groups?
     ilo  = 0
     do k=1, nsym
        ihi = ilo + nbpsy(k)-1
        if(j3.ge.ilo .and. j3.le.ihi) isym3=k-1
        if(j4.ge.ilo .and. j4.le.ihi) isym4=k-1
        ilo = ilo + nbpsy(k)
     enddo

     irx = ieor( irrep(isym3), irrep(isym4) )			! BW

     if( irprod .ne. irx) goto 30

     !        detach electrons on selected orbitals
     icij(ispin,n1,ilength)=ibclr(icij(ispin,n1,ilength),j1shift)	! BW
     icij(ispin,n2,ilength)=ibclr(icij(ispin,n2,ilength),j2shift)	! BW

     if(test) then
        kold=0
        do k=0, nbft-1
           n=k/int_bits+1
           kshift=k-(n-1)*int_bits
           if( btest(icij(ispin,n,ilength),kshift) ) kold=kold+1	! BW
        enddo
        if(kold.ne.nhalf-2) then
           write(*,*)'j1 j2',j1,j2
           STOP 'Orbits not deleted correctly in doubles'
        endif
     endif

     !        attach electrons on selected orbitals
     icij(ispin,n3,ilength)=ibset(icij(ispin,n3,ilength),j3shift) 	! BW
     icij(ispin,n4,ilength)=ibset(icij(ispin,n4,ilength),j4shift) 	! BW

     if(test) then  
        lnew = 0
        do l=0, nbft-1
           n = l/int_bits + 1
           lshift = l- (n-1)*int_bits
           if(btest(icij(ispin,n,ilength),lshift))lnew=lnew+1	! BW
        enddo
        if(lnew .ne. nhalf) then
           write(*,*)'New nhalf ',lnew,' old nhalf ',nhalf
           write(*,*)'Substituted to orbitals ',j3,j4
           STOP 'Changed occupation in doubles: parallel subst'
        endif
     endif

     !        substitution complete, let's get out of here
     return

  elseif(ispin.eq.3) then

     if(test) then
        iold1=0
        iold2=0
        do i=0, nbft-1
           n = i/int_bits+1
           ishift = i-(n-1)*int_bits
           if( btest(icij(1,n,ilength),ishift) ) iold1=iold1+1	! BW
           if( btest(icij(2,n,ilength),ishift) ) iold2=iold2+1	! BW
        enddo
     endif

     !        randomly pick two orbitals from the current configurations
     it=0
50   call random(seed,rand)
     iorb1 = int(rand*(nhalf1-nfreeze)) + 1

60   call random(seed,rand)
     iorb2 = int(rand*(nhalf2-nfreeze)) + 1

     if(iorb1.le.0.or.iorb1.gt.(nhalf1-nfreeze)) then
        write(*,*) 'iorb1', iorb1
        STOP 'Error doubles: random 2'
     endif
     if(iorb2.le.0.or.iorb2.gt.(nhalf2-nfreeze)) then
        write(*,*) 'iorb2', iorb2
        STOP 'Error doubles: random 2'
     endif

     icount1= 0
     icount2= 0 
     do j=1, nactive
        ja=iactive(j)-1
        n = ja/int_bits + 1
        jshift = ja - (n-1)*int_bits
        if( btest(icij(1,n,ilength),jshift) ) then
           icount1=icount1+1
           if (icount1.eq.iorb1) then
              j1 = ja
              j1shift = jshift
              n1 = n
           endif
        endif
        if( btest(icij(2,n,ilength),jshift) ) then
           icount2=icount2+1
           if (icount2.eq.iorb2) then
              j2 = ja
              j2shift = jshift
              n2 = n
           endif
        endif
     enddo

     !        orbitals j1 and j2 belong to which irrep?
     ilo  = 0
     do k=1, nsym
        ihi = ilo + nbpsy(k)-1
        if(j1.ge.ilo .and. j1.le.ihi) isym1=k-1
        if(j2.ge.ilo .and. j2.le.ihi) isym2=k-1
        ilo = ilo + nbpsy(k)
     enddo

     irprod = ieor( irrep(isym1), irrep(isym2) )		! BW

     !        randomly generate two orbitals with the same symmetry product 
70   call random(seed,rand)
     j3 = int(rand*nactive)+1
     if(j3.le.0 .or. j3.gt. nactive) STOP 'Error in doubles selection'
     j3=iactive(j3)-1

     !        is/was this orbital occupied in the previous configuration?
     n3 = j3/int_bits+1
     j3shift = j3 -(n3-1)*int_bits

     it=it+1
     if(it.gt.500) return ! got stuck here

     if(btest(icij(1,n3,iprev),j3shift)) goto 70		! BW

80   call random(seed,rand)
     j4 = int(rand*nactive)+1
     if(j4.le.0 .or. j4.gt. nactive) STOP 'Error in doubles selection'
     j4=iactive(j4)-1

     !        is/was this orbital occupied in the previous configuration?
     n4 = j4/int_bits+1
     j4shift = j4 -(n4-1)*int_bits

     it=it+1
     if(it.gt.500) return ! got stuck here

     if(btest(icij(2,n4,iprev),j4shift)) goto 80		! BW

     !        orbitals j3 and j4 belong to which irrep groups?
     ilo  = 0
     do k=1, nsym
        ihi = ilo + nbpsy(k)-1
        if(j3.ge.ilo .and. j3.le.ihi) isym3=k-1
        if(j4.ge.ilo .and. j4.le.ihi) isym4=k-1
        ilo = ilo + nbpsy(k)
     enddo

     irx = ieor( irrep(isym3), irrep(isym4) )		! BW

     if( irprod .ne. irx) goto 70

     !        detach electrons on selected orbitals
     icij(1,n1,ilength)=ibclr(icij(1,n1,ilength),j1shift)	! BW
     icij(2,n2,ilength)=ibclr(icij(2,n2,ilength),j2shift)	! BW

     if(test) then
        kold1=0
        kold2=0
        do k=0, nbft-1
           n = k/int_bits + 1
           kshift = k - (n-1)*int_bits
           if( btest(icij(1,n,ilength),kshift) ) kold1=kold1+1	! BW
           if( btest(icij(2,n,ilength),kshift) ) kold2=kold2+1	! BW
        enddo
        if(kold1.ne.nhalf1-1.or.kold2.ne.nhalf2-1) then
           write(*,*)'k1 k2 nhalf1 nhalf2',kold1,kold2,nhalf1,nhalf2
           write(*,*)'n1 j1 j1shift',n1,j1,j1shift
           write(*,*)'n2 j2 j2shift',n2,j2,j2shift
           STOP 'Orbits not deleted correctly in doubles'
        endif
     endif

     !        attach electrons on selected orbitals
     icij(1,n3,ilength)=ibset(icij(1,n3,ilength),j3shift) 		! BW
     icij(2,n4,ilength)=ibset(icij(2,n4,ilength),j4shift) 		! BW

     if(test) then  
        lnew1 = 0
        lnew2 = 0
        do l=0, nbft-1
           n=l/int_bits + 1
           lshift = l - (n-1)*int_bits
           if(btest(icij(1,n,ilength),lshift)) lnew1=lnew1+1	! BW
           if(btest(icij(2,n,ilength),lshift)) lnew2=lnew2+1	! BW
        enddo
        if(lnew1 .ne. nhalf1 .or. lnew2 .ne. nhalf2) then
           write(*,*)'New nhalf ',lnew1,lnew2,' old nhalf ',iold1,iold2
           write(*,*)'Substituted to orbitals ',j3,j4
           STOP 'Changed occupation in doubles: anti-parallel subst'
        endif
     endif

     !       substitution complete, let's get out of here
     return

  endif

  STOP 'Did not make a double substitution'

end subroutine doubles
