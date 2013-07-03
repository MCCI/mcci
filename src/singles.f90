subroutine singles(ilength,ispin,seed)
  use commonarrays, only: icij, n_alpha, n_beta, nsym, nbft, nbpsy, nfreeze, nactive, iactive
  use dyn_par
  implicit real*8   (a-h,o-z)
  logical            btest
  !common  /config/  icij(2,iword,maxc), nword
  !common  /occupy/  ntotal, n_alpha, n_beta, i_sx2
  !common  /aodat/   nsym, nbft, nbpsy(irmax)
  !common  /froze/   nfreeze, ifreeze(maxocc)
  !common  /active/  nactive, iactive(maxocc)

  if(ispin.eq.1) nhalf = n_alpha
  if(ispin.eq.2) nhalf = n_beta

  if(test) then
     iold = 0
     do j=0, nbft-1
        n = j/int_bits+1
        jshift = j-(n-1)*int_bits
        if( btest(icij(ispin,n,ilength),jshift) ) iold=iold+1	! BW
     enddo
     if(iold.ne.nhalf) then
        write(*,*)'Wrong no. of particles upon entry', iold
        STOP 'in singles'
     endif
  endif

  !     randomly pick an orbital from the current configuration
  it = 0
10 call random(seed,rand)
  iorb = int(rand*(nhalf-nfreeze)) + 1
  if(iorb.lt.1.or.iorb.gt.(nhalf-nfreeze)) STOP 'Error singles: random'

  icount = 0
  do j=1, nactive

     ja=iactive(j)-1

     n1 = ja/int_bits+1
     j1 = ja-(n1-1)*int_bits
     if( btest(icij(ispin,n1,ilength),j1) ) icount=icount+1	! BW

     if (icount.eq.iorb) then 

        !           orbital ja belongs to which symmetry group?
        ilo  = 0
        do k=1, nsym
           ihi = ilo + nbpsy(k)-1
           if(ja.ge.ilo .and. ja.le.ihi) isym=k
           ilo = ilo + nbpsy(k)
        enddo
        ioffset = 0
        do k=1, isym-1
           ioffset = ioffset + nbpsy(k)
        enddo

        !     check to see if orbital excitations to this irrep are possible

        is=0
        icountsymk=0
        do k=1,nactive
           !        orbital iactive(k) belongs to which symmetry group?
           ilo = 0
           do i=1, nsym
              ihi = ilo + nbpsy(i)-1
              if((iactive(k)-1).ge.ilo .and.(iactive(k)-1).le.ihi) isymk=i
              ilo = ilo + nbpsy(i)
           enddo
           if(isymk.eq.isym) then
              icountsymk=icountsymk+1
              !           orbital iactive(k) is occupied?
              n = (iactive(k)-1)/int_bits+1
              jshift =(iactive(k)-1)-(n-1)*int_bits
              if( btest(icij(ispin,n,ilength),jshift) ) is=is+1
           endif
        enddo

        it = it + 1
        if(it.gt.500) return

        if(is.eq.icountsymk) then
           goto 10
        endif

        !     randomly generate new orbital of same symmetry

        it=0
20      call random(seed,rand)
        jj = int(rand*nactive)+1
        if(jj.lt.1.or.jj.gt.nactive) STOP 'Error singles: New orbit outside of range'
        jj = iactive(jj)-1

        !     orbital jj belongs to which symmetry group?

        ilo = ioffset
        ihi = ilo + nbpsy(isym)-1
        if(jj.ge.ilo.and.jj.le.ihi) then

           !        already occupied? compare with old configuration to avoid making
           !        substitution to orbital just deleted or still occupied orbital

           n2 = jj/int_bits + 1
           j2 = jj - (n2-1)*int_bits
           if(btest(icij(ispin,n2,ilength),j2)) then           ! BW
              it = it + 1
              if(it.gt.500) return
              goto 20
           endif

           icij(ispin,n1,ilength)=ibclr(icij(ispin,n1,ilength),j1)     ! BW

           icij(ispin,n2,ilength)=ibset(icij(ispin,n2,ilength),j2)     ! BW

        else
           it = it + 1
           if(it.gt.500) return
           goto 20
        endif


        if(test) then
           inew = 0
           do k=0, nbft-1
              n = k/int_bits + 1
              kshift = k -(n-1)*int_bits
              if(btest(icij(ispin,n,ilength),kshift))inew=inew+1	! BW
           enddo
           if(inew.ne.nhalf) then
              write(*,*)'New nhalf ',inew,' old nhalf ',iold
              write(*,*)'Substitution to orbital ',jj
              STOP 'Changed occupation in singles'
           endif
	    endif

!           substitution complete, let's get out of here
	    return
     
          endif

      enddo

      STOP 'Did not make a single substitution'

      end
