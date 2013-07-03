subroutine get_int(ecore)
  use commonarrays, only: nsym, nbft, nbpsy, e1ints, e2ints, ipoint, nword
  use dyn_par
  implicit real*8  (a-h,o-z)
  !common  /aodat/   nsym, nbft, nbpsy(irmax)
  !common  /ints/    e1ints(max1), e2ints(max2)
  !common  /para/    me, nproc !don't seem to be used here
  !common  /pointer/ ipoint(max2)
  !common  /config/  icij(2,iword,maxc), nword

  !     many variable names have been taken from the COLUMBUS
  !     documentation

  !     store values for pointer array for indexing symmetric matrices
  ipoint(1) = 0
  do i = 2,  max2 
     ipoint(i) = ipoint(i-1) + i - 1
  enddo

  open(10, file='moints.ascii', form='formatted')
  read(10,'(6i9)') ntitle, nsym, nbft, ninfo, nenrgy, nmap
  !     write(*,*) ntitle, nsym, nbft
  nword = (nbft-1)/int_bits + 1
  !     W.Gy I've changed the nominator in the previous line (nbft-1)
  !     that makes the following couple of lines unnecessary. The problem of
  !     'nword' and 'iword' has been solved very simple. The function for 'nword'
  !     above is tested and proved.

  !     if (nword > iword) then
  !        write(50,*)
  !        write(50,*) 'Problem in get_int.f'
  !        write(50,*) 'nword is bigger than iword: nword,iword=',nword,iword
  !        write(50,*) 'Change iword in params. Stopping.'
  !        stop
  !     end if

  !     skip titles
  do i=1, ntitle
     read(10,*)
  enddo
  read(10,*) (nbpsy(i), i=1, nsym )
  !     write(*,*) nsym, i
  !     write(*,*) nbpsy

  !     skip representation labels
  read(10,*)
  !     skip fsplit and record lengths
  read(10,*)
  !     skip orbital labels
  nskip = nbft/8 
  !     write(*,*) nskip
  if(nskip*8 .ne. nbft) nskip = nskip + 1
  !     write(*,*) nskip, nbft
  do i = 1, nskip
     read(10,*)
  enddo
  !     skip energy type (assume it to be nuclear repulsion)
  read(10,*)

  read(10,*) ecore
  !     write(*,*) ecore
  !
  !  New control structure added to skip 1e property integrals
  !  and only read in KE, Coulomb and ECP integrals where necessary
  !  Based on the fact that the above integrals have the following
  !  convention:
  !
  !   itypea = 0
  !			itypeb = 0 overlap 
  !			itypeb = 1 KE
  !			itypeb = 2 1e Coulomb
  !			itypeb = 3 ECP
  !
  !   itypea = 1 or 2 ---> property integrals
  !
  !   itypea = 3 ---> two-electron integrals
  !
  !   mn 27 May 2002
  !
  !  Some F90 contructs follow.
  !
  e1ints(:)=0.0d0  ! initialise eints to zero.  F90

  do                 ! PD June 7/02 - gets rid of labels and gotos

     read(10,'(7i8)') num, lrecl, idum1, itypea, itypeb, idum2, last
     !       write(*,*) num, lrecl, itypea, itypeb
     !
     !     itypea is zero, check if itypeb = 0 --> overlap integrals, and if so skip 
     !     else
     !     read in 1e integrals
     !
     if(itypea.eq.0)then

        if(itypeb.eq.0)then
           iskip = lrecl - 1
           do i=1,iskip
              read(10,*)
           enddo
        else
           iskip = lrecl - num - 1
           do i= 1, iskip
              read(10,*)
           enddo
           do i = 1, num
              read(10,*) eint, i1, j1
              if(j1.gt.i1) STOP 'Error get_int: 1e ints #2 '
              ij = ipoint(i1)+j1
              e1ints(ij) = e1ints(ij)  + eint
           enddo
        endif
        !
        !   itypea = 1 or 2, skip property integrals
        !
     elseif(itypea.eq.1.or.itypea.eq.2) then
        iskip = lrecl - 1
        do i= 1, iskip
           read(10,*) 
        enddo

        !
        !     itypea = 3, read in 2e integrals
        !
     elseif(itypea.eq.3)then

        do i=1, num
           read(10,*) eint, i2, j2, k2, l2
           if(j2.gt.i2 .or. l2.gt.k2) STOP 'get_int: 2e ints #1'
           ij   = ipoint(i2)+j2
           kl   = ipoint(k2)+l2
           if( kl .gt. ij) then
              itmp = ij
              ij =kl
              kl = itmp
           endif
           ijkl = ipoint(ij)+kl
           e2ints(ijkl)=eint
        enddo
        if( last .eq. 2) exit    ! We have reached the end of the file. F90
     endif

  end do
  return
end subroutine get_int
