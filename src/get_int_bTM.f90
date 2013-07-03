subroutine get_int_bTM(ecore)
      use commonarrays
      use dyn_par
!     Read in the one- and two-electron integrals from disk.
!     Do this in binary format, from the modified dscf/mpgrad turbomole calculation

!     TURBOMOLE has written out all the oneints (occupied and virtual) from the
!     dscf run. MPGRAD then restricted the oneints and twoints to the set specified
!     in the MPGRAD control file. These are the ones we read in here.

      implicit real*8  (a-h,o-z)


!     output variables 
      real*8 :: ecore
!      common  /pointer/ ipoint(max2)
!      common  /config/  icij(2,iword,maxc), nword    ! we set nword here
!      common  /aodat/   nsym, nbft, nbpsy(irmax)     ! all these are set here
!      common  /ints/    e1ints(max1), e2ints(max2)   ! and these


!     local variables
      integer :: ierr
      integer :: iline
      integer :: isym
      integer :: isym_temp

      character (len=4) :: ityp  (irmax)             ! irrep labels
      integer           :: idim  (irmax)             ! dimensions of real irreps
      integer           :: nlamda(irmax)             ! number of occurences of each irrep
      integer           :: norbs (irmax)             ! highest occupied occurrence of each irrep
                                                     ! - this is from occupancies set in the
                                                     ! - MPGRAD control file.

      integer :: num_oneints
      integer :: num_twoints
      integer :: i,j
      real*8  :: eint
      integer :: k,l
      integer :: ij, jl, ijkl


!     many variable names have been taken from the COLUMBUS documentation
!     and some from TURBOMOLE (nlamda, nsym, idim, ityp, norbs)

!     Start of executable

!     store values for pointer array for indexing symmetric matrices
      ipoint(1) = 0
      do i = 2,  max2 
         ipoint(i) = ipoint(i-1) + i - 1
      enddo


      open(unit=10,file="moints.bTM",iostat=ierr, status='old',form='unformatted')

      if (ierr/=0) then
         write(50,*)
         write(50,*) "Can not open file 'moints.bTM' in get_int_bTM.f"
         write(50,*) "Stopping..."
         stop
      end if

      rewind(10)

      read(10) nsym
      
      if (nsym > irmax) then
         write(50,*)
         write(50,*) 'Problem in get_int_bTM.f'
         write(50,*) 'The number of real irreps from the Turbomole calculation is ',nsym
         write(50,*) 'but the maximum allowed number is ',irmax
         write(50,*) 'Increase irmax.'
         write(50,*) 'Stopping...'
         stop
      end if

      ityp  (:) = '    '
      idim  (:) = 0
      nlamda(:) = 0
      norbs (:) = 0

      nbft = 0
      do isym = 1,nsym

         read(10) isym_temp, ityp(isym), idim(isym),& 
                  nlamda(isym), norbs(isym)

         nbpsy(isym) = idim(isym)*norbs(isym)

         nbft = nbft + idim(isym)*norbs(isym)

      end do

      if (nbft > maxbfs) then
         write(50,*)
         write(50,*) 'Problem in get_int_bTM.f'
         write(50,*) 'nbft     =',nbft
         write(50,*) 'maxbfs   =',maxbfs
         write(50,*) 'so maxbfs is too small: increase it.'
         write(50,*) 'Stopping...'
         stop
      end if
         


!     Now we can work out nword

      nword = (nbft-1)/int_bits + 1

      if (nword > iword) then
         write(50,*)
         write(50,*) 'Problem in get_int_bTM.f'
         write(50,*) 'nword is bigger than iword: nword,iword=',&
                      nword,iword
         write(50,*) 'Change iword in params. Stopping.'
         stop
      end if

!     Read nuclear repulsion energy

      read(10) ecore

!     calculate the number of lines containing one-electron integrals

      num_oneints = 0
      do isym=1,nsym
         num_oneints = num_oneints + &
                       idim(isym)*(( norbs(isym)*(norbs(isym)+1) )/2)
      end do

      if (num_oneints > max1) then
         write(50,*)
         write(50,*) 'Problem in get_int_bTM.f'
         write(50,*) 'We are told there are ',num_oneints
         write(50,*) 'one-electron integrals in the moints.bTM file'
         write(50,*) 'But max1 has been set at ',max1
         write(50,*)
         write(50,*) 'Increase max1.'
         write(50,*) 'Stopping...'
         stop
      end if         
         

      e1ints(:)=0.0d0  ! initialise e1ints to zero.  

      do iline=1,num_oneints

         read(10) i,j, eint

         if (j>i) then
            write(50,*) 
            write(50,*) 'Problem in get_int_bTM.f'
            write(50,*) 'Found j>i while reading in the oneints.'
            write(50,*) 'i,j =',i,j
            write(50,*) 'Stopping...'
            stop
         end if

         ij = ipoint(i)+j
         e1ints(ij) =  eint

      end do

!     Now for the two-electron integrals

      e2ints(:)=0.0d0  ! initialise e2ints to zero.  
      num_twoints = 0

!     Read till the end of the file

      do
         read(10,end=3000) i,j,k,l,eint

!     LOTS of checks to see if i,j,k,l are sensible

         if ( j>i .or. l>k .or. k>i .or. (i==k .and. l>j) ) then
            write(50,*) 
            write(50,*) 'Problem in get_int_bTM.f'
            write(50,*) 'Found i,j,k,l not in canonical ordering.'
            write(50,*) 'i,j =', i,j
            write(50,*) 'k,l =', k,l
            write(50,*) 'Stopping...'
            stop            
         end if
         
         if ((i>max2) .or. (k>max2)) then
            write(50,*) 
            write(50,*) 'Problem in get_int_bTM.f'
            write(50,*) 'Found i or k greater that max2.'
            write(50,*) 'i, k =', i, k
            write(50,*) 'max2 =', max2
            write(50,*) 'Stopping...'
            stop            
         end if

         ij   = ipoint(i)+j
         kl   = ipoint(k)+l

         if (ij>max2)  then
            write(50,*) 
            write(50,*) 'Problem in get_int_bTM.f'
            write(50,*) 'Found ij greater that max2.'
            write(50,*) 'i, j =', i, j
            write(50,*) 'ij   =', ij
            write(50,*) 'max2 =', max2

            write(50,*) 
            write(50,*) 'This means there are too many two-electron integrals'
            write(50,*) 'Stopping...'
            stop            
         end if


         ijkl = ipoint(ij)+kl

         if (ijkl>max2)  then
            write(50,*) 
            write(50,*) 'Problem in get_int_bTM.f'
            write(50,*) 'Found ijkl greater that max2.'
            write(50,*) 'i,j  =', i,j
            write(50,*) 'j,l  =', k,l
            write(50,*) 'ij   =', ij
            write(50,*) 'kl   =', kl
            write(50,*) 'ijkl =', ijkl
            write(50,*) 'max2 =', max2

            write(50,*) 
            write(50,*) 'This means there are too many two-electron integrals'
            write(50,*) 'Stopping...'
            stop            
         end if

         e2ints(ijkl)=eint

         num_twoints = num_twoints+1

      end do

 3000 close(10)

      write(50,*) 
      write(50,*) 'read ', num_twoints, ' two electron integrals'

      end
