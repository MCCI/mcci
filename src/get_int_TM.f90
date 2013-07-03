subroutine get_int_TM(ecore)
  use commonarrays, only:ipoint,icij,nword,nsym,nbft,nbpsy,e1ints,e2ints
  use dyn_par
  !     get one- and two-electron integrals from Turbomole, plus some other data

  implicit real*8  (a-h,o-z)

  !     output variables 
  real*8 :: ecore
  !common  /pointer/ ipoint(max2)
  !common  /config/  icij(2,iword,maxc), nword    ! we set nword here
  !common  /aodat/   nsym, nbft, nbpsy(irmax)     ! all these are set here
  !common  /ints/    e1ints(max1), e2ints(max2)   ! and these


  !     local variables
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
  integer :: i,j
  real*8  :: eint
  integer :: num_twoints
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

  open(10, file='moints.TM', form='formatted')
  rewind(10)

  !     Skip the text only lines at the beginning, up to and including 
  !     "number of real irred. reps..."

  do iline = 1,40
     read(10,*)
  end do

  read(10,*) nsym


  if (nsym > irmax) then
     write(50,*)
     write(50,*) 'Problem in get_int_TM.f'
     write(50,*) 'The number of real irreps from the Turbomole calculation is ',nsym
     write(50,*) 'but the maximum allowed number is ',irmax
     write(50,*) 'Increase irmax.'
     write(50,*) 'Stopping...'
     stop
  end if



  do i=1,6  ! Skip blank, "for each irr", "----", blank, "number :", blank
     read (10,*)
  end do

  ityp  (:) = '    '
  idim  (:) = 0
  nlamda(:) = 0
  norbs (:) = 0

  nbft = 0
  do isym = 1,nsym

     read(10,1000) isym_temp, ityp(isym), idim(isym),nlamda(isym), norbs(isym)

     nbpsy(isym) = idim(isym)*norbs(isym)

     nbft = nbft + idim(isym)*norbs(isym)

     if (isym /= isym_temp) then
        write(50,*)
        write(50,*) 'Problem in get_int_TM.f'
        write(50,*) 'isym     =',isym
        write(50,*) 'isym_temp=',isym_temp
        write(50,*) 'These should be equal.'
        write(50,*) 'Stopping...'
        stop
     end if
  end do

1000 format(2x,i3,7x,a4,6x,i3,8x,i5,5x,i5)


  if (nbft > maxbfs) then
     write(50,*)
     write(50,*) 'Problem in get_int_TM.f'
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
     write(50,*) 'Problem in get_int_TM.f'
     write(50,*) 'nword is bigger than iword: nword,iword=',nword,iword
     write(50,*) 'Change iword in params. Stopping.'
     stop
  end if

  !     Read nuclear repulsion energy

  read(10,*)       ! blank
  read(10,*)       ! "Nuclear rep.."
  read(10,*) ecore

  !     Skip "mos" checksum

  read(10,*)       ! blank
  read(10,*)       ! "Checksum of the ..."
  read(10,*)       ! the checksum itself

  !     Skip heading for T+V_nuc+ECP matrix element

  do iline=1,4      ! blank, blank, "T+V_nuc+...", "----"
     read(10,*)
  end do

  !     Echo heading for number of lines

  do iline=1,2     ! blank, "Number of data lines..."
     read(10,*)
  end do


  read(10,*) num_oneints

  if (num_oneints > max1) then
     write(50,*)
     write(50,*) 'Problem in get_int_TM.f'
     write(50,*) 'We are told there are ',num_oneints
     write(50,*) 'one-electron integrals in the moints.TM file'
     write(50,*) 'But max1 has been set at ',max1
     write(50,*)
     write(50,*) 'Increase max1.'
     write(50,*) 'Stopping...'
     stop
  end if


  !     Skip over description of data line format

  do iline=1,3     ! blank, "Format of data line:",  "i_mo, ..."
     read(10,*)
  end do


  e1ints(:)=0.0d0  ! initialise e1ints to zero.  

  do iline=1,num_oneints

     read(10,*) i,j, eint

     if (j>i) then
        write(50,*) 
        write(50,*) 'Problem in get_int_TM.f'
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

  !     Skip over header

  do iline=1,4     ! blank, blank, "Two-electron...",  "---..."
     read(10,*)
  end do

  !     Skip over description of data line format

  do iline=1,3    ! blank, "Format of ...", " i_mo, ..."
     read(10,*)
  end do

  !     Now read till the end of the file

  do
     read(10,*,end=3000) i,j,k,l,eint

     !     LOTS of checks to see if i,j,k,l are sensible

     if ( j>i .or. l>k .or. k>i .or. (i==k .and. l>j) ) then
        write(50,*) 
        write(50,*) 'Problem in get_int_TM.f'
        write(50,*) 'Found i,j,k,l not in canonical ordering.'
        write(50,*) 'i,j =', i,j
        write(50,*) 'k,l =', k,l
        write(50,*) 'Stopping...'
        stop            
     end if

     if ((i>max2) .or. (k>max2)) then
        write(50,*) 
        write(50,*) 'Problem in get_int_TM.f'
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
        write(50,*) 'Problem in get_int_TM.f'
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
        write(50,*) 'Problem in get_int_TM.f'
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

  end do

3000 close(10)

end subroutine get_int_TM
