subroutine read_mcci_in(iword,maxc,maxocc,maxbfs,               &
                        icij,inflg,n_alpha,n_beta,ntotal,i_sx2, &
                        SCF_integral_filename,ieig,nfreeze,     &
                        ifreeze,nactive,iactive)

  use commonarrays, only: me
  use dyn_par,      only: maxtry, lmin,nbyte,int_bits,npfull,   &
                          lref,lkeep,conv_average,conv_history, &
                          cmin,hmin,davidson_stop,bmin,bmax,    &
                          cref,frac,conv_thresh_e,conv_thresh_l,&
                          test,time,time_all,generate_cfgs,     &
                          nobrnch_first,nodiag,i_want_conv,     &
                          npfull_conv,auto_cmin,cmin_finish

  ! written by Paul Delaney 25 October 2005.
  
  use precision 

  implicit none                   ! just say no

  !common  /para/    me, nproc

  ! input variables
  integer, intent(in)  :: iword
  integer, intent(in)  :: maxc
  integer, intent(in)  :: maxocc
  integer, intent(in)  :: maxbfs

  ! parameters which used to be in params, whether with defaults or not

  ! parameters for which there are no defaults.

  integer,        intent(out) :: inflg
  integer,        intent(out) :: n_alpha
  integer,        intent(out) :: n_beta
  integer,        intent(out) :: i_sx2
  ! also maxtry,cmin

  ! parameters for which there are defaults.

  character (len=12), intent(out) :: SCF_integral_filename  ! standard one- and two-electron integrals of the mos.
  integer,            intent(out) :: ieig
  integer,            intent(out) :: nfreeze
  integer,            intent(out) :: ifreeze(maxocc)
  integer,            intent(out) :: nactive
  integer,            intent(out) :: iactive(maxocc)
  ! also lmin,nbyte,int_bits,npfull,lref,lkeep,hmin,davidson_stop,bmin,bmax,cref,frac,test,time,time_all,generate_cfgs,nobrnch_first,nodiag
  ! conv_average,conv_history,conv_thresh_e,conv_thresh_l,i_want_conv,npfull_conv
  ! other output variables
  integer,            intent(out) :: ntotal
  integer,            intent(out) :: icij(2,iword,maxc)   ! mo_up and mo_dn tell us this


  ! local variables
  integer, parameter  :: line_length           = 79             ! the length of the line we read in from the control file
  integer, parameter  :: num_nondefault_params =  6             ! the number of program parameters for which there are no default values


  character (len=6)            :: format_string
  character (len=line_length)  :: line
  character (len=line_length)  :: keyword_string
  character (len=line_length)  :: value_string

  logical                      :: nondefault_param_set (num_nondefault_params)
  character (len=line_length)  :: nondefault_param_name(num_nondefault_params)

  integer             :: i
  logical             :: restart
  logical             :: syntax_error
  integer             :: n_alpha_list
  integer             :: n_beta_list
  real (kind=pr)      :: s
  integer             :: mo_up(maxocc)
  integer             :: mo_dn(maxocc)
  logical             :: is_up_occupied(maxbfs)
  logical             :: is_dn_occupied(maxbfs)
  integer             :: j, jshift, n
  !integer             :: me, nproc   !declaration moved to commonarrays.f90

  ! start of executable

  ! values of parameters for which there are defaults; nactive and iactive have defaults depending on nbft, nfreeze
  ! so work out these default values properly outside this subroutine, once we know nbft.

  SCF_integral_filename = "moints.ascii"  ! standard one- and two-electron integrals of the mos.
  ieig              =   1
  nfreeze           =   0
  ifreeze(:)        =  -1
  mo_up(:)          =   0
  mo_dn(:)          =   0
  nactive           =   0
  iactive(:)        =  -1
  lmin              = 100          ! min vector length
  nbyte             =   4          ! number of bytes / integer word
  int_bits          =  8*nbyte     ! bits / integer word
  npfull            =  10          ! perform full pruning after npfull steps
  lref              =  0           ! length of ref vector; <=length; if =0 then :=length
  lkeep             =  0           ! cnfgs #1 to #lkeep never touched in pruning
  conv_average      = 7            ! 7 fullprune steps are being averaged
  conv_history      = 7            ! 7 averaged Dvalues are being tracked
  hmin              = 1.0e-18_pr   ! h matrix threshold 
  davidson_stop     = 1.0e-07_pr   ! stop for convergence in davidson
  bmin              = 0.0_pr       ! min vector boost
  bmax              = 1.0_pr       ! max vector boost
  cref              = 0.1_pr       ! c>cref is treated as a ref. for branching
  frac              = 0.5_pr       ! rand < frac for branching
  conv_thresh_e     = 3.67e-06_pr  ! = (0.0001 eV)/(27.2113845 eV/a.u.) 
  conv_thresh_l     = 0.002_pr     ! changes in Dlength/length < 0.002
  test              = .FALSE.      ! error checking
  time              = .FALSE.      ! timing information, general
  time_all          = .FALSE.      ! timing information, detailed
  generate_cfgs     = .TRUE.       ! running with/without generating new cfgs
  nobrnch_first     = .FALSE.      ! no branching in the first step
  nodiag            = .FALSE.      ! only for collecting CSFs
  i_want_conv       = .FALSE.      ! stop mcci by using convergence crit.
  npfull_conv       = .TRUE.       ! convergence test only in npfull steps
  auto_cmin         = .FALSE.      ! Auto-adjust cmin
  cmin_finish       = 1e-4_pr
  open(10, file='mcci.in', form='formatted', status='old')
  rewind(10)

  nondefault_param_name(1)  = 'restart'
  nondefault_param_name(2)  = 'n_up'
  nondefault_param_name(3)  = 'n_dn'
  nondefault_param_name(4)  = 's'
  nondefault_param_name(5)  = 'cmin'
  nondefault_param_name(6)  = 'maxtry'

  nondefault_param_set = .FALSE.

  n_alpha_list = 0
  n_beta_list  = 0  

  select case (line_length)
  case (10:99)
     write(format_string, '(a2,i2,a1)' ) '(a',line_length,')'            ! (a79),  for example
  case (100:999)
     write(format_string, '(a2,i3,a1)' ) '(a',line_length,')'            ! (a132), for example
  case default
     if (me.eq.0) then
        write(50,*)
        write(50,*) 'Problem in subroutine read_mcci_in.f90'
        write(50,*) 'The parameter line_length is ',line_length
        write(50,*) 'and we are only set up to deal with line_length in the range 10 ... 999'
        write(50,*) 'Stopping...'
     endif
     stop
  end select

  do
     read(10, format_string, end=7000) line

     call parse_line(line, keyword_string, value_string)

     keyword_string = to_lower_case(keyword_string)

     select case (trim(adjustl(keyword_string)))

     case ('skip')  ! This means a blank line, or one with only a comment.

        ! First, list the parameters for which there are no defaults, and so must be set in the control file

     case ('restart')
        read(value_string,*) restart
        select case (restart)
        case (.TRUE.)
           inflg = 1
        case (.FALSE.)
           inflg = 0
        end select
        nondefault_param_set (1) = .TRUE.

     case ('n_up')
        read(value_string,*) n_alpha
        call non_negative_integer(n_alpha, 'n_dn')
        nondefault_param_set (2) = .TRUE.

     case ('n_dn')
        read(value_string,*) n_beta
        call non_negative_integer(n_beta, 'n_dn')
        nondefault_param_set (3) = .TRUE.

     case ('s')
        read(value_string,*) s
        call non_negative_real(s, 's')

        ! convert to spin times 2
        i_sx2 = int(2*s + 0.1_pr)    ! plus .1 in case of round-off
        nondefault_param_set (4) = .TRUE.

     case ('cmin')
        read(value_string,*) cmin
        call non_negative_real(cmin, 'cmin')
        nondefault_param_set (5) = .TRUE.

     case ('maxtry')
        read(value_string,*) maxtry
        call non_negative_integer(maxtry, 'maxtry')
        nondefault_param_set (6) = .TRUE.

        ! Now for the parameters for which there are defaults, and so need not be set in the control file

     case ('mo_up')
        call read_orbital_list(value_string, mo_up, n_alpha_list, syntax_error)
        if (syntax_error) then
           if (me.eq.0) then
              write(50,*) 'Could not read the list of up orbitals in mcci.in'
              write(50,*) 'the list is ', value_string
           endif
           stop
        end if

     case ('mo_dn')
        call read_orbital_list(value_string, mo_dn, n_beta_list, syntax_error)
        if (syntax_error) then
           if (me.eq.0) then
              write(50,*) 'Could not read the list of down orbitals in mcci.in'
              write(50,*) 'the list is ', value_string
           endif
           stop
        end if

     case ('scf_integral_filename')
        read(value_string,*) SCF_integral_filename
        SCF_integral_filename = trim(adjustl(SCF_integral_filename))
        select case (trim(SCF_integral_filename))
        case ('moints.bTM')
        case ('moints.TM')
        case ('moints.ascii')
        case ('FCIDUMP')
        case default
           if (me.eq.0) then
              write(50,*) 
              write(50,*) 'SCF_integral_filename must be one of moints.bTM, moints.TM, moints.ascii or FCIDUMP.'
              write(50,*) 'Stopping...'
           endif
           stop
        end select

     case ('ieig')
        read(value_string,*) ieig
        call non_negative_integer(ieig,'ieig')

     case ('frozen_doubly')
        call read_orbital_list(value_string, ifreeze, nfreeze, syntax_error)
        if (syntax_error) then
           if (me.eq.0) then
              write(50,*) 'Could not read the list of frozen doubly occupied orbitals in mcci.in'
              write(50,*) 'the list is ', value_string
           endif
           stop
        end if

     case ('mo_active')
        call read_orbital_list(value_string, iactive, nactive, syntax_error)
        if (syntax_error) then
           if (me.eq.0) then
              write(50,*) 'Could not read the list of active orbitals in mcci.in'
              write(50,*) 'the list is ', value_string
           endif
           stop
        end if
        iactive(nactive+1:) = -1  ! restore Werner's default, overwritten by read_orbital_list

     case ('generate_cfgs')
        read(value_string,*) generate_cfgs 
        select case (generate_cfgs)
        case (.FALSE.)
           cref = 10.0_pr
           lmin = 0   
           bmax = 0.0_pr
           frac = 0.0_pr
           npfull = 1
        case (.TRUE.)
        end select

     case ('lmin')
        read(value_string,*) lmin
        call non_negative_integer(lmin,'lmin')

     case ('nbyte')
        read(value_string,*) nbyte
        call non_negative_integer(nbyte,'nbyte')

     case ('int_bits')
        read(value_string,*) int_bits
        call non_negative_integer(int_bits,'int_bits')

     case ('npfull')
        read(value_string,*) npfull
        call non_negative_integer(npfull,'npfull')

     case ('lref')
        read(value_string,*) lref
        call non_negative_integer(lref,'lref')

     case ('lkeep')
        read(value_string,*) lkeep
        call non_negative_integer(lkeep,'lkeep')

     case ('conv_average')
        read(value_string,*) conv_average
        call non_negative_integer(conv_average,'conv_average')

     case ('conv_history')
        read(value_string,*) conv_history
        call non_negative_integer(conv_history,'conv_history')

     case ('hmin')
        read(value_string,*) hmin
        call non_negative_real(hmin, 'hmin')

     case ('davidson_stop')
        read(value_string,*) davidson_stop
        call non_negative_real(davidson_stop, 'davidson_stop')

     case ('bmin')
        read(value_string,*) bmin
        call non_negative_real(bmin, 'bmin')

     case ('bmax')
        read(value_string,*) bmax
        call non_negative_real(bmax, 'bmax')

     case ('cref')
        read(value_string,*) cref
        call non_negative_real(cref, 'cref')

     case ('frac')
        read(value_string,*) frac
        call non_negative_real(frac, 'frac')

     case ('conv_thresh_e')
        read(value_string,*) conv_thresh_e
        call non_negative_real(conv_thresh_e, 'conv_thresh_e')

     case ('conv_thresh_l')
        read(value_string,*) conv_thresh_l
        call non_negative_real(conv_thresh_l, 'conv_thresh_l')

     case ('auto_cmin')
        read(value_string,*) auto_cmin

     case ('cmin_finish')
        read(value_string,*) cmin_finish
        call non_negative_real(cmin_finish,'cmin_finish')

     case ('test')
        read(value_string,*) test

     case ('time')
        read(value_string,*) time

     case ('time_all')
        read(value_string,*) time_all

     case ('nobrnch_first')
        read(value_string,*) nobrnch_first

     case ('nodiag')
        read(value_string,*) nodiag

     case ('i_want_conv')
        read(value_string,*) i_want_conv

     case ('npfull_conv')
        read(value_string,*) npfull_conv

     !The following three keywords are related to memory parameters
     !and are ignored by read_mcci_in
     case ('maxc')
     case ('kmax')
     case ('maxh')

     case default
        if (me.eq.0) then
           write(50,*)
           write(50,*) 'Problem in subroutine read_mcci_in.f90'
           write(50,*) 'The line we''re processing from the control file is:'
           write(50,*) line
           write(50,*) 'and this has been parsed to yield'
           write(50,*) 'Keyword :', keyword_string
           write(50,*) 'Value   :', value_string
           write(50,*)
           write(50,*) 'This Keyword is not recognised.'
           write(50,*) 'Stopping...'
        endif
        stop

     end select

  end do

7000 continue
    
  if ( any( .NOT. nondefault_param_set) ) then
     if (me.eq.0) then
        write(50,*)
        write(50,*) 'Problem in subroutine read_mcci_in.f90'
        write(50,*)
        write(50,*) 'In reading in the control file "mcci.in", we could not'
        write(50,*) 'find values for the following parameter(s):'
        write(50,*)
     endif
     
     do i=1,num_nondefault_params
        if (me.eq.0) then
           if (.NOT. nondefault_param_set(i)) write(50,*) nondefault_param_name(i)
        endif
     end do
     
     if (me.eq.0) then
        write(50,*) 'These parameter(s) do not have default values, and so must be set explicitly.'
        write(50,*) 'Stopping...'
     endif
     stop
  end if

  close(10)
      
  if (n_alpha /= n_alpha_list) then
     if (me.eq.0) then
        write(50,*)
        write(50,*) 'Problem in subroutine read_mcci_in.f90'
        write(50,*) 'n_alpha = ',n_alpha,' but ',n_alpha_list,' up orbitals are listed'
        write(50,*) 'Stopping...'
     endif
     stop
  end if

  if (n_beta /= n_beta_list) then
     if (me.eq.0) then
        write(50,*)
        write(50,*) 'Problem in subroutine read_mcci_in.f90'
        write(50,*) 'n_beta = ',n_beta,' but ',n_beta_list,' up orbitals are listed'
        write(50,*) 'Stopping...'
     endif
     stop
  end if

  ntotal = n_alpha + n_beta

  icij(:,:,1) = 0

  is_up_occupied(:) = .FALSE.

  ! set bits in reference
  do j=1,n_alpha
     jshift      = mo_up(j) - 1
     n           = jshift/int_bits +1
     jshift      = jshift -(n-1)*int_bits
     icij(1,n,1) = ibset(icij(1,n,1),jshift)                       ! BW
     is_up_occupied(mo_up(j)) = .TRUE.
  enddo

  is_dn_occupied(:) = .FALSE.

  do j=1,n_beta
     jshift      = mo_dn(j) - 1
     n           = jshift/int_bits +1
     jshift      = jshift -(n-1)*int_bits
     icij(2,n,1) = ibset(icij(2,n,1),jshift)                       ! BW
     is_dn_occupied(mo_dn(j)) = .TRUE.
  enddo

  ! check to see if the frozen orbitals are really doubly occupied in the reference

  do i=1,nfreeze
     
     if (.not.(is_up_occupied(ifreeze(i)).and.is_dn_occupied(ifreeze(i)))) then
     if (me.eq.0) then
        write(50,*)
        write(50,*) 'Frozen orbital isn''t doubly occupied in the reference. Stopping.'
     endif
     stop
     end if
  end do
  
  
  if (me.eq.0) then
     write(50,*)
     write(50,*) 'The following parameters were read from mcci.in'
     write(50,*)
     write(50,*) 'parameters with no default values:'
     write(50,*) '----------------------------------'
  
     write(50,*)
     write(50,8010) 'inflg                        = ', inflg
     write(50,8010) 'n_alpha                      = ', n_alpha
     write(50,8010) 'n_beta                       = ', n_beta
     write(50,8000) 's                            = ', s
     write(50,8002) 'cmin                         = ', cmin
     write(50,8010) 'maxtry                       = ', maxtry

     write(50,*)
     write(50,*)
     write(50,*) 'parameters with default values:'
     write(50,*) '-------------------------------'

     write(50,*)
     write(50,*)    'SCF_integral_filename        = ', SCF_integral_filename

     write(50,*)  
     write(50,8010) 'ieig                         = ', ieig

     !write(50,'(1x,a,100i4)')    'mo_up(1:n_alpha)             =   ', mo_up(1:n_alpha)
     !write(50,'(1x,a,100i4)')    'mo_dn(1:n_beta)              =   ', mo_dn(1:n_beta)

     write(50,*)  
     write(50,8010) 'nfreeze                      = ', nfreeze
     if (nfreeze /= 0) &
        !write(50,'(1x,a,100i4)')    'ifreeze(1:nfreeze)           =   ', ifreeze(1:nfreeze)

     write(50,*)  
     write(50,8010) 'nactive                      = ', nactive
     if (nactive /= 0) &
        write(50,'(1x,a,100i4)')    !'iactive(1:nactive)           =   ', iactive(1:nactive)

     write(50,8010) 'lmin                         = ', lmin
     write(50,8010) 'nbyte                        = ', nbyte
     write(50,8010) 'int_bits                     = ', int_bits
     write(50,8010) 'npfull                       = ', npfull
     write(50,8010) 'lref                         = ', lref
     write(50,8010) 'lkeep                        = ', lkeep
     write(50,8010) 'conv_average                 = ', conv_average
     write(50,8010) 'conv_history                 = ', conv_history
     write(50,8002) 'hmin                         = ', hmin
     write(50,8002) 'davidson_stop                = ', davidson_stop
     write(50,8000) 'bmin                         = ', bmin
     write(50,8000) 'bmax                         = ', bmax
     write(50,8000) 'cref                         = ', cref
     write(50,8000) 'total frac                   = ', frac
     write(50,8000) 'conv_thresh_e                = ', conv_thresh_e
     write(50,8000) 'conv_thresh_l                = ', conv_thresh_l
     write(50,8020) 'test                         = ', test
     write(50,8020) 'time                         = ', time
     write(50,8020) 'time_all                     = ', time_all
     write(50,8020) 'generate_cfgs                = ', generate_cfgs
     write(50,8020) 'nobrnch_first                = ', nobrnch_first
     write(50,8020) 'nodiag                       = ', nodiag       
     write(50,8020) 'i_want_conv                  = ', i_want_conv  
     write(50,8020) 'npfull_conv                  = ', npfull_conv  
  endif

8000 format (1x, a, 2x, f9.4)       ! normal reals
8002 format (1x, a, 4x, es9.2)      ! reals that are very big or small
8010 format (1x, a, i6)             ! integers
8020 format (1x, a, 5x, l1)         ! logicals
    
contains
    
  subroutine parse_line(line, keyword_string, value_string)

    implicit none                          ! just say no

    ! input variables
    character (len=*), intent(in) :: line

    ! output variables
    character (len=len(line)), intent(out) :: keyword_string
    character (len=len(line)), intent(out) :: value_string

    ! local variables
    integer                    :: comment_position
    character (len=len(line))  :: working_line
    integer                    :: i
    integer                    :: equals_leftmost_position
    integer                    :: equals_rightmost_position
    integer                    :: equals_position

    ! Start of executable

    working_line = line


    ! First, look for "!", the comment sign, and convert it and everything rightwards of it to spaces.

    comment_position = index(working_line,'!')

    if (comment_position /= 0) working_line = working_line(1:comment_position-1)


    ! Remove all funny control characters, like CTRL-x's and TABS, and convert them all to spaces.

    do i=1,len(line)
       if (iachar(working_line(i:i)) <= 31) working_line(i:i) = ' ' 
    end do



    ! Now see if we have a purely blank line.        
    if (len_trim(working_line) == 0) then
       keyword_string = 'skip'
       return
    end if

    ! From here on, we assume we have a line with a KEYWORD = VALUE statement in it
    ! The keyword is defined to be the connected component before the = sign, and the value that after.

    ! First, check we have one "=" present

    equals_leftmost_position  = index(working_line, '=', back=.FALSE.)
    equals_rightmost_position = index(working_line, '=', back=.TRUE.)

    if (equals_leftmost_position == 0) then
       if(me.eq.0) then
          write(50,*)
          write(50,*) 'Problem in subroutine read_mcci_in.f90'
          write(50,*) 'The line we are processing from the control file is:'
          write(50,*) line
          write(50,*) 'This seems to be a KEYWORD=VALUE type line, but we find no equal-to sign "=" '
          write(50,*) 'Stopping...'
       endif
       stop
    end if

    if (equals_leftmost_position /= equals_rightmost_position) then
       if(me.eq.0) then
          write(50,*)
          write(50,*) 'Problem in subroutine read_mcci_in.f90'
          write(50,*) 'The line we are processing from the control file is:'
          write(50,*) line
          write(50,*) 'This seems to be a KEYWORD=VALUE type line, but we find more than one equal-to sign "=" '
          write(50,*) 'Stopping...'
       endif
       stop
    end if

    ! Now we know we have precisely one equal-to sign
    equals_position = equals_leftmost_position
    
    keyword_string  = working_line(1:equals_position-1)
    value_string    = working_line(equals_position+1:)
    
  end subroutine parse_line

  
  subroutine read_orbital_list(line_in, list_out, num_orbs, syntax_error)

    ! given a string line_in, read a set of molecular orbitals from it, using the Turbomole format.
    ! Return the (not necc in any order) list of orbitals in list_out 
    ! The number of orbitals so specified is returned in num_orbs
    
    character (len=*), intent(in)  :: line_in
    integer,           intent(out) :: list_out(:)
    integer,           intent(out) :: num_orbs
    logical,           intent(out) :: syntax_error

    ! local variables
    integer :: line_length
    integer :: list_length
    integer :: scan_start
    integer :: int_end
    integer :: sep_end
    logical :: find_int_success
    logical :: comma
    logical :: minus
    logical :: find_sep_success
    logical :: need_int
    integer :: start_int
    integer :: end_int
    integer :: counter


    ! start of executable
    line_length = len(line_in)
    list_length = size(list_out)

    num_orbs    = 0
    list_out(:) = 0
    syntax_error = .FALSE.

    if (line_length == 0) return

    scan_start = 1
    need_int     = .FALSE.

    do
       call find_int(line_in,scan_start,int_end,find_int_success)

       if (.not. find_int_success) then
          if (need_int) syntax_error = .TRUE.
          return
       end if

       num_orbs = num_orbs + 1
       if (num_orbs > list_length) then
          if (me.eq.0) then
             write(50,*)
             write(50,*) 'Too many orbitals in read_orbital_list. Stopping.'
          endif
          stop
       end if

       read(line_in(scan_start:int_end),*) list_out(num_orbs)

       if (int_end == line_length) return
       
       call find_separator(line_in,int_end+1,sep_end,comma,minus,find_sep_success)
       if (.not. find_sep_success) then
          syntax_error = .TRUE.
          return
       end if

       ! set status of need_int for next iteration
       
       if (comma) then
          need_int = .TRUE.
       else
          need_int = .FALSE.
       end if

       if (.not. minus) then
          scan_start = sep_end+1
          cycle
       end if

       if (minus) then
          call find_int(line_in,sep_end+1,int_end,find_int_success)
          
          if (.not. find_int_success) then
             syntax_error = .TRUE.
             return
          end if

          start_int = list_out(num_orbs)
          read(line_in(sep_end+1:int_end),*) end_int
          if (end_int < start_int) then
             syntax_error = .TRUE.
             return
          end if

          if (num_orbs + (end_int-start_int) > list_length) then
          write(*,*) num_orbs + (end_int-start_int),list_length
             if (me.eq.0) then
                write(50,*)
                write(50,*) 'Too many orbitals in read_orbital_list. Stopping.'
             endif
             stop
          end if

          do counter = 1, (end_int-start_int)
             list_out(num_orbs+counter) = start_int + counter
          end do
          
          num_orbs = num_orbs+(end_int-start_int)
          
          if (int_end == line_length) exit

          call find_separator(line_in,int_end+1,sep_end,comma,minus,find_sep_success)

          if (.not. find_sep_success) then
             syntax_error = .TRUE.
             return
          end if
          
          if (minus) then
             syntax_error = .TRUE.
             return
          end if

          if (comma) then
             need_int = .TRUE.
          else
             need_int = .FALSE.
          end if

          scan_start = sep_end+1

       end if
       
    end do
    
  end subroutine read_orbital_list


  subroutine find_int(line_in,scan_start,int_end,find_int_success)

    character (len=*), intent(in)  :: line_in
    integer,           intent(in)  :: scan_start
    integer,           intent(out) :: int_end
    logical,           intent(out) :: find_int_success

    ! local variables 
    integer :: posn

   if ((scan_start < 1) .or. (scan_start>len(line_in)) ) then
       find_int_success = .FALSE.
       return
    end if

    ! try to find the start character of the integer

    posn = scan_start

    do 
       select case (line_in(posn:posn))
       case (' ')
          posn = posn + 1
          if (posn > len(line_in)) then
             find_int_success = .FALSE.
             return
          end if

       case ('1','2','3','4','5','6','7','8','9','0')
          exit

       case default
          find_int_success = .FALSE.
          return
       end select

    end do
          

    ! try to find the last character of the integer

    do 
       select case (line_in(posn:posn))
       case ('1','2','3','4','5','6','7','8','9','0')
          if (posn /= len(line_in)) then
             posn = posn + 1
          else
             find_int_success = .TRUE.
             int_end = posn
             return
          end if

       case (' ',',','-')
          find_int_success = .TRUE.
          int_end = posn-1
          return
  
       case default
          find_int_success = .FALSE.
          return
       end select

    end do

  end subroutine find_int


  subroutine find_separator(line_in,scan_start,sep_end,comma,minus,find_sep_success)

    character (len=*), intent(in)  :: line_in
    integer,           intent(in)  :: scan_start
    integer,           intent(out) :: sep_end
    logical,           intent(out) :: comma
    logical,           intent(out) :: minus
    logical,           intent(out) :: find_sep_success

    ! local variables 
    integer :: posn
    logical :: found_a_space

   if ((scan_start < 1) .or. (scan_start>len(line_in)) ) then
       find_sep_success = .FALSE.
       return
    end if

    ! separators are :  1 one or more blanks and then (EOL OR a character not a blank or comma or minus)
    !                   2 any number of blanks and then a comma
    !                   3 any number of blanks and then a minus

    posn          = scan_start
    found_a_space = .FALSE.
    comma         = .FALSE.
    minus         = .FALSE.

    do 
       select case (line_in(posn:posn))
       case (' ')
          found_a_space = .TRUE.
          if (posn /= len(line_in)) then
             posn = posn + 1
          else
             sep_end = posn
             find_sep_success = .TRUE.
             return
          end if
       case (',')
          sep_end = posn
          comma = .TRUE.
          find_sep_success = .TRUE.
          return
       case ('-')
          sep_end = posn
          minus = .TRUE.
          find_sep_success = .TRUE.
          return
       case default
          if (found_a_space) then
             find_sep_success = .FALSE.
          else
             find_sep_success = .TRUE.
          end if
          return
       end select
    end do

  end subroutine find_separator


  subroutine non_negative_integer(check_int, check_int_name)

    ! sanity check on integers which are not allowed to be negative

    ! input variables        
    integer,           intent(in) :: check_int
    character (len=*), intent(in) :: check_int_name
     
    ! start of executable
    
    if (check_int < 0) then
       if(me.eq.0) then
          write(50,*)
          write(50,*) 'Control file error:'
          write(50,*) '-------------------'
       
          write(50,*) 'In the control file "mcci.in"'
          write(50,*) 'the parameter ', check_int_name, ' has been set to the value ', check_int
          write(50,*) 'This parameter must have a non-negative value.'
          write(50,*) 'Stopping...'
       endif
       stop
    end if

  end subroutine non_negative_integer


  subroutine non_negative_real(check_real, check_real_name)

    ! sanity check on reals which are not allowed to be negative
    
    ! input variables        
    real      (kind=pr), intent(in) :: check_real
    character (len=*),   intent(in) :: check_real_name
    
    ! start of executable
    
    if (check_real < 0.0_pr) then
       if(me.eq.0) then
          write(50,*)
          write(50,*) 'Control file error:'
          write(50,*) '-------------------'
       
          write(50,*) 'In the control file "mcci.in"'
          write(50,*) 'the parameter ', check_real_name, ' has been set to the value ', check_real
          write(50,*) 'This parameter must have a non-negative value.'
          write(50,*) 'Stopping...'
       endif
       stop
    end if

  end subroutine non_negative_real


  function to_lower_case(string)
     
    ! replace all upper case characters in "string" with lower case ones
    
    ! input variables
    character (len=*), intent(in) :: string
    
    ! output variables
    character (len=len(string))              :: to_lower_case
    
    ! local variables
    
    character (len=26), parameter :: capitals = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    character (len=26), parameter :: lowers   = "abcdefghijklmnopqrstuvwxyz"
    integer   :: posn
    character :: letter
    integer   :: number
    
    to_lower_case = string
     
    do 
       posn = scan(to_lower_case, capitals)    ! position of leftmost character of out_string which is in capitals
       if (posn == 0) exit
       letter  = to_lower_case(posn:posn)
       
       number = scan(capitals,letter)       ! number should be between 1 and 26
       to_lower_case(posn:posn) = lowers(number:number)
    end do

  end function to_lower_case
  
end subroutine read_mcci_in

