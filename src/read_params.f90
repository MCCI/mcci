subroutine read_params()
  !This subroutine is similar to read_mcci_in.f90.  mcci.in is
  !scanned for parameter keywords.
  use commonarrays
  use dyn_par
  use precision 
  ! March 2010

  implicit none   
  integer  :: i

  ! local variables
  integer, parameter   :: line_length           = 79             ! the length of the line we read in from the control file
  integer, parameter   :: num_nondefault_params =  4             ! the number of program parameters for which there are no default values
  integer              :: rubbish(4)

  ! variables related to reading info from turbomole
  integer, allocatable :: ityp(:),idim(:),nlamda(:),norbs(:)
  integer              :: isym, isym_temp

  character (len=6)            :: format_string
  character (len=line_length)  :: line
  character (len=line_length)  :: keyword_string
  character (len=line_length)  :: value_string
  character (len=12)  :: SCF_integral_filename

  logical                      :: nondefault_param_set (num_nondefault_params)
  character (len=line_length)  :: nondefault_param_name(num_nondefault_params)

  open(10, file='mcci.in', form='formatted', status='old')
  rewind(10)

  nondefault_param_name(1)  = 'maxh'
  nondefault_param_name(2)  = 'kmax'
  nondefault_param_name(3)  = 'maxc'
  nondefault_param_name(4)  = 'SCF_integral_filename'

  nondefault_param_set = .FALSE.
  int_bits = 32

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

     case ('maxh')
        read(value_string,*) maxh
        call non_negative_integer(maxh, 'maxh')
        nondefault_param_set (1) = .TRUE.

     case ('kmax')
        read(value_string,*) kmax
        call non_negative_integer(kmax, 'kmax')
        nondefault_param_set (2) = .TRUE.

     case ('maxc')
        read(value_string,*) maxc
        call non_negative_integer(maxc, 'maxc')
        nondefault_param_set (3) = .TRUE.

     case ('scf_integral_filename')
        read(value_string,*) SCF_integral_filename
        SCF_integral_filename = trim(adjustl(SCF_integral_filename))
        select case (trim(SCF_integral_filename))
        case ('moints.TM')
        case ('moints.ascii')
        case ('FCIDUMP')
        case default
           if (me.eq.0) then
              write(50,*)
              write(50,*) 'SCF_integral_filename must be moints.TM, moints.ascii or FCIDUMP.'
              write(50,*) 'Stopping...'
           endif
           stop
        end select
        nondefault_param_set(4) = .TRUE.

     case default

     end select

  end do
     
7000 continue

  if ( any( .NOT. nondefault_param_set) ) then
     if (me.eq.0) then
        write(50,*)
        write(50,*) 'Problem in subroutine read_params.f90'
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

  select case (trim(SCF_integral_filename))

  case('moints.ascii')
     open(20, file='moints.ascii', form='formatted')
     read(20,'(6i9)') rubbish(1), irmax, maxbfs, rubbish(2), rubbish(3), rubbish(4)
     close(20)
  
  case ('FCIDUMP')
  call get_Molproparams(maxbfs,irmax) 

  case('moints.TM')
     open(20,file='moints.TM',form='formatted')
     rewind(20)
     do i=1,40
        read(20,*)
     enddo

     read(20,*) irmax
     allocate(ityp   (irmax))
     allocate(idim   (irmax))
     allocate(nlamda(irmax))
     allocate(norbs  (irmax))
     do i=1,6  ! Skip blank, "for each irr", "----", blank, "number :", blank
          read (20,*)
     end do
     ityp  (:) = 0 ! Edited from ='     '
     idim  (:) = 0
     nlamda(:) = 0
     norbs (:) = 0

     maxbfs = 0

     do isym = 1,irmax

        read(20,1000) isym_temp, ityp(isym), idim(isym),& 
                       nlamda(isym), norbs(isym)


        maxbfs = maxbfs + idim(isym)*norbs(isym)

        if (isym /= isym_temp) then
           write(50,*)
           write(50,*) 'Problem in get_int_skip_TM.f'
           write(50,*) 'isym     =',isym
           write(50,*) 'isym_temp=',isym_temp
           write(50,*) 'These should be equal.'
           write(50,*) 'Stopping...'
           stop
        end if
     end do
     close(20) 
  end select

  iword  = maxbfs/int_bits + 1
  max1   = maxbfs*(maxbfs+1)/2
  max2   =  max1*(max1+1)/2
  maxocc = 2*maxbfs
  maxs   = maxh/50

      
8010 format (1x, a, i6)             ! integers
1000 format(2x,i3,7x,a4,6x,i3,8x,i5,5x,i5) !moints.TM
    
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
          write(50,*) 'Problem in subroutine read_params.f90'
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
          write(50,*) 'Problem in subroutine read_params.f90'
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
  
end subroutine read_params

