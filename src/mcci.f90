program mcci
  use commonarrays
  use dyn_par

  implicit real*8   (a-h,o-z)
  include          'mpif.h'
  real*4            twb,twg1,twp0,twp,twi,twf,dwsum
  real*4            tub,tug1,tup0,tup,tui,tuf,dusum
  real*4            tsb,tsg1,tsp0,tsp,tsi,tsf,dssum
  real*4            t_cumulative
  character*24      date
  character*255     cmin_char,civ_command
  logical           prune_all
  real*8            maxde,maxdn
  real*8,allocatable :: w(:)
  real*8            state(7),stateave(2),de(7)
  real*8            statnave(2),dn(7)
  integer           statn(7)
  integer           ialloc(22)
  logical           ecnvrgd,ncnvrgd,crit,nobranch_flag,cmin_flag
  integer           irefractor

  ! call pbeginf                            ! TCGMSG
  call mpi_init(ierr)                       ! MPI

  ! who am i and how many others are there?
  ! me = nodeid()                             ! TCGMSG
  ! nproc = nnodes()                          ! TCGMSG
  call mpi_comm_rank(MPI_COMM_WORLD, me, ierr)    ! MPI
  call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr) ! MPI

  irefractor  = 0
  cmin_flag = .false.
  if(me.eq.0) then
     open(50,file='e_summary')
     write(50,'(24a)') fdate()
     call header
  end if

  call read_params()
  call allocate_memory
  call init(seed,ecore,inflg,ieig)
  ! initialize variables for diagonalization steps
  c(1)   = 1.0
  length = 1
  dnorm  = 1.0
  if (inflg.eq.0) then
     call h_s_sparse(1,0)
     ! call dump('hdump',h,ijh)
     ! call dump('sdump',s,ijs)
     eref   = h(1)/s(1)
     vnorm  = dsqrt(s(1))
     c(1)   = c(1)/vnorm
  else
     call restart(length)
     call chk_list(length,0)
     call h_s_sparse(length,0)
     eref   = h(1)/s(1)
  endif

  write(*,*) '***initialization has been done***'

  lb4prune = length
  last_ok =  length

  !Write info to e_summary
  if(me.eq.0) call initial_info

  !for branching with nproc>1, how many extra configs before an exchange?
  f_boost = bmin

  t_cumulative = 0.0

  call init_bk(ieig,length,inflg)

  !Divide frac among processors

  open(51,file='fracValue')
  write(51,*) 'frac before', frac
  frac = frac/nproc
  temp_frac = frac
  write(51,*) 'frac after',temp_frac
  write(51,*) 'nproc', nproc
  close(51)
  nobranch_flag = .false.

  !This is the loop responsible for iterating the cycles in mcci.  In each cycle,
  !new configurations are appended onto the current CI vector.  The resultant
  !hamiltonian is diagonalised, and poor configurations are pruned from the vector.
  i_try = 1
  do while (i_try.le.maxtry)

     !Timing data is will be written to e_summary
     if(me.eq.0.and.time) call timer(twi,tui,tsi)

     if(me.eq.0) then
        open(50,file='e_summary',status='old',position='append')
        write(50,*)
        write(50,*)
        write(50,*)'Diagonalization',i_try
        write(50,*)
     endif

     if(i_try .lt. maxtry) then

        !The number of new configurations to be appended
        nu_configs = int(frac*dble(length))
        if(nu_configs.lt.1) nu_configs=0
        !The number of new configurations must be large enough to
        !satisfy lmin
        if(nu_configs.lt.lmin-length) nu_configs=lmin-length

        if(nproc .gt. 1) then
           ltarget    = length + nu_configs    ! desired length
           inxcss     = int(f_boost*dble(length))
           nu_configs = nu_configs + inxcss    ! configs+excess
        else
           ltarget    = length + nu_configs    ! desired length 
        endif

        call gen_seed(seed)

        if(nobrnch_first.and.(i_try.eq.1)) then
           nu_configs=0
           ltarget=length
        endif

        call branch(nu_configs,length,lb4brnch,seed)

        call genealogy(length,lb4brnch,seed)

        call chk_list(length,lb4brnch)

        if(nproc.gt.1) then

           if(me.eq.0) then
              write(50,*)'Boosted CI vector length =',length
              write(50,*)'f boost= ',f_boost
           endif

           compare = dble(ltarget-length)/dble(lb4brnch)

           ! adjust boost keeping within bounds
           f_boost = min(bmax,max(bmin,f_boost+compare))

           if(length.gt.ltarget) length = ltarget

        endif

        if(me.eq.0 .and. time_all) call timer(twb,tub,tsb)

     endif

     call h_move(length,lb4prune)
     call s_move(length,lb4prune)

     if(length.gt.last_ok) then 
        call h_s_sparse(length,last_ok)
     else
        call h_s_sparse(length,0)
     endif
     ! call dump('hdump',h,ijh)
     ! call dump('sdump',s,ijs)

     if(me.eq.0 .and. time_all) call timer(twg1,tug1,tsg1)

     dwsum  = 0.0
     dusum  = 0.0
     dssum  = 0.0
     if(nodiag) goto 10

     !The davidson routine is responsible for diagonalising the hamiltonian
     call davidson(length,ieig,idiag)

     ! write(50,*)'edavidson', e(ieig) + ecore
     eval = e(ieig)

     if(me.eq.0) then
        write(50,*)'Iterations to convergence',idiag
        write(50,1002) eval+ecore
        write(50,*)'Branched  CI vector length=',length
     endif

     entot = eval+ecore

     call energy(length,eval,dnorm)

     call mxv_sparse(length,h,ijh,c,w)

10   continue

     if(me.eq.0) call energy_eval 

     if(i_try.lt.maxtry)then
        prune_all = .false. 
        if(me.eq.0 .and. time_all) call timer(twp0,tup0,tsp0)
        if( mod(i_try,npfull).eq.0.or.i_try.eq.(maxtry-1).or.i_try.eq.1) prune_all = .true.

        lb4prune = length

        call prune(length,lb4brnch,dnorm,eval,i_got_hit,prune_all,ieig)

        if(me.eq.0) then
           if(prune_all) then
              write(50,*)'Pruned    CI vector length=',length,&
                   '   *** Full pruning ***'              
              lconv = length
           else
              write(50,*)'Pruned    CI vector length=',length
           endif
        endif


        if(i_got_hit.ne.maxc) then
           last_ok = i_got_hit-1
        else
           last_ok = length                  ! no pruning took place
        endif

        if(nproc.gt.1) then
           lb4exc = length
           do jj=0, nproc-1
              call exc(lb4exc,lb4brnch,length,jj)    !parallel
           enddo
           if(me.eq.0) write(50,*)'Exchanged CI vector length=',length
           call chk_list(length,lb4exc)
        endif


        if(me.eq.0 .and. time_all) call timer(twp,tup,tsp)

        if(me.eq.0) then
           if(time_all) then 
              write(50,*)'BRANCHING                                    total'
              write(50,2004) twb - twi
              write(50,2005) tub - tui
              write(50,2006) tsb - tsi
              write(50,*)'H GENERATION FOR BRANCHED CONFIGURATIONS     total'
              write(50,2004) twg1 - twb
              write(50,2005) tug1 - tub
              write(50,2006) tsg1 - tsb
              write(50,*)'H DIAGONALIZATION average                    total'
              write(50,2001) dwsum/idiag,dwsum
              write(50,2002) dusum/idiag,dusum
              write(50,2003) dssum/idiag,dssum
              write(50,*)'PRUNING                                      total'
              write(50,2004) twp - twp0
              write(50,2005) tup - tup0
              write(50,2006) tsp - tsp0
           endif
        endif
     endif

     if(me.eq.0) then
        if(time) then 
           call timer(twf,tuf,tsf)
           twf = twf - twi
           tuf = tuf - tui
           tsf = tsf - tsi
           write(50,2007) twf
           write(50,2008) tuf
           write(50,2009) tsf
           t_cumulative = t_cumulative + twf
           write(50,2010) t_cumulative
           write(50,'(24a)') fdate()
        endif
        close(50)
     endif

     !The root processor checks for convergence
     if(me .eq. 0 .and. maxtry-i_try .gt. 2) then
        if(frac == 0.0) then
           open(50,file='e_summary',status='old',position='append')
           write(50,*) '========================='
           write(50,*) 'Converged for Cmin:',cmin
           write(50,*) '========================='
           close(50)

           write(cmin_char,'(f8.5)') cmin
           civ_command = "cp civ_out civ_out_" // trim(adjustl(cmin_char))
           call system(civ_command)          

           open(110,file='cmin_dat',status='old',position='append')
           write(110,*) i_try,cmin,entot,lconv
           close(110)

           if(.not.auto_cmin .or. cmin <= cmin_finish) then
              maxtry = i_try+2
           else if(cmin > 0.1) then
              cmin = cmin-0.1
           else if(cmin > 0.01) then
              cmin = cmin-0.001
           else if(cmin > 0.001) then
              cmin = cmin-0.0001 
           else if(cmin > 0.0001) then 
              cmin = cmin-0.00001
           else
              maxtry = i_try + 2
           end if


        end if
     end if
     frac = temp_frac
     if(me.eq.0 .and. i_want_conv) call convergence_test

     if(nproc.gt.1) then
        call MPI_BCAST(maxtry,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(cmin,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(frac,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
     endif
     i_try = i_try + 1

  enddo

1000 format(1x,'Eref= ',f20.14)
1001 format(1x,'Ecor= ',f20.14,' from Davidson diag.')
1002 format(1x,'E   = ',f20.14,' from iter. diag.')
1003 format(1x,'Egnd= ',f20.14,' end of iteration ')
1004 format(1x,'Ecor= ',f20.14,' best guess on a node from cHc',' after pruning') 
1005 format(1x,'Ecor= ',f20.14,' starting guess from cHc')

2001 format(1x,'wall time',f12.3,' s      wall time',f12.3,' s')
2002 format(1x,'user time',f12.3,' s      user time',f12.3,' s')
2003 format(1x,'sys  time',f12.3,' s      sys  time',f12.3,' s')

2004 format(1x,'         ', 12x ,'        wall time',f12.3,' s')
2005 format(1x,'         ', 12x ,'        user time',f12.3,' s')
2006 format(1x,'         ', 12x ,'        sys  time',f12.3,' s')

2007 format(1x,'Cycle WALL time',f12.3,' s')
2008 format(1x,'      USER time',f12.3,' s')
2009 format(1x,'      SYS  time',f12.3,' s')

2010 format(1x,'Cumul. run time',f12.3,' s')

  !  call pend                              ! TCGMSG
  call mpi_finalize(ierr)                   ! MPI


  !End of mcci.f90


  contains

  subroutine allocate_memory
     allocate(e(kmax),             stat=ialloc(1))
     allocate(icij(2,iword,maxc),  stat=ialloc(2))
     allocate(h(maxh),             stat=ialloc(3))
     allocate(s(maxs),             stat=ialloc(4))
     allocate(ijh(maxh),           stat=ialloc(5))
     allocate(ijs(maxs),           stat=ialloc(6))
     allocate(c(maxc),             stat=ialloc(7))
     allocate(ifreeze(maxocc),     stat=ialloc(8))
     allocate(iactive(maxocc),     stat=ialloc(9))
     allocate(hf(kmax*(kmax+1)/2), stat=ialloc(10))
     allocate(sf(kmax*(kmax+1)/2), stat=ialloc(11))
     allocate(b(maxc,kmax),        stat=ialloc(12))
     allocate(irrep(0:irmax-1),    stat=ialloc(13))
     allocate(list(2,maxocc),      stat=ialloc(14))
     allocate(my_pair(2,maxocc),   stat=ialloc(15))
     allocate(icase(16),           stat=ialloc(16))
     allocate(ipoint(max2),        stat=ialloc(17))
     allocate(e1ints(max1),        stat=ialloc(18))
     allocate(e2ints(max2),        stat=ialloc(19))
     allocate(nbpsy(irmax),        stat=ialloc(20))
     allocate(cnorm(maxc),         stat=ialloc(21))
     allocate(w(maxc),             stat=ialloc(22))
     if(any(ialloc /=0)) STOP "Error allocating memory in mcci.f90"
  end subroutine allocate_memory

  !Initial information regarding the branching factor etc. is
  !written to e_summary
  subroutine initial_info
     write(50,*)'Calculating',ieig,'th state in this irrep'
     write(50,*)'Running on',nproc,' nodes'
     write(50,*)'Branching factor      f =',frac
     write(50,*)'Davidson tolerance stop =',davidson_stop
     write(50,*)'Coef.    tolerance cmin =',cmin
     write(50,*)'H        tolerance hmin =',hmin
     if(inflg.ne.0)write(50,*)'RESTARTED'
     close(50)
  end subroutine initial_info

  subroutine energy_eval
     open(40,file='civ_out',form='formatted')
        open(60,file='weight',form='formatted')
        write(60,*) 'Energy = ', eval+ecore
        eval = 0.0
        do i=1,length
           do n=1,nword
              if(n.eq.1) then
                 write(40,'(i6,2x,e24.17,2x,i11,2x,i11)')&
                      i,c(i)/dsqrt(dnorm),icij(1,n,i),icij(2,n,i)
              else
                 write(40,'(33x,i11,2x,i11)')&
                      icij(1,n,i),icij(2,n,i)
              endif
           enddo
           x = c(i)*w(i)/dnorm
           write(60,'(i6,2x,e24.17,2x,e24.17)') i, x, x/e(ieig)
           eval = eval + x
        enddo
        write(60,*) 'Energy sum',eval, 'Core energy',ecore
        write(60,*) 'Energy check', eval+ecore
        close(40,err=111)
111     continue
        close(60)
  end subroutine energy_eval

  subroutine convergence_test

     if(i_try.eq.1) then
        do i=1,conv_average
           state(i)=0.0
           statn(i)=0
        enddo
        do i=1,2
           stateave(i)=0.0
           statnave(i)=0.0
        enddo
        do i=1,conv_history
           de(i)=0.0
           dn(i)=0.0
        enddo
        ecnvrgd  = .false.
        ncnvrgd = .false.
        open(70,file='convergence',form='formatted')
        write(70,*) '************CONVERGENCE TEST IS RUNNING***********'
        close(70)
        open(110,file='cmin_dat',form='formatted')
        write(110,*)'****CMIN DATA****'
        close(110)
     endif

     if(npfull_conv) then
        crit = (mod(i_try,npfull).eq.1).and.(.not.(ecnvrgd.and.ncnvrgd).and.i_try.gt.1)
        nnn = npfull
     else
        crit = (.not.(ecnvrgd.and.ncnvrgd)).and.(i_try.gt.1)
        nnn = 1
     endif

     if(i_try .eq. 1) then
        open(70,file='convergence',status='old',position='append')
        write(70,*) 'Convergence checking will begin after ',&
                    (conv_average+conv_history)*nnn, ' cycles.'
        write(70,*) '**************************************************'
        close(70)
     endif
        
     if (crit) then

        do i=1,conv_average-1
           state(i)=state(i+1)
        enddo
        state(conv_average)= entot
        stateave(1)=stateave(2)
        stateave(2)=0.0
        do i=1,conv_average
           stateave(2)=stateave(2)+state(i)
        enddo
        stateave(2)=stateave(2)/dble(conv_average)
        do i=1,conv_history-1
           de(i)=de(i+1)
        enddo
        de(conv_history)=stateave(2)-stateave(1)
        maxde=0.0
        do i=1,conv_history
           if(abs(de(i)).gt.maxde) maxde=abs(de(i))
        enddo

        do i=1,conv_average-1
           statn(i)=statn(i+1)
        enddo
        statn(conv_average)=lconv
        statnave(1)=statnave(2)
        statnave(2)=0.0
        do i=1,conv_average
           statnave(2)=statnave(2)+dble(statn(i))
        enddo
        statnave(2)=statnave(2)/dble(conv_average)
        do i=1,conv_history-1
           dn(i)=dn(i+1)
        enddo
        if (statnave(1).ne.0) then
           dn(conv_history)=statnave(2)/statnave(1) - 1.0
        else
           dn(conv_history)=999.
        endif
        maxdn=0.0
        do i=1,conv_history
           if(abs(dn(i)).gt.maxdn) maxdn=abs(dn(i))
        enddo

        if (i_try.ge.((conv_average+conv_history)*nnn)) then 
           ecnvrgd = .false.
           if(maxde.lt.conv_thresh_e) ecnvrgd=.true.
           ncnvrgd = .false.
           if(maxdn.lt.conv_thresh_l) ncnvrgd=.true.
           open(70,file='convergence',status='old',position='append')
           write(70,*) 'i_try:',i_try,'ecnvrgd:',ecnvrgd,'ncnvrgd:',ncnvrgd
           write(70,*) 'de:',de(conv_history),'dn:',dn(conv_history)
           ! write(70,*) 'maxde:',maxde,'maxdn:',maxdn
           close(70,err=112)
112        continue

           if(i_try < irefractor+10) then
              ecnvrgd = .false.
              ncnvrgd = .false.
           endif
           if((ecnvrgd.and.ncnvrgd)) then
              ecnvrgd = .false.
              ncnvrgd = .false.

              irefractor = i_try
              frac = 0.0
           endif
        endif
     endif

  end subroutine convergence_test

end program mcci

subroutine header
  write(50,*)
  write(50,*)'                   m c c i  3.0'
  write(50,*) 
  write(50,*)'              written by J.C. Greer'
  write(50,*)
  write(50,*)'   ============================================='
  write(50,*)'   J.C. Greer, J. Chem. Phys. 103 (1995) p. 1821'
  write(50,*)'   J.C. Greer, J. Comp. Phys. 146 (1998) p. 181 '
  write(50,*)'   ============================================='
  write(50,*)
  return
end subroutine header

