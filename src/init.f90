subroutine init(seed,ecore,inflg,ieig)

  use commonarrays 
  use dyn_par
  implicit real*8   (a-h,o-z)

  real*8       seed,ecore
  integer      inflg, ieig
  character*12 SCF_integral_filename

!  common  /aodat/   nsym, nbft, nbpsy(irmax)
!  common  /occupy/  ntotal, n_alpha, n_beta, i_sx2
!  common  /config/  icij(2,iword,maxc), nword
!  common  /froze/   nfreeze, ifreeze(maxocc)
!  common  /active/  nactive, iactive(maxocc)
!  common  /para/    me, nproc

  logical is_frozen(maxbfs)


  call sym_init

  call gen_seed(seed)

  call read_mcci_in(iword,maxc,maxocc,maxbfs,icij,&
       inflg,n_alpha,n_beta,ntotal,i_sx2,SCF_integral_filename,&
       ieig,nfreeze,ifreeze,nactive,iactive)

  if (me.eq.0) write(50,*)

  if (SCF_integral_filename .eq. 'moints.ascii') then
     if (me.eq.0) write(50,*) 'Getting the one- and two-electron integrals from moints.ascii'
     call get_int(ecore)
  end if

  if (SCF_integral_filename .eq. 'moints.TM') then
     if (me.eq.0) write(50,*) 'Getting the one- and two-electron integrals from moints.TM'
     call get_int_TM(ecore)
  end if

  if (SCF_integral_filename .eq. 'moints.bTM') then
     if (me.eq.0) write(50,*) 'Getting the one- and two-electron integrals from moints.bTM'
     call get_int_bTM(ecore)
  end if

  if (SCF_integral_filename .eq. 'FCIDUMP') then
     if (me.eq.0) write(50,*) 'Getting the one- and two-electron integrals from Molpro FCIDUMP'
     call get_intMolpro(ecore)
  end if

  !     now we know nbft, set up the defaults for the complete active space part,
  !     where if nactive has been set to 0 or left unset, all nonfrozen orbitals are active

  if (nactive == 0) then

     nactive = nbft-nfreeze

     do i = 1,maxbfs
        is_frozen(i) = .FALSE.
     end do

     do i=1,nfreeze
        is_frozen(ifreeze(i)) = .TRUE.
     end do

     icount = 0
     do i = 1,nbft
        if (.not. is_frozen(i)) then
           icount=icount+1
           iactive(icount) = i
        end if
     end do

  end if

  return  
end subroutine init
