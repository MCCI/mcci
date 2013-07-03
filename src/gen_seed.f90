subroutine gen_seed(seed)
  use commonarrays, only: me
  integer         nseed  !, me, nproc
  real*8          seed
  !common /para/   me,nproc

  !     generate seed for random number generator
  call system_clock(nseed)                      ! f90 intrinsic
  !     nseed= 18800925
  !     if(me.eq.0) write(50,*) 'Warning: seed set by hand'
  nseed = nseed / ( me + 1)
  nfrac = nseed/199017
  nseed = nseed - nfrac*199017
  seed  = float(nseed)

  return
end subroutine gen_seed
