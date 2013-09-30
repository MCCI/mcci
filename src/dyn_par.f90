module dyn_par

!     As we move parameters from hardcoding (params) to dynamic setting, we place them
!     in a mdoule here.
!     This makes them easy to pass to the parts of mcci.f90 that need them - just follow

      integer ::  maxbfs,irmax,max1,max2,maxocc,maxc,kmax,maxh,maxs,iword

      integer maxtry          ! the number of diagonalisations
      integer lmin            ! min vector length
      integer nbyte           ! number of bytes / integer word
      integer int_bits        ! bits / integer word
      integer npfull          ! perform full pruning after npfull steps
      integer lref            ! length of ref vector; <=length; if =0 then =length
      integer lkeep           ! cnfgs #1 to #lkeep never touched in pruning
      integer conv_average    ! for making E(iter) and length(iter) func. smoother
      integer conv_history    ! how many DE amd Dlength values being tracked

      real*8  cmin
      real*8  cmin_finish
      real*8  hmin            ! h matrix threshold
      real*8  davidson_stop   ! stop for convergence in davidson
      real*8  bmin            ! min vector boost
      real*8  bmax            ! max vector boost
      real*8  cref            ! c>cref is treated as a ref. for branching
      real*8  frac            ! rand < frac for branching
      real*8  conv_thresh_e   ! <DE> < conv_thresh_e to stop
      real*8  conv_thresh_l   ! <Dlength> < conv_thresh_l to stop

      logical auto_cmin
      logical test            ! error checking
      logical time            ! timing information, general
      logical time_all        ! timing information, detailed
      logical generate_cfgs   ! running with/without generating new cfgs
      logical nobrnch_first   ! no branching in the first step
      logical nodiag          ! only for collecting CSFs
      logical i_want_conv     ! stop mcci by using convergence crit.
      logical npfull_conv     ! convergence test only in npfull steps
      logical run_pt2             ! Run mccipt2 after mcci calculation

end module dyn_par
