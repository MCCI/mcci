subroutine h_move(length,llast) 
  use commonarrays, only: h, ijh
  use dyn_par
  implicit  real*8    (a-h,o-z)
  !common   /hands/   h(maxh), ijh(maxh), s(maxs), ijs(maxs)

  !     this subroutine adjusts the h and ijh arrays to allow
  !     new configurations to be inserted during the next
  !     call to gen_h_m

  if( llast .ne. ijh(1)-2 ) STOP 'error h_move: llast doesn''t match stored length'

  !     number of configurations to increase/decrease
  ldiff = length-(ijh(1)-2)

  mv_begin = ijh(1)             ! begin of off diagonals
  mv_end   = ijh(ijh(1)-1)-1    ! end   of off diagonals

  if(ldiff.gt.0) then
     !        move off-diagonal information up for h
     !        move row information up for ijh
     do k=mv_end, mv_begin, -1
        ijh(k+ldiff) = ijh(k) 
        h(k+ldiff)   =   h(k)
     enddo
  elseif(ldiff.lt.0) then
     !        move off-diagonal information down for h
     !        move row information down for ijh
     do k=mv_begin,mv_end
        ijh(k+ldiff) = ijh(k) 
        h(k+ldiff)   =   h(k)
     enddo
  endif

  !     correct index information
  if(ldiff.gt.0) then
     do k=1,llast+1
        ijh(k) = ijh(k) + ldiff
     enddo
  elseif(ldiff.lt.0) then
     do k=1,llast+ldiff+1
        ijh(k) = ijh(k) + ldiff
     enddo
  endif

  if(ijh(1) .ne. length + 2) STOP 'h_move: troubles'

  return
end subroutine h_move
