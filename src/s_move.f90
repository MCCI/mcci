subroutine s_move(length,llast) 
  use commonarrays, only: s, ijs
  use dyn_par
  implicit  real*8    (a-h,o-z)
  !common   /hands/   h(maxh), ijh(maxh), s(maxs), ijs(maxs)

  !     this subroutine adjusts the s and ijs arrays to allow
  !     new configurations to be inserted during the next
  !     call to gen_h_s

  if( llast .ne. ijs(1)-2 ) STOP 'error h_move: llast doesn''t match stored length'

  !     number of configurations to increase/decrease
  ldiff = length-(ijs(1)-2)

  mv_begin = ijs(1)             ! begin of off diagonals
  mv_end   = ijs(ijs(1)-1)-1    ! end   of off diagonals

  if(ldiff.gt.0) then
     !        move off-diagonal information up for s
     !        move row information up for ijs
     do k=mv_end, mv_begin, -1
        ijs(k+ldiff) = ijs(k) 
        s(k+ldiff)   =   s(k)
     enddo
  elseif(ldiff.lt.0) then
     !        move off-diagonal information down for s
     !        move row information down for ijs
     do k=mv_begin,mv_end
        ijs(k+ldiff) = ijs(k) 
        s(k+ldiff)   =   s(k)
     enddo
  endif

  !     correct index information
  if(ldiff.gt.0) then
     do k=1,llast+1
        ijs(k) = ijs(k) + ldiff
     enddo
  elseif(ldiff.lt.0) then
     do k=1,llast+ldiff+1
        ijs(k) = ijs(k) + ldiff
     enddo
  endif

  if(ijs(1) .ne. length + 2) STOP 's_move: troubles'

  return
end subroutine s_move
