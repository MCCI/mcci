subroutine init_bk(ieig,length,inflg)
  use commonarrays, only: c, b
  use dyn_par
  implicit     real*8 (a-h,o-z)
  !common /bk/   b(maxc,kmax)
  !common /coef/ c(maxc)  
  integer       i, j, k, ieig, length, inflg

  do j=1, kmax
     do i= 1, maxc
        b(i,j) = 0.0
        if (i.eq.j) b(i,j) = 1.0
     enddo
  enddo

  if (inflg .ne. 0) then
     do i= 1, length
        b(i,ieig) = c(i)     ! copy input civ_in vector to b(:,ieig)
        !            b(i,1) = c(i)        ! copy input civ_in vector to b(:,1)
     enddo
  end if

  !     orthogonalize b_k vectors  
  do kk = 2, kmax
     do k=kk-1,1, -1
        dot = 0.0
        do i = 1, maxc
           dot = dot + b(i,k)*b(i,kk)
        enddo
        do i = 1, maxc
           b(i,kk) = b(i,kk) - dot*b(i,k)
        enddo
     enddo
     vnorm = 0.0
     do i = 1, maxc
        vnorm = vnorm + b(i,kk)*b(i,kk)
     enddo
     vnorm = dsqrt(vnorm)
     do i= 1,length
        b(i,kk) = b(i,kk)/vnorm
     enddo
  enddo

  return
end subroutine init_bk
