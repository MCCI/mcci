subroutine davidson(length,ieig,idiag)
  !use commonarrays, only: c
  use commonarrays, only: hf, sf, b, h, ijh, s, ijs, c, e
  use dyn_par
  !implicit none
  implicit   real*8  (a-h,o-z)
  parameter           (eps=1.0d-18,maxit=100)
  !common  /reduced/   hf(kmax*(kmax+1)/2), sf(kmax*(kmax+1)/2)
  !common  /bk/        b(maxc,kmax)
  !common  /hands/     h(maxh), ijh(maxh), s(maxs), ijs(maxs)
  !common  /coef/      c(maxc)
  !common  /eigen/     e(kmax)
  integer             length, info, idiag
  integer             i, ici, k, kk, kkk, klim
  real*8              z(kmax,kmax), btemp(maxc,kmax)
  real*8              vnorm, rnorm
  real*8              a(maxc), d(maxc), r(maxc)
  real*8              work(3*kmax)

  !     write(34,*)
  !     write(34,*) 'entering davidson...'

  idiag = 1

111 continue

  klim = min(kmax,length)
  do kk=ieig, klim
     !     generate reduced h and s matrices need to fix lower limit
     call h_s_reduced(length,1,kk)

     !     solve nonorthogonal eigenvalue problem using routine from LAPACK
     !     hred * z = e * sred * z
     call dspgv(1 , 'V' ,'U', kk, hf, sf, e, z, kmax, work, info)
     if(info.ne.0) THEN
     PRINT *,kk,kmax
     PRINT *, 'Problems with dspgv '
     END IF
     do ici=1, maxc
        r(ici) = 0.0
     enddo

     do k=1, kk
        !     note that a and d are calculated in h_s_reduced,
        !     but are currently not being saved for use 
        !     a = H * b
        call mxv_sparse(length,h,ijh,b(1,k),a)
        !     d = S * b
        call mxv_sparse(length,s,ijs,b(1,k),d)
        !     calculate residual 
        do ici=1, length
           r(ici) = r(ici)+z(k,ieig)*(a(ici)-e(ieig)*d(ici))
        enddo
     enddo

     !     Undo the premultiplication by S: S_inv might be hard to calculate, so...
     do ici=1,length
        r(ici)=r(ici)/s(ici) 
     end do
     !     approximate correction to get closer to the true residual: r is now the
     !     exact residual if S is diagonal which should be nearly true
     !     (S diagonally dominant) - PD. Does not seem to make much difference to
     !     convergence anyway.

     !     residual norm
     vnorm = 0.0
     do ici = 1, length
        vnorm = vnorm + s(ici)*r(ici)*r(ici) ! approximate norm:diag. dom.
     enddo
         rnorm = dsqrt(vnorm)
!        write(34,*) 'idiag,kk,rnorm=',idiag,kk,rnorm

         if(rnorm.lt.davidson_stop) goto 222 ! finished?
         
!     form correction vector and orthogonalize to b_ks
         if(kk.lt.klim) then
            do ici=1, length
               if( dabs(r(ici)) .lt. eps) then
                  b(ici,kk+1) = 0.0
               else
                  b(ici,kk+1) = r(ici)/(e(ieig)*s(ici)-h(ici))
               endif
            enddo
            vnorm = 0.0
            do ici = 1, length
               vnorm = vnorm + b(ici,kk+1)*b(ici,kk+1)
            enddo
            vnorm = dsqrt(vnorm)
            do ici= 1,length
               b(ici,kk+1) = b(ici,kk+1)/vnorm
            enddo
!     orthogonalize b_k vectors               
            do k=kk,1, -1
               dot = 0.0
               do ici = 1, length
                  dot = dot + b(ici,k)*b(ici,kk+1)
               enddo
               do ici = 1, length
                  b(ici,kk+1) = b(ici,kk+1) - dot*b(ici,k)
               enddo
            enddo
            vnorm = 0.0
            do ici = 1, length
               vnorm = vnorm + b(ici,kk+1)*b(ici,kk+1)
            enddo
            vnorm = dsqrt(vnorm)
            do ici= 1,length
               b(ici,kk+1) = b(ici,kk+1)/vnorm
            enddo
         endif
         
      enddo                     ! if allowed, go back again with one more
                                ! b_k vector and rediagonalise in this new
                                ! space.
      
!     If we get to here, we've brought in as many new b_k vectors as we
!     can, and we're still not converged. So, form the first ieig 
!     eigenvectors from the last LAPACK diagonalisation, make these the
!     new b_k for k=1...ieig, and forget the other b_k vectors.
      
!     generate new b_k vectors from eigenvectors of hred
      do i=1, maxc
         do k=1, kmax
            btemp(i,k) = 0.0
         enddo
      enddo
      do i= 1, ieig
         do ici=1, length
            do k=1, klim        ! we've just expanded in this many b_k's
               btemp(ici,i) = btemp(ici,i) + z(k,i)*b(ici,k)
            enddo
         enddo
      enddo
!     We guess these eigenvectors are linearly independent: true if they
!     come from different eigenvalues, and likely to be true other wise
!     if the LAPACK routine is doing its job.
      
!     normalize b_k vectors ! don't know why, or the orthogonalisation either
!     - no difference in space spanned. In fact it means the z matrix, the
!     next times around will not be as close to the identity as it could be
!     (on that part where k1,k2 <= ieig).
      do k=1, ieig
         vnorm = 0.0
         do ici=1, length
            b(ici,k) = btemp(ici,k)
            vnorm  = vnorm + b(ici,k)*b(ici,k)
         enddo
         vnorm = dsqrt(vnorm)
         do ici=1, length
            b(ici,k) = b(ici,k)/ vnorm
         enddo
         do ici=length+1, maxc
            b(ici,k) = 0.0
         enddo
      enddo
!     orthogonalize the new b_k vectors
      do kkk=2, ieig
         do k=kkk-1, 1, -1
            dot = 0.0
            do ici = 1, length
               dot = dot + b(ici,k)*b(ici,kkk)
            enddo
            do ici = 1, length
               b(ici,kkk) = b(ici,kkk) - dot*b(ici,k)
            enddo
         enddo
         vnorm = 0.0
         do ici=1, length
            vnorm = vnorm + b(ici,kkk)*b(ici,kkk)
         enddo
         vnorm= dsqrt(vnorm)
         do ici = 1, length
            b(ici,kkk) = b(ici,kkk)/vnorm
         enddo
      enddo
      
!     We know we're not converged.
      idiag = idiag + 1 
      if(idiag.gt.maxit) STOP 'TOO MANY DAVIDSON ITERATIONS'
      goto 111                  ! start afresh with this older, wiser set
                                ! of 1..ieig b_k's. 
      
 222  continue                  ! jump to here when converged.
      
!     We guess that the next call to davidson will be with a hamiltonian
!     not too different from this one: so the final b_k, k=1..ieig will
!     still be good approximate eigenvectors of the new projected hamiltonian.
!     So, save these to b(ici,k),k=1..ieig, so the first LAPACK
!     diagonalisation next Hamiltonian will (if the Hamiltonian were exactly
!     the same) would give us the same approximate Psi, and the same
!     residual as we have now (ie. ~0 residual). This stops us redoing work
!     next call that's done already. 
!     As it is, if we've run over klim LAPACK calls, we've already saved
!     the b_k, k=1..ieig, to memory at this point, and not changed them
!     since, so we have a not too bad space already, and any subsequent calls
!     will start with this space, so with the residual set to the value at
!     the point where we saved the (Gram-Schmidt orthonormalised) b_k. We
!     could do better, though, and save the b_k (k=1..ieig) at the very
!     end, and so start next time round with this residual (modulo the
!     Hamiltonian change) which should be smaller.

!     generate new b_k vectors from eigenvectors of hred
      do i=1, maxc
         do k=1, kmax
            btemp(i,k) = 0.0
         enddo
      enddo
      do i= 1, ieig
         do ici=1, length
            do k=1, kk        ! we've just expanded in this many b_k's
               btemp(ici,i) = btemp(ici,i) + z(k,i)*b(ici,k)
            enddo
         enddo
      enddo
!     We guess these eigenvectors are linearly independent: true if they
!     come from different eigenvalues, and likely to be true otherwise
!     if the LAPACK routine is doing its job.

      do k=1, ieig
         do ici=1, maxc         ! better than length, for branching
            b(ici,k) = btemp(ici,k)
         enddo
      enddo
 
!     update c vector            
      do ici=1, maxc
         c(ici) = 0.0
      enddo

      do ici= 1, length
         c(ici) = b(ici,ieig)
      enddo
      
      return
      end

