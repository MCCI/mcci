      subroutine perturb(lb4,length,ea)
      use commonarrays, only: c, h, ijh, s, ijs

      implicit  real*8  (a-h,o-z)
      include           'params'
      !common /hands/     h(maxh), ijh(maxh), s(maxs), ijs(maxs)
      !common /coef/      c(maxc)

      if(ijh(1) .ne. length + 2) STOP 'energy: h array mismatch'
      if(ijs(1) .ne. length + 2) STOP 'energy: s array mismatch'

      eval  = 0.0
      dnorm = 0.0
      do k=lb4branch+1, length
         dnom=h(k)-ea
         nom1=0.0
         do jci = ijs(k), ijs(k)+lb4-1
            nom1  = nom1  + ea*s(jci)*c(ijs(jci))
         enddo
         nom2=0.0
         do jci = ijh(k), ijh(k)+lb4-1
            nom2  = nom2  + h(jci)*c(ijh(jci))
         enddo
         c(k)=(nom1-nom2)/dnom
      enddo

      return
      end
