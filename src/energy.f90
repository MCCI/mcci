subroutine energy(length,eval,dnorm)
  use commonarrays, only: c, h, ijh, s, ijs
  use dyn_par
  implicit  real*8  (a-h,o-z)       
  !common /hands/     h(maxh), ijh(maxh), s(maxs), ijs(maxs)
  !common /coef/      c(maxc) 

  if(ijh(1) .ne. length + 2) STOP 'energy: h array mismatch'
  if(ijs(1) .ne. length + 2) STOP 'energy: s array mismatch'

  eval  = 0.0
  dnorm = 0.0
  do ici=1, length
     eval  = eval  + c(ici)*h(ici)*c(ici)
     dnorm = dnorm + c(ici)*s(ici)*c(ici)
     do jci = ijh(ici), ijh(ici+1) - 1
        eval  =  eval  + 2.0*c(ici)*h(jci)*c(ijh(jci))
     enddo
     do jci = ijs(ici), ijs(ici+1) - 1
        dnorm =  dnorm + 2.0*c(ici)*s(jci)*c(ijs(jci))
     enddo
  enddo

  eval = eval/dnorm

  return
end subroutine energy
