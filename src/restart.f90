subroutine restart(length)
  use commonarrays, only: c, icij, nword
  use dyn_par
  integer           length, idummy !icij !declaration moved to commonarrays.f90 
  integer*8  bigtemp1,bigtemp2 !to prevent gfortran problem if reading in most negative integer.
  !common  /config/  icij(2,iword,maxc), nword
  !common  /coef/    c(maxc)

  open(unit=60,file='civ_in')
  i = 0
11 continue
  i = i + 1	
  do n=1, nword
     !        read(60,*,end=22) idummy,icij(1,n,i),icij(2,n,i)
     if(n.eq.1) then
        !           read(60,*,end=22) c(i),icij(1,n,i),icij(2,n,i)
        read(60,*,end=22) idummy,c(i),bigtemp1,bigtemp2
        icij(1,n,i)=bigtemp1
        icij(2,n,i)=bigtemp2
     else
         read(60,*,end=22) bigtemp1,bigtemp2
        icij(1,n,i)=bigtemp1
        icij(2,n,i)=bigtemp2 
     endif
  enddo
  goto 11
22 continue
  length = i-1

  return
end subroutine restart
