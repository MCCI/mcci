subroutine exc(l_end,l_begin,l_to,iexc)
  use commonarrays, only: c, icij, me
  use dyn_par
  implicit real*8   (a-h,o-z) 
  !integer           me, nproc !declaration moved to commonarrays.f90
  integer           ntype, ibyte, ntime, lbuf, iexc
  integer           lpts, lptr, MAXBUF
  !integer           icij !declaration moved to commonarrays.f90
  integer           l_end, l_begin, l_to, ls, lr
  include          'mpif.h' 
  parameter         (MAXBUF = 1048576) ! must be multiple of 8

  !common  /coef/    c(maxc) 
  !common  /config/  icij(2,iword,maxc), nword
  !common  /para/    me, nproc


  !     l_end       end of buffer to broadcast from node iexc
  !     l_begin+1   beginning of buffer to broadcast from node iexc
  !     l_to        end of buffer on receiving nodes, upon return l_to=l_to+lr

  ntype = 1

  ls = l_end - l_begin      ! buffer to send   / 8
  if(ls.lt.0) ls = 0

  ntype = ntype + 1
  if(me.eq.iexc) then
     !        call brdcst(ntype,ls,4,iexc)                          ! TCGMSG
     call MPI_BCAST(ls,4,MPI_BYTE,iexc,MPI_COMM_WORLD,ierr)! MPI
     ibyte = 2*nbyte*iword*ls
  else
     !        call brdcst(ntype,lr,4,iexc)		               ! TCGMSG
     call MPI_BCAST(lr,4,MPI_BYTE,iexc,MPI_COMM_WORLD,ierr)! MPI
     ibyte = 2*nbyte*iword*lr
  endif

  ntime = ibyte/MAXBUF + 1

  idiff = ibyte - (ntime-1)*MAXBUF

  lpts = l_begin
  lptr = l_to

  do k = 1, ntime
     lbuf = MAXBUF
     if(k .eq. ntime) lbuf = idiff
     ntype = ntype + 1
     if(me.eq.iexc) then
        !           call brdcst(ntype,icij(1,1,lpts+1),lbuf,iexc)         ! TCGMSG
        call MPI_BCAST(icij(1,1,lpts+1),lbuf,MPI_BYTE,iexc,MPI_COMM_WORLD,ierr)     ! MPI
     else
        !           call brdcst(ntype,icij(1,1,lptr+1),lbuf,iexc)         ! TCGMSG
        call MPI_BCAST(icij(1,1,lptr+1),lbuf,MPI_BYTE,iexc,MPI_COMM_WORLD,ierr)   ! MPI
     endif
     lpts = lpts + MAXBUF/2/nbyte/iword
     lptr = lptr + MAXBUF/2/nbyte/iword
  enddo

  if(me.ne.iexc) then 
     do ici = l_to+1, l_to+lr+1
        c(ici) = 0.0
     enddo
     l_to = l_to + lr
  endif

  if(l_to.gt.maxc-1) STOP 'error exc: exceeded max. array length'

  return
end subroutine exc
