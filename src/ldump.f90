subroutine ldump(i_am_nu,i_am_mu,nu_doubly)
  use commonarrays
  use dyn_par

  implicit real*8    (a-h,o-z)
!  common  /mxcoinc/  list(2,maxocc)
!  common  /aodat/    nsym, nbft, nbpsy(irmax)
!  common  /occupy/   ntotal, n_alpha, n_beta, i_sx2
!  common  /pairs/    my_pair(2,maxocc)

  !     dump maximum coincidence lists to standard output

  do nn=1,ntotal

     if(   list(i_am_nu,nn).le.nbft.and.list(i_am_mu,nn).le.nbft)then

        write(*,'(i4,a,i4,a,i5,i5)')&
             &list(i_am_mu,nn),'u',list(i_am_nu,nn),'u'&
             &,my_pair(i_am_mu,nn),my_pair(i_am_nu,nn)

     elseif(   list(i_am_nu,nn).gt.nbft.and.list(i_am_mu,nn).le.nbft)then


        write(*,'(i4,a,i4,a,i5,i5)')&
             &list(i_am_mu,nn),'u',list(i_am_nu,nn)-nbft,'d'&
             &,my_pair(i_am_mu,nn),my_pair(i_am_nu,nn)


     elseif(list(i_am_nu,nn).le.nbft.and.list(i_am_mu,nn).gt.nbft)then


        write(*,'(i4,a,i4,a,i5,i5)')&
             &list(i_am_mu,nn)-nbft,'d',list(i_am_nu,nn),'u'&
             &,my_pair(i_am_mu,nn),my_pair(i_am_nu,nn)


     else

        write(*,'(i4,a,i4,a,i5,i5)')&
             &list(i_am_mu,nn)-nbft,'d',list(i_am_nu,nn)-nbft,'d'&
             &,my_pair(i_am_mu,nn),my_pair(i_am_nu,nn)

     endif

     if(nn.eq.nu_doubly) write(*,*)
  enddo

  return
end subroutine ldump
