subroutine chk_list(l2,l1)
  use commonarrays, only: icij,nword
  use dyn_par
  integer            l1, l2, lshift, jshift, i, j
  logical            test1
  !common  /config/   icij(2,iword,maxc), nword

  ! check the buffer from l1+1 to l2 against all
  ! previous items in the list. Returns the new list length
  ! as l2.

  lshift = 0
  do i = l1+1, l2

     do n=1,nword
        icij(1,n,i-lshift) = icij(1,n,i)
        icij(2,n,i-lshift) = icij(2,n,i)
     enddo

     jshift = lshift
     do j = 1, i-jshift-1

        test1= .TRUE.

        ! same config.
        do n=1,nword
           if(  (icij(1,n,i).ne.icij(1,n,j))&
                .or.  (icij(2,n,i).ne.icij(2,n,j)) ) test1 = .FALSE.
        enddo

        if(test1) lshift = lshift + 1

     enddo
  enddo

  l2 = l2 - lshift

  return
end subroutine chk_list
