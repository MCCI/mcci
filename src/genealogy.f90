subroutine genealogy(length,llast,seed)
  use commonarrays, only: nbft, icij, nword, ntotal, n_alpha, n_beta, i_sx2
  use dyn_par
  implicit real*8   (a-h,o-z)
  integer           ijdoubles(iword)
  logical           l1, l2
  !common  /aodat/   nsym, nbft, nbpsy(irmax)
  !common  /config/  icij(2,iword,maxc), nword
  !common  /occupy/  ntotal, n_alpha, n_beta, i_sx2

  !     random walk through branching diagram, generates linearly independent
  !     spin eigenstates after projection from primitives

  do ici = llast+1, length

     do n=1,nword
        ijdoubles(n) = iand(icij(1,n,ici),icij(2,n,ici))
     enddo

     n_2p = 0
     do j=0,nbft-1
        n = j/int_bits+1
        jshift = j -(n-1)*int_bits
        if(btest(ijdoubles(n),jshift)) n_2p = n_2p+1
     enddo
     n_2p = ntotal-2*n_2p

     !        n_up - n_down = 2s    (principle case m_s = s)
     !        n_up + n_down = n_2p, n_up and n_down are for unpaired electrons, 
     !                              hence they are not n_alpha and n_beta

     n_up_max2 = n_2p+i_sx2      ! max up   steps x 2
     n_dn_max2 = n_2p-i_sx2      ! max down steps x 2

     m_sx2 =   (n_alpha - n_beta) 
     if(m_sx2.ne.i_sx2) STOP 'geneaology: not principle spin case'

     msum   = 0
     n_up   = 0
     n_down = 0
     do j=0, nbft-1

        n = j/int_bits+1
        jshift = j-(n-1)*int_bits

        l1 = btest(icij(1,n,ici),jshift)
        l2 = btest(icij(2,n,ici),jshift)

        i1 = 0
        i2 = 0
        if(l1) i1=1
        if(l2) i2=1

        if(ieor(i1,i2).ne.0) then

           if(msum.eq.0.and.n_up.eq.0) then

              !                 first spin is always alpha
              msum = 1
              n_up = 1
              if(l2) then
                 icij(2,n,ici)=ibclr(icij(2,n,ici),jshift)
                 icij(1,n,ici)=ibset(icij(1,n,ici),jshift)
              endif

           else

              if(msum.eq.0) then
                 !                    the only way is up
                 msum = msum+1
                 n_up = n_up+1
                 if(l2) then
                    icij(2,n,ici)=ibclr(icij(2,n,ici),jshift)
                    icij(1,n,ici)=ibset(icij(1,n,ici),jshift)
                 endif

              elseif(2*n_up .lt. n_up_max2) then

                 call random(seed,rand)

                 if(rand.le.0.5.or.2*n_down.eq.n_dn_max2) then
                    n_up = n_up +1
                    msum = msum +1
                    if(l2) then
                       icij(2,n,ici)=ibclr(icij(2,n,ici),jshift)
                       icij(1,n,ici)=ibset(icij(1,n,ici),jshift)
                    endif
                 else
                    n_down = n_down+1
                    msum   = msum - 1
                    if(l1) then
                       icij(1,n,ici)=ibclr(icij(1,n,ici),jshift)
                       icij(2,n,ici)=ibset(icij(2,n,ici),jshift)
                    endif
                 endif

              else

                 !                    exceeded upward steps, must go down
                 msum   = msum - 1
                 n_down = n_down+1
                 if(l1) then
                    icij(1,n,ici)=ibclr(icij(1,n,ici),jshift)
                    icij(2,n,ici)=ibset(icij(2,n,ici),jshift)
                 endif
              endif

           endif
        endif
     enddo

     if(n_2p.ne.n_up+n_down) STOP'genealogy: unpaired electrons nos. have changed'

     if(i_sx2.ne.   (n_up-n_down)) then
        write(*,*) 'is nu nd',i_sx2,n_up,n_down
        write(*,*) 'genealogy: wrong spin projection m_s'
     endif

  enddo

  return
end subroutine genealogy
