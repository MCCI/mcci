subroutine sch0(e,i_am_mu,i_am_nu,nu_doubly,kck,n_2p,ep)
  use commonarrays, only: nbft, ntotal, i_sx2, e1ints, e2ints, ipoint, list
  use dyn_par
  implicit real*8    (a-h,o-z)
  !common  /aodat/   nsym, nbft, nbpsy(irmax)
  !common  /occupy/  ntotal, n_alpha, n_beta, i_sx2
  !common  /ints/    e1ints(max1), e2ints(max2)
  !common  /pointer/ ipoint(max2)
  !common  /mxcoinc/ list(2,maxocc)

  !     diagonal Slater Condon Harris rules

  !     zero energy
  e = 0.0

  !     sum one particle contributions for doubly and singly occupied orbitals
  do n = 1, ntotal
     korb = list(i_am_nu,n) 
     if( korb .gt. nbft) korb = korb - nbft
     kk = ipoint(korb) + korb
     e = e + e1ints(kk)
  enddo

  !     sum two particle contributions for doubly occupied orbitals
  do n = 2, nu_doubly, 2
     do m = 2, nu_doubly, 2

        !           get the Coulomb terms
        korb = list(i_am_nu,n)
        if(korb.gt.nbft) korb = korb - nbft

        lorb = list(i_am_nu,m)
        if(lorb.gt.nbft) lorb = lorb - nbft

        kk = ipoint(korb) + korb
        ll = ipoint(lorb) + lorb

        if(kk.ge.ll) then
           kkll = ipoint(kk) + ll
        else	
           kkll = ipoint(ll) + kk
        endif

        !           (ij|v|ij)
        e = e + 2.0*e2ints(kkll)

        !           get the exchange terms
        if(korb.ge.lorb) then
           kl = ipoint(korb) + lorb    
        else	
           kl = ipoint(lorb) + korb
        endif

        klkl = ipoint(kl) + kl       

        !           (ii|v|jj) = (ij|v|ji)
        e = e - e2ints(klkl)

     enddo
  enddo

  !     sum two particle contributions between doubly and singly 
  !     occupied orbitals
  do n = 2, nu_doubly, 2
     do m = nu_doubly+1, ntotal

        !           get the Coulomb terms
        korb = list(i_am_nu,n) 
        if(korb.gt.nbft) korb = korb - nbft

        lorb = list(i_am_nu,m)
        if(lorb.gt.nbft) lorb = lorb - nbft

        if(korb .ge. lorb) then
           kk = ipoint(korb) + korb
           ll = ipoint(lorb) + lorb
        else	
           kk = ipoint(lorb) + lorb
           ll = ipoint(korb) + korb
        endif

        kkll = ipoint(kk) + ll

        !           (ij|v|ij)
        e = e + 2.0*e2ints(kkll)

        !           get the exchange terms
        if(korb .ge. lorb) then 
           kl   = ipoint(korb) + lorb   
        else	
           kl   = ipoint(lorb) + korb
        endif

        klkl = ipoint(kl) + kl      

        !           (ii|v|jj) = (ij|v|ji)
        e = e - e2ints(klkl)

     enddo
  enddo

  sc1 = ck(i_sx2,n_2p,kck)

  e = sc1*e

  !     sum two particle contributions between singly occupied orbitals
  do n = nu_doubly+2, ntotal
     do m = nu_doubly+1, n-1

        ispin = 1
        if(list(i_am_mu,n).gt.nbft) ispin = 0

        jspin = 1
        if(list(i_am_mu,m).gt.nbft) jspin = 0

        !           get the Coulomb terms
        korb = list(i_am_nu,n)
        kspin = 1
        if(korb.gt.nbft) then
           korb = korb - nbft
           kspin = 0
        endif

        lorb = list(i_am_nu,m)
        lspin = 1
        if(lorb.gt.nbft) then
           lorb = lorb - nbft
           lspin = 0 
        endif

        kk = ipoint(korb) + korb
        ll = ipoint(lorb) + lorb

        if(kk.ge.ll) then
           kkll = ipoint(kk) + ll
        else	
           kkll = ipoint(ll) + kk
        endif

        !           (ij|v|ij)
        e = e + sc1*e2ints(kkll)

        if(korb.ge.lorb) then 
           kl = ipoint(korb) + lorb   
        else 
           kl = ipoint(lorb) + korb
        endif

        klkl = ipoint(kl) + kl        

        !           must calculate a new Sanibel coefficient for the exchange terms
        if(ispin.eq.jspin .or. kspin.eq.lspin) then
           sc2 = sc1                     	! k(ml|ml) = k
        elseif(ispin.eq.kspin .and. jspin.eq.lspin) then
           sc2 = ck(i_sx2,n_2p,kck+1)  	! k(ml|ml) = k+1
        elseif(ispin.eq.lspin .and. jspin.eq.kspin) then
           sc2 = ck(i_sx2,n_2p,kck-1)  	! k(ml|ml) = k-1
        endif

        !           get the exchange terms
        !           (ii|v|jj) = (ij|v|ji)
        e = e - sc2*e2ints(klkl)

     enddo
  enddo

  e = e*ep

  return
end subroutine sch0
