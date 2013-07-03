subroutine branch(nu_configs,length,llast,seed)
  use commonarrays, only: icij, nword, c, n_alpha, n_beta, nfreeze
  use dyn_par
  !implicit none
  implicit real*8   (a-h,o-z)
  !common  /config/  icij(2,iword,maxc), nword
  !common  /coef/    c(maxc)
  !common  /occupy/  ntotal,n_alpha,n_beta,i_sx2
  !common  /froze/   nfreeze, ifreeze(maxocc)
  !common  /active/  nactive, iactive(maxocc)

  parameter (one3= 1.0/3.0   )
  parameter (two3= 2.0/3.0   )
  ilength =length

  if (lref.eq.0) then
     lref_in_br=length
  else
     lref_in_br=lref  
  endif

  do while (ilength-length .lt. nu_configs)

     call random(seed,rand) 
     i=int(dble(lref_in_br)*rand + 1)

     ilength = ilength + 1
     if(ilength.gt.maxc) STOP 'Error branch: exceeded maximum configs'

     ! setup new configurations based on old configs
111  continue
     do n=1,nword
        icij(1,n,ilength) = icij(1,n,i)
        icij(2,n,ilength) = icij(2,n,i)
     enddo

     ! make a double or single substitution?
     isubst = 1
     call random(seed,rand)
     if( rand .ge. 0.5 ) isubst = 2

     !        single 
     if (isubst .eq. 1) then
        ispin = 1 ! alpha
        call random(seed,rand)
        if ( rand .ge. 0.5 ) ispin = 2 ! beta
        if(ispin.eq.1 .and. n_alpha-nfreeze.lt. 1) then
           if(n_beta-nfreeze .ge.1) then
              ispin = 2       
           else
              STOP' branch: no orbitals to substitute'
           endif
        endif
        if(ispin.eq.2 .and.n_beta -nfreeze.lt. 1) then
           if(n_alpha-nfreeze.ge.1) then
              ispin = 1       
           else
              STOP' branch: no orbitals to substitute'
           endif
        endif
        call singles(ilength,ispin,seed)

     else
        !        double
        ispin = 1 ! alpha, alpha
        call random(seed,rand)
        if(one3 .le.rand .and.two3 .gt.rand) ispin = 2 ! beta,  beta
        if ( two3 .le.rand) ispin = 3 ! alpha, beta
        if(ispin.eq.1 .and.n_alpha-nfreeze.lt. 2) then
           if(n_beta-nfreeze .ge.2) then
              ispin = 2       
           elseif(n_alpha-nfreeze.ge.1 .and.n_beta -nfreeze.ge.1) then
              ispin = 3
           else
              return ! better luck with singles?
           endif
        endif
        if(ispin.eq.2 .and.n_beta -nfreeze.lt. 2) then
           if(n_alpha-nfreeze.ge.2) then
              ispin = 1       
           elseif(n_alpha-nfreeze.ge.1 .and.n_beta -nfreeze.ge.1) then
              ispin = 3
           else
              return ! better luck with singles?
           endif
        endif
        if(ispin.eq.3 .and.(n_alpha-nfreeze.lt. 1 .or.n_beta -nfreeze.lt. 1)) then
           if(n_alpha-nfreeze.ge.2) then
              ispin = 1       
           elseif(n_beta-nfreeze.ge.2) then
              ispin = 2
           else  
              return ! better luck with singles?
           endif
        endif
        call doubles(i,ilength,ispin,seed)
     endif
     if(ilength-length .ge. nu_configs) goto 222
  enddo

222 llast  = length
  length = ilength

  do ici=llast+1, length
     c(ici) = 0.0
  enddo

  return
end subroutine branch
