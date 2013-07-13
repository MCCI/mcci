      subroutine get_intMolpro(ecore)
   use commonarrays, only:nsym,nword,nbft,nbpsy,e1ints,e2ints,ipoint
       use dyn_par
       implicit real*8  (a-h,o-z)
       integer norb,nelec,ms2,orbsym(500),isym
       NAMELIST /FCI/ norb,nelec,ms2,orbsym,isym


         open(UNIT=14,FILE='FCIDUMP')


          READ(14,NML=FCI)
          nbft=norb
          nsym=orbsym(nbft)
         nbpsy=0 
         !!create nbpsy from orbsym
         do i=1,nbft
         nbpsy(orbsym(i))=nbpsy(orbsym(i))+1
         end do       
        
!         PRINT *,nbpsy
          

    


!   store values for pointer array for indexing symmetric matrices
      ipoint(1) = 0
      do i = 2,  max2 
         ipoint(i) = ipoint(i-1) + i - 1
      enddo

!cccccccccccccccccc Read in integrals


      nword = ((nbft-1)/int_bits) + 1

      i2=1

      do while (i2/=0)

      READ(14,*) eint, i2, j2, k2, l2
      
      if (i2/=0) THEN 
      if (k2==0) THEN  !1 body
          if(j2.gt.i2) STOP 'Error get_int: 1e ints #2 ' !molpro 1e ints reverse labels
	 ij = ipoint(i2)+j2 
         e1ints(ij) =  eint


      ELSE  !two body      
          

           if(j2.gt.i2 .or. l2.gt.k2) STOP 'Error get_int: 2e ints #1'
	 ij   = ipoint(i2)+j2
	 kl   = ipoint(k2)+l2
 	 if( kl .gt. ij) then
	     itmp = ij
	     ij =kl
	     kl = itmp
         endif
	 ijkl = ipoint(ij)+kl
	 e2ints(ijkl)=eint
                   
      END IF
      END IF
      end do

      !last eint has i2=0 so ecore
      ecore=eint
      CLOSE(14) 
      
       
      return
      end
