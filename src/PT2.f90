
!This version can use up to 128 orbitals and compiles with gfortran or ifort

subroutine PT2(ilength,llast,deltaE)
  use commonarrays, only: icij, n_alpha, n_beta, nsym, nbft, nbpsy,&
 nfreeze, nactive, iactive,ntotal,ifreeze,irrep,Ipt,nword,c
  use dyn_par
  implicit none
  integer ilength,ispin,i,il1,isym1,isym2,j,j1,j2
  integer jfroze,jj,k,n1,n2,n_spin,ik
  integer n_2p,irx,irprod,mo2,isym3,mok,isym4,orb2,orbk
  double precision seed,rand,eval,dnorm
  integer   molist(nbft,2),delmo(2),totalI,delorb(2)
   integer orb1,mo1,orblist(ntotal,2),dup,n,maxpt
    integer*8,allocatable:: list(:)
    integer i1,ibit,idup,left,icyc,firstwrd
    integer pttoadd(2000),sumptadd,llast
    double precision contribthresh,deltaE,e,s_overlap
    double precision Ek,ptnum,contrib
    integer orbsym2(nbft)
    logical            btest,isfrozen
         
         PRINT *,'***Running MCCIPT2***'

         if(nbft.gt.128) THEN
          PRINT *,'Too many basis functions for this version of MCCIPT2.'
          return
          end if
        
          maxpt=1000*maxc !size of array for candidate configurations in PT2
         
         contribthresh=cmin
!Use nbpsy to create orbsym2
         orbsym2=0
         j1=1
         do i=1,nsym
          
         do j=1,nbpsy(i)
         orbsym2(j1)=i
         j1=j1+1
         end do
         end do
!!!!!!!!!!!!!!!!!

          
           totalI=0
       
            
            ALLOCATE(Ipt(2,iword,maxpt))
              Ipt=0
              do i=1,ilength
              molist=0 ! zero all MOs then add back active space
     


! All possible singles first
  
          do ispin=1,2
!         Create list of MOs 

!Only include active space and then don't try to delete frozen orbitals from this later as they are deleted here
          
          do j=1,nbft

          
          do jj=1,nactive
          
          if(iactive(jj).eq.j) THEN
          molist(j,ispin)=j
          exit
          end if
          end do

          end do


          delmo(ispin)=nbft-nactive
!
          delorb(ispin)=0 !number of deleted occupied mos as they are frozen 
!

!         Create lists of occupied and unoccupied orbitals

 		orb1=1
        do j=0, nbft-1

         n1 = (j/int_bits)+1
	 j1 = j-((n1-1)*int_bits)
   
         if( btest(icij(ispin,n1,i),j1) ) THEN
         orblist(orb1,ispin)=j+1

! remove frozen here
         isfrozen=.false.
         do jfroze=1,nfreeze
         if(ifreeze(jfroze).eq.orblist(orb1,ispin)) THEN
          orb1=orb1-1
          delorb(ispin)=delorb(ispin)+1
          isfrozen=.true.
         exit
         END IF
         end do
          orb1=orb1+1
        


         if(isfrozen) cycle ! don't try to delete frozen orbitals as they are not in the active space.
        
         molist(j+1,ispin)=0 !Delete those MOs that are occupied 
         delmo(ispin)=delmo(ispin)+1
         END IF
         end do
! put deleted MOs to end
! last delmo(ispin) should be zero
         do j=0,delmo(ispin)-1
! if it is not find first zero and swap
         if (molist(nbft-j,ispin).ne.0) THEN
         
         do jj=1,(nbft-j)-1

         if (molist(jj,ispin).eq.0) THEN
         molist(jj,ispin)=molist(nbft-j,ispin)
         EXIT
         END IF
         end do
          
         end if

         end do


        

! for each spin check swap with each unoccupied mo
        
         if(ispin.eq.1) then 
         n_spin=n_alpha
         else
           n_spin=n_beta    
         end if
         do orb1=1,n_spin-delorb(ispin)
         do mo1=1,nbft-delmo(ispin)

! Check that the replacement will give the required symmetry
          isym1=orbsym2(molist(mo1,ispin))-1
          isym2=orbsym2(orblist(orb1,ispin))-1
           if((isym1).ne.isym2) cycle
           





!add new configuration to PT2 list 
    

         k=orblist(orb1,ispin)-1
     
         n1 = k/int_bits+1
	 j1 = k-(n1-1)*int_bits
         k=molist(mo1,ispin)-1
      
         n2 = k/int_bits+1
	 j2 = k-(n2-1)*int_bits
         
         
         totalI=totalI+1
         
         
   

         if (totalI.gt.maxpt) THEN
         PRINT *, 'Number of PT2 configs exceeds maxpt' 
          PRINT *, 'Configs=',totalI,'maxpt=1000*maxc'
         DEALLOCATE(Ipt)
         return
         end if
         
         Ipt(:,:,totalI)=icij(:,:,i)
           
 	 Ipt(ispin,n1,totalI)=ibclr(Ipt(ispin,n1,totalI),j1)	! BW

         Ipt(ispin,n2,totalI)=ibset(Ipt(ispin,n2,totalI),j2)	! BW

         end do
             
         end do
          end do



!doubles
! Double substitutions: one alpha, one beta
              do orb1=1,n_alpha-delorb(1)
              do orb2=1,n_beta-delorb(2)
         do mo1=1,nbft-delmo(1)
         do mo2=1,nbft-delmo(2)


! check symmetry product candidate two added orbitals equals that of  deleted
! irrep is a c style array
         isym1=orbsym2(orblist(orb1,1))-1
         isym2=orbsym2(orblist(orb2,2))-1
         isym3=orbsym2(molist(mo1,1))-1
         isym4=orbsym2(molist(mo2,2))-1

         irprod = ieor( irrep(isym1), irrep(isym2) )	
         irx = ieor( irrep(isym3), irrep(isym4) )	
         if( irprod .ne. irx) cycle
!!!!!!!!!!!!!!!!!!
                  totalI=totalI+1

                 if (totalI.gt.maxpt) THEN
         PRINT *, 'Number of PT2 configs exceeds maxpt' 
          PRINT *, 'Configs=',totalI,'maxpt=1000*maxc'
         DEALLOCATE(Ipt)
         return
         end if
         
         Ipt(:,:,totalI)=icij(:,:,i)
         
  
           
          do ispin=1,2

          orbk=orb1
          mok=mo1
          if(ispin.eq.2) THEN
          orbk=orb2
          mok=mo2
          END IF
          k=orblist(orbk,ispin)-1
         n1 = k/int_bits+1
	 j1 = k-(n1-1)*int_bits
             k=molist(mok,ispin)-1        
         n2 = k/int_bits+1
	 j2 = k-(n2-1)*int_bits   
         Ipt(ispin,n1,totalI)=ibclr(Ipt(ispin,n1,totalI),j1)	! BW
         Ipt(ispin,n2,totalI)=ibset(Ipt(ispin,n2,totalI),j2)	! BW  
         end do




        end do
        end do         
        end do
        end do


! alpha alpha and beta beta
! Double alpha then Double beta
         do ispin=1,2

         if(ispin.eq.1) then 
         n_spin=n_alpha
         else
           n_spin=n_beta    
         end if

        do orb1=1,n_spin-1-delorb(ispin)
               
              do orb2=orb1+1,n_spin-delorb(ispin)
                do mo1=1,nbft-delmo(ispin)-1
                  do mo2=mo1+1,nbft-delmo(ispin)


! check symmetry product candidate two added orbitals equals that of  deleted
         isym1=orbsym2(orblist(orb1,ispin))-1
         isym2=orbsym2(orblist(orb2,ispin))-1
         isym3=orbsym2(molist(mo1,ispin))-1
         isym4=orbsym2(molist(mo2,ispin))-1

         irprod = ieor( irrep(isym1), irrep(isym2) )	
         irx = ieor( irrep(isym3), irrep(isym4) )	
         if( irprod .ne. irx) cycle
!!!!!!!!!!!!
         
                   totalI=totalI+1
         if (totalI.gt.maxpt) THEN
         PRINT *, 'Number of PT2 configs exceeds maxpt' 
          PRINT *, 'Configs=',totalI,'maxpt=1000*maxc'
         DEALLOCATE(Ipt)
         return
         end if
         
         Ipt(:,:,totalI)=icij(:,:,i)
         
  
           
          do ik=1,2

          orbk=orb1
          mok=mo1
          if(ik.eq.2) THEN
          orbk=orb2
          mok=mo2
          END IF
          k=orblist(orbk,ispin)-1
         n1 = (k/int_bits)+1
	 j1 = k-((n1-1)*int_bits)
             k=molist(mok,ispin)-1        
         n2 = (k/int_bits)+1
	 j2 = k-((n2-1)*int_bits)
         

            
         Ipt(ispin,n1,totalI)=ibclr(Ipt(ispin,n1,totalI),j1)	! BW
         Ipt(ispin,n2,totalI)=ibset(Ipt(ispin,n2,totalI),j2)	! BW    
         end do
          

        end do
        end do         
        end do
        end do
        end do






        


      
        end do ! loop over configs

WRITE(50,*)  'PT2 configurations without removing duplicates', totalI
 call flush(50)

 call  genealogyPT(totalI,0)



!!!!!!!!!!!!!! Sorting
         PRINT *, 'Removing duplicates using quicksort method.'
        
       
         j=0
         ALLOCATE(list(totalI))
         list=0
! switch to integer*8 to allow nword=2 to be sorted in one go
         firstwrd=nword  ! put up to first 2 in the integer*8 array
         if(firstwrd.gt.2) firstwrd=2
      
         do i=1,totalI

! 1 to firstwrd first

         do n=1,firstwrd
       
         do ibit=0,31
         if(btest(Ipt(1,n,i),ibit))& 
list(i)=ibset(list(i),ibit+(32*(n-1)))
         end do
         end do         

         end do
          
           PRINT *, 'start alpha1 quick sort'
         call sort(1,totalI,list)
        

           if(nword.gt.2) THEN
           PRINT *, 'start alpha2 quick sort'
           left=1
           
         do i=2,totalI

	if (list(i).ne.list(left)) THEN
          if(left.ne.i-1) THEN !copy alpha2 to appropriate part of list
         do i1=left,i-1
         list(i1)=0
         do n=3,nword
       
         do ibit=0,31
         if(btest(Ipt(1,n,i1),ibit))& 
list(i1)=ibset(list(i1),ibit+(32*(n-3)))
         end do
         end do       
         end do
          call sort(left,i-1,list)
    	  end if 
        left=i
	END IF
	end do
        END IF
            

! all alphas should now be sorted 
! list is complicated by alpha1 and alpha 2 now so use actual Ipt
!  if Ipt for alphas is the same for all nword then sort beta1
          

              PRINT *, 'start beta1 quick sort'
          left=1
          do i=2,totalI
          idup=1
          do n=1,nword
         
          if (Ipt(1,n,i).ne.Ipt(1,n,left)) THEN
          idup=0
          EXIT
          END IF
          end do
           
          if(idup.eq.0) THEN
          if(left.ne.i-1) THEN !copy beta1 to appropriate part of list
         do i1=left,i-1
         list(i1)=0
         do n=1,firstwrd 
       
         do ibit=0,31
         if(btest(Ipt(2,n,i1),ibit))& 
list(i1)=ibset(list(i1),ibit+(32*(n-1)))
         end do
         end do       
         end do
          call sort(left,i-1,list)
    	  end if 


          left=i
          END IF
          end do

          IF (nword.gt.2) THEN
! now check if alpha1,alpha2, and beta1 are the same and if so sort using beta2
          PRINT *, 'start beta2 quick sort'
          left=1
          do i=2,totalI
          idup=1
          do n=1,nword
         
          if (Ipt(1,n,i).ne.Ipt(1,n,left)) THEN
          idup=0
          EXIT
          END IF
          end do


          if(idup.eq.1) THEN  ! if alphas are the same then check beta1

          do n=1,firstwrd 
         
          if (Ipt(2,n,i).ne.Ipt(2,n,left)) THEN
          idup=0
          EXIT
          END IF
          end do
          END IF

           
          if(idup.eq.0) THEN
          if(left.ne.i-1) THEN !copy beta1 to appropriate part of list
         do i1=left,i-1
         list(i1)=0
         do n=3,nword
       
         do ibit=0,31
         if(btest(Ipt(2,n,i1),ibit))& 
list(i1)=ibset(list(i1),ibit+(32*(n-3)))
         end do
         end do       
         end do
          call sort(left,i-1,list)
    	  end if 


          left=i
          END IF
          end do
          END IF


         DEALLOCATE(list)


!!!!!!!!!!!!!!!!!!
!Delete duplicates

        left=1
        do i=2,totalI

        do n=1,nword
        if (Ipt(1,n,i).ne.Ipt(1,n,left)) THEN
   	 left=i
         EXIT
	 END IF
         
         if (Ipt(2,n,i).ne.Ipt(2,n,left)) THEN
         left=i
         EXIT
         END IF

     
        end do

        if(left.eq.i) cycle

         
   	 Ipt(:,:,i)=0
         j=j+1
     
	

	end do        

        WRITE(50,*)  'Deleted due to replication in PT2 space',j
       
          j=0
   
!  Check for duplicates from reference space - 
! Slower but ref space is small in comparision to PT space
         do j1=1,totalI
         
! need to check that at least one integer is not zero
         icyc=1 
         do n=1,nword  
         IF (Ipt(1,n,j1).ne.0) THEN
         icyc=0
         EXIT
         END IF

         IF (Ipt(2,n,j1).ne.0) THEN
         icyc=0
         EXIT
         END IF
         
         end do

         if(icyc.eq.1) cycle
         
         do j2=1,ilength
         icyc=0
         do n=1,nword
         if ( Ipt(1,n,j1).ne.icij(1,n,j2)) then
         icyc=1
         Exit
         END IF

          if ( Ipt(2,n,j1).ne.icij(2,n,j2)) then
         icyc=1
         Exit
         END IF
         
         end do

         if(icyc.eq.1) cycle

         
      
         Ipt(:,:,j1)=0
         j=j+1
         exit


         end do
          
         end do
         
           WRITE(50,*) 'Deleted due to replication in reference space',j



!!!!!!!! Calculate DeltaE
          call energy(ilength,eval,dnorm)

            Ek=eval

            deltaE=0.0D0
           

            sumptadd=0
            do I=1,totalI
! only non-duplicates i.e. non-zero entries check alpha
            idup=1
            do n=1,nword
            if(Ipt(1,n,I).ne.0) THEN
            idup=0
            exit
            end if
            end do
! If all alpha are zero then check beta in case we have all spin beta
            if(idup.eq.1) THEN
            do n=1,nword
            if(Ipt(2,n,I).ne.0) THEN
            idup=0
            exit
            end if
            end do
            end if 

            if(idup.eq.1) cycle
!


! copy member of I to end of icij      
          icij(:,:,ilength+1)= Ipt(:,:,I)
     
         
       
! get <I|H|K> K is mcci vector
          ptnum=0.0D0

          do j=1,ilength
          
          call  muHnu(ilength+1,j,e,n_2p,s_overlap)
!c          ptnum=ptnum+(e*c(j))

!ccccccccccc
           ptnum=ptnum+((e-(s_overlap*Ek))*c(j))
!ccccccccccccccc          

     
          end do
      call  muHnu(ilength+1,ilength+1,e,n_2p,s_overlap)
! mcci vector is divided by SQRT dnorm single trial PT function by SQRT s_overlap                   
         contrib=((ptnum**2)) 
          

       
           
          if(ABS(contrib).gt.&  
ABS(dnorm*contribthresh*((Ek*s_overlap)-e))) THEN  !if removing division then swap gt to lt
           
          sumptadd=sumptadd+1
          pttoadd(sumptadd)=I 
          if(sumptadd.gt.2000) THEN
          PRINT *, 'Exiting PT2: more than 2000 new configurations found.'
          return
          END IF
          END IF

         if(sumptadd.eq.0) THEN
         contrib=contrib/(dnorm*((s_overlap*Ek)-e))  
                 

             

           deltaE=deltaE+contrib
         END IF


          end do

      do j=1,sumptadd
       
       I=pttoadd(j)

       icij(:,:,ilength+j)= Ipt(:,:,I)
       end do
       llast=ilength
       ilength=ilength+sumptadd 


  DEALLOCATE(Ipt)
  return

      end



subroutine genealogyPT(length,llast)
  use commonarrays, only: nbft, Ipt, nword, ntotal, n_alpha, n_beta, i_sx2 
  use dyn_par
  implicit real*8   (a-h,o-z)
  integer           ijdoubles(iword)
  logical           l1, l2
  !common  /aodat/   nsym, nbft, nbpsy(irmax)
  !common  /config/  Ipt(2,iword,maxc), nword
  !common  /occupy/  ntotal, n_alpha, n_beta, i_sx2

  !     random walk through branching diagram, generates linearly independent
  !     spin eigenstates after projection from primitives
  seed=1.0D0
 
  do ici = llast+1, length

     do n=1,nword
        ijdoubles(n) = iand(Ipt(1,n,ici),Ipt(2,n,ici))
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

        l1 = btest(Ipt(1,n,ici),jshift)
        l2 = btest(Ipt(2,n,ici),jshift)

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
                 Ipt(2,n,ici)=ibclr(Ipt(2,n,ici),jshift)
                 Ipt(1,n,ici)=ibset(Ipt(1,n,ici),jshift)
              endif

           else

              if(msum.eq.0) then
                 !                    the only way is up
                 msum = msum+1
                 n_up = n_up+1
                 if(l2) then
                    Ipt(2,n,ici)=ibclr(Ipt(2,n,ici),jshift)
                    Ipt(1,n,ici)=ibset(Ipt(1,n,ici),jshift)
                 endif

              elseif(2*n_up .lt. n_up_max2) then
            
                 call random(seed,rand)
           
                 if(rand.le.0.5.or.2*n_down.eq.n_dn_max2) then
                    n_up = n_up +1
                    msum = msum +1
                    if(l2) then
                       Ipt(2,n,ici)=ibclr(Ipt(2,n,ici),jshift)
                       Ipt(1,n,ici)=ibset(Ipt(1,n,ici),jshift)
                    endif
                 else
                    n_down = n_down+1
                    msum   = msum - 1
                    if(l1) then
                       Ipt(1,n,ici)=ibclr(Ipt(1,n,ici),jshift)
                       Ipt(2,n,ici)=ibset(Ipt(2,n,ici),jshift)
                    endif
                 endif

              else

                 !                    exceeded upward steps, must go down
                 msum   = msum - 1
                 n_down = n_down+1
                 if(l1) then
                    Ipt(1,n,ici)=ibclr(Ipt(1,n,ici),jshift)
                    Ipt(2,n,ici)=ibset(Ipt(2,n,ici),jshift)
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
end subroutine genealogyPT

!  swaps corresponding alpha and beta when it does list
	recursive subroutine sort(left,right,list)
          use commonarrays, only: Ipt,nword
	implicit  none
	integer*8 list(*),temp,pivotval
	integer i,left,right
	integer pivot
	integer newL,newR


	pivot=(left+right)/2
	pivotval=list(pivot)


	newL=left
	newR=right

	do while (newL<newR)
	do while(list(newL)<pivotval)
	 newL=newL+1
	end do

	do while(list(newR)>pivotval)
	 newR=newR-1
	end do

	if (newR>newL) THEN
	temp=list(newR)
	list(newR)=list(newL)
	list(newL)=temp
        ! move all alphas
        do i=1,nword
	temp=Ipt(1,i,newR)
	Ipt(1,i,newR)=Ipt(1,i,newL)
	Ipt(1,i,newL)=temp
        end do


	! move betas too 
        do i=1,nword
	temp=Ipt(2,i,newR)
	Ipt(2,i,newR)=Ipt(2,i,newL)
	Ipt(2,i,newL)=temp
        end do

	!if on pivot but not the same
	!also need to change pivot
	if (newL==pivot) THEN

	if (list(newL)==list(newR)) THEN !This sorts problem of duplicate of pivotvalue

	 newL=newL-1
	ELSE

 	 pivot=newR
	 newR=newR+1
	END IF
	END IF

	if (newR==pivot) THEN

	if (list(newL)==list(newR)) THEN

 	newR=newR+1
	ELSE

 	pivot=newL
 	newL=newL-1
	END IF
	END IF

	newL=newL+1
	  newR=newR-1
	END IF
  





	end do





	if (pivot+1<right) THEN
 	call sort(pivot+1,right,list)

	END IF
	if (pivot-1>left) THEN
 	call sort(left,pivot-1,list)
	END IF
	END subroutine sort

