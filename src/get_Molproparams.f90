      subroutine get_Molproparams(nbft,nsym)
       integer nbft,nsym
       integer norb,nelec,ms2,orbsym(500),isym ! 500 MO limit here
       NAMELIST /FCI/ norb,nelec,ms2,orbsym,isym

        
         open(UNIT=14,FILE='FCIDUMP')


          READ(14,NML=FCI)
          nbft=norb
          nsym=orbsym(nbft)
          
        
          
      CLOSE(14) 
      
       
      return
      end
