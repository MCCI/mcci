subroutine sym_init
  use commonarrays, only: irrep
  use dyn_par
  !common /sym/  irrep(0:irmax-1)

  !     since the irrep product table looks like
  !                \  0  1  2  3  4  5  6  7 
  !              0    0
  !              1       0
  !              2          0
  !              3             0
  !              4                0
  !              5                   0
  !              6                      0
  !              7                         0
  !
  !
  !     the only totally symmetric ground state wave
  !     functions have the symmetry product for the alpha orbitals
  !     equal to the symmetry product for the beta orbitals

  do i=0, irmax-1
     irrep(i) = i
  enddo

  return
end subroutine sym_init
