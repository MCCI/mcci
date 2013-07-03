subroutine wxmxv_sparse(length,m,ijm,w,v,wmv)
  integer  ijm(*), length, ici, jci, k
  real*8   m(*), w(*), v(*), mv(length), wmv

  if( ijm(1) .ne. length + 2) STOP 'wxmxv: length mismatch'  

  do ici=1, length
     !        diagonal
     mv(ici) = m(ici)*v(ici)
  enddo

  do ici=1, length
     do k=ijm(ici), ijm(ici+1)-1
        !           lower off diagonals
        mv(ici) = mv(ici) + m(k)*v(ijm(k))
        !           upper off diagonals
        jci =ijm(k)
        mv(jci) = mv(jci) + m(k)*v(ici)
     enddo
  enddo

  wmv = 0.0d0
  do ici = 1, length
     wmv = wmv + w(ici)*mv(ici)
  enddo

  return
end subroutine wxmxv_sparse
