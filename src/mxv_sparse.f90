subroutine mxv_sparse(length,m,ijm,v,mv)
integer  ijm(*), length, ici, jci, k
real*8   m(*),v(*),mv(*)

if( ijm(1) .ne. length + 2) STOP 'gen_mxv: length mismatch'  

do ici=1, length
   !diagonal
   mv(ici) = m(ici)*v(ici)
enddo

do ici=1, length
   do k=ijm(ici), ijm(ici+1)-1
      !lower off diagonals
      mv(ici) = mv(ici) + m(k)*v(ijm(k))
      !upper off diagonals
      jci =ijm(k)
      mv(jci) = mv(jci) + m(k)*v(ici)
   enddo
enddo

return
end
