subroutine h_s_reduced(length,kl,ku)
  use commonarrays, only: h, ijh, s, ijs, hf, sf, b
  use dyn_par
  implicit   real*8  (a-h,o-z)
  !common  /hands/     h(maxh), ijh(maxh), s(maxs), ijs(maxs)
  !common  /reduced/   hf(kmax*(kmax+1)/2), sf(kmax*(kmax+1)/2) 
  !common  /bk/        b(maxc,kmax)
  real*8              ai(maxc), di(maxc)  !, h, s, hf, sf, b
  integer             i, j, k, kl, ku, length
  !integer             ijh, ijs !declaration has moved to commonarrays module

  !     matrices are UPPER packed
  if(kl.le.0) kl=l 
  do j = kl, ku
     do i=1, j

        call mxv_sparse(length,h,ijh,b(1,j),ai)
        call mxv_sparse(length,s,ijs,b(1,j),di)

        hf(i+(j-1)*j/2) = 0.0d0
        sf(i+(j-1)*j/2) = 0.0d0
        do k = 1, length   
           hf(i+(j-1)*j/2) = hf(i+(j-1)*j/2) + b(k,i)*ai(k)
           sf(i+(j-1)*j/2) = sf(i+(j-1)*j/2) + b(k,i)*di(k)
        enddo

     enddo
  enddo

  return
end subroutine h_s_reduced
