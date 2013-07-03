module commonarrays
  use              dyn_par
  
  integer ,allocatable     ::   icij(:,:,:)
  integer ,allocatable     ::   ijh(:), ijs(:)
  integer ,allocatable     ::   ifreeze(:)
  integer ,allocatable     ::   iactive(:)
  integer ,allocatable     ::   nbpsy(:)
  integer ,allocatable     ::   irrep(:)
  integer ,allocatable     ::   list(:,:)
  integer ,allocatable     ::   my_pair(:,:)
  integer ,allocatable     ::   icase(:)
  integer ,allocatable     ::   ipoint(:)
  real*8    ,allocatable     ::   e1ints(:), e2ints(:)
  real*8    ,allocatable     ::   cnorm(:)
  real*8    ,allocatable     ::   e(:)
  real*8    ,allocatable     ::   hf(:), sf(:)
  real*8    ,allocatable     ::   b(:,:)
  real*8    ,allocatable     ::   c(:)
  real*8    ,allocatable     ::   h(:), s(:)
  integer                  ::   nword,nfreeze,nactive
  integer                  ::   me, nproc
  integer                  ::   ntotal, n_alpha, n_beta, i_sx2
  integer                  ::   nsym, nbft

end module commonarrays
