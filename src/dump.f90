subroutine dump(name,a,ija)
  implicit   real*8    (a-h,o-z)
  character*10        name
  dimension  a(*), ija(*)

  open(70,file=name)

  !     diagonals
  do ici = 1, ija(1)-2
     write(70,'(d18.10,i6)') a(ici) ,ija(ici)
  enddo

  write(70,*)'location n+1'
  write(70,*)'---------------------------'
  write(70,'(d18.10,i6)') a(ija(1)-1) ,ija(ija(1)-1)

  write(70,*)'off diagonals'
  write(70,*)'---------------------------'

  do ici = ija(1), ija(ija(1)-1)-1
     write(70,'(d18.10,i6)') a(ici) ,ija(ici)
  enddo

  close(70)

  return
end subroutine dump
