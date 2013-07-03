subroutine timer(tw,tu,ts)
  real*4   tarray(2),etime
  real*4   tw,tu,ts

  tw = etime(tarray)
  tu = tarray(1)
  ts = tarray(2)

  return
end subroutine timer
