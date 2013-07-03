function ck(i_sx2,n_2p,kck)
  integer i_sx2,n_2p,kck,m1,m2,m3
  real*8  x, y, xm1, xm2, xmult, binomial, ck

  !Sanibel coefficients: principle case

  m1 = 2*i_sx2+2
  m1 = (-1)**kck*m1 
  m2 = n_2p+i_sx2+2

  xm1 = float(m1)
  xm2 = float(m2)
  xmult =  xm1 / xm2

  y = float(n_2p)
  y = y + float(i_sx2)
  m3 =nint( y/2.0d0 )

  x=nint(binomial(m3,kck))

  ck = xmult/x

end function ck

function binomial(n,k)
  integer k,n
  real*8  binomial, factorial
  binomial = dexp(factorial(n)-factorial(k)-factorial(n-k))
  return
end function binomial

function factorial(n)
  integer n
  real*8  a(100),gammaCI,factorial,xx
  save a
  data a/100*-1.0/
  if(n.lt.0) STOP'error: factorial of negative integer'
  xx =float(n+1)
  if(n.le.99) then
     if(a(n+1).lt.0.0) a(n+1) = gammaCI(xx)
     factorial=a(n+1)
  else
     factorial=gammaCI(xx)
  endif
  return
end function factorial

function gammaCI(xx)
  real*8  gammaCI,xx
  real*8  ser,stp,tmp,x,y,cof(6)
  integer j
  save    cof,stp
  data    cof,stp /76.18009172947146d0,&
       -86.50532032941677d0,&
       24.01409824083091d0,&
       -1.231739572450155d0,&
       .1208650973866179d-2,&
       -.5395239384953d-5,&
       2.5066282746310005d0/
  x=xx
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*dlog(tmp)-tmp 
  ser=1.000000000190015d0
  do j=1,6
     y=y+1.d0
     ser=ser+cof(j)/y
  enddo
  gammaCI = tmp+dlog(stp*ser/x)
  return 
end function gammaCI
