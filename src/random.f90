subroutine random(seed,rand)
  real*8     seed, rand, m, a, c, k
  parameter (m = 199017.0, a = 24298.0, c = 1.0)
  seed=a*seed+c
  k    = dint(seed/m)
  seed = seed-m*k
  rand = seed/(m+1.0)
  return
end subroutine random
