# if n is prime
import random
from qsieve import rand_bits, gen_prime, gen_semiprime

# (x+a)^n == x^n+a mod n

def aks_test(n):
  print("testing", n)
  for _ in range(10):
    a = rand_bits(20)%n
    x = rand_bits(20)%n
    lhs = pow(x+a, n, n)
    rhs = (pow(x, n, n)+a)%n
    print(f"aks {a=} {x=} {lhs=} {rhs=}")
    if lhs != rhs:
      return False
  return lhs == rhs

if __name__ == "__main__":
  t1 = gen_prime(50)
  print(aks_test(t1))
  t2 = gen_semiprime(50)
  print(aks_test(t2))
