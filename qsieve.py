# quadratic sieve to start
import math
import random
def is_prime(n): return n > 1 and all(n%d for d in range(2, math.isqrt(n)+1))
def gen_prime(bits):
  ret = random.randint((1 << (bits-1))+1, 1 << bits)
  while not is_prime(ret): ret += 1
  return ret

# generate primes up to B
B = 17
PRIMES = [n for n in range(B+1) if is_prime(n)]

def b_smooth_factorize(num):
  factors = [0]*len(PRIMES)
  for i,p in enumerate(PRIMES):
    while num%p == 0:
      factors[i] += 1
      num //= p
      if num == 1: return factors
  return None

if __name__ == "__main__":
  # generate number for factoring
  while 1:
    p,q = gen_prime(20), gen_prime(20)
    if p==q: continue
    semiprime_factor_me = p*q
    print(f"factoring {semiprime_factor_me} into {p} {q}")
    break

  # first we need to find B-smooth numbers that are perfect squares
  start = int(math.sqrt(semiprime_factor_me)+1)
  relations = []
  while len(relations) < 5:
    relation = b_smooth_factorize((start*start)%semiprime_factor_me)
    if relation:
      # start^2 === relation
      print(start, relation)
      relations.append((start, relation))
    start += 1

  # then we need to solve to make a perfect square from the relations
  # TODO: real algorithm here
  for i in range(1, 1<<len(relations)):
    factors = [0]*len(PRIMES)
    for j,(_,relation) in enumerate(relations):
      if i&(1<<j):
        for k,n in enumerate(relation):
          factors[k] += n
    if all(x%2 == 0 for x in factors):
      nums = []
      for j,(num,_) in enumerate(relations):
        if i&(1<<j):
          nums.append(num)
      break
  else:
    raise RuntimeError("failed to find congruence")

  # then we solve with GCD
  print("FOUND", nums, factors)
  rhs = math.prod((p**(f//2)) for p,f in zip(PRIMES, factors))
  factor = math.gcd(math.prod(nums) - rhs, semiprime_factor_me)
  other_factor = semiprime_factor_me//factor
  print(factor, other_factor)
  assert semiprime_factor_me == factor*other_factor





