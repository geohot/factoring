# quadratic sieve to start
import math

semiprime_factor_me = 539873

B = 17

# generate primes up to B
PRIMES = [2, 3, 5, 7, 11, 13, 17]

def b_smooth_factorize(num):
  factors = [0]*len(PRIMES)
  for i,p in enumerate(PRIMES):
    while num%p == 0:
      factors[i] += 1
      num //= p
      if num == 1: return factors
  return None

if __name__ == "__main__":
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





