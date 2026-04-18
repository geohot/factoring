# quadratic sieve to start
import math
import random

# SLOW
#def is_prime(n): return n > 1 and all(n%d for d in range(2, math.isqrt(n)+1))
# FAST (can implement Miller-Rabin but meh)
from sympy import isprime as is_prime

def gen_prime(bits):
  ret = random.randint((1 << (bits-1))+1, 1 << bits)
  while not is_prime(ret): ret += 1
  return ret

# generate
SEMIPRIME_BITS = 18

# params for the solver
B = 50

# generate primes up to B
PRIMES = [n for n in range(B+1) if is_prime(n)]
NUM_RELATIONS = int(len(PRIMES)*1.2)+2
print(f"num primes {len(PRIMES)}, num relations {NUM_RELATIONS}")

def b_smooth_factorize(num):
  factors = [0]*len(PRIMES)
  for i,p in enumerate(PRIMES):
    while num%p == 0:
      factors[i] += 1
      num //= p
      if num == 1: return factors
  return None

def process_congruence(semiprime_factor_me, nums, factors):
  # then we solve with GCD
  rhs = math.prod((p**(f//2)) for p,f in zip(PRIMES, factors))
  neg_factor = math.gcd(math.prod(nums) - rhs, semiprime_factor_me)
  pos_factor = math.gcd(math.prod(nums) + rhs, semiprime_factor_me)

  # try neg and pos factors
  for factor in [neg_factor, pos_factor]:
    if factor == 1 or factor == semiprime_factor_me: continue
    other_factor = semiprime_factor_me//factor
    print("FOUND", nums, factors, "factors into", factor, other_factor)
    assert semiprime_factor_me == factor*other_factor
    return True
  return False

def qsieve(semiprime_factor_me):
  # first we need to find B-smooth numbers that are perfect squares
  # TODO: real qsieve doesn't check all of these, it finds likely candidates
  start = math.isqrt(semiprime_factor_me)+1
  relations = []
  while len(relations) < NUM_RELATIONS:
    # TODO: chatgpt says it should be - and not % here, but that's slower
    qx = start*start % semiprime_factor_me
    assert qx > 0
    relation = b_smooth_factorize(qx)
    if relation:
      # start^2 === relation
      print(start, relation)
      relations.append((start, relation))
    start += 1

  # then we need to solve to make a perfect square from the relations
  print(f"collected {len(relations)=}")
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
      if process_congruence(semiprime_factor_me, nums, factors): break
  else:
    raise RuntimeError("failed to find congruence")

if __name__ == "__main__":
  # generate number for factoring
  while 1:
    p,q = gen_prime(SEMIPRIME_BITS), gen_prime(SEMIPRIME_BITS)
    if p==q: continue
    semiprime_factor_me = p*q
    print(f"factoring {semiprime_factor_me} into {p} {q}")
    break

  qsieve(semiprime_factor_me)

