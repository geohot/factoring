# quadratic sieve to start
import math
import random
def is_prime(n): return n > 1 and all(n%d for d in range(2, math.isqrt(n)+1))
def gen_prime(bits):
  ret = random.randint((1 << (bits-1))+1, 1 << bits)
  while not is_prime(ret): ret += 1
  return ret

# params for the solver
B = 29
NUM_RELATIONS = 15

# generate primes up to B
PRIMES = [n for n in range(B+1) if is_prime(n)]
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
  start = int(math.sqrt(semiprime_factor_me)+1)
  relations = []
  while len(relations) < NUM_RELATIONS:
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
      if process_congruence(semiprime_factor_me, nums, factors): break
  else:
    raise RuntimeError("failed to find congruence")

if __name__ == "__main__":
  # generate number for factoring
  while 1:
    p,q = gen_prime(18), gen_prime(18)
    if p==q: continue
    semiprime_factor_me = p*q
    print(f"factoring {semiprime_factor_me} into {p} {q}")
    break

  qsieve(semiprime_factor_me)

