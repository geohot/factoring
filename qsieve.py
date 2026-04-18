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
# generate number for factoring (N)
while 1:
  p,q = gen_prime(SEMIPRIME_BITS), gen_prime(SEMIPRIME_BITS)
  if p==q: continue
  N = p*q
  print(f"factoring {N} into {p} {q}")
  break
# no cheating
del p,q

# params for the solver
B = 200

# generate primes up to B filtered by quadratic residue
# https://en.wikipedia.org/wiki/Euler%27s_criterion
FACTOR_BASE = [2] + [p for p in range(3, B+1, 2) if is_prime(p) and pow(N, (p-1)//2, p) == 1]
NUM_RELATIONS = int(len(FACTOR_BASE)*1.2)+2
print(f"{len(FACTOR_BASE)=} {NUM_RELATIONS=}")

def b_smooth_factorize(num):
  factors = [0]*len(FACTOR_BASE)
  for i,p in enumerate(FACTOR_BASE):
    while num%p == 0:
      factors[i] += 1
      num //= p
      if num == 1: return factors
  return None

def process_congruence(N, nums, factors):
  # then we solve with GCD
  rhs = math.prod((p**(f//2)) for p,f in zip(FACTOR_BASE, factors))
  neg_factor = math.gcd(math.prod(nums) - rhs, N)
  pos_factor = math.gcd(math.prod(nums) + rhs, N)

  # try neg and pos factors
  for factor in [neg_factor, pos_factor]:
    if factor == 1 or factor == N: continue
    other_factor = N//factor
    print("FOUND")
    print(nums)
    print([(p,f) for p,f in zip(FACTOR_BASE, factors) if f > 0])
    print("factors into", factor, other_factor)
    assert N == factor*other_factor
    return True
  return False

def qsieve(N):
  # first we need to find B-smooth Q(x) values
  # TODO: real qsieve doesn't check all of these, it finds likely candidates
  start = math.isqrt(N)+1
  end = math.isqrt(2*N)  # after this it gets dumb
  relations = []
  while len(relations) < NUM_RELATIONS and start < end:
    # TODO: chatgpt says it should be - and not % here, but that's slower
    # added end bound for the search to fix
    qx = start*start - N
    assert qx > 0
    relation = b_smooth_factorize(qx)
    if relation:
      # start^2 === relation
      print(start, relation)
      relations.append((start, relation))
    start += 1

  # then we need to solve to make a perfect square from the relations
  print(f"collected {len(relations)=}")

  # we need to find a basis among the parity masks
  basis = {}
  for i,(_,relation) in enumerate(relations):
    mask = sum(1 << i for i,n in enumerate(relation) if n&1)
    combo = 1 << i
    while mask:
      # what's the largest number in mask
      # did you know bit_length was python built-in on int?
      pivot = mask.bit_length() - 1

      # if we don't have this number, we add it
      if pivot not in basis:
        basis[pivot] = (mask, combo)
        break

      # if we do have this number, we remove it
      mask ^= basis[pivot][0]
      combo ^= basis[pivot][1]

    # if we get a hit, reconstruct nums and factors from combo
    if mask == 0:
      nums = []
      factors = [0]*len(FACTOR_BASE)
      for j,(num,relation) in enumerate(relations):
        if combo&(1<<j):
          nums.append(num)
          for k,n in enumerate(relation):
            factors[k] += n
      if process_congruence(N, nums, factors): break
  else:
    raise RuntimeError("failed to find congruence")

if __name__ == "__main__":
  qsieve(N)

