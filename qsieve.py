# quadratic sieve to start
import math
import random
import tqdm
import time

# SLOW
#def is_prime(n): return n > 1 and all(n%d for d in range(2, math.isqrt(n)+1))
# FAST (can implement Miller-Rabin but meh)
from sympy import isprime as is_prime

def gen_prime(bits):
  ret = random.randint((1 << (bits-1))+1, 1 << bits)
  while not is_prime(ret): ret += 1
  return ret

# generate
SEMIPRIME_BITS = 40
# generate number for factoring (N)
while 1:
  p,q = gen_prime(SEMIPRIME_BITS), gen_prime(SEMIPRIME_BITS)
  if p==q: continue
  N = p*q
  print(f"factoring {N} into {p} {q} with {N.bit_length()} bits")
  break
# no cheating
del p,q

# params for the solver
B = 5000

# params for log sieve
BLOCK_SIZE = 4096
LOG_SIEVE_THRESHOLD = math.log(B)

# generate primes up to B filtered by quadratic residue
# https://en.wikipedia.org/wiki/Euler%27s_criterion
FACTOR_BASE = [2] + [p for p in range(3, B+1, 2) if is_prime(p) and pow(N, (p-1)//2, p) == 1]
NUM_RELATIONS = len(FACTOR_BASE) + max(8, len(FACTOR_BASE) // 10)
print(f"{B=} {max(FACTOR_BASE)=} {len(FACTOR_BASE)=} {NUM_RELATIONS=} {LOG_SIEVE_THRESHOLD=:.2f}")

def format_factors(factors):
  return ' * '.join([f"{p}^{f}" if f > 1 else f"{p}" for p,f in zip(FACTOR_BASE, factors) if f > 0])

def process_congruence(N, nums, factors):
  # then we solve with GCD
  rhs = math.prod((p**(f//2)) for p,f in zip(FACTOR_BASE, factors))
  neg_factor = math.gcd(math.prod(nums) - rhs, N)
  pos_factor = math.gcd(math.prod(nums) + rhs, N)

  # try neg and pos factors
  for factor in [neg_factor, pos_factor]:
    if factor == 1 or factor == N: continue
    other_factor = N//factor
    print(f"found square, {len(nums)} used relations with {sum(x>0 for x in factors)} factors")
    print("factors into", factor, other_factor)
    assert N == factor*other_factor
    return True
  return False

def b_smooth_factorize(num):
  factors = [0]*len(FACTOR_BASE)
  for i,p in enumerate(FACTOR_BASE):
    while num%p == 0:
      factors[i] += 1
      num //= p
      if num == 1: return factors
  return None

def qsieve(N):
  a = math.isqrt(N)
  delta = a*a - N
  # this is the magic polynomial of Quadratic Sieve
  def Q(x):
    # TODO: MPQS and SIQS use multiple polynomials here, not just one
    #return (a+x)*(a+x) - N
    # FOIL to x*x + 2*a*x + a*a - N
    return x*x + 2*a*x + delta

  # for log(Q(x)) approx
  a_f, delta_f = float(a), float(delta)
  def Qf(x): return x*x + 2*a_f*x + delta_f

  st = time.perf_counter()
  # compute the ROOTS for each prime in FACTOR_BASE
  # Q(x) == 0 mod p
  ROOTS = {p:list(set([x for x in range(p) if Q(x) % p == 0])) for p in FACTOR_BASE}

  # precompute the logs of the primes and put roots in a list
  ROOTS_LIST = []
  for p in FACTOR_BASE:
    log_p = math.log(p)
    for root in ROOTS[p]:
      ROOTS_LIST.append((root,p,log_p))
  print(f"{len(ROOTS_LIST)=} in {time.perf_counter()-st:.2f} s")

  # first we need to find B-smooth Q(x) values
  st = time.perf_counter()
  relations = []
  progress = tqdm.tqdm(total=NUM_RELATIONS)
  log_sieve_false_positive = 0
  searched = 0

  # after the max here it gets dumb
  # TODO: we can explore both negative and positive x
  for x_block in range(1, math.isqrt(2*N)-a, BLOCK_SIZE):
    # we approximate log(Q(x)) in float
    scores = [math.log(Qf(x_block + j)) for j in range(BLOCK_SIZE)]

    # sieve with the dividing roots
    for root,p,log_p in ROOTS_LIST:
      j = (root - x_block) % p
      while j < BLOCK_SIZE:
        scores[j] -= log_p
        j += p

    # check for success
    for j in range(BLOCK_SIZE):
      # if we fully divided it, it's a good relation
      if scores[j] < LOG_SIEVE_THRESHOLD:
        x = x_block+j
        if relation:=b_smooth_factorize(Q(x)):
          progress.set_description(f"{a+x}")
          relations.append((a+x, relation))
        else:
          # NOTE: this can miss if the threshold is high
          # the log doesn't include multiple
          log_sieve_false_positive += 1

    # are we done after this block?
    searched += BLOCK_SIZE
    progress.update(min(NUM_RELATIONS, len(relations))-progress.n)
    if len(relations) >= NUM_RELATIONS: break
  progress.close()

  # then we need to solve to make a perfect square from the relations
  print(f"collected {len(relations)=} in {time.perf_counter()-st:.2f} s")
  print(f"searched {searched} values of Q")
  print(f"got {log_sieve_false_positive} false positives from log sieve")

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

