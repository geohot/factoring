# quadratic sieve to start
import math, random, tqdm, time, itertools, collections

# SLOW
#def isprime(n): return n > 1 and all(n%d for d in range(2, math.isqrt(n)+1))
# FAST (can implement Miller-Rabin but meh)
from sympy import isprime

# find x such that (x*x) % p == n
# SLOW
#def sqrt_mod(n, p):
#  n %= p
#  for x in range(p):
#    if (x * x) % p == n: return x
#  raise ValueError(f"no sqrt_mod for {n=} and {p=}")
# FAST (can implement with Tonelli-Shanks but meh)
from sympy import sqrt_mod

# param for the generator
BITS = 100

# use this heuristic for B
B = int(math.exp(0.5 * math.sqrt((BITS * math.log(2)) * math.log(BITS * math.log(2)))))

# params for log sieve
BLOCK_SIZE = 4096
BLOCKS_PER_SUPERBLOCK = 64
SUPERBLOCK_SIZE = BLOCK_SIZE*BLOCKS_PER_SUPERBLOCK
LOG_SIEVE_THRESHOLD = 2*math.log(B) + 1

def gen_prime(bits):
  ret = random.randint((1 << (bits-1))+1, 1 << bits)
  while not isprime(ret): ret += 1
  return ret

# generate number for factoring (N)
while 1:
  p,q = gen_prime(BITS//2), gen_prime(BITS//2)
  if p==q: continue
  N = p*q
  print(f"factoring {N} into {p} {q} with {N.bit_length()} bits")
  break
# no cheating
del p,q

# generate primes up to B filtered by quadratic residue
# https://en.wikipedia.org/wiki/Euler%27s_criterion
FACTOR_BASE = [2] + [p for p in range(3, B+1, 2) if isprime(p) and pow(N, (p-1)//2, p) == 1]
NUM_RELATIONS = len(FACTOR_BASE) + max(8, len(FACTOR_BASE) // 10)
FACTOR_BASE_SMALL = [x for x in FACTOR_BASE if x < BLOCK_SIZE]
print(f"{B=} {max(FACTOR_BASE)=} {len(FACTOR_BASE)=} {len(FACTOR_BASE_SMALL)=}")
print(f"{NUM_RELATIONS=} {LOG_SIEVE_THRESHOLD=:.2f} {BLOCK_SIZE=} {SUPERBLOCK_SIZE=}")

def format_factors(factors):
  return ' * '.join([f"{p}^{f}" if f > 1 else f"{p}" for p,f in zip(FACTOR_BASE, factors) if f > 0])

def process_congruence(N, relations, combo):
  nums = []
  factors = [0]*len(FACTOR_BASE)
  all_extra = 1
  for j, (ax, relation, neg, extra) in enumerate(relations):
    if combo&(1<<j):
      nums.append(ax)
      for k,n in enumerate(relation):
        factors[k] += n
      all_extra *= extra

  # then we solve with GCD
  lhs = math.prod(nums)
  rhs = all_extra * math.prod((p**(f//2)) for p,f in zip(FACTOR_BASE, factors))
  neg_factor = math.gcd(lhs - rhs, N)
  pos_factor = math.gcd(lhs + rhs, N)

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
      if num == 1: return factors, num
  return factors, num

# this is the magic polynomial of Quadratic Sieve
def Q(x, a, delta):
  # TODO: MPQS and SIQS use multiple polynomials here, not just one
  #return (x+a)*(x+a) - N
  # FOIL to x*x + 2*a*x + a*a - N
  return x*x + 2*a*x + delta

def do_superblock_sieve(x_superblock, ROOTS_LIST, a_f, delta_f):
  ret = []
  for x_block in range(x_superblock, x_superblock+SUPERBLOCK_SIZE, BLOCK_SIZE):
    # we approximate log(Q(x)) in float
    scores = [math.log(abs(Q(x_block + j, a_f, delta_f))) for j in range(BLOCK_SIZE)]

    # sieve with the dividing roots
    # NOTE: we tried bucket sieve, it wasn't much faster and it was more complex
    for root,p,log_p in ROOTS_LIST:
      j = (root - x_block) % p
      while j < BLOCK_SIZE:
        scores[j] -= log_p
        j += p

    # check for success
    for j in range(BLOCK_SIZE):
      if scores[j] < LOG_SIEVE_THRESHOLD: ret.append(x_block+j)
  return ret

def qsieve(N):
  a = math.isqrt(N)
  delta = a*a - N

  # for log(Q(x)) approx
  a_f, delta_f = float(a), float(delta)

  st = time.perf_counter()
  # compute the ROOTS for each prime in FACTOR_BASE
  # Q(x) % p == 0 (solve for x)
  assert FACTOR_BASE[0] == 2, "two isn't used?"
  ROOTS = {2:[(1-a)%2]} # 0 == (x*x + 2*a*x + a*a - N) % 2 == x + a - 1 -> 1-a == x
  for p in FACTOR_BASE[1:]:
    # odd primes have two roots
    # 0 == ((x+a)^2 - N) % p
    # N == (x+a)^2 % p
    s = sqrt_mod(N, p)  # s*s == N % p
    # so x+a == +/- s % p
    r1 = ( s - a) % p
    r2 = (-s - a) % p
    assert r1 != r2
    ROOTS[p] = [r1, r2]
  print(f"computed {len(ROOTS)=} in {time.perf_counter()-st:.2f} s")

  # precompute the logs of the primes and put roots in a list
  ROOTS_LIST = []
  for p in FACTOR_BASE:
    log_p = math.log(p)
    for root in ROOTS[p]:
      ROOTS_LIST.append((root,p,log_p))

  # first we need to find B-smooth Q(x) values
  # we use a log sieve to prefilter everything
  st = time.perf_counter()
  progress = tqdm.tqdm(total=NUM_RELATIONS)
  searched = 0
  relations = []
  partials = {}
  log_sieve_false_positive = 0
  sieve_time_s = 0.0
  # NOTE: we explore both negative and positive x
  BOUND = math.isqrt(2*N)-a  # after the max here it gets dumb
  superblock_schedule_order = itertools.chain.from_iterable(zip(
    range(1, BOUND, SUPERBLOCK_SIZE),
    range(-SUPERBLOCK_SIZE, -BOUND, -SUPERBLOCK_SIZE)))

  for x_superblock in superblock_schedule_order:
    # here we extract the real relations from the log_sieve
    st = time.perf_counter()
    sieved = do_superblock_sieve(x_superblock, ROOTS_LIST, a_f, delta_f)
    sieve_time_s = time.perf_counter() - st
    for x in sieved:
      relation, num = b_smooth_factorize(abs(Q(x, a, delta)))
      if num == 1: relations.append((a+x, relation, x<0, 1))
      elif isprime(num):
        if num in partials:
          # got a match!
          lhs0, relation0, neg0 = partials[num]
          combined_relation = [u + v for u, v in zip(relation0, relation)]
          relations.append(((a+x)*lhs0, combined_relation, neg0 ^ (x<0), num))
          del partials[num]
        else:
          partials[num] = (a+x, relation, x<0)
      else: log_sieve_false_positive += 1
    searched += 1
    progress.set_description(f"{searched:5d} superblocks, {len(relations)/searched:.2f} relations/superblock")
    progress.update(min(NUM_RELATIONS, len(relations))-progress.n)
    if len(relations) >= NUM_RELATIONS: break
  progress.close()
  print(f"collected {len(relations)=} across {searched} blocks in {time.perf_counter()-st:.2f} s")
  print(f"got {log_sieve_false_positive=} with {sieve_time_s:.2f} s in the sieve")

  # then we need to solve to make a perfect square from the relations
  # we need to find a basis among the parity masks
  print(f"matrix size is {len(relations)} x {len(FACTOR_BASE)+1}")
  st = time.perf_counter()
  congruence_false_positive = 0
  basis = {}
  for i, (_, relation, neg, _) in enumerate(relations):
    mask = sum(1 << i for i,n in enumerate(relation) if n&1)
    # NOTE: it works fine without this line? because process_congruence is retried?
    if neg: mask |= 1 << len(FACTOR_BASE)
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
      if process_congruence(N, relations, combo): break
      else: congruence_false_positive += 1
  else:
    raise RuntimeError("failed to find congruence")
  print(f"solved with {congruence_false_positive=} in {time.perf_counter()-st:.2f} s")
  assert congruence_false_positive < 10

if __name__ == "__main__":
  qsieve(N)

