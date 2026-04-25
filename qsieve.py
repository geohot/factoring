# quadratic sieve to start
import math, random, tqdm, time, itertools, collections, functools
from dataclasses import dataclass, field

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

# tinygrad for speed
from tinygrad import Tensor, getenv

# param for the generator
BITS = 150

# use this heuristic for B
B = int(math.exp(0.5 * math.sqrt((BITS * math.log(2)) * math.log(BITS * math.log(2)))))

# hack to deal with inefficencies in this implementation (makes it way faster)
#B *= 2
GPU = getenv("GPU", 1)

# large prime bound
LP_BOUND = 32*B

# params for log sieve
BLOCK_SIZE = 4096*256
BLOCKS_PER_SUPERBLOCK = 8 # TODO: on just the GPU, this should be a lot bigger
SUPERBLOCK_SIZE = BLOCK_SIZE*BLOCKS_PER_SUPERBLOCK
LOG_SIEVE_THRESHOLD = math.log(LP_BOUND) + 0.5

def gen_prime(bits):
  ret = random.randint((1 << (bits-1))+1, 1 << bits)
  while not isprime(ret): ret += 1
  return ret

# generate number for factoring (N)
while 1:
  p,q = gen_prime(BITS//2), gen_prime(BITS//2)
  if p==q: continue
  N = p*q
  break
print(f"factoring {N} into {p} {q} with {N.bit_length()} bits")
del p,q # no cheating

# generate primes up to B filtered by quadratic residue
# https://en.wikipedia.org/wiki/Euler%27s_criterion
FACTOR_BASE = [2] + [p for p in range(3, B+1, 2) if isprime(p) and pow(N, (p-1)//2, p) == 1]
NUM_RELATIONS = len(FACTOR_BASE) + max(8, len(FACTOR_BASE) // 10)
LS_SCALE = 1 << 32
print(f"{B=} {max(FACTOR_BASE)=} {len(FACTOR_BASE)=}")
print(f"{NUM_RELATIONS=} {LOG_SIEVE_THRESHOLD=:.2f} {BLOCK_SIZE=} {SUPERBLOCK_SIZE=} {math.log(LS_SCALE):f=}")

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

class QFunction:
  def __init__(self, N, A, B):
    # required
    assert (B*B - N) % A == 0
    self.N = N
    C = (B*B - N) // A
    self.A, self.B, self.C = A, B, C
  def __call__(self, x):
    return self.A*x*x + 2*self.B*x + self.C

def do_block_sieve_gpu(Q:QFunction, gpu_roots, gpu_p, x_block):
  x = Tensor([x_block]) + Tensor.arange(0, BLOCK_SIZE)
  #sieve = Q(rng).abs().log()
  sieve = ((float(Q.A)/LS_SCALE)*x*x + 2*(float(Q.B)/LS_SCALE)*x + (float(Q.C)/LS_SCALE)).abs().log() + math.log(LS_SCALE)
  sieve = sieve - ((((gpu_roots - x.reshape(-1, 1)) % gpu_p) == 0).float() * gpu_p.log()).sum(1)
  hit = (sieve<LOG_SIEVE_THRESHOLD).tolist()

  ret = []
  # check for success
  for j in range(BLOCK_SIZE):
    if hit[j]: ret.append(x_block+j)
  return ret

def do_block_sieve(Q, ROOTS_LIST, x_block):
  ret = []
  # we approximate log(Q(x)) in float
  scores = [math.log(abs(Q(x))) for x in range(x_block, x_block+BLOCK_SIZE)]

  # sieve with the roots
  for root,p,log_p in ROOTS_LIST:
    j = (root - x_block) % p
    while j < BLOCK_SIZE:
      scores[j] -= log_p
      j += p

  # check for success
  for j in range(BLOCK_SIZE):
    if scores[j] < LOG_SIEVE_THRESHOLD: ret.append(x_block+j)
  return ret

def roots_for_factor_base(Q:QFunction, FACTOR_BASE):
  ROOTS = {}
  for p in FACTOR_BASE:
    if p == 2:
      # NOTE: this works for everything
      ROOTS[p] = [x for x in range(p) if Q(x)%p == 0]
    elif Q.A % p == 0:
      ROOTS[p] = [(-Q.C * pow(2*Q.B, -1, p)) % p]
    else:
      # odd primes have two roots
      # 0 == ((x+a)^2 - N) % p
      # N == (x+a)^2 % p
      s = sqrt_mod(N, p)  # s*s == N % p
      inv_A = pow(Q.A, -1, p)
      # so x+a == +/- s % p
      r1 = (( s - Q.B) * inv_A) % p
      r2 = ((-s - Q.B) * inv_A) % p
      assert r1 != r2
      ROOTS[p] = [r1, r2]
  return ROOTS

@dataclass
class RelationState:
  relations: list = field(default_factory=list)
  partials: dict = field(default_factory=dict)
  seen: set = field(default_factory=set)
  # tracking
  searched: int = 0
  log_sieve_false_positive:int = 0
  log_sieve_duplicate:int = 0
  progress: tqdm.tqdm = field(default_factory=lambda: tqdm.tqdm(total=NUM_RELATIONS))

def qsieve_get_relations(N, A, B, superblock_schedule_order, rs:RelationState):
  # A must fully factorize
  A_relation, rem = b_smooth_factorize(A)
  assert rem == 1
  Q = QFunction(N, A, B)

  st = time.perf_counter()
  # compute the ROOTS for each prime in FACTOR_BASE
  # Q(x) % p == 0 (solve for x)
  ROOTS = roots_for_factor_base(Q, FACTOR_BASE)
  #print(f"computed {len(ROOTS)=} in {time.perf_counter()-st:.2f} s")

  # precompute the logs of the primes and put roots in a list
  ROOTS_LIST = []
  for p in FACTOR_BASE:
    log_p = math.log(p)
    for root in ROOTS[p]:
      ROOTS_LIST.append((root,p,log_p))

  if GPU:
    my_superblock_sieve = functools.partial(do_block_sieve_gpu, Q,
                                            Tensor([x[0] for x in ROOTS_LIST]),
                                            Tensor([x[1] for x in ROOTS_LIST]))
  else:
    my_superblock_sieve = functools.partial(do_block_sieve, Q, ROOTS_LIST)

  # first we need to find B-smooth Q(x) values
  # we use a log sieve to prefilter everything
  st = time.perf_counter()
  partial_match = 0
  sieve_time_s = 0.0

  begin_relations = len(rs.relations)
  for x_superblock in superblock_schedule_order:
    start_relations = len(rs.relations)
    # here we extract the real relations from the log_sieve
    inner_st = time.perf_counter()
    sieved = my_superblock_sieve(x_superblock)
    sieve_time_s += time.perf_counter() - inner_st
    for x in sieved:
      lhs = A*x + B
      key = lhs%N
      if key in rs.seen:
        rs.log_sieve_duplicate += 1
        continue
      q = Q(x)
      relation, num = b_smooth_factorize(abs(q))
      # the A relation needs to be included in all relations for that polynomial
      relation = [u + v for u, v in zip(A_relation, relation)]
      rs.seen.add(key)
      neg = q < 0
      if num == 1: rs.relations.append((lhs, relation, neg, 1))
      elif num <= LP_BOUND and isprime(num):
        if num in rs.partials:
          # got a match!
          lhs0, relation0, neg0 = rs.partials[num]
          combined_relation = [u + v for u, v in zip(relation0, relation)]
          rs.relations.append((lhs*lhs0, combined_relation, neg^neg0, num))
          del rs.partials[num]
          partial_match += 1
        else:
          rs.partials[num] = (lhs, relation, neg)
      else: rs.log_sieve_false_positive += 1
    rs.searched += 1
    rs.progress.set_description(f"{rs.searched:5d} b, {(len(rs.relations)-start_relations):3d} r/b, "
                             f"{(len(rs.relations)-begin_relations)-partial_match:3d}+{partial_match:3d} p, "
                             f"{rs.log_sieve_false_positive*100./(len(rs.relations)+rs.log_sieve_false_positive):.1f}% fp, "
                             f"{rs.log_sieve_duplicate} dup, A={Q.A}")
    rs.progress.update(min(NUM_RELATIONS, len(rs.relations))-rs.progress.n)
    if len(rs.relations) >= NUM_RELATIONS: break
  et = time.perf_counter()-st
  #print(f"collected {(len(rs.relations)-start_relations)} relations across {searched} blocks in {et:.2f} s, {et*1000./searched:.2f} ms/superblock")
  #print(f"{len(rs.partials)=} unmatched partial")
  #print(f"got {log_sieve_false_positive=} {log_sieve_duplicate=} with {sieve_time_s:.2f} s in the sieve")

def qsieve(N):
  import random
  rs = RelationState()
  while 1:
    #p1 = random.choice(FACTOR_BASE)
    #p2 = random.choice(FACTOR_BASE)
    #p3 = random.choice(FACTOR_BASE)
    #if p1 == p2 or p2 == p3 or p1 == p3: continue
    #A = p1*p2*p3
    A = random.choice([1]+FACTOR_BASE)
    for fudge in range(-10000, 10000):
      B = math.isqrt(N)+fudge
      if (B*B - N) % A == 0: break
    else:
      #print(f"A={A} is bad")
      continue

    # NOTE: we explore both negative and positive x
    superblock_schedule_order = itertools.chain.from_iterable(zip(
      range(1, SUPERBLOCK_SIZE, BLOCK_SIZE),
      range(-BLOCK_SIZE, -SUPERBLOCK_SIZE, -BLOCK_SIZE)))
    qsieve_get_relations(N, A, B, superblock_schedule_order, rs)
    if len(rs.relations) >= NUM_RELATIONS: break
  rs.progress.close()

  # then we need to solve to make a perfect square from the relations
  # we need to find a basis among the parity masks
  print(f"matrix size is {len(rs.relations)} x {len(FACTOR_BASE)+1}")
  st = time.perf_counter()
  congruence_false_positive = 0
  basis = {}
  for i, (_, relation, neg, _) in enumerate(rs.relations):
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
      if process_congruence(N, rs.relations, combo): break
      else: congruence_false_positive += 1
  else:
    raise RuntimeError("failed to find congruence")
  print(f"solved with {congruence_false_positive=} in {time.perf_counter()-st:.2f} s")
  assert congruence_false_positive < 10

if __name__ == "__main__":
  qsieve(N)

