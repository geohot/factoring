from qsieve import gen_semiprime

N = gen_semiprime(16)

a = 2

hits = 0
for r in range(2**16):
  if pow(2, r, N) == 1:
    print("HIT", r)
    hits += 1
  if hits > 10: break
