import numpy as np
import time

start_time = time.time()
s = 0.0
for counter in range(100000):
  s += np.random.random()
print s, time.time() - start_time

start_time = time.time()
rands = list(np.random.random(100000))
for counter in range(100000):
  s += rands[counter]
print s, time.time() - start_time
