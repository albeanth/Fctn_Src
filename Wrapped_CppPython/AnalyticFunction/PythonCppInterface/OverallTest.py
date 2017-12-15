import TestInt
from math import pi
from numpy import zeros
import time

go = TestInt.GaussianIntegration()

## 3D TEST
start = time.perf_counter()
test_integral_serial = go.GaussInt_3D_Serial(0.0,pi,100, 0.0,pi,100, 0.0,pi,100)
end = time.perf_counter()
SerialTime = end-start
print('Serial Calulation took {0:.4f}'.format(SerialTime))

Threads = [2,4,6,8,12,16,20,24,32]
ParallelTime = zeros(len(Threads))

for idx,NUMTHREADS in enumerate(Threads):
  start = time.perf_counter()
  test_integral_parallel = go.GaussInt_3D_Parallel(NUMTHREADS, 0.0,pi,100, 0.0,pi,100, 0.0,pi,100)
  end = time.perf_counter()
  ParallelTime[idx] = end-start
  #print('Parallel Calulation took {0:.4f}'.format(ParallelTime))

print('\n\n----- Summary Table -----')
print('# Threads    Time       Speedup   F_parallel')
print('     {0: d}      {1:.4f}       --         --'.format(1, SerialTime))
for idx in range(0,len(Threads)):
  print('     {0: d}      {1:.4f}     {2:.4f}     {3:.4f}'.format(Threads[idx], ParallelTime[idx], SerialTime/ParallelTime[idx], ((ParallelTime[idx]/SerialTime) - 1)*(Threads[idx]/(1-Threads[idx]))))

  
  
