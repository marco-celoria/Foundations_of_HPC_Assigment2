import numpy as np
from matplotlib import pyplot as plt

st1_omp=36.8572
st2_omp=21.6287
st4_omp=13.3102
st8_omp=9.72101
st16_omp=7.26727

st1_mpi=35.735761
st2_mpi=20.225521
st4_mpi=11.679689
st8_mpi=6.830646
st16_mpi=3.889892


data0 = np.array([
 [1.000,  st1_omp],
 [2.000,  st2_omp],
 [4.000,  st4_omp],
 [8.000,  st8_omp],
 [16.000, st16_omp],
])
x0, y0 = data0.T
plt.plot(x0,y0,color='tab:blue',linestyle='dotted',marker='o',label="OpenMP")

data1 = np.array([
[1.000,   st1_mpi],
 [2.000,  st2_mpi],
 [4.000,  st4_mpi],
 [8.000,  st8_mpi],
 [16.000, st16_mpi],
])
x1, y1 = data1.T
plt.plot(x1,y1,color='tab:red',linestyle='dashed', marker='v',label="MPI")

plt.xlabel("P")
plt.ylabel("T [s]")
plt.legend(loc="upper right")
plt.title("Time of KdTree build function")
plt.savefig("Strong_time.png")
plt.clf()

data0 = np.array([
 [1.000,  st1_omp/ st1_omp],
 [2.000,  st1_omp/ st2_omp],
 [4.000,  st1_omp/ st4_omp],
 [8.000,  st1_omp/ st8_omp],
 [16.000, st1_omp/ st16_omp],
])
x0, y0 = data0.T
plt.plot(x0,y0,color='tab:blue',linestyle='dotted',marker='o',label="OpenMP")

data1 = np.array([
[1.000,   st1_mpi/ st1_mpi],
 [2.000,  st1_mpi/ st2_mpi],
 [4.000,  st1_mpi/ st4_mpi],
 [8.000,  st1_mpi/ st8_mpi],
 [16.000, st1_mpi/ st16_mpi],
])
x1, y1 = data1.T
plt.plot(x1,y1,color='tab:red',linestyle='dashed', marker='v',label="MPI")

#plt.xscale("log")
#plt.yscale("log")
plt.xlabel("P")
plt.ylabel(r"$S_p=T(1)/T(P)$")
plt.legend(loc="upper left")
plt.title("Strong Speedup")
plt.savefig("Strong_speedup.png")
plt.clf()



wt1_omp=36.5936
wt2_omp=44.6332
wt4_omp=58.5614
wt8_omp=78.8836
wt16_omp=118.857


wt1_mpi=35.882094
wt2_mpi=41.431630
wt4_mpi=49.488500
wt8_mpi=58.178536
wt16_mpi=67.924564


data0 = np.array([
 [1.000, wt1_omp],
 [2.000, wt2_omp],
 [4.000, wt4_omp],
 [8.000, wt8_omp ],
 [16.000, wt16_omp],
])
x0, y0 = data0.T
plt.plot(x0,y0,color='tab:blue',linestyle='dotted',marker='o',label="OpenMP")

data1 = np.array([
  [1.000, wt1_mpi],
 [2.000, wt2_mpi],
 [4.000, wt4_mpi],
 [8.000, wt8_mpi ],
 [16.000, wt16_mpi],
])
x1, y1 = data1.T
plt.plot(x1,y1,color='tab:red',linestyle='dashed', marker='v',label="MPI")

plt.xlabel("P")
plt.ylabel(r"$T_N(P)$"+" [s]")
plt.legend(loc="upper left")
plt.title("Time of KdTree build function")
plt.legend(loc="upper left")
plt.savefig("Weak_time.png")
plt.clf()


data0 = np.array([
  [1.000, wt1_omp/wt1_omp],
 [2.000, wt1_omp/wt2_omp],
 [4.000, wt1_omp/wt4_omp],
 [8.000, wt1_omp/wt8_omp ],
 [16.000, wt1_omp/wt16_omp],
])
x0, y0 = data0.T
plt.plot(x0,y0,color='tab:blue',linestyle='dotted',marker='o',label="OpenMP")

data1 = np.array([
  [1.000, wt1_mpi/wt1_mpi],
 [2.000, wt1_mpi/wt2_mpi],
 [4.000, wt1_mpi/wt4_mpi],
 [8.000, wt1_mpi/wt8_mpi ],
 [16.000, wt1_mpi/wt16_mpi],
])
x1, y1 = data1.T
plt.plot(x1,y1,color='tab:red',linestyle='dashed', marker='v',label="MPI")
plt.ylabel(r"$E_w=T(1)/T_N(P)$")
plt.xlabel("P")
plt.legend(loc="upper right")
plt.title("Weak Efficiency")
plt.savefig("Weak_speedup.png")
plt.clf()



