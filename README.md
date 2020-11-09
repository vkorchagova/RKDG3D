# RKDG3D v.0.1

MFEM-based RKDG compressible flow solver for 3D flows

EoS: perfect gas, Covolume

Meshes: unstructured (supports triangular and quadrangle elements) supported by MFEM 

Riemann solvers: Rusanov, Local Lax --- Fridriechs, HLL, HLLC

Averaging for Riemann solvers: Einfeldt, TRRS

Limiters: FinDiff, BJ, Venkatakrishnan, Michalak

Boundary conditions: Constant (fixed), Slip, Open

Time step: static, dynamic

Parallelisation: MPI (by MFEM)
