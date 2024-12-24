# MGndeAdjointFoam
A neutron space-time kinetics solver based on OpenFOAM 10

The neutron diffusion solver includes 1 steady-state solver and 5 transient-state solvers.

The power iteration method is used to calculate keff in steady-state solver.
The transient solver adopts Driect Method(DM), Adiabatic Method, Improved Quasi-Static(IQS) method, Predictor-Corrector Quasi-Static(PCQS) method and Point-Kinetics Equation(PKE-$\psi_0$) method.

The results of 2D-TWIGL benchmark, 3D-IAEA benchmark, 3D-LMW benchmark, 1D-SMR problem, 3D-SMR problem are summarized here.
