
3 Options for decomposition:
- 1st: first order: V * C_P
- 2nd: second order: V * C_P * C_E
- MM: Michaelis-Menten: V * C_P * C_E / (K + C_P)

3 Options for uptake:
- 1st: equal to diffusion: diff
- 2nd: second order: V * diff * C_M
- MM: Michaelis-Menten: V * diff * C_M / (K + diff)

submit with:
bsub -q mpi -W 6:00 -o RUN-%J.out -n 1 Rscript --slave Optim.R S1P1
