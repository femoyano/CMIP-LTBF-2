
submit with:
bsub -q mpi -W 6:00 -o RUN-%J.out -n 1 Rscript --slave Optim.R 1
