plot(obs$soc.t.ha ~ obs$time, ylim=c(0, NA))
lines(run_fiteq$soc ~ run_fiteq$time, col =2)
lines(run_fittr$soc ~ run_fittr$time, col =2)
