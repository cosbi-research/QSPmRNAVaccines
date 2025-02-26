function result = logmeanandci(vect)

n = length(vect);
log_vect = log(vect);
log_mean = mean(log_vect);
log_stdev = std(log_vect);
t_crit = tinv(1-0.05/2, n-1);
log_ci_low = log_mean - t_crit*log_stdev/sqrt(n); 
log_ci_upp = log_mean + t_crit*log_stdev/sqrt(n);

mean_v = exp(log_mean);
low95ci = exp(log_ci_low);
upp95ci = exp(log_ci_upp);
result = [mean_v, low95ci, upp95ci];
end

