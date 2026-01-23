devtools::load_all()
set.seed(123)
tmpl <- tmpl_stsp_full(1,1,2,sigma_L='chol')
model <- r_model(tmpl, bpoles=1, sd=0.5)
data <- sim(model, n.obs=50, a1=NA)
y <- data$y
pf_result <- pfilter(model, y, N_particles=500, method='sir')
cat('Class:', class(pf_result), '\n')
cat('Names:', names(pf_result), '\n')
# Test each plot type
try({
  pdf(file=NULL)
  plot(pf_result, type='states')
  plot(pf_result, type='ess')
  plot(pf_result, type='weights')
  plot(pf_result, type='likelihood')
  dev.off()
  cat('All plots succeeded\n')
})