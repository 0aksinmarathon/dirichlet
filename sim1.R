library(mvtnorm)

DPsim = function(sigma,alpha, n, t, v_n, var_variance, error_variance, fix_eff, epsilon = 0.01){
  v=1
  k=1
  W = c()
  V = c()
  B = t(diag(dim(sigma)[1])[1,1:dim(sigma)[1]])
  
  while (sum(W) < 1 - epsilon){
    v = rbeta(1, 1, alpha)
    b = rmvnorm(1, mean = rep(0, dim(sigma)[1]), sigma = sigma)
    
    if (k == 1){
      B = t(data.frame(c(b)))
      W[k] = v
    }else{
      B = rbind(B, b)
      W[k] = v * prod(V)
    }
    V[k] = 1-v
    k = k + 1  
  }
  
  print(W)
  print(B)
  print(sum(W))
  
  b = rmvnorm(1, mean = rep(0, dim(sigma)[1]), sigma = sigma)
  v = rbeta(1, 1, alpha)
  W[k] = 1 - sum(W)
  B = rbind(B, b)
  print(W)
  print(B)
  print(sum(W))
  
  int = rep(1, n*t)
  time = rep(1:t,n)
  time2 = time**2
  time2
  df = cbind(int, time, time2)
  df
  
  vnames = c()
  for (i in 1:v_n){
    gen = rnorm(n, 0, var_variance)
    var = rep(gen, each = t)
    df = cbind(df, var)
    varname = paste("var", i, sep = "")
    vnames = c(vnames, varname)
  }
  
  colnames(df) = c("int", "time", "time2", vnames)
  
  rand_assign = rmultinom(n, prob = W, size = 1)
  rand_eff = t(rand_assign)%*%B
  rand_effs = apply(rand_eff, 2, rep, each = t)
  df[,1:dim(sigma)[1]] = df[,1:dim(sigma)[1]] + rand_effs
  df
  
  pure_error = rnorm(n*t, 0, error_variance)
  
  
  resp = (df %*% fix_eff) + pure_error
  ID = rep(1:n, each = t)
  df = cbind(ID, resp, df, pure_error, rand_effs)
  colnames(df)[2] = "resp"
  return(list(as.data.frame(df), rand_eff))
}

sigma = diag(c(3,2,1))
alpha = 3
epsilon = 0.01
v = 1
n = 500
t = 10
v_n = 3
var_variance = 1
error_variance = 3
fix_eff = c(1, 1, 1, 1, 1,1)

df = DPsim(sigma=sigma,alpha=alpha, n=n, t=t, v_n=v_n, var_variance= var_variance, error_variance=error_variance, fix_eff=fix_eff, epsilon = epsilon)


mcmc <- list(nburn=3000,nsave=10000,nskip=1,ndisplay=1000)

prior_dp1_sim<-list(a0 = 3, b0 = 3, nu0=4.01,tau1=0.01,tau2=0.01, beta0 = rep(0,v_n+(3-dim(sigma)[1])), Sbeta0 = diag(1000,v_n+(3-dim(sigma)[1])),
                    tinv=diag(10,dim(sigma)[1]),mub=rep(0,dim(sigma)[1]),Sb=diag(1000,dim(sigma)[1]))

prior_dpm1_sim<-list(a0 = 3, b0 = 3, beta0 = rep(0,v_n+(3-dim(sigma)[1])), Sbeta0 = diag(1000,v_n+(3 -dim(sigma)[1])), tau1=0.01,tau2=0.01, nu0=4.01, 
                     tinv=diag(10,dim(sigma)[1]), nub=4.01, tbinv=diag(10,dim(sigma)[1]), mb=rep(0,dim(sigma)[1]), Sb=diag(1000,dim(sigma)[1]))

fit1a_dp_sim <- DPlmm(fixed=resp~var1+var2+var3,
                      random=~time+time2|ID, prior=prior_dp1_sim,mcmc=mcmc,
                      state=state,status=TRUE, data = df)

fit1a_dpm_sim <- DPMlmm(fixed=resp~var1+var2+var3,
                        random=~time+time2|ID, prior=prior_dpm1_sim,mcmc=mcmc,
                        state=state,status=TRUE, data = df)

fit1_sim<-MCMCglmm(resp~var1+var2+var3+time+time2, random=~us(1+time+time2):ID, data=df, family = "gaussian", verbose=FALSE,
                   nitt=10000, burnin=3000, thin=1)
summary(fit1a_dp_sim)
summary(fit1a_dpm_sim)
summary(fit1_sim)

sigma = diag(c(3,2,1))
alpha = 1
epsilon = 0.01
v = 1
n = 500
t = 10
v_n = 3
var_variance = 1
error_variance = 3
fix_eff = c(1, 1,1,1,1,1)

index = seq(1, 5000, 10)

df = DPsim(sigma=sigma,alpha=alpha, n=n, t=t, v_n=v_n, var_variance= var_variance, error_variance=error_variance, fix_eff=fix_eff, epsilon = epsilon)
df 

v10 = as.data.frame(df[1])$V10 # first random effect: v10
hist(v10, breaks = 1000) #histrogram of random effects 
table(v10[index]) # frequency table

uni_sort = sort(unique(v10)) # unique list + sorted
length(unique(v10)) # number of unique rand effects

grid = c(seq(-7,7, 0.01))

d = c()
for (i in 1: length(uni_sort)){
  d = rbind(dnorm(x = grid , mean = uni_sort[i]) * table(v10[index])[i],  d)
}
# d: for each row, density of norm dist with mean = random effect and weight = frequency

plot(apply(X = d, MARGIN = 2, FUN = sum), cex = 0.1)
# for each column (for each point in the grid), sum up the weighted densities

d2 = apply(X = d, MARGIN = 2, FUN = sum)
d2 = as.data.frame(t(rbind(d2, grid)))

colnames(d2) = c("density", "x")
ggplot(d2,aes(x=x,y=density))+
  geom_line()
