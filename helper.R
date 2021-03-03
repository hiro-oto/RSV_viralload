logP = function(Datas, pars)
{
  pars$mu = makemu(Datas, pars)
  dens1 = dnorm(Datas$VLs, pars$mu, pars$sigma, log = T)*Datas$SW
  dens2 = pnorm(Datas$VLs, pars$mu, pars$sigma, log = T)*Datas$CS
  dens3 = dlnorm(Datas$Onsetday - pars$Tinf[,1], pars$logmu, pars$logsd, log = T)*Datas$SWonset
  dens4 = dnorm(pars$logmu,1.48,0.0042,log=TRUE)
  dens5 = dnorm(pars$logsd,0.22,0.056,log=TRUE)
  pars$LL = sum(dens1) + sum(dens2) + dens4 + dens5 + sum(dens3, na.rm = T)
  return(pars) 
}


mh = function(oldpars, newpars, Datas)
{
  reject = FALSE
  if(newpars$sigma <= 0) reject = TRUE
  if(newpars$theta1 <= 0) reject = TRUE
  if(newpars$theta2 <= 0) reject = TRUE
  if(newpars$theta3 <= 0) reject = TRUE
  if(max(newpars$Tinf) >= 14) reject = TRUE
  if(min(newpars$Tinf) < -14) reject = TRUE
  
  i = which(Datas$SWonset==0)# Without onset
  if(max(pars$DifminPosTinf[i]) >10) reject = TRUE
  if(min(pars$DifminPosTinf[i]) < 2) reject = TRUE
  
  j = which(Datas$SWonset==1)# With onset
  if(max(newpars$DifTinfOnset[j]) > 7) reject = TRUE
  
  if(newpars$logmu < 0) reject = TRUE
  if(newpars$logsd < 0) reject = TRUE
  
  
  if(!reject)
  {
    newpars = logP(Datas, newpars)
    if(is.na(newpars$LL))reject=TRUE
  }
  if(!reject)
  {
    logaccept = newpars$LL - oldpars$LL
    ln = -rexp(1)  
    if(ln > logaccept) reject = TRUE
  }
  if(reject) return(oldpars)
  return(newpars)
}


find_abc = function(pars)
{
  temp   = pars$theta1*pars$theta2/pars$theta3
  optim  = optimize(function(a) abs(exp(a*log(a-1)+(1-a)-lgamma(a))-temp), interval=c(1,2000000))
  pars$a = optim$minimum
  pars$b = (pars$a-1)/pars$theta1
  pars$c = pars$theta3
  return(pars)
}

initialise_pars = function(file,data=datas)
{
  loadpars = TRUE
  if(!file.exists(file))loadpars=FALSE
  if(!loadpars) # no previous run?! make something up then!
  {
    Tinf = datas$Onsetday -6
    Tinf[asympt] = initialTinf+1
    
    DifTinfOnset = datas$Onsetday - Tinf
    
    Q = datas$TS
    for(k in 1:dim(Q)[2]){Q[,k] = Tinf}
    pars = list(mu = 1, sigma = 0.1, Tinf = Q, DifTinfOnset = DifTinfOnset,
                logmu = 1.5, logsd = 0.2,
                theta1 = 3, theta2 = 4, theta3 = 50)
  }
  if(loadpars) # use the last iteration of the last run
  {
    load('storage/mcmc.rdata')
    k = length(storage$sigma)
    pars = list()
    pars$sigma = storage$sigma[k]
    pars$logmu = storage$logmu[k]
    pars$logsd = storage$logsd[k]
    
    pars$theta1 = storage$theta1[k]
    pars$theta2 = storage$theta2[k]
    pars$theta3 = storage$theta3[k]
    pars$Tinf = matrix(0, datas$Npax, dim(datas$VLs)[2])
    for(j in 1:dim(datas$VLs)[2])pars$Tinf[,j] = storage$Tinf[k,]
  }
  pars$DifTinfOnset = datas$Onsetday - pars$Tinf[,1] 
  return(pars)
}


makemu = function(Datas, pars)
{
  mu = matrix(0, Datas$Npax, dim(datas$VLs)[2])
  for(j in 1:dim(Datas$VLs)[2])
  {
    i = which(Datas$TS[,j] <= pars$Tinf[,j])
    mu[i,j] = 0
    i = which(Datas$TS[,j] > pars$Tinf[,j])
    deltat = Datas$TS[i,j] - pars$Tinf[i,j]
    mu[i,j] = pars$c*(pars$b^pars$a)*deltat^(pars$a-1)*exp(-pars$b*deltat)/gamma(pars$a)
  }
  return(mu)
}



