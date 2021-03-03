# RSV-viral-load

set.seed(123)

##---- Load functions
source('code/helper.r')

##---- Load data
load('combined_data/data_viralload_2020-12-21.Rdata')

datas$SW = 0*datas$CS
for(i in 1:dim(datas$SW)[1]){for(j in 1:dim(datas$SW)[2]){if(datas$VLs[i,j]>1)datas$SW[i,j]=1}}
rm(i,j)
datas$SWonset = ifelse(datas$Onsetday == -100, 0, 1)
datas$Onsetday = ifelse(datas$Onsetday == -100, NA, datas$Onsetday)

#-----------Initial run only-----------
asympt = which(is.na(datas$Onsetday) == T)

temp = c()
for(pax in 1:datas$Npax)
{
  t = which(datas$SW[pax,]==1)
  temp[pax] = datas$TS[pax, min(t)]
}
datas$minPositive = temp

initialTinf = rep(0, length(asympt))
initialTinf = datas$minPositive[asympt]-6 #Initial run only
#-----------Initial run only END---------

##---- Initialise
pars = initialise_pars('storage/mcmc.rdata',datas)
pars$DifminPosTinf = datas$minPositive - pars$Tinf[,1] 
#--------------
pars = find_abc(pars)
pars = logP(datas,pars)

##---- MCMC
MCMCITS = 10000
THINNING= 10
storage = list(
  sigma = rep(0, MCMCITS),
  theta1 = rep(0, MCMCITS),
  theta2 = rep(0, MCMCITS),
  theta3 = rep(0, MCMCITS),
  a = rep(0, MCMCITS),
  b = rep(0, MCMCITS),
  c = rep(0, MCMCITS),
  Tinf = matrix(0, MCMCITS, datas$Npax),
  logmu = rep(0, MCMCITS),
  logsd = rep(0, MCMCITS),
  LL = rep(0, MCMCITS))

for (iteration in 1:MCMCITS)
{
  if(iteration %% 100 == 0) cat(iteration, ' in ', MCMCITS,' (logL = ',round(pars$LL),')\n',sep='')
  
  for(thin in 1:THINNING)
  {
    
    oldpars = pars; pars$theta1 = rnorm(1, pars$theta1, 10); pars = find_abc(pars); pars = mh(oldpars, pars, datas)
    oldpars = pars; pars$theta2 = rnorm(1, pars$theta2, 10); pars = find_abc(pars); pars = mh(oldpars, pars, datas)
    oldpars = pars; pars$theta3 = rnorm(1, pars$theta3, 10); pars = find_abc(pars); pars = mh(oldpars, pars, datas)
    
    oldpars = pars; pars$sigma = rnorm(1, pars$sigma, 0.1); pars = mh(oldpars, pars, datas)
    oldpars = pars; pars$logmu = rnorm(1, pars$logmu, 0.1); pars = mh(oldpars, pars, datas)
    oldpars = pars; pars$logsd = rnorm(1, pars$logsd, 0.1); pars = mh(oldpars, pars, datas)
    
    for (pax in 1:datas$Npax)
    {
      oldpars = pars
      pars$Tinf[pax,] = rnorm(1, pars$Tinf[pax,], 1)
      pars$DifTinfOnset[pax] = datas$Onsetday[pax] - pars$Tinf[pax,1]
      pars$DifminPosTinf[pax] = datas$minPositive[pax] - pars$Tinf[pax,1]
      pars = mh(oldpars, pars, datas)
    }
  }  
  storage$sigma[iteration] = pars$sigma
  storage$theta1[iteration] = pars$theta1 
  storage$theta2[iteration] = pars$theta2 
  storage$theta3[iteration] = pars$theta3 
  storage$a[iteration] = pars$a
  storage$b[iteration] = pars$b
  storage$c[iteration] = pars$c
  storage$Tinf[iteration,] = pars$Tinf[,1]
  storage$logmu[iteration] = pars$logmu
  storage$logsd[iteration] = pars$logsd
  storage$LL[iteration] = pars$LL
}