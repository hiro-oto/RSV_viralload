library(mvtnorm)
library(matrixcalc)
abcde = c(1.482,0.213,0.00000813,0,0.00299)
set.seed(874)
F1 = rnorm(10000,0,1)
f2 = rnorm(10000,0,1)

errors = function(abcde)
{
  if(abcde[3]<0)return(99999999) 
  if(abcde[5]<0)return(99999999)
  
  #set.seed(9543)
  #SIGMA = matrix(c(abcde[3],abcde[4],abcde[4],abcde[5]),2,2)
  #if(!is.positive.semi.definite(SIGMA))return(99999999)
  
  mi = F1*sqrt(abcde[3])+abcde[1]
  si = F1*sqrt(abcde[5])+abcde[2]
  #si = temp[,2]
  si2 = abs(si)
  T1 = qlnorm(0.05,mi,si2)
  T2 = qlnorm(0.25,mi,si2)
  T3 = qlnorm(0.50,mi,si2)
  T4 = qlnorm(0.75,mi,si2)
  T5 = qlnorm(0.95,mi,si2)
  q1 = quantile(T1,c(0.025,0.5,0.975))
  q2 = quantile(T2,c(0.025,0.5,0.975))
  q3 = quantile(T3,c(0.025,0.5,0.975))
  q4 = quantile(T4,c(0.025,0.5,0.975))
  q5 = quantile(T5,c(0.025,0.5,0.975))
  e1 = c(2.5,3.1,3.8)
  e2 = c(3.3,3.8,4.4)
  e3 = c(3.9,4.4,4.9)
  e4 = c(4.5,5.1,5.7)
  e5 = c(5.2,6.3,7.3)
  
  output = sum(abs(e1-q1)+abs(e2-q2)+abs(e3-q3)+abs(e4-q4)+abs(e5-q5))
  output = output + 100*sum(si<0)
  return(output)
}

obj = errors(abcde)

SIMS = 10000
storage = matrix(0,SIMS,6)
for(sim in 1:SIMS)
{
  if(sim%%100==0)cat(sim,'of',SIMS,abcde,obj,'\n')
  abcde0 = abcde
  obj0 = obj
  
  abcde = rnorm(5,abcde,c(0.001,0.001,0.000001,1,0.0001))
  obj = errors(abcde)
  if(obj>obj0){abcde=abcde0;obj=obj0}
  storage[sim,1:5]=abcde
  storage[sim,6] = obj
}

abcde <-
  c(1.4816758100025733, 0.21751664240190205, 1.7336505088157274e-05, 
    1.6637171054639046, 0.0031156839503659749)

mi = F1*sqrt(abcde[3])+abcde[1]
si = F1*sqrt(abcde[5])+abcde[2]

T1 = qlnorm(0.05,mi,si)
T2 = qlnorm(0.25,mi,si)
T3 = qlnorm(0.50,mi,si)
T4 = qlnorm(0.75,mi,si)
T5 = qlnorm(0.95,mi,si)


