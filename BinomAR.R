###########################################
## Accept Reject algorithm to draw from
## Binomial(n,p)
###########################################
# setting the seed makes it so that the same sets of
# random variables are realized.
set.seed(1)  

# Function draws one value from Binom(n,p)
# n = number of trials
# p = probability of success
draw_binom <- function(n, p)
{
  accept <- 0 # Will track the acceptance
  try <- 0 # Will track the number of proposals
  
  # upper bound calculated in the notes
  x <- 0:n
  #all_c <- choose(n,x) * (1-p)^(n - 2*x) * p^(x-1)  # from notes
  p.star <- 1/(n*p + 1)
  all_c = choose(n,x)*((p)^x) * ((1-p)^(n-x))
  all_c = all_c /(p.star * (1 -p.star)^x)
  c <- max(all_c)  + 0.001 # what is the value of c ?
                            #as c inreases , loops increase. For accuracy, add 0.001 
                            #0.001 is for machine tolerance 
  while(accept == 0)
  {
    try <- try + 1
    
    U <- runif(1) 
    prop <- rgeom(1, prob = p.star) 
    
    ratio= choose(n,prop)*((p)^prop) * ((1-p)^(n-prop))
    ratio = ratio/(p.star * (1 -p.star)^prop)
    ratio = ratio/c
    
    if(U < ratio)
    {
      accept <- 1
      rtn <- prop
    }
  }
  return(c(rtn, try ,c))
}
draw_binom(n = 10, p = .25 )  #look for numerical instability 


###
# If we want X1, ..., Xn ~ Binom(n.p)
# we need to call the function multiple times

# sample size
N <- 1e3
samp <- numeric(N)
n.try <- numeric(N)
for(t in 1:N)
{
  # I use as a dummy variable often
  foo <- draw_binom(n = 10, p = .25)
  samp[t] <- foo[1]
  n.try[t] <- foo[2]
}
mean(samp) #should be n*p = 2.5
mean(n.try)


###########################################
## A closer look at Binomial and Geometric
###########################################
# Turns out, this choice of Binomial and Geometric
# can work, but not always. In the code below, 
# increase n to see what happens

p <- .25
n <- 1000
x <- 0:(n)
mass.geom <- dgeom(x, p)
mass.bin <- dbinom(x, size = n, prob = p)

all_c <- choose(n,x) * (1-p)^(n - 2*x) * p^(x-1)
(c <- max(all_c))

#mass of geometric on the same support as binomial in red 
plot(x, mass.geom, pch = 16, col = "red", type= "n")

#mass of binomial in blue
points(mass.bin, pch = 16, col = "red", type= "h")

#mass of geom which does not have common support in dashed blue
points(mass.geom, pch = 16, col = "blue", type = "h", lty = 2)


# Matching the means:
# choosing p* for rgeom so that np = (1-p*)/p*
p.star <- 1/(n*p + 1)
mass.geom <- dgeom(x, p.star)
mass.bin <- dbinom(x, size = n, prob = p)

all_c <- choose(n,x) * (1-p.star)^(n - 2*x) * p.star^(x-1)
(c <- max(all_c))


plot(mass.geom, pch = 16, col = "red", type= "n")
points(mass.bin, pch = 16, col = "red", type= "h")
points(mass.geom, pch = 16, col = "blue", type = "h", lty = 2)





##############################
## Truncated poisson 
##############################
#m = 30 , lamda = 20 
#taking poisson as target distribution 

draw_tpois <- function(m, lambda)
{
  accept <- 0 # Will track the acceptance
  try <- 0 # Will track the number of proposals
  x <- 0:m
  
  
  #calculating the pmf of tpois
  foo = -lambda + x*log(lambda) - lfactorial(x) 
  foo = exp(foo)
  c = 1/sum(foo)  #calculate the upper bound experrsio  with poisson(usual poisson)with same lambda
  #as proporsal) 
  tpmf = c*foo
  
  
  while(accept == 0)
  {
    try <- try + 1
    
    U <- runif(1) 
    prop <- rpois(1, lambda = lambda) #draw proposal from pois 
    
    p_j = ifelse(prop<= m ,tpmf[prop] , 0 ) #proposal value can be greater than m 
    q_j = ppois(prop,lambda = lambda)- ppois(prop-1 , lambda = lambda)
    
    ratio = (p_j/(c*q_j))
    
    
    if(U < ratio)
    {
      accept <- 1
      rtn <- prop
    }
  }
  return(c(rtn, try))
}
draw_tpois(m = 30, lambda  = 20)




#simulating 1000 times 
# sample size
N <- 1e3
samp <- numeric(N)
n.try <- numeric(N)
for(t in 1:N)
{
  # I use as a dummy variable often
  foo <- draw_tpois(m = 30, lambda = 30)
  samp[t] <- foo[1]
  n.try[t] <- foo[2]
}
mean(samp)
mean(n.try)  #mean is 2 tries ! Isn't that a fair thing to do ?  

## trying to see if it makes sense or not  
m = 30 
lambda = 20
x <- 0:m


#calculating the pmf of tpois
foo = -lambda + x*log(lambda) - lfactorial(x) 
foo = exp(foo)
c = 1/sum(foo)  #calculate the upper bound experrsio  with poisson(usual poisson)with same lambda
#as proporsal) 
tpmf = c*foo
plot(x , tpmf) 
points(x ,dpois(x = x , lambda =20))  #it closely follows trancated pois! 





##Q.4 

find_c = function(p , lambda , m)
{
  x = 1:m
  pmf_geom = p*(1-p)^x 
  lpmf_pois = -lambda + x * log(lambda) - lfactorial(x)
  pmf_pois = exp(lpmf_pois)
  
  ratio  = pmf_geom / pmf_pois
  c = max(ratio)
  return(c)
}


y = sapply(10:30 , function(i) find_c(0.5,0.5,i))

plot(10:30 , y , pch = 16 , main = "value of c ")  #the values shoot up! 
#upper bound is not convergent! 










