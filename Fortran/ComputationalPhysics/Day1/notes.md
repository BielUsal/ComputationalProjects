# Random Number Generation
Fortran provides us with only one distribution:
Rand(0) $~ U(0,1)$
We want to sample from 'arbitrary' pdfs. For that, we will use inverse-CDF sampling, i.e. use the fact that if $X$ is a random variable that we want to sample and $y~U(0,1)$ is the uniform, by usual pdf transform results:
$$x = F^-1_x(y)$$ 
