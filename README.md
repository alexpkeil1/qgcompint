### qgcompint: quantile g-computation with effect measure modification 

### Quick start
    devtools::install_github("alexpkeil1/qgcompint")
    library(qgcomp)
    library(qgcompint)
	 set.seed(40)
    dat <- data.frame(y=runif(50),
	                  x1=runif(50),
	                  x2=runif(50),
	                  z=rbinom(50,1,0.5),
	                  r=rbinom(50,1,0.5))
	 
	 # quantile g-computation without effect measure modification
	 qfit <- qgcomp.noboot(f=y ~ z + x1 + x2, 
	           expnms = c('x1', 'x2'), 
	           data=dat, q=2, 
	           family=gaussian())
	 # no output given here          
	 
	 # with effect measure modification by Z
	 (qfitemm <- qgcomp.emm.noboot(f=y ~ z + x1 + x2,
	           emmvar="z", 
	           expnms = c('x1', 'x2'), 
	           data=dat, q=2, 
	           family=gaussian()))

   

    > ## Qgcomp weights/partial effects at z = 0
	> Scaled effect size (positive direction, sum of positive coefficients = 0)
	> None

	> Scaled effect size (negative direction, sum of negative coefficients = -0.278)
	>    x2    x1 
	> 0.662 0.338 

    > ## Qgcomp weights/partial effects at z = 1
	> Scaled effect size (positive direction, sum of positive effects = 0.0028)
	> x1 
	>  1 

	> Scaled effect size (negative direction, sum of negative effects = -0.0128)
	> x2 
	>  1 

	> Mixture slope parameters (Delta method CI):

	> 			  Estimate Std. Error Lower CI Upper CI t value  Pr(>|t|)
	> (Intercept)  0.58062    0.11142  0.36224  0.79900  5.2112 4.787e-06
	> psi1        -0.27807    0.20757 -0.68490  0.12876 -1.3397    0.1872
	> z           -0.10410    0.15683 -0.41148  0.20329 -0.6637    0.5103
	> z:mixture    0.26811    0.26854 -0.25822  0.79444  0.9984    0.3235

	> Estimate (CI), z=1: 
	> -0.0099575 (-0.34389, 0.32398)
	
### Current package capabilities/limitations
- Single modifiers only (e.g. interaction terms between a modifier and the mixture can only be estimated for a single modifier at a time). This also implies that no interaction terms between the modifier and other covariates can be considered.
- linear only specifications (i.e. linear effects of a mixture, but not covariates) for Cox model

### Interpretation
- coefficients
  - psi1 coefficient: effect of the mixture in the referent category of the effect measure modifier (`z` in this case)
  - z coefficient: main effect of the effect measure modifier (will change based on name of effect measure modifier)
  - z:mixture coefficient: interaction term between the effect measure modifier and the entire mixture
- weights/partial effects (only given if effect measure modifier is binary)
  - first set of weights/sum of negative/positive coefficients: 
    	- weights: proportion of positive/negative partial effect in the reference stratum effect measure modifier (here, `z=0`)
  	  - partial effect: sum of main term coefficients for exposure in a given direction (these will sum to psi1, the effect estimate in the reference stratum effect measure modifier (here, `z=0`))
  - second set:
  	  - weights: proportion of positive/negative partial effect in the index stratum effect measure modifier (here, `z=1`)
  	  - partial effect: sum of main term coefficients + interaction term coefficients for exposure in a given direction (these will sum to the mixture effect estimate in the index stratum effect measure modifier (here, `z=1`))
- Additional estimates (only given if effect measure modifier is binary): The effect of the mixture in the index category of the effect measure modifier (here, `z=1`)
