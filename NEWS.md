
# qgcompint: quantile g-computation with statistical interaction

# v0.10.2
## Major changes
- added non-linearity for estimating equation methods
- added plotting and bounds for estimating equation methods
- new dependency: MASS

## Minor changes
- Documentation improvements

## Bug fixes
- major, rare: fixed an issue with formula interpretation for complex formulas with nested interaction terms


# v0.9.0
## Major changes
- added new dependencies: numDeriv and rootSolve


## Minor changes
- added a number of generic functions, including `anova` to help make global chi-squared tests of interaction

## Bug fixes
- n/a


# v0.8.0
## Major changes
- Added estimating functions approach `qgcomp.emm.glm.ee` which allows a variety of new approaches (see ?qgcomp::qgcomp.glm.ee for more info)

## Minor changes
- n/a

## Bug fixes
- `qgcomp.emm.glm.boot` had a rare, fatal bug when setting q=NULL


# v0.7.0
## Minor changes
- Added `getjointeffects` function to get mixture effects at values of an interacting variable versus a common referent level of exposure and the interacting variable

## Bug fixes
- n/a

# v0.6.6
## Minor changes
- Documentation changes
- Handling of parallel processing through futures: bringing in line with package recommendations and planned deprecations.
- Reduce help example times

## Bug fixes
- correcting calculations in 'modelbound' function
- correcting starting values in bootstrap fits for binary outcomes when rr = TRUE


# v0.6.2
## Minor changes
- implementing data output for survival models
- Improved vignette
- added examples to `plot` documentation for manipulating plots

## Bug fixes
- Typos in plotting
- Doc unicode scrubbing

# v0.5.0
## Minor changes
- Better error handling for non-linear fits
- Improved vignette


# v0.4.0
## Major changes
- Better error handling in likelihood functions

## Bug fixes
- Fixed bug in log-Jacobian calculation for Gaussian model
