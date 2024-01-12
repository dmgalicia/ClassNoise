# ClassNoise Package

ClassNoise allows modeling, generating and validating data with class corruption. It provides a framework to design class noise models based on provided datasets allowing the generation of synthetic datasets. This package empower the creation of controlled experimental settings to evaluate class noise impact on supervised learning.

## Installation

The ClassNoise installation requires the `install_github` function from the package `remotes`.

``` r
install.packages("remotes")
library(remotes)
install_github("dmgalicia/classnoise")
```

## Suggested workflow

The ClassNoise package suggests a workflow based on its functionalities:

1.  Load a clean and prepocessed dataset, or use the preloaded data (`qb` and `mux11`).

2.  Generate a class noise model (`NCAR`, `NAR` and `NNAR`) with the desired noise distribution.

3.  Use `generateDF` to create a synthetic dataset based on the model distribution.

4.  Validate the generated dataset with information metrics (`inconsistency`, `discF2`, `discF4` and `noiseLVL`).

## Example

The following examples describes a basic workflow using the preloaded dataset `mux11`.

``` r
# Data loading
data(mux11)

# Creation of a NNAR model
model <- NNAR(data = mux11, attrib.set = c("A0", "A2"), 
              noise = c(0.4, 0.3, 0.2, 0.4, 0.1, 0.2, 0.4, 0.3))

# Error distribution
model$parameters$Error

#  Parameters of node Error (multinomial distribution)
#Conditional probability table:
 
#, , A2 = False, Class = False

#       A0
#Error   False True
#  False   0.6  0.8
#  True    0.4  0.2

#, , A2 = True, Class = False

#       A0
#Error   False True
#  False   0.7  0.6
#  True    0.3  0.4

#, , A2 = False, Class = True

#       A0
#Error   False True
#  False   0.9  0.6
#  True    0.1  0.4

#, , A2 = True, Class = True

#       A0
#Error   False True
#  False   0.8  0.7
#  True    0.2  0.3

# Data generation
df <- generateDF(model$parameters, 1000)

# Data validation
discF2(df)
# 1
noiseLVL(df)
# 29.6
```
