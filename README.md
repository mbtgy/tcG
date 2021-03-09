# tcG
A R package for fitting meta-Gaussian models to precipitation data

## Introduction 

This package contains all the functions needed to fit the models presented in [Boutigny et al. (to be published)](link).

The main functions of the package are `tcG.fit`, `extGP.fit` that are used to fit the models and `res.plot` can be used to draw plots. You can check the documentation pages.

Two types of models can be fitted with this package, all aiming at modeling precipitation which we will note $Y$. 

### Meta-Gaussian models

A meta-Gaussian model can be written as

<img src="https://render.githubusercontent.com/render/math?math=%5Clarge%20Y%20%3D%200*%5Cmathbb%7BI%7D_%7BX%3C0%7D%20%2B%20%5Cpsi(X)*%5Cmathbb%7BI%7D_%7BX%5Cgeq0%7D%20%5Ctext%7B%2C%20%20%20%20%20%20%20with%7D%20%5C%3B%20X%5Csim%20%5Cmathcal%7BN%7D(%5Cmu%2C1)%20">

where <img src="https://render.githubusercontent.com/render/math?math=%5Cmathbb%7BI%7D"> is the indicator function equals to 1 if the condition is true and 0 else. <img src="https://render.githubusercontent.com/render/math?math=%5Cmu"> is the parameter that controls the probability of dry measurement.
The transformation <img src="https://render.githubusercontent.com/render/math?math=%5Cpsi"> is called the anamorphosis, and 4 options are available in the package:

- `gp` (for Generalized Pareto): <img src="https://render.githubusercontent.com/render/math?math=%5Clarge%20%5Cpsi(x)%20%3D%20y_m%2B%5Csigma%20x%5E%7B%5Cfrac%7B1%7D%7B%5Calpha%7D%7D%5Cexp%7B%5Cfrac%7B%5Cxi%20x%5E2%7D%7B2%7D%7D%20">
- `power`: <img src="https://render.githubusercontent.com/render/math?math=%5Clarge%20%5Cpsi(x)%20%3D%20y_m%2B%5Csigma%20x%5E%7B%5Cfrac%7B1%7D%7B%5Calpha%7D%7D">
- `quadratic-power`: <img src="https://render.githubusercontent.com/render/math?math=%5Clarge%20%5Cpsi(x)%20%3D%20y_m%2B%5Csigma_1%20x%5E%7B%5Cfrac%7B1%7D%7B%5Calpha%7D%7D%2B%5Csigma_2%20x%5E%7B%5Cfrac%7B2%7D%7B%5Calpha%7D%7D">
- `power-exp`: <img src="https://render.githubusercontent.com/render/math?math=%5Clarge%20%5Cpsi(x)%20%3D%20%5Csigma_2%20(%5Cexp(%5Csigma_1x%5E%7B1%2F%5Calpha%7D)-1)">

<img src="https://render.githubusercontent.com/render/math?math=y_m"> is the minimal value that can be observed. The choice of the anamorphosis will always be controlled by the argument `name`.


### Extended Generalized Pareto model

This model is the one of [Naveau et al. (2016)](https://agupubs.onlinelibrary.wiley.com/doi/pdfdirect/10.1002/2015WR018552). It  models only the positive part of the distribution and can be written as

<img src="https://render.githubusercontent.com/render/math?math=%5Clarge%20Y_%2B%20%3D%20%20y_m%20%2B%20%5Csigma%20H%5E%7B-1%7D_%5Cxi(U%5E%7B1%2F%5Calpha%7D)">

where <img src="https://render.githubusercontent.com/render/math?math=U%20%5Csim%20Unif(0%2C1)">, <img src="https://render.githubusercontent.com/render/math?math=H_%5Cxi"> is the cdf of a GPD and <img src="https://render.githubusercontent.com/render/math?math=y_m"> is the minimal value that can be observed.

## Data

`guip` is the data set provided with the package. It's the 2006-2017 Automn rainfall series at a 6 minutes time step recorded in Guipavas (France) provided by [Météo France](https://donneespubliques.meteofrance.fr).

## Examples

Examples showing how to use the package are shown in the RMarkdown available in `/vignette`.



# Contact

Please feel free to ask any questions at: `marie.boutigny1[at]gmail.com`

