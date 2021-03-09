# tcG
A R package for fitting meta-Gaussian models to precipitation data

## Introduction 

This package contains all the functions needed to fit the models presented in [Boutigny et al. (to be published)](link).

The main functions of the package are `tcG.fit`, `extGP.fit` that are used to fit the models and `res.plot` can be used to draw plots. You can check the documentation pages.

Two types of models can be fitted with this package, all aiming at modeling precipitation which we will note $Y$. 

### Meta-Gaussian models

A meta-Gaussian model can be written as

$$ Y = 0*\mathbb{I}_{X<0} + \psi(X)*\mathbb{I}_{X\geq0} \text{,       with } X\sim \mathcal{N}(\mu,1) $$
where $\mathbb{I}$ is the indicator function equals to 1 if the condition is true and 0 else. $\mu$ is the parameter that controls the probability of dry measurement.
The transformation $\psi$ is called the anamorphosis, and 4 options are available in the package:

- `gp` (for Generalized Pareto): $\psi(x) = y_m+\sigma x^{\frac{1}{\alpha}}\exp{\frac{\xi x^2}{2}} $
- `power`: $\psi(x) = y_m+\sigma x^{\frac{1}{\alpha}}$
- `quadratic-power`: $\psi(x) = y_m+\sigma_1 x^{\frac{1}{\alpha}}+\sigma_2 x^{\frac{2}{\alpha}}$ 
- `power-exp`: $\psi(x) = \sigma_2 (\exp(\sigma_1x^{1/\alpha})-1)$

$y_m$ is the minimal value that can be observed. The choice of the anamorphosis will always be controlled by the argument `name`.


### Extended Generalized Pareto model

This model is the one of [Naveau et al. (2016)](https://agupubs.onlinelibrary.wiley.com/doi/pdfdirect/10.1002/2015WR018552). It  models only the positive part of the distribution and can be written as

$$ Y_+ =  y_m + \sigma H^{-1}_\xi(U^{1/\alpha})$$
where $U \sim Unif(0,1)$, $H_\xi$ is the cdf of a GPD and $y_m$ is the minimal value that can be observed.

## Data

`guip` is the data set provided with the package. It's the 2006-2017 Automn rainfall series at a 6 minutes time step recorded in Guipavas (France) provided by [Météo France](https://donneespubliques.meteofrance.fr).

