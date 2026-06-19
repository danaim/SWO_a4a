---
toc: false
header-includes:
  - \usepackage{xcolor}
---

# REVIEW OF THE ASSESSMENT FOR ALL INITIATIVE (a4a) STATISTICAL CATCH-AT-AGE FRAMEWORK FOR THE ICCAT SOFTWARE CATALOGUE


**Authors:** \textcolor{red}{Ai, Tserpes etc.  - ask Ai}


## SUMMARY

The Assessment for All Initiative (a4a) is a statistical catch-at-age (SCA) stock assessment framework originally developed by the European Commission Joint Research Centre (JRC) to simplify, standardize, and scale up fish stock assessment. The framework is built on the Fisheries Library for R (FLR) and uses an age-structured population dynamics model with a flexible, submodel-based formulation that allows the analyst to express fishing mortality, catchability, recruitment, variance, and initial population abundance as linear models using standard R formula syntax—including generalized additive models (GAMs) via the `mgcv` package. The model has been applied to more than 200 stock assessments in the Mediterranean Sea through the General Fisheries Commission for the Mediterranean (GFCM) and the Scientific, Technical and Economic Committee for Fisheries (STECF). It supports full uncertainty propagation through maximum likelihood estimation (MLE), MCMC, and simulation, as well as short-term projections and Management Strategy Evaluation (MSE) conditioning. The a4a R-package (`FLa4a`) and its comprehensive handbook are openly available via the FLR r-universe and GitHub repositories. The authors recommend ICCAT consider the registration of the a4a methodology in the ICCAT software catalogue.


**KEYWORDS:** ICCAT Software Catalogue, FLa4a, R-package, FLR, catch-at-age, stock assessment, MSE


---

\newpage

## Introduction

The SCRS in 2026 recommended that the a4a model (assessment for all initiative) will be considered for inclusion in the ICCAT software catalogue in 202X (Anon., 202X). Following this recommendation, this document provides a review of the a4a for the ICCAT software catalogue in 202X.

## The a4a Model

The Assessment for All (a4a) framework is a statistical catch-at-age (SCA) stock assessment model integrated within the Fisheries Library for R (FLR). It provides a flexible environment that enables scientists to construct and diagnose age-structured assessments by utilizing standardized FLR data structures and statistical submodels. This framework is designed to handle varying levels of data complexity while maintaining a rapid operational timeframe.

### Background and development

The Assessment for All Initiative (a4a) was launched by the European Commission Joint Research Centre (JRC) in 2012 with the primary goal of developing, testing, and distributing methods capable of assessing a large number of stocks within operational time frames. The initiative arose from the recognition that the implementation of the EU Data Collection Framework generates large volumes of data—including for stocks not regularly assessed—and that traditional stock assessment tools take too long to set up and run for such a broad data landscape (Jardim et al., 2014).

The a4a framework is built upon the Fisheries Library for R (FLR; Kell et al., 2007), a widely used open-source software framework. By integrating flexible statistical submodels with the established FLR data structures, a4a enables scientists with varying statistical backgrounds to rapidly construct, fit, and diagnose age-structured stock assessments. The framework has been applied in both European and international contexts, and this document reviews the a4a methodology for potential inclusion in the ICCAT software catalogue.

### Population Dynamics and Statistical Framework

The a4a framework implements a statistical catch-at-age (SCA) model. It follows standard equations to track the number of individuals across ages and years:

$$
N_{a+1,y+1} = N_{a,y} \exp(-F_{a,y} - M_{a,y})
$$

where $N_{a,y}$ is the number of individuals at age $a$ in year $y$, $F_{a,y}$ is fishing mortality, and $M_{a,y}$ is natural mortality. The model is fitted to two types of observations: catches $C_{a,y}$ and abundance indices $I_{a,y,s}$, each assumed to be log-normally distributed about the model predictions:

$$
C_{a,y} \sim \text{Normal}(\hat{C}_{a,y},\ \sigma^2_{a,y})
$$

$$
I_{a,y,s} \sim \text{Normal}(\hat{I}_{a,y,s},\ \sigma^2_{a,y,s})
$$

The total log-likelihood combines the catch $(\ell_C)$ and index $(\ell_I)$ components, and optionally a stock–recruitment $(\ell_{SR})$ component:

$$
\ell = \ell_C + \ell_I + \ell_{SR}
$$

### Submodel Structure

A defining feature of the a4a framework is its **submodel formulation**, where each key quantity is specified through flexible R formulas, allowing fishing mortality, catchability, variance, initial population abundance, and recruitment to be parameterized as linear or additive models. Five submodels are in operation:

- **Fishing mortality** (`fmodel`): model for $F_{a,y}$
- **Catchability** (`qmodel`): model for abundance index catchability $Q_{a,y,s}$
- **Recruitment** (`srmodel`): stock–recruitment relationship (Ricker, Beverton-Holt, smooth hockey-stick, or geometric mean)
- **Observation variance** (`vmodel`): models for $\sigma^2_{a,y}$ and $\sigma^2_{a,y,s}$
- **Initial population** (`n1model`): model for $N_{a,y=1}$

Using this approach, log-transformed quantities are parameterized as linear combinations of age and year covariates, including spline smoothers from the `mgcv` package. For example, fishing mortality on the log scale is expressed as:

$$
\log F_{a,y} = \sum_k \beta_k x_{a,y,k}
$$

This reduces the parameter space while allowing structured, interpretable models.

### Multi-Stage Modelling Approach

The a4a framework adopts a **multi-stage modelling approach**. Prior to fitting the stock assessment model, individual growth (via the `a4aGr` class and von Bertalanffy or other models), natural mortality (via the `a4aM` class), and length-to-age conversion (via the `l2a` method using stochastic slicing) are handled independently. This reduces parameter confounding, decreases computational demands, and provides an intuitive pathway to propagate biological uncertainty into the overall stock assessment. Uncertainty in growth and natural mortality parameters is incorporated using multivariate normal distributions, triangle distributions, or statistical copulas.

---

## Main Advantages of a4a

The main advantages of the a4a framework for application to ICCAT stocks are:

1. **Flexibility through formula-based submodels**: all key quantities are specified using standard R formula syntax (including GAMs), enabling a wide range of model structures without requiring hard-coded options.
2. **Uncertainty quantification**: parameter uncertainty is fully propagated via MLE with variance-covariance inversion, MCMC (via ADMB), and simulation methods.
3. **MSE and forecasting ready**: the framework includes projection tools, harvest control rules (HCR), and a dedicated MP fitting mode optimized for use in Management Strategy Evaluations.
4. **Comprehensive diagnostics**: residuals (standardized, Pearson, raw deviances), retrospective analysis (Mohn's rho), hindcast (MASE score), predictive skill plots, and aggregated catch diagnostics are all implemented.
5. **Scalability**: a4a was designed to assess large numbers of stocks rapidly, an important feature for ICCAT given the diversity of migratory species under its mandate.
6. **Integration with FLR**: full interoperability with the established FLR ecosystem, including `FLCore`, `FLBRP`, `FLasher`, and `a4adiags`.

---

## Application to Fisheries Stocks

The a4a framework has been widely applied in European and international contexts:

- **Mediterranean Sea**: More than 200 stock assessments as of 2024, conducted through the GFCM and STECF assessment working groups.
- **ICES stocks**: Used in the Strategic Initiative on Stock Assessment Methods (2013 World Conference on Stock Assessment Methods) and evaluated comparatively against other SCA methods (Cadrin, 2025).
- **Case studies documented in the a4a handbook**: North Sea plaice (*Pleuronectes platessa*, ICES area IV), European hake (*Merluccius merluccius*, GFCM GSAs 1, 5, 6, 7), red mullet (*Mullus barbatus*, GFCM GSAs 1 and 9), and North Sea herring (*Clupea harengus*).
- **ICCAT context**: The a4a framework has been applied to Mediterranean Swordfish assessment in 2020 (Mantopoulou-Palouka et al., 2020).

---

## Software Developments

The ICES, JRC, and the broader FLR community have organized manuals, software, and training resources around a4a.

### FLa4a R-Package

*Developed by Millar, Jardim, Mosqueira, and collaborators*

This package is the primary software component for inclusion in the ICCAT software catalogue.

- **GitHub**: https://github.com/flr/FLa4a  
- **FLR r-universe**: https://flr.r-universe.dev/FLa4a
- **CRAN/FLR install**: `install.packages(c("FLCore","FLa4a","FLBRP","FLasher",
"a4adiags","ggplotFL"), 
repos=c(FLR="https://flr.r-universe.dev", CRAN="https://cloud.r-project.org"))`

The FLa4a package fits statistical catch-at-age models to fisheries catch data and abundance indices (both scientific survey and commercial CPUE). It supports growth modelling, natural mortality modelling, stochastic length-to-age conversion, and projection under harvest control rules.

### a4adiags Diagnostics Package

*Developed by Mosqueira and Winker*

A companion R-package providing extended diagnostics including hindcast evaluation (MASE), following the recommendations of Carvalho et al. (2021) and Kell et al. (2016).

- **GitHub**: https://github.com/flr/a4adiags

---

## Conclusions

The a4a framework offers a flexible, well-documented, and widely applied age-structured stock assessment methodology with full uncertainty propagation and MSE compatibility. Its formula-based submodel structure distinguishes it from other SCA approaches by providing unparalleled flexibility without sacrificing interpretability, making it suitable for the diverse range of species and data situations encountered in ICCAT. The authors recommend ICCAT consider the registration of the a4a methodology (`FLa4a`) in the ICCAT software catalogue.

---

## References

Carvalho, F., Winker, H., Courtney, D., Kapur, M., Kell, L.T., Cardinale, M., Schirripa, M., Kitakado, T., and Methot, R. 2021. A cookbook for using model diagnostics in integrated stock assessments. *Fisheries Research*, 240: 105959.

Jardim, E., Millar, C.P., Mosqueira, I., Scott, F., Osio, G.C., Ferretti, M., Alzorriz, N., and Orio, A. 2015. What if stock assessment is as simple as a linear model? The a4a initiative. *ICES Journal of Marine Science*, 72(1): 232–236.

Kell, L.T., Mosqueira, I., Grosjean, P., Fromentin, J.-M., Garcia, D., Hillary, R., Jardim, E., Mardle, S., Pastoors, M.A., Poos, J.J., Scott, F., and Scott, R.D. 2007. FLR: an open-source framework for the evaluation and development of management strategies. *ICES Journal of Marine Science*, 64(4): 640–646.

Kell, L.T., Kimoto, A., and Kitakado, T. 2016. Evaluation of the prediction skill of stock assessment using hindcasting. *Fisheries Research*, 183: 119–127.

Mantopoulou-Palouka, D., & Tserpes, G. (2020). Assessment of the Mediterranean swordfish stock by means of Assessment for All. Collective Volume of Scientific Papers ICCAT, 77(3), 482–507.

Mosqueira, I. and Winker, H. 2025. a4adiags: Diagnostics for FLa4a stock assessment models. R package. https://github.com/flr/a4adiags

Punt, A.E., Butterworth, D.S., de Moor, C.L., De Oliveira, J.A.A., and Haddon, M. 2016. Management strategy evaluation: best practices. *Fish and Fisheries*, 17(2): 303–334.
