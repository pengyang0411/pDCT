# Robust Bayesian Estimation of Treatment Effects on Longitudinal Outcomes in Partially Decentralized Clinical Trials
## Project overview

This project contains code for simulation studies and model fitting for partially decentralized clinical trials (pDCTs) with longitudinal outcomes. The main goal is to evaluate the performance of the proposed **Digital Twins** method against several comparator models under different data-generating scenarios.

The simulation framework allows for:\
- repeated longitudinal outcomes,\
- a mixture of onsite and offsite visits,\
- treatment-specific offsite bias,\
- baseline covariate effects,\
- estimation of the average treatment effect (ATE) over time.

---

## File structure

### `run_all_scen.R`
Main script used to generate simulation datasets, fit models under each scenario, and save simulation results.

### `source/source.R`
Contains core functions used by the simulation pipeline, including:\
- data generation,\
- tsBART training,\
- Bayesian functional mixed-effects model fitting.

---

## Main functions in `source.R`

### `sim_data()`
Simulate longitudinal pDCT data.

#### Inputs
- `n`: number of subjects\
- `n_T`: number of visits per subject\
- `n_p`: number of baseline covariates\
- `d_s`: magnitude of offsite bias shift\
- `sig_1`: standard deviation for onsite measurements\
- `sig_2`: standard deviation for offsite measurements\
- `if_null`: whether to simulate under the null treatment effect\
- `if_plot`: whether to display simulated data

#### Output
A list containing:\
- `Z`: baseline covariates\
- `Y`: observed outcomes\
- `X`: treatment assignment\
- `S`: visit-level site assignment schedule\
- `Y0`, `Y1`: potential outcomes under control and treatment\
- `phi`: subject-specific random-effect parameters\
- `S_all`, `X_all`, `t_all`, `ID_all`: expanded visit-level variables\

---

---

### `train_tsBART()`
Fits separate tsBART models to generate digital twin predictions.


---

### `run_LDT()`
Fits the Bayesian longitudinal digital twins model or its reduced comparator versions.

#### Key options
- `if_site = TRUE/FALSE`: whether to include site-assignment effects\
- `if_digital_twins = TRUE/FALSE`: whether to include digital twins\
- `if_covariates = TRUE/FALSE`: whether to include baseline covariates\
- `if_model_diagnostic = TRUE/FALSE`: whether to save diagnostic plots\

#### Output
A list containing posterior draws and summaries for:\
- mean trajectory,\
- treatment effect trajectory,\
- treatment-specific offsite bias trajectory,\
- fixed-effect coefficients,\
- residual variances,\
- posterior ATE and posterior probabilities.

---

## Simulation scenarios in `run_all_scen.R`

The script evaluates the following scenarios:

1. **Baseline**
   - No treatment effect
   - No treatment-specific offsite bias
   - Time-invariant nonlinear covariate effects

2. **Linear**
   - No treatment effect
   - Offsite bias follows a quadratic-plus-linear time trend
   - Linear covariate effects

3. **L_Increased_ATE**
   - Linearly increasing treatment effect over time

4. **NL_Increased_ATE**
   - Nonlinearly increasing treatment effect over time

5. **Time_varying**
   - Nonlinear treatment effect
   - Time-varying covariate effects

---

## Methods compared

The following methods are fitted in each simulation scenario:

- `f(t) + X`  
  Functional mixed-effects model with treatment only

- `f(t) + X + S`  
  Functional mixed-effects model with treatment and site assignment

- `f(t) + X + S + Z`  
  Functional mixed-effects model with treatment, site assignment, and baseline covariates

- `Digital Twins`  
  Proposed two-stage method using tsBART-generated digital twins

---

## Output

For each simulation replicate, the script saves:

- posterior mean ATE at each visit,\
- posterior probability of positive treatment effect at each visit,\
- results under each scenario and method,

to:

```r
Res/Sim_ID_<i>.RData
```
