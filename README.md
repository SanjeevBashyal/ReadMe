## mHM Fortran Code Overview

This folder contains the core Fortran implementation of the mesoscale Hydrologic Model (mHM). The files are organized as process modules (e.g., interception, snow, soil moisture, runoff), orchestration interfaces, I/O utilities, and model configuration. This README summarizes each file, the principal equations implemented, and how modules interact during a time step.

### How the model runs (high-level flow)
- `mhm_driver.f90` launches mHM, parses CLI options, initializes, runs a simulation or optimization, and finalizes.
- `mo_mhm_interface.F90` initializes configuration and data, runs the model loop via `mo_mhm_interface_run.f90`, and handles finalization (restart/output cleanup).
- `mo_mhm_interface_run.f90` drives the domain/time-step loop. For each time step and domain, it:
  - Prepares meteorological forcings (`meteo_handler`), PET correction, and state slices.
  - Calls the pure hydrologic core `mo_mhm.mHM` to compute L1-level fluxes and states.
  - Optionally triggers routing in mRM, writes outputs per output cadence.
- `mo_mhm.f90` is the hydrological core that executes, per grid cell:
  1) canopy interception, 2) snow accumulation/melt, 3) soil moisture and ET reduction, 4) runoff generation (unsat/sat), 5) total runoff aggregation, 6) optional neutrons state.

The process order and options are configured via the process matrix (`processMatrix`) read from `mhm.nml`.

---

### Notation (symbols used in equations)
- Precipitation and canopy:
  - \(P\): precipitation (per time step), \(C\): canopy water storage, \(C_{\max}\): max canopy storage, \(F\): throughfall
  - \(E_p\): potential evapotranspiration, \(E_\text{canopy}\): canopy evaporation
- Snow/temperature:
  - \(T\): air temperature, \(T_\mathrm{th}\): threshold temperature for snow/rain
  - \(d_d\): degree-day factor; \(d_d^\mathrm{incr}\): increment; \(d_d^\mathrm{max}\): maximum; \(d_d^0\): no-precip value
- Soil moisture and ET:
  - \(\theta\): layer soil water storage; \(\theta_\mathrm{sat}\): saturation; \(\theta_\mathrm{fc}\): field capacity; \(\theta_\mathrm{pwp}\): wilting point
  - \(I_k\): infiltration into layer \(k\); \(\beta\): nonlinear exponent for infiltration attenuation
  - \(f_\mathrm{roots}\): fraction of roots; \(c_1\): Jarvis threshold; \(f_{SM}\): soil moisture stress factor
- Reservoirs and runoff:
  - \(S_u\): unsaturated store, \(S_g\): groundwater store, \(S_\mathrm{th}\): threshold for fast interflow
  - \(k_0, k_1, k_2\): recession coefficients (fast, slow, baseflow), \(k_p\): percolation coefficient, \(\alpha\): nonlinearity, \(k_\mathrm{karst}\): karst factor
  - \(q_\mathrm{fast}, q_\mathrm{slow}, q_\mathrm{base}\): runoff components, \(q_D\): direct (sealed) runoff, \(q_T\): total runoff
  - \(f_\mathrm{sealed}\): sealed area fraction
- PET and radiation:
  - \(R_a\): extraterrestrial radiation, \(R_n\): net radiation, \(\Delta\): slope of saturation vapour pressure curve
  - \(\gamma\): psychrometric constant, \(\lambda\): latent heat of vaporization
  - \(r_a\): aerodynamic resistance, \(r_s\): surface resistance, \(\rho\): air density, \(c_p\): heat capacity of air
  - \(e_s-e\): vapour pressure deficit, \(a_{sh}, a_s\): partitioning constants (see PM notes)
- Neutrons:
  - \(N\): neutron counts, \(N_0\): dry neutron counts (site specific)
  - \(a_0,a_1,a_2\): empirical constants (Desilets), \(\bar{\theta}\): depth-weighted grav. soil water
  - \(A_\text{fast}\): angular integral (COSMIC), \(\Lambda_\text{fast}\): attenuation argument

---

### Orchestration and global state
- `mhm_driver.f90`
  - Main entry. Initializes MPI (if compiled), parses command line (`mo_mhm_cli`), calls `mhm_interface_init`, then `mhm_interface_run` or `mhm_interface_run_optimization`, and finally `mhm_interface_finalize`.
- `mo_mhm_interface.F90`
  - Provides the public API: `mhm_interface_init/run/finalize`, optimization wrapper, and parameter getters.
  - Reads and checks namelists, configures coupling and routing, reads input data, initializes domains/states, and orchestrates evaluation or optimization calls.
- `mo_mhm_interface_run.f90`
  - Encapsulates the simulation loop. On each step: update meteo forcings, compute corrected PET, call `mo_mhm.mHM`, optionally execute mRM routing, and write outputs. Handles warming period and optional Baseflow Index (BFI) aggregation.
- `mo_global_variables.f90`
  - Declares global variables used in reading/writing/startup and shared run-time arrays for L1 states/fluxes (e.g., `L1_soilMoist`, `L1_total_runoff`) and configuration data (e.g., `meteo_handler`, `evap_coeff`, `neutron_integral_AFast`).

---

### Hydrologic core (process modules)
- `mo_mhm.f90` (hydrologic kernel)
  - Per cell sequence:
    1) `mo_canopy_interc.canopy_interc`: canopy storage update, throughfall, canopy evaporation.
    2) `mo_snow_accum_melt.snow_accum_melt`: degree-day melt, split rain/snow, snowpack update, effective precipitation.
    3) `mo_soil_moisture.soil_moisture`: multi-layer infiltration, ET reduction (Feddes/Jarvis variants), sealed fraction processes.
    4) `mo_runoff.runoff_unsat_zone`: fast and slow interflow, percolation, karst loss; update unsat/sat storages.
    5) `mo_runoff.runoff_sat_zone`: baseflow from groundwater storage.
    6) `mo_runoff.L1_total_runoff`: mosaic combination of permeable and sealed fractions.
    7) Optional neutrons: `mo_neutrons.DesiletsN0` or `mo_neutrons.COSMIC`.
  - Sets initial soil moisture to 0.5 of field capacity when states are not provided at the first step.

- `mo_canopy_interc.f90`
  - Throughfall F and canopy evaporation E:
    - Canopy update and throughfall: F = max(P + C − C_max, 0), with new canopy storage C' = min(P + C, C_max).
    - Canopy evaporation: E = PET · (C/C_max)^(2/3).
  - Outputs: `throughfall`, updated `interc` (canopy storage), `aet_canopy`.

- `mo_snow_accum_melt.f90`
  - Split throughfall into rain/snow using temperature threshold T_thresh.
  - Degree-day factor:
    - If P ≤ (DD_max − DD_0)/DD_incr: DD = DD_0 + DD_incr · P, else DD = DD_max.
  - Melt when T > T_thresh: melt = min(DD · (T − T_thresh), snowpack).
  - Effective precipitation: `prec_effect = melt + rain`; updates `snow_pack`, `melt`, `rain`, `snow`.

- `mo_soil_moisture.f90`
  - Impervious (sealed) fraction:
    - Storage update with threshold runoff; sealed ET scales with relative storage: aET_sealed = (PET/evap_coeff − aET_canopy) · (storage_sealed / water_thresh_sealed), bounded to ≥ 0.
  - Multi-layer soil (pervious):
    - Layer-wise infiltration from above: I_k = I_(k−1) · (θ_k/θ_sat_k)^β_k (implemented via exp/log for stability).
    - Storage update: θ_t = θ_(t−1) + I − ET.
    - ET demand allocated top-down; reduced by soil moisture stress:
      - Feddes (processCase 1 or 4): f = f_roots if θ ≥ θ_fc; f = f_roots · (θ − θ_pwp)/(θ_fc − θ_pwp) if θ_pwp < θ < θ_fc; else 0.
      - Jarvis (processCase 2 or 3): θ_norm = (θ − θ_pwp)/(θ_sat − θ_pwp); f = f_roots if θ_norm ≥ c1; else f = f_roots · θ_norm / c1.
    - Outputs: per-layer `aet`, `infiltration`, updated `soil_moist`, `aet_sealed`, `runoff_sealed`, `storage_sealed`.

- `mo_runoff.f90`
  - Unsaturated zone runoff (`runoff_unsat_zone`):
    - Fast interflow (threshold behavior): q_fast = max(k0·(S_u − S_thresh), 0) bounded by storage.
    - Slow interflow: q_slow = min(k1 · S_u^(1+α), S_u − ε).
    - Percolation: perc = k_p · S_u; groundwater storage increases by `karst_loss`-weighted percolation.
  - Saturated zone runoff (`runoff_sat_zone`):
    - Baseflow: q_base = k2 · S_gw; storage reduced accordingly.
  - Total runoff (`L1_total_runoff`):
    - q_T = (q_fast + q_slow + q_base) · (1 − f_sealed) + q_direct · f_sealed.

- `mo_pet.f90`
  - Hargreaves-Samani:
    - PET = c · R_a · (T_avg + b) · sqrt(T_max − T_min), with extraterrestrial radiation R_a from Duffie & Beckman formula.
  - Priestley-Taylor:
    - PET = α · Δ/(γ + Δ) · R_n · (86400/λ), with slope Δ from Tetens relation.
  - Penman-Monteith (transpiration form with a_sh/a_s terms):
    - PET = (Δ·R_n + ρ·c_p·(e_s − e)·a_sh/r_a) / (Δ + γ·a_sh/a_s·(1 + r_s/r_a)) · (86400/λ).
  - Supporting functions: `extraterr_rad_approx`, `slope_satpressure`, `sat_vap_pressure` (Tetens equation).

- `mo_neutrons.f90` (optional nested model)
  - Desilets relation (semi-empirical): neutrons = N0 · (a1 + a0 / (θ̄ + a2)), with θ̄ depth-weighted soil water content (uses horizon depths, bulk density, lattice water).
  - COSMIC (simplified neutron transport): integrates layer contributions using pretabulated integral `neutron_integral_AFast` to approximate angular attenuation.
  - Used when `processMatrix(10,1)` indicates neutron modeling mode.

- `mo_mhm_bfi.f90`
  - Calculates Baseflow Index (BFI) for domains from observed discharge using Eckhardt filter; can be used during optimization or evaluation.

---

### Initialization, I/O, utilities
- `mo_startup.f90`
  - Contains startup routines for allocating variables and preparing domain-level data (invoked from `mhm_interface_init`).
- `mo_init_states.f90`
  - Sets default/initial states and fluxes when not read from restart files.
- `mo_restart.f90`
  - Reading/writing model restart states for warm starts and continuity.
- `mo_write_fluxes_states.f90`
  - NetCDF (or similar) dataset creation and incremental updates for hydrology outputs at L1; used by `mo_mhm_interface_run`.
- `mo_write_ascii.f90`
  - Writes out configuration snapshots and optimization results as ASCII files.
- `mo_file.F90`
  - Paths to namelists and file descriptors used during initialization.
- `mo_mhm_read_config.f90`
  - Reads mHM-specific configuration (complements common config readers).
- `mo_read_optional_data.f90`
  - Loads optional observation datasets (soil moisture, ET, TWS, neutrons) used for optimization/objectives.
- `mo_objective_function.F90`
  - Objective function definitions for calibration/optimization (runoff-only and multi-variable cases).
- `mo_mhm_eval.f90`
  - The evaluation wrapper used by the interface to execute a parameterized run (called by run/optimization).
- `mo_mhm_cli.f90`
  - Command line parsing for mHM driver.
- `mo_mhm_constants.f90`
  - Model constants (physical/empirical) used across PET, neutrons, and other process equations.
- `mo_mhm_messages.F90`
  - User messages/logging helpers (startup, directory checks, finalization).
- `mo_write_ascii.f90`, `mo_file.F90`
  - ASCII output helpers and file path definitions.

---

### Interconnections during one hydrologic time step
Within `mo_mhm_interface_run.mhm_interface_run_do_time_step`:
1) `meteo_handler` provides corrected PET, temperature, and precipitation for L1 cells.
2) `mo_mhm.mHM` is called with:
   - Configuration: `processMatrix`, time, horizons.
   - Forcings: PET/Temp/Prec per cell.
   - Effective parameters (from MPR): α, k0, k1, k2, k_p, θ_sat, θ_FC, θ_exponent, wilting point, etc.
   - Previous states: canopy storage, snowpack, soil moisture, unsat/sat storages.
   - Outputs (in/out): fluxes (ET components, infiltration, interflows, percolation, baseflow, total runoff), and updated states.
3) If routing is active, mRM is invoked to route `L1_total_runoff` to gauges/outlets at L11 scale.
4) Writing of outputs is done per cadence; restart files can be written at finalize.

---

### Model parameters: sources, mapping, and usage per process
This section explains how parameters are provided to the model, how they are transformed into effective fields, and where they are used in each process module.

Sources and mapping
- Global parameters (γ): Read from namelists (e.g., `mhm_parameters.nml`) and accessible via `global_parameters(:,3)`. They control process options and regionalization (MPR) and are the calibration targets for optimization.
- MPR effective parameters: The Multiscale Parameter Regionalization (MPR) maps global parameters and physiographic data to effective spatial fields used at runtime. After `mpr_eval`, the effective parameters are stored mostly as `L1_*` arrays (and `L11_*` for routing). Examples:
  - Soil and ET: `L1_soilMoistSat` (\(\theta_\mathrm{sat}\)), `L1_soilMoistFC` (\(\theta_\mathrm{fc}\)), `L1_soilMoistExp` (\(\beta\)), `L1_wiltingPoint` (\(\theta_\mathrm{pwp}\)), `L1_fRoots` (\(f_\mathrm{roots}\)), `L1_alpha` (\(\alpha\)).
  - Runoff reservoirs: `L1_kFastFlow` (\(k_0\)), `L1_kSlowFlow` (\(k_1\)), `L1_kBaseFlow` (\(k_2\)), `L1_kPerco` (\(k_p\)), `L1_karstLoss` (\(k_\mathrm{karst}\)), thresholds `L1_unsatThresh` (\(S_\mathrm{th}\)), `L1_sealedThresh` (\(W_\mathrm{sealed}\)).
  - PET: `L1_HarSamCoeff` (Hargreaves \(c\)), `L1_PrieTayAlpha` (Priestley–Taylor \(\alpha\)), `L1_aeroResist` (\(r_a\)), `L1_surfResist` (\(r_s\)), `L1_petLAIcorFactor`, aspect factor `L1_fAsp`.
  - Land cover/sealing: `L1_fSealed`, `L1_maxInter` (canopy interception capacity).
  - Neutrons: `L1_No_Count` (N0), `L1_bulkDens`, `L1_latticeWater`, `L1_COSMICL3`.
- Forcings and weights: Provided by `meteo_handler` per time step:
  - `L1_pet_calc`, `L1_temp_calc`, `L1_prec_calc`; monthly `evap_coeff` for free-water surfaces; PET/Temp/Prec weighting fields may be applied depending on configuration.
- Configuration/process selection: `processMatrix` selects process variants and indexes into parameter subsets where applicable (e.g., PET option, soil moisture stress function, routing mode and parameters for mRM).

Parameter usage by process
- Canopy interception (`mo_canopy_interc.canopy_interc`)
  - Inputs: `pet_calc(k)`, `interc_max(k)` = `L1_maxInter`, `prec_calc(k)`.
  - Role: `interc_max` bounds canopy storage; PET enters canopy evaporation formula E = PET · (C/C_max)^(2/3).

- Snow accumulation and melt (`mo_snow_accum_melt.snow_accum_melt`)
  - Inputs (per cell): `deg_day_incr` = `L1_degDayInc`, `deg_day_max` = `L1_degDayMax`, `deg_day_noprec` = `L1_degDayNoPre`, `temperature_thresh` = `L1_tempThresh`, `thrfall` from interception, `temperature` = `L1_temp_calc`, `prec` = `L1_prec_calc`.
  - Role: Degree-day factor depends on precipitation via (DD_0 + DD_incr·P) capped by `deg_day_max`. `temp_thresh` splits rain/snow and triggers melt when exceeded.

- Soil moisture and ET reduction (`mo_soil_moisture.soil_moisture`)
  - Impervious sub-process:
    - Inputs: `frac_sealed` = `L1_fSealed` (\(f_\mathrm{sealed}\)), `water_thresh_sealed` = `L1_sealedThresh` (\(W_\mathrm{sealed}\)), `pet` = `L1_pet_calc` (\(\mathrm{PET}\)), monthly `evap_coeff` (\(k_\mathrm{evap}\)), `aet_canopy` from interception (\(E_\text{canopy}\)).
    - Role: Threshold storage partitions direct runoff vs storage; ET over sealed area scales with storage and `evap_coeff` (water surfaces correction).
  - Pervious multi-layer sub-process:
    - Per-layer inputs: `soil_moist_sat` = `L1_soilMoistSat(:,gg,lc)` (\(\theta_\mathrm{sat}\)), `soil_moist_FC` = `L1_soilMoistFC` (\(\theta_\mathrm{fc}\)), `soil_moist_exponen` = `L1_soilMoistExp` (\(\beta\)), `wilting_point` = `L1_wiltingPoint` (\(\theta_\mathrm{pwp}\)), `frac_roots` = `L1_fRoots` (\(f_\mathrm{roots}\)), `jarvis_thresh_c1` = `L1_jarvis_thresh_c1` (\(c_1\)), PET and canopy ET to set residual PET demand.
    - Role: Infiltration attenuation exponent \(\beta\) = `L1_soilMoistExp`; storage caps by \(\theta_\mathrm{sat}\). ET reduction depends on process case:
      - Feddes: thresholds \(\theta_\mathrm{pwp}\), \(\theta_\mathrm{fc}\) and \(f_\mathrm{roots}\) define linear reduction.
      - Jarvis: normalized \(\theta\) uses \(\theta_\mathrm{sat}\), \(\theta_\mathrm{pwp}\), with threshold \(c_1\) and scaling by \(f_\mathrm{roots}\).

- Runoff generation (`mo_runoff.runoff_unsat_zone` and `runoff_sat_zone`)
  - Unsaturated reservoir:
    - Inputs: `k0` = `L1_kFastFlow` (\(k_0\), upper outlet), `k1` = `L1_kSlowFlow` (\(k_1\), lower outlet), `alpha` = `L1_alpha` (\(\alpha\)), `kp` = `L1_kPerco` (\(k_p\)), `karst_loss` = `L1_karstLoss` (\(k_\mathrm{karst}\)), `unsat_thresh` = `L1_unsatThresh` (\(S_\mathrm{th}\)), and input from last soil horizon `infiltration(:,nHorizons)`.
    - Role: Compute fast (thresholded) and slow (nonlinear via \(\alpha\)) interflow, then percolation via \(k_p\); groundwater gain scaled by \(k_\mathrm{karst}\).
  - Saturated reservoir:
    - Inputs: `k2` = `L1_kBaseFlow` (\(k_2\)).
    - Role: Baseflow \(= k_2 \cdot S_g\); storage reduced accordingly.
  - Total runoff (`L1_total_runoff`):
    - Inputs: `fSealed` = `L1_fSealed` (\(f_\mathrm{sealed}\)), `runoff_sealed` (\(q_D\)), `fast_interflow` (\(q_\mathrm{fast}\)), `slow_interflow` (\(q_\mathrm{slow}\)), `baseflow` (\(q_\mathrm{base}\)).
    - Role: Mosaic combination of pervious and impervious contributions.

- PET calculation (`mo_pet.f90`) — if PET is computed (not provided as input)
  - Hargreaves-Samani (processMatrix(5,1)=1):
    - Inputs: `HarSamCoeff` = `L1_HarSamCoeff` (\(c\)), `tavg` `tmax` `tmin` (from meteo), latitude (grid), constant offset \(b\) from `mo_mhm_constants`. `extraterr_rad_approx` uses day-of-year and latitude to compute \(R_a\).
  - Priestley-Taylor (processMatrix(5,1)=2):
    - Inputs: `PrieTayAlpha` = `L1_PrieTayAlpha` (\(\alpha\)), `Rn` (net radiation from forcings), `tavg`. Uses `slope_satpressure(tavg)` (Tetens) and constants (\(\gamma\), \(\lambda\)).
  - Penman-Monteith (processMatrix(5,1)=3):
    - Inputs: `netrad` ( \(R_n\) ), `tavg`, vapour pressure and aerodynamic resistance from forcings, `aeroResist` = `L1_aeroResist` (\(r_a\)), `surfResist` = `L1_surfResist` (\(r_s\)), leaf partitioning \(a_{sh}\), \(a_s\) (fixed assumptions). All physical constants from `mo_constants`/`mo_mhm_constants`.
  - PET corrections (any PET case):
    - LAI-based correction factor `L1_petLAIcorFactor` and aspect factor `L1_fAsp` applied by `meteo_handler%get_corrected_pet` to produce `L1_pet_calc`.

- Neutrons (optional; `mo_neutrons.f90`)
  - Desilets:
    - Inputs: soil moisture profile (from `L1_soilMoist`), horizon depths (`HorizonDepth_mHM`), `bulkDens` = `L1_bulkDens`, `latticeWater` = `L1_latticeWater`, `N0` = `L1_No_Count` (\(N_0\)).
    - Role: Depth-weighted θ̄ via exponential kernel (~D86) and semi-empirical neutron relation with constants from `mo_mhm_constants`.
  - COSMIC:
    - Inputs: soil moisture profile, horizon depths, `neutron_integral_AFast` (pretabulated, for \(A_\text{fast}\)), interception and snowpack (surface layer water), `L1_No_Count` (\(N_0\)), `L1_bulkDens`, `L1_latticeWater`, `L1_COSMICL3`.
    - Role: Layer-wise attenuation and upward fast neutron flux integration; tabular integral avoids costly per-step integration.

Where parameters are set in code
- After initialization (`mhm_interface_init`), `mo_mpr_eval.mpr_eval` computes effective parameters at L1 (and L11 if routing); they are stored in `mo_mpr_global_variables` arrays with names `L1_*`, then referenced in `mo_mhm_interface_run` when calling `mHM`.
- PET corrections call in `mhm_interface_run_do_time_step` passes:
  - `petLAIcorFactorL1 = L1_petLAIcorFactor(:, iLAI, yId)`
  - `fAsp = L1_fAsp(:,1,1)`, `HarSamCoeff`, `PrieTayAlpha`, `aeroResist`, `surfResist` with appropriate time indices.
- The hydrologic kernel call (`mo_mhm.mHM`) wires each `L1_*` effective parameter into the matching argument (e.g., `alpha => L1_alpha`, `k0 => L1_kFastFlow`, thresholds, soil properties, etc.).
- The `processMatrix` determines which process variants are active (e.g., soil moisture stress case, PET option, routing mode) and which subset of the global parameter vector is used for routing calibration.

---

### Process equations (LaTeX)

Canopy interception and evaporation
- Canopy storage update and throughfall for precipitation \(P\), canopy storage \(C\), and maximum storage \(C_{\max}\):
\[
F = \max(P + C - C_{\max},\ 0), \quad C' = \min(P + C,\ C_{\max})
\]
- Canopy evaporation with potential evapotranspiration \(E_p\):
\[
E_\text{canopy} = E_p \left(\frac{C}{C_{\max}}\right)^{2/3}
\]

Snow accumulation and melt (degree-day)
- Split throughfall into snow/rain using threshold temperature \(T_\mathrm{th}\):
\[
\text{snow} = \begin{cases}
F, & T \le T_\mathrm{th} \\[2pt]
0, & T > T_\mathrm{th}
\end{cases}
\quad
\text{rain} = \begin{cases}
0, & T \le T_\mathrm{th} \\[2pt]
F, & T > T_\mathrm{th}
\end{cases}
\]
- Degree-day factor with precipitation \(P\), increment \(d_d^\mathrm{incr}\), max \(d_d^\mathrm{max}\), and no-precip value \(d_d^0\):
\[
d_d = \min\!\left(d_d^\mathrm{max},\ d_d^0 + d_d^\mathrm{incr}\, P\right)
\]
- Melt and snowpack update; effective precipitation:
\[
\text{melt} = \min\!\left(d_d \,(T - T_\mathrm{th})_+,\ \text{snowpack}\right), \quad
\text{prec}_\mathrm{eff} = \text{melt} + \text{rain}
\]

Soil moisture: pervious layers (infiltration and ET)
- Infiltration cascade from layer \(k-1\) to \(k\) with storage \(\theta_k\), saturation \(\theta_{k,\mathrm{sat}}\), exponent \(\beta_k\), input \(I_{k-1}\):
\[
I_k = I_{k-1}\left(\frac{\theta_k}{\theta_{k,\mathrm{sat}}}\right)^{\beta_k}
\]
- Storage update in a layer (omitting indices \((k,t)\) for brevity):
\[
\theta_t = \theta_{t-1} + I - \mathrm{ET}
\]
- ET demand allocation: from top to bottom, residual PET after canopy is successively satisfied.
- Soil moisture stress factor \(f_{SM}\) (two options):
  - Feddes (with field capacity \(\theta_\mathrm{fc}\), wilting point \(\theta_\mathrm{pwp}\), root fraction \(f_\mathrm{roots}\)):
\[
f_{SM} =
\begin{cases}
f_\mathrm{roots}, & \theta \ge \theta_\mathrm{fc} \\[4pt]
f_\mathrm{roots}\dfrac{\theta - \theta_\mathrm{pwp}}{\theta_\mathrm{fc} - \theta_\mathrm{pwp}}, & \theta_\mathrm{pwp} < \theta < \theta_\mathrm{fc} \\[10pt]
0, & \theta \le \theta_\mathrm{pwp}
\end{cases}
\]
  - Jarvis (with saturation \(\theta_\mathrm{sat}\), \(c_1\) threshold):
\[
\theta_\mathrm{norm} = \frac{\theta - \theta_\mathrm{pwp}}{\theta_\mathrm{sat} - \theta_\mathrm{pwp}}, \quad
f_{SM} =
\begin{cases}
f_\mathrm{roots}, & \theta_\mathrm{norm} \ge c_1 \\[4pt]
f_\mathrm{roots}\dfrac{\theta_\mathrm{norm}}{c_1}, & \theta_\mathrm{norm} < c_1
\end{cases}
\]
- Actual ET in layer:
\[
\mathrm{ET} = f_{SM}\cdot \mathrm{PET}_\text{residual}
\]

Impervious (sealed) fraction
- Storage thresholding with threshold \(W_\mathrm{sealed}\) and storage \(S_\mathrm{sealed}\):
\[
Q_\mathrm{direct} = \max\!\left(S_\mathrm{sealed} + \text{prec}_\mathrm{eff} - W_\mathrm{sealed},\ 0\right),\quad
S_\mathrm{sealed}' = \min\!\left(S_\mathrm{sealed} + \text{prec}_\mathrm{eff},\ W_\mathrm{sealed}\right)
\]
- ET over sealed surfaces with monthly coefficient \(k_\mathrm{evap}\) and canopy ET \(E_\mathrm{canopy}\):
\[
E_\mathrm{sealed} = \max\!\left( \left(\frac{\mathrm{PET}}{k_\mathrm{evap}} - E_\mathrm{canopy}\right)\frac{S_\mathrm{sealed}}{W_\mathrm{sealed}},\ 0 \right)
\]

Runoff generation (reservoirs)
- Unsaturated store \(S_u\) with threshold \(S_\mathrm{th}\), fast outlet \(k_0\), slow outlet \(k_1\) with nonlinearity \(\alpha\), percolation \(k_p\), karst factor \(k_\mathrm{karst}\):
\[
q_\mathrm{fast} = \max\!\big(k_0(S_u - S_\mathrm{th}),\, 0\big),\quad
S_u \leftarrow S_u - q_\mathrm{fast}
\]
\[
q_\mathrm{slow} = \min\!\big(k_1\, S_u^{\,1+\alpha},\, S_u\big),\quad
S_u \leftarrow S_u - q_\mathrm{slow}
\]
\[
\text{perc} = k_p\, S_u,\quad
S_u \leftarrow S_u - \text{perc},\quad
S_g \leftarrow S_g + k_\mathrm{karst}\,\text{perc}
\]
- Saturated store \(S_g\) (baseflow with \(k_2\)):
\[
q_\mathrm{base} = k_2\, S_g,\quad S_g \leftarrow S_g - q_\mathrm{base}
\]
- Total runoff at L1 with sealed area fraction \(f_\mathrm{sealed}\) and direct runoff \(q_D\):
\[
q_T = (q_\mathrm{fast} + q_\mathrm{slow} + q_\mathrm{base})\,(1 - f_\mathrm{sealed}) + q_D\, f_\mathrm{sealed}
\]

PET formulations
- Hargreaves–Samani (with extraterrestrial radiation \(R_a\), constants \(c\), \(b\)):
\[
\mathrm{PET} = c \, R_a \, (T_\mathrm{avg} + b)\, \sqrt{T_\mathrm{max} - T_\mathrm{min}}
\]
- Priestley–Taylor (with net radiation \(R_n\), slope \(\Delta\), psychrometric constant \(\gamma\), latent heat \(\lambda\)):
\[
\mathrm{PET} = \alpha \frac{\Delta}{\gamma + \Delta}\, \frac{R_n\, 86400}{\lambda}
\]
- Penman–Monteith (transpiration form; aerodynamic resistance \(r_a\), surface resistance \(r_s\), air density \(\rho\), heat capacity \(c_p\), VPD \(e_s-e\), partition \(a_{sh}, a_s\)):
\[
\mathrm{PET} = \frac{86400}{\lambda}\,
\frac{\Delta R_n + \rho c_p (e_s - e)\, a_{sh}/r_a}{\Delta + \gamma \, a_{sh}/a_s \left(1 + r_s/r_a\right)}
\]

Neutrons (optional)
- Desilets relation with depth-weighted gravimetric soil water \(\bar{\theta}\) and site-specific constants \(a_0,a_1,a_2\) and \(N_0\):
\[
N = N_0 \left(a_1 + \frac{a_0}{\bar{\theta} + a_2}\right)
\]
- COSMIC (concept): upward fast neutron flux attenuated by integrated mass of soil/soil-water columns; evaluated via a pretabulated angular integral \(A_\text{fast}(\Lambda_\text{fast})\).

---

### Execution flow and loops (domains, timesteps, cells)

Initialization
1) `mhm_driver` calls `mhm_interface_init`, which:
   - Reads namelists, sets coupling, reads common and mHM configs, prepares the process matrix.
   - Initializes domain metadata (`domainMeta`), levels (L0/L1/L11 where applicable), and timers.
   - Runs `mhm_initialize` to allocate/prepare variables; `meteo_handler%init_level2`.
   - If routing is on, initializes mRM configuration.
   - Optionally reads optimization datasets (SM/ET/TWS/neutrons) based on objective function selection.

Run preparation
2) `mhm_interface_run_prepare`:
   - Stores parameter set (global or provided), selects domain indices (all or optimization subset).
   - Reads restart states if configured; otherwise calls `variables_default_init` and `mpr_eval` to derive effective parameters (`L1_*`, `L11_*`). For routing, initializes routing state if needed.

Domain loop and time loop
3) For each selected domain:
   - `mhm_interface_run_prepare_domain` fixes slices/indices (`s1:e1`, masks), prepares routing mapping (`s11:e11`) if needed, and initializes datetime.
   - Time loop: for `tt = 1 .. nTimeSteps(domain)`
     - Update LAI timestep and meteo (`meteo_handler%update_timestep`).
     - Compute corrected PET (`get_corrected_pet`) using effective fields (`L1_petLAIcorFactor`, `L1_fAsp`, `L1_HarSamCoeff`, `L1_PrieTayAlpha`, `L1_aeroResist`, `L1_surfResist`) and latitude/weights; read temperature and precipitation for current step.
     - Call hydrologic core `mHM(...)` with:
       - Configuration: `processMatrix`, conversion factors, horizon depths, counts.
       - Forcings: `L1_pet_calc`, `L1_temp_calc`, `L1_prec_calc`.
       - States (INOUT): `L1_inter`, `L1_snowPack`, `L1_sealSTW`, `L1_soilMoist(:, :)`, `L1_unsatSTW`, `L1_satSTW`, `L1_neutrons`.
       - Fluxes (INOUT): `L1_aET*`, `L1_infilSoil`, `L1_fastRunoff`, `L1_slowRunoff`, `L1_percol`, `L1_baseflow`, `L1_runoffSeal`, `L1_total_runoff`, `L1_preEffect`, `L1_melt`, `L1_rain`, `L1_snow`, `L1_Throughfall`.
       - Effective parameters: `L1_alpha`, `L1_kFastFlow`, `L1_kSlowFlow`, `L1_kBaseFlow`, `L1_kPerco`, `L1_karstLoss`, `L1_unsatThresh`, `L1_sealedThresh`, soil properties (`L1_soilMoistSat/FC/Exp`, `L1_wiltingPoint`, `L1_fRoots`, `L1_jarvis_thresh_c1`), neutron parameters (if enabled).
     - Inside `mHM`, an OpenMP-parallel cell loop runs for `k = 1 .. nCells` executing the process modules in order: interception → snow → soil moisture/ET → runoff reservoirs → total runoff → optional neutrons.
     - If routing is active:
       - Determine routing time-step factor relative to hydrologic time step (adaptive vs. fixed Muskingum).
       - Aggregate or disaggregate runoff to routing cadence, add gauge inflows, and call `mRM_routing` with L11 mapping and parameters; optionally compute river temperature coupling energy terms; manage sub-stepping accumulation and reset as needed.
     - Update optimization accumulators (SM/ET/TWS/neutrons) respecting warming period and aggregation cadence.
     - Write outputs when `timeStep_model_outputs` cadence triggers, creating/updating datasets; similarly for mRM outputs and monthly GW head if coupled.
     - Increment datetime and update year-dependent indices (e.g., land cover id).
   - `mhm_interface_run_finalize_domain` clears routing buffers (and river temperature lateral fields if active).

Finalize
4) `mhm_interface_run_finalize` optionally returns routed runoff timeseries and BFI; deallocates transient arrays.
5) `mhm_interface_finalize` optionally writes restart files (if enabled and not optimizing), triggers mRM final output, prints finish messages, and deallocates global variables.

Parallelization and MPI
- The domain/time-step loops may run under MPI with master/subprocess coordination; reading and some operations occur only on domain masters or specific ranks, as guarded by `rank` and `domainMeta%isMasterInComLocal` checks.
- The per-cell hydrologic computations are OpenMP-parallelized in `mo_mhm.mHM` for computational efficiency.

---

### Notes for readers and students
- The process options are controlled via the `processMatrix`, e.g. PET method (Hargreaves/Priestley-Taylor/Penman-Monteith) and soil moisture stress function (Feddes/Jarvis).
- Many empirical constants (Tetens, Duffie-Beckman radiation, Penman-Monteith parameters) are consolidated in `mo_mhm_constants.f90` and `mo_constants` (external/common).
- mHM relies on pre-processing and parameter regionalization (MPR) to derive spatially distributed effective parameters used at runtime; this happens before calling the time step hydrological kernel.

---

### File index (quick reference)
- Model driver and interfaces: `mhm_driver.f90`, `mo_mhm_interface.F90`, `mo_mhm_interface_run.f90`, `mo_mhm.f90`
- Core processes: `mo_canopy_interc.f90`, `mo_snow_accum_melt.f90`, `mo_soil_moisture.f90`, `mo_runoff.f90`, `mo_pet.f90`, `mo_neutrons.f90`, `mo_mhm_bfi.f90`
- Initialization/States/Restart/Output: `mo_startup.f90`, `mo_init_states.f90`, `mo_restart.f90`, `mo_write_fluxes_states.f90`, `mo_write_ascii.f90`
- Config and globals: `mo_mhm_read_config.f90`, `mo_global_variables.f90`, `mo_mhm_constants.f90`, `mo_mhm_messages.F90`, `mo_file.F90`, `mo_objective_function.F90`, `mo_mhm_eval.f90`, `mo_mhm_cli.f90`


