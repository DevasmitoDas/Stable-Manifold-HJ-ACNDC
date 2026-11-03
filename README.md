# Stable Manifold Method for a Hysteretic Soft Pneumatic Actuator (MATLAB)

This MATLAB code reproduces the computations and the figures of

> N. Sakamoto, D. Das, *Saturation-Aware Stable-Manifold Optimal Control of a Hysteretic Soft Pneumatic Actuator*, ACNDC 2027 (`Documents/NS_DD_ACNDC_2027.pdf`).

It follows the structure of the stable-manifold toolbox `Stable-manifold-HJ-master`, which was originally distributed as the NonlinearControlToolbox by N. Sakamoto. The package has the same 8 top-level folders. The work is split into the same five steps, `DirectApp.m` → `P1` … `P5`, and the files use the same prefixes (`P#_`, `SP_`, `i_`, `f_`). The control problem lives in `inputData/5D_Pneumatic_nu_sat_symb/`.

---

## 1. Where the paper's numbers come from (provenance check)

The four supplied folders were compared with the paper's numbers and figures.

| Folder | Used for the paper? | Evidence |
|---|---|---|
| `Sakamoto-Das/Stable-manifold-HJ-master/soft_pneumatic_hj_stable_manifold.m` (Revision 5), with `inputData/5D_Pneumatic_nu_sat_symb/` | **Yes. This code produced every result in the paper.** | Its output folder `soft_pneumatic_hj_nu_saturated_results/` matches the paper. The paper's four PNGs are pixel-identical crops of that folder's PNGs (mean absolute pixel difference 0). Its `numerical_summary.txt` and MAT file match every number in the paper: poles, T_m = 672.5 s, the K₅S⁻¹ gain, Table 1, the fit 7363 → 284, held-out 0.099/0.377/1 with a 30.6 % pass rate, Table 2 (11.340 / 8.979, HJ use 76.6 %) and Table 3 |
| `Sakamoto-Das/soft_pneumatic_hj_stable_manifold.m` (Revision 4) and `soft_pneumatic_hj_nu_results/` | No (predecessor) | Hard-clipped command u ∈ [0,1] with no σ(u) in the drift. Its poles (−6.723 ± 11.770i, …), costs (11.169 / 7.641) and figures differ from the paper |
| `Sakamoto-Das/soft_pneumatic_ladder_handoff_*.m`, `soft_pneumatic_direct_input_*.m` | No | Earlier direct-input studies with other figure sets |
| `ACNDC 2027/soft_pneumatic_bifurcation_control*.m` | **No** | Different model: fractional (Caputo) viscoelasticity, Duffing stiffness, bifurcation maps and a reference governor. It has no Bouc–Wen hysteresis, no tanh saturation and no HJ/Riccati design |
| `Stable-manifold-HJ-master` | Structure only | Supplies the layout, file naming and the P1–P5 pipeline that this package follows |

The present code reorganises Revision 5 into the toolbox structure, and its numerical algorithm is unchanged. MATLAB was not available where this package was assembled. Its algorithm was therefore cross-checked with an independent NumPy re-implementation, and the toolbox-free formulas of `numericProblem.m` were compared with the `matlabFunction` output of Revision 5 (agreement ≈ 1e-13 for every λ_h). The re-implementation reproduced the native Revision-5 result file to 1e-12 or better. The match covers the design matrices, all 18 manifold trajectories (iterations 5–8, anchors to 5e-15), the 55 fit coefficients (1e-13) and the LQR, guarded-HJ, disturbance, relaxation and open-loop trajectories (≤ 2e-13, with no HJ/LQR switching decision different).

---

## 2. Quick start

1. Open MATLAB in this folder (the repository root).
2. Run
   ```matlab
   >> DirectApp
   ```
   The run is non-interactive. Set `askFlag = 1` in `DirectApp.m` to get the prompts of the reference toolbox.
3. Results are written to `result/`:
   - `figures/open_loop_loops.png`, `figures/relaxation_hysteresis.png`, `figures/nominal_physical_states.png`, `figures/nominal_commands_and_hj_selection.png` are the paper figures, with the same file names as in the paper sources. Compare them with `Documents/NS_DD_ACNDC_2027_figures/`.
   - `numerical_summary.txt` lists all the reports next to the values printed in the paper.
   - `closed_loop_metrics.csv`, `relaxation_metrics.csv` and `open_loop_metrics.csv` hold the case metrics.
   - `<date>_..._data.mat` is saved after P4 and `<date>_..._results.mat` at the end, as in the reference toolbox.
   - `symbolic_generated/lambda_*_Rnu_*/` holds the exported symbolic derivation (symbolic route only).

The full run consists of the manifold computation, the fit, 2 nominal runs, 8 disturbance runs, 4 relaxation runs and the open-loop run. Each closed-loop run has 24 000 samples, so the full run takes a while. To produce only the nominal figures quickly, set `cfg.runDisturbanceSweep = false` in `systemData.m`.

### Requirements

| Component | Needed for |
|---|---|
| MATLAB R2016b or newer (implicit expansion) | Everything. Revision 5 ran on R2024b |
| Symbolic Math Toolbox | The symbolic route (`symbFlag = 1`, default and used for the paper). Without it the code falls back to `numericProblem.m`, which uses the same formulas |
| Control System Toolbox | **Not needed.** CARE and Lyapunov are solved by `f_care.m` and `f_lyap.m` |
| YALMIP + SDPT3 (bundled in `LMI/`) | Only the optional LMI certificate (`cfg.lmiCheckFlag = 1`) |

---

## 3. Repository layout

```
Stable-manifold-HJ-SoftPneumatic/
├── DirectApp.m              Main entry point: P1 → P5, saves result/*.mat
├── P1_data_setting.m        systemData → CARE, F, L, T (Diagonalize), Trs, self-checks, Report_P1
├── P2_hamiltonian_check.m   Successive approximation along eigen-directions of F: |H(0)| and residual vs k
├── P3_decide_xi.m           Seeds xi (Decide_xi.m): ±e_i and 8 Halton directions, radius 2e-4
├── P4_iterative_calcu.m     Integral iteration per seed, acceptance tests, anchors (Table 1)
├── P5_solution_check.m      Value fit, guarded HJ vs LQR, experiments, figures, report
├── SP_appNth.m              One successive approximation k → k+1 (eq. (iter))
├── SP_manifold_check.m      Acceptance tests of a trajectory (Table 1)
├── SP_ham_calcu.m, SP_shooting_plot.m       Hamiltonian on anchors / anchors on (y_i, y_j) planes
├── SP_value_fit.m, SP_model_evaluation.m    Collocation set + fit / gradient errors, held-out rho_H
├── SP_response_calcu.m, SP_norm_calcu.m     Nominal closed loops / Table 2 and cost comparison
├── SP_graphplot_paper.m, SP_report.m        Paper figures / numerical summary and CSV
├── i_*.m                    canonical equations, nonlinear terms n_s, n_u, forward/backward sweeps, plant RHS
├── f_*.m                    design, CARE, Lyapunov, kernels, Hamiltonian, rate laws, simulation, metrics
├── inputData/5D_Pneumatic_nu_sat_symb/      Control problem (symbolic and analytic routes)
├── polybase/                Positive-definite value ansatz Vhat = (y'P5y)exp(g(y)) and its Gauss-Newton fit
├── VariousFunction/         sigma, sigma', sigma^-1, smooth |.|, Halton, percentile, settling time, ordblk, reig
├── SP_iterData/             Paper experiments: open_loop_simulation, relaxation_sweep, disturbance_sweep
├── LMI/                     Bundled YALMIP + SDPT3, optional LMI_riccati_check.m
├── NeuralNetwork/           Not used by this problem (README only, kept for the layout)
├── shift&rotation/          Not used by this problem (README only, kept for the layout)
├── Documents/               Paper PDF and original figures, HJ guide, Beeler–Tran–Banks paper
└── result/                  Created at run time (git-ignored)
```

As in the reference toolbox, `DirectApp.m`, the `P*` files, the `SP_*` files, `systemData.m`, `Decide_xi.m` and the `SP_iterData/*` files are **scripts** that share one workspace. The numerical work is done in functions (`f_*`, `i_*`, `polybase/*`), which take explicit structures (`cfg`, `par`, `des`, `prob`), so no globals are needed.

**Coordinates.** All HJ quantities use the scaled error y = S⁻¹(z − z*) and the scaled costate μ = Sλ. The seeds and the iteration use the diagonalised coordinates (s, χ) = T⁻¹(y, μ). The fifth state is the **raw** command u. The applied command σ(u) is always in (0, 1).

**Toolbox names.** In the toolbox convention p = Px (V = ½xᵀPx), so `P = des.Gamma = 2*P5` and `P5 = des.P` is the paper's value matrix (V ≈ yᵀP₅y). The other toolbox names map as `Rbar = D`, `Diagonalize = T`, `Diag_inv = T⁻¹` and `Ham = [A −D; −2Q −Aᵀ]`.

---

## 4. Equations and where they are implemented

| Paper | Content | File |
|---|---|---|
| (mech)–(hyst) | Actuator: angle, pressure lag, Bouc–Wen hysteresis with smooth ‖·‖ | `symbSetting.m`, `numericProblem.m`, `VariousFunction/f_smooth_abs.m` |
| (sat) | σ(u) = u_m + u_h·tanh((u − u_m)/u_h), σ′, σ⁻¹ | `VariousFunction/f_SatTanh*.m`, `f_DSatTanh.m` |
| (aug) | ż = F_σ(z) + Gν + Ed, exosystem for d | `i_state_eq_for_control.m` |
| (cost), (preH) | Q_z, R_ν, scaling S, pre-Hamiltonian | `systemData.m`, `symbSetting.m` |
| Prop. 1, (nustar), (Hbar), (lam5) | ν* = −μ₅/(2R_ν s_u), H̄, σ′(u) in λ̇₅ | `symbSetting.m`/`symbCompile.m`, `f_hamcalcu.m`, `f_inputcalcu.m`, `f_self_checks.m` |
| (canon) | Canonical equations | `i_canonical_eq.m` |
| Lemma 1, T, L | CARE, F_c, Lyapunov L, block diagonalisation | `f_design.m`, `f_care.m`, `f_lyap.m`, `P1_data_setting.m` |
| (iter) | Integral iteration, exact exponential weights, relaxation 0.8, non-uniform mesh, T_m | `SP_appNth.m`, `i_directforx.m`, `i_directforp.m`, `i_nonlinear_term.m`, `f_manifold_kernels.m`, `f_semigroup_weights.m` |
| Table 1 | Acceptance tests | `SP_manifold_check.m`, `P4_iterative_calcu.m` |
| (ansatz), Prop. 2 | V̂ = (yᵀP₅y)e^{g(y)}, 55 coefficients | `polybase/f_value_eval.m`, `f_poly_*.m` |
| (fit) | Weighted HJB collocation, anchors, ridge, Gauss–Newton | `SP_value_fit.m`, `polybase/f_value_fit.m`, `f_value_colloc.m`, `f_value_anchor_res.m` |
| (rho) | Held-out HJB residual ρ_H | `SP_model_evaluation.m` |
| Sec. 4.2, (guard) | Guarded HJ/LQR selection, rate bound 20 | `f_feedback_rate.m` |
| Sec. 5 | Ts = 0.01 s, 240 s, 2 RK4 substeps, J₂₄₀, tail RMSE, settling | `f_simulate_case.m`, `f_case_metrics.m` |
| Fig. 1 | Open-loop loops | `SP_iterData/open_loop_simulation.m` |
| Fig. 3, Remark 2 | Relaxation sweep with retuned LQR | `SP_iterData/relaxation_sweep.m` |
| Figs. 4–5, Table 2 | Nominal regulation | `SP_response_calcu.m`, `SP_norm_calcu.m` |
| Table 3 | Disturbance sweep | `SP_iterData/disturbance_sweep.m` |
| All figures | Same names and layouts as the paper | `f_paper_figures.m` |

---

## 5. Settings (`inputData/5D_Pneumatic_nu_sat_symb/systemData.m`)

| Group | Values used in the paper |
|---|---|
| Flags | `satFlag = 0` (σ is in the drift), `meshFlag = 0`, `kFlag = xiRFlag = xiTHFlag = iniFlag = 1`, `symbFlag = 1` |
| Plant | J = 0.02, k = 0.6, c = 0.08, b_p = 0.35, k_h = 0.25, τ_p = 0.2, k_u = 1, A_h = 1, γ = 0.2, β = 0.8, λ_h = 0.02, ε = 1e-4 |
| Problem | θ* = 0.35, Q_z = diag(80, 5, 2, 1, 0.2), R_ν = 0.02, S = diag(0.35, 1, 0.6, 0.4, 0.6), σ: ℝ → (0, 1) |
| Manifold | radius 2e-4, 8 Halton directions, T_m = max(40, 8/\|Re λ_max\|), mesh 0.01/0.1/0.5 s, tol 1e-7, relaxation 0.8, ≤ 80 iterations |
| Fit | 1000 Halton points (scales 0.1–1), 9 points per axis, 250 LQR states, ‖y‖ ≤ 1.6, ω_a = 0.05, λ_c = 0.005, floor 0.03, \|g\| ≤ 3 |
| Selector | ‖y‖ ≤ 1.6, ρ_H ≤ 0.05, decrease fraction 0.02, \|ν\| ≤ 20 |
| Simulation | Ts = 0.01 s, 2 RK4 substeps, 240 s from rest (u(0) = 0, σ(0) = 0.1192), tail window 30 s, band 0.007 rad |
| Experiments | B_d ∈ {0.15, 0.5, 1, 2} with forcing scale 20; λ_h ∈ {0.005, 0.02, 0.05, 0.1}; open loop σ(u) = 0.5 + 0.35 sin(πt/2), 40 s, B_d = 0.15 |

---

## 6. Expected results

`result/numerical_summary.txt` prints every value next to the paper's value. Revision 5 produced the following:

- **P1.** Poles −6.620 ± 11.619i, −12.493, −4.158 and −0.01190. T_m = 672.5 s. K₅S⁻¹ = [−1.674, 11.236, 17.590, 26.795, 20.883]. Riccati and block-diagonalisation residuals are about 2e-13 and 5e-13.
- **P4 (Table 1).** 18 of 18 seeds accepted after 5–8 iterations, giving 1 050 anchors. Worst values: relative |H̄| 1.28e-5, defect 2.13e-7, fixed-point residual 1.34e-8, terminal norm 7.24e-8, boundary residual 6.8e-21.
- **Fit.** The objective drops from 7 363 to 284 in 9 Gauss–Newton steps on 1 294 states. Anchor gradient error: median 1.84e-5, p90 3.58e-5, max 1.03e-4. Held-out ρ_H: 0.099 / 0.377 / ≈1. HJ-eligible: 153 of 500 states (30.6 %).
- **Table 2.** J₂₄₀ = 11.340 (LQR) and 8.979 (guarded HJ). HJ is used on 18 382 of 24 000 intervals, with 10 domain and 5 608 residual fallbacks. Settling time 1.0 s vs 20.5 s.
- **Table 3.** ΔJ = −7.7 %, −0.5 %, −2.1 % and −7.4 % for B_d = 0.15, 0.5, 1 and 2.
- **Relaxation.** h(240) = 0.140, 0.0159, 2.2e-4 and 1.8e-7.

The guarded selector switches between HJ and LQR many times around ρ_H = 0.05 (roughly 14–68 s). With the symbolic route the code evaluates the same compiled functions as Revision 5. A different MATLAB release or BLAS, or the analytic route, can change single switching decisions at the 1e-15 level. Such changes leave LQR, the manifold and the fit unchanged, but the HJ use count and the guarded-HJ costs can then differ in the last digits.

---

## 7. Differences from the reference toolbox

- **No input saturation function in the HJB problem.** The actuator saturation σ(u) is part of the drift, and the input ν = u̇ is unconstrained in the design. Hence `satFlag = 0`, and the reference `f_Sat`/`f_SatS` functions are replaced by `f_SatTanh*`.
- **Integral iteration.** The forward and backward integrals of the successive approximation are evaluated with exact exponential weights on a non-uniform mesh (`f_manifold_kernels.m`) instead of `ode45`. Seeds are relaxed, accepted only after the Table 1 tests, and there is no negative-time extension.
- **Value representation.** One scalar, positive-definite value function replaces the component-wise polynomial or neural-network fit of p(x).
- **Feedback.** The controller is sampled, with a held rate and a guarded HJ/LQR selector, instead of the continuous `ode45` closed loop.
- **Interaction.** The run is non-interactive by default (flags in `systemData.m`), the toolbox-free CARE and Lyapunov solvers are used, and the costs follow the paper's quadrature (`f_costopt.m` is kept as a trapezoidal cross-check).

---

## 8. Limitations (as stated in the paper)

- Exact costates exist only for ‖y‖ ≲ 2.6e-4, while the controller acts up to ‖y‖ = 1.6.
- The guard is a pointwise check at the samples for d = 0. It does not certify the inter-sample behaviour, the switching or disturbance rejection.
- The rate bound is applied after an unconstrained-rate HJB design.
- The parameters are normalised demonstration values. They have not been identified on a device.

## References

1. N. Sakamoto, A. J. van der Schaft, "Analytical approximation methods for the stabilizing solution of the Hamilton–Jacobi equation," *IEEE TAC* 53(10), 2008.
2. N. Sakamoto, "Case studies on the application of the stable manifold approach for nonlinear optimal control design," *Automatica* 49(2), 2013.
3. `Documents/HJ_Guide 20231110.pdf` (user guide of the toolbox) and `Documents/NS_DD_ACNDC_2027.pdf` (the paper).

Third-party code: YALMIP (J. Löfberg) and SDPT3 (Toh, Todd, Tütüncü, GPL v2) in `LMI/`. `ordblk.m` and `reig.m` (R. Y. Chiang, M. G. Safonov) in `VariousFunction/`.
