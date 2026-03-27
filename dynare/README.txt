Dynare Folder - Quick File Guide

pf_cbam.mod: Main two-country CBAM model, perfect foresight with introduction of CBAM. {0,1} switch to have Taylor rule or Ramsey monetary policy. Calls a ss_pf.m.
ss_pf.m: Steady-state solver used by pf_cbam.mod for the transition model block. Big system to account for the fact that disutility for labor is fixed at time 0, takes into account all the GE effects of CBAM.

compare_pf_nominal_rigidities.m: Calls pf_cbam.mod and a generated flex-price version (kappa=0) and compares perfect-foresight paths.
pf_cbam_transition_flex.mod: Auto-generated flexible-price counterfactual used by compare_pf_nominal_rigidities.m.


model.mod: Model version with CBAM used to observe the economy reaction to a shock. Calls ss_foreign.m.
ss_foreign.m: small fixed-point steady-state helper for the foreign decision of pricing and abatement.

model_no_cbam.mod: Same model CBAM, used as the counterfactual baseline. Calls ss_no_cbam.m.
ss_no_cbam.m: Steady-state helper for the no-CBAM model, to compute the pricing decision.


compare_models.m: Calls model.mod and model_no_cbam.mod and compares IRFs/shock responses with and without CBAM.




