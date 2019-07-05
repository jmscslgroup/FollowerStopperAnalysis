/*
 * carPlant_smoothUp_followerStopper_data.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "carPlant_smoothUp_followerStopper".
 *
 * Model version              : 1.31
 * Simulink Coder version : 9.0 (R2018b) 24-May-2018
 * C source code generated on : Fri Jul  5 13:50:39 2019
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "carPlant_smoothUp_followerStopper.h"
#include "carPlant_smoothUp_followerStopper_private.h"

/* Block parameters (default storage) */
P_carPlant_smoothUp_followerS_T carPlant_smoothUp_followerSto_P = {
  /* Mask Parameter: FilteredDerivativeDiscreteorCon
   * Referenced by: '<S5>/[A,B]'
   */
  -20.0,

  /* Mask Parameter: FilteredDerivativeDiscreteorC_f
   * Referenced by: '<S6>/[A,B]'
   */
  -20.0,

  /* Mask Parameter: FollowerVehicle_Af
   * Referenced by: '<S45>/.5.*A.*Pabs.//R.//T'
   */
  2.0,

  /* Mask Parameter: LeadVehicle_Af
   * Referenced by: '<S241>/.5.*A.*Pabs.//R.//T'
   */
  2.0,

  /* Mask Parameter: FilteredDerivativeDiscreteorC_p
   * Referenced by: '<S5>/[A,B]'
   */
  20.0,

  /* Mask Parameter: FilteredDerivativeDiscreteorC_m
   * Referenced by: '<S6>/[A,B]'
   */
  20.0,

  /* Mask Parameter: FollowerVehicle_Cd
   * Referenced by: '<S45>/Constant'
   */
  0.3,

  /* Mask Parameter: LeadVehicle_Cd
   * Referenced by: '<S241>/Constant'
   */
  0.3,

  /* Mask Parameter: FollowerVehicle_Cl
   * Referenced by: '<S45>/Constant1'
   */
  0.1,

  /* Mask Parameter: LeadVehicle_Cl
   * Referenced by: '<S241>/Constant1'
   */
  0.1,

  /* Mask Parameter: FollowerVehicle_Cpm
   * Referenced by: '<S45>/Constant2'
   */
  0.1,

  /* Mask Parameter: LeadVehicle_Cpm
   * Referenced by: '<S241>/Constant2'
   */
  0.1,

  /* Mask Parameter: FollowerVehicle_Cs
   * Referenced by: '<S45>/Cs'
   */
  { 0.0, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.3,
    0.32999999999999996, 0.36, 0.39, 0.42, 0.45, 0.48000000000000004, 0.51, 0.54,
    0.57000000000000006, 0.60000000000000009, 0.63, 0.66, 0.69000000000000006,
    0.72, 0.75, 0.78, 0.81, 0.84000000000000008, 0.87, 0.9 },

  /* Mask Parameter: LeadVehicle_Cs
   * Referenced by: '<S241>/Cs'
   */
  { 0.0, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.3,
    0.32999999999999996, 0.36, 0.39, 0.42, 0.45, 0.48000000000000004, 0.51, 0.54,
    0.57000000000000006, 0.60000000000000009, 0.63, 0.66, 0.69000000000000006,
    0.72, 0.75, 0.78, 0.81, 0.84000000000000008, 0.87, 0.9 },

  /* Mask Parameter: FollowerVehicle_Cy_f
   * Referenced by: '<S41>/Cyf'
   */
  12000.0,

  /* Mask Parameter: LeadVehicle_Cy_f
   * Referenced by: '<S237>/Cyf'
   */
  12000.0,

  /* Mask Parameter: FollowerVehicle_Cy_r
   * Referenced by: '<S41>/Cyr'
   */
  11000.0,

  /* Mask Parameter: LeadVehicle_Cy_r
   * Referenced by: '<S237>/Cyr'
   */
  11000.0,

  /* Mask Parameter: FollowerVehicle_Cym
   * Referenced by: '<S45>/Cym'
   */
  { 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12,
    0.13, 0.14, 0.15, 0.15999999999999998, 0.16999999999999998, 0.18, 0.19,
    0.19999999999999998, 0.21, 0.21999999999999997, 0.22999999999999998, 0.24,
    0.25, 0.26, 0.27, 0.27999999999999997, 0.29, 0.3 },

  /* Mask Parameter: LeadVehicle_Cym
   * Referenced by: '<S241>/Cym'
   */
  { 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12,
    0.13, 0.14, 0.15, 0.15999999999999998, 0.16999999999999998, 0.18, 0.19,
    0.19999999999999998, 0.21, 0.21999999999999997, 0.22999999999999998, 0.24,
    0.25, 0.26, 0.27, 0.27999999999999997, 0.29, 0.3 },

  /* Mask Parameter: KinematicSteering1_Db
   * Referenced by: '<S205>/Backlash'
   */
  0.0,

  /* Mask Parameter: KinematicSteering_Db
   * Referenced by: '<S186>/Backlash'
   */
  0.0,

  /* Mask Parameter: FollowerVehicle_Fznom
   * Referenced by: '<S8>/vehicle model'
   */
  5000.0,

  /* Mask Parameter: LeadVehicle_Fznom
   * Referenced by: '<S11>/vehicle model'
   */
  5000.0,

  /* Mask Parameter: FollowerVehicle_Izz
   * Referenced by: '<S8>/vehicle model'
   */
  4000.0,

  /* Mask Parameter: LeadVehicle_Izz
   * Referenced by: '<S11>/vehicle model'
   */
  4000.0,

  /* Mask Parameter: FilteredDerivativeDiscreteor_me
   * Referenced by: '<S5>/Gain'
   */
  1.0,

  /* Mask Parameter: FilteredDerivativeDiscreteorC_h
   * Referenced by: '<S6>/Gain'
   */
  1.0,

  /* Mask Parameter: FollowerVehicle_NF
   * Referenced by: '<S8>/vehicle model'
   */
  2.0,

  /* Mask Parameter: LeadVehicle_NF
   * Referenced by: '<S11>/vehicle model'
   */
  2.0,

  /* Mask Parameter: FollowerVehicle_NR
   * Referenced by: '<S8>/vehicle model'
   */
  2.0,

  /* Mask Parameter: LeadVehicle_NR
   * Referenced by: '<S11>/vehicle model'
   */
  2.0,

  /* Mask Parameter: FollowerVehicle_Pabs
   * Referenced by: '<S45>/.5.*A.*Pabs.//R.//T'
   */
  101325.0,

  /* Mask Parameter: LeadVehicle_Pabs
   * Referenced by: '<S241>/.5.*A.*Pabs.//R.//T'
   */
  101325.0,

  /* Mask Parameter: DragForce_R
   * Referenced by: '<S45>/.5.*A.*Pabs.//R.//T'
   */
  287.058,

  /* Mask Parameter: DragForce_R_e
   * Referenced by: '<S241>/.5.*A.*Pabs.//R.//T'
   */
  287.058,

  /* Mask Parameter: KinematicSteering1_StrgRatio
   * Referenced by: '<S206>/Gain'
   */
  100.0,

  /* Mask Parameter: KinematicSteering_StrgRatio
   * Referenced by: '<S187>/Gain'
   */
  100.0,

  /* Mask Parameter: KinematicSteering1_StrgRng
   * Referenced by: '<S205>/Saturation'
   */
  3.9269908169872414,

  /* Mask Parameter: KinematicSteering_StrgRng
   * Referenced by: '<S186>/Saturation'
   */
  3.9269908169872414,

  /* Mask Parameter: FilteredDerivativeDiscreteor_hf
   * Referenced by: '<S19>/Time constant'
   */
  0.1,

  /* Mask Parameter: FilteredDerivativeDiscreteor_md
   * Referenced by: '<S23>/Time constant'
   */
  0.1,

  /* Mask Parameter: FollowerVehicle_Tair
   * Referenced by: '<S45>/.5.*A.*Pabs.//R.//T'
   */
  273.0,

  /* Mask Parameter: LeadVehicle_Tair
   * Referenced by: '<S241>/.5.*A.*Pabs.//R.//T'
   */
  273.0,

  /* Mask Parameter: KinematicSteering1_TrckWdth
   * Referenced by: '<S207>/Constant1'
   */
  1.0,

  /* Mask Parameter: KinematicSteering_TrckWdth
   * Referenced by: '<S188>/Constant1'
   */
  1.0,

  /* Mask Parameter: KinematicSteering1_WhlBase
   * Referenced by: '<S207>/Constant'
   */
  1.524,

  /* Mask Parameter: KinematicSteering_WhlBase
   * Referenced by: '<S188>/Constant'
   */
  1.524,

  /* Mask Parameter: LeadVehicle_X_o
   * Referenced by: '<S290>/X_o'
   */
  10.0,

  /* Mask Parameter: FollowerVehicle_X_o
   * Referenced by: '<S94>/X_o'
   */
  0.0,

  /* Mask Parameter: LeadVehicle_Y_o
   * Referenced by: '<S290>/Y_o'
   */
  0.0,

  /* Mask Parameter: FollowerVehicle_Y_o
   * Referenced by: '<S94>/Y_o'
   */
  0.0,

  /* Mask Parameter: FollowerVehicle_a
   * Referenced by:
   *   '<S8>/vehicle model'
   *   '<S45>/Constant3'
   *   '<S95>/a'
   *   '<S96>/a'
   */
  1.4,

  /* Mask Parameter: LeadVehicle_a
   * Referenced by:
   *   '<S11>/vehicle model'
   *   '<S241>/Constant3'
   *   '<S291>/a'
   *   '<S292>/a'
   */
  1.4,

  /* Mask Parameter: FollowerVehicle_b
   * Referenced by:
   *   '<S8>/vehicle model'
   *   '<S45>/Constant3'
   *   '<S98>/b'
   *   '<S99>/b'
   */
  1.6,

  /* Mask Parameter: LeadVehicle_b
   * Referenced by:
   *   '<S11>/vehicle model'
   *   '<S241>/Constant3'
   *   '<S294>/b'
   *   '<S295>/b'
   */
  1.6,

  /* Mask Parameter: FollowerVehicle_beta_w
   * Referenced by:
   *   '<S45>/Cs'
   *   '<S45>/Cym'
   */
  { 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12,
    0.13, 0.14, 0.15, 0.15999999999999998, 0.16999999999999998, 0.18, 0.19,
    0.19999999999999998, 0.21, 0.21999999999999997, 0.22999999999999998, 0.24,
    0.25, 0.26, 0.27, 0.27999999999999997, 0.29, 0.3 },

  /* Mask Parameter: LeadVehicle_beta_w
   * Referenced by:
   *   '<S241>/Cs'
   *   '<S241>/Cym'
   */
  { 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12,
    0.13, 0.14, 0.15, 0.15999999999999998, 0.16999999999999998, 0.18, 0.19,
    0.19999999999999998, 0.21, 0.21999999999999997, 0.22999999999999998, 0.24,
    0.25, 0.26, 0.27, 0.27999999999999997, 0.29, 0.3 },

  /* Mask Parameter: FollowerVehicle_d
   * Referenced by:
   *   '<S8>/vehicle model'
   *   '<S95>/d'
   *   '<S96>/d'
   *   '<S98>/d'
   *   '<S99>/d'
   */
  0.0,

  /* Mask Parameter: LeadVehicle_d
   * Referenced by:
   *   '<S11>/vehicle model'
   *   '<S291>/d'
   *   '<S292>/d'
   *   '<S294>/d'
   *   '<S295>/d'
   */
  0.0,

  /* Mask Parameter: FollowerVehicle_g
   * Referenced by: '<S8>/vehicle model'
   */
  9.81,

  /* Mask Parameter: LeadVehicle_g
   * Referenced by: '<S11>/vehicle model'
   */
  9.81,

  /* Mask Parameter: FollowerVehicle_h
   * Referenced by:
   *   '<S8>/vehicle model'
   *   '<S95>/h'
   *   '<S96>/h'
   *   '<S98>/h'
   *   '<S99>/h'
   */
  0.35,

  /* Mask Parameter: LeadVehicle_h
   * Referenced by:
   *   '<S11>/vehicle model'
   *   '<S291>/h'
   *   '<S292>/h'
   *   '<S294>/h'
   *   '<S295>/h'
   */
  0.35,

  /* Mask Parameter: LeadVehicle_latOff
   * Referenced by: '<S293>/latOff'
   */
  0.0,

  /* Mask Parameter: FollowerVehicle_latOff
   * Referenced by: '<S97>/latOff'
   */
  0.0,

  /* Mask Parameter: LeadVehicle_longOff
   * Referenced by: '<S293>/longOff'
   */
  0.0,

  /* Mask Parameter: FollowerVehicle_longOff
   * Referenced by: '<S97>/longOff'
   */
  0.0,

  /* Mask Parameter: FollowerVehicle_m
   * Referenced by: '<S8>/vehicle model'
   */
  2000.0,

  /* Mask Parameter: LeadVehicle_m
   * Referenced by: '<S11>/vehicle model'
   */
  2000.0,

  /* Mask Parameter: FilteredDerivativeDiscreteorC_g
   * Referenced by: '<S19>/Minimum sampling to time constant ratio'
   */
  10.0,

  /* Mask Parameter: FilteredDerivativeDiscreteorC_d
   * Referenced by: '<S23>/Minimum sampling to time constant ratio'
   */
  10.0,

  /* Mask Parameter: LeadVehicle_psi_o
   * Referenced by: '<S378>/psi_o'
   */
  0.0,

  /* Mask Parameter: FollowerVehicle_psi_o
   * Referenced by: '<S182>/psi_o'
   */
  0.0,

  /* Mask Parameter: LeadVehicle_r_o
   * Referenced by: '<S378>/r_o'
   */
  0.0,

  /* Mask Parameter: FollowerVehicle_r_o
   * Referenced by: '<S182>/r_o'
   */
  0.0,

  /* Mask Parameter: FollowerVehicle_sigma_f
   * Referenced by: '<S178>/Constant1'
   */
  0.1,

  /* Mask Parameter: LeadVehicle_sigma_f
   * Referenced by: '<S374>/Constant1'
   */
  0.1,

  /* Mask Parameter: FollowerVehicle_sigma_r
   * Referenced by: '<S178>/Constant2'
   */
  0.1,

  /* Mask Parameter: LeadVehicle_sigma_r
   * Referenced by: '<S374>/Constant2'
   */
  0.1,

  /* Mask Parameter: DeadMansSwitch_stepSize
   * Referenced by: '<S4>/Simulate step size'
   */
  0.05,

  /* Mask Parameter: div0protectpoly_thresh
   * Referenced by:
   *   '<S383>/Constant'
   *   '<S384>/Constant'
   */
  0.001,

  /* Mask Parameter: div0protectabspoly3_thresh
   * Referenced by:
   *   '<S219>/Constant'
   *   '<S220>/Constant'
   */
  1.4901161193847656E-8,

  /* Mask Parameter: div0protectabspoly_thresh
   * Referenced by:
   *   '<S213>/Constant'
   *   '<S214>/Constant'
   */
  1.4901161193847656E-8,

  /* Mask Parameter: div0protectabspoly3_thresh_g
   * Referenced by:
   *   '<S200>/Constant'
   *   '<S201>/Constant'
   */
  1.4901161193847656E-8,

  /* Mask Parameter: div0protectabspoly_thresh_b
   * Referenced by:
   *   '<S194>/Constant'
   *   '<S195>/Constant'
   */
  1.4901161193847656E-8,

  /* Mask Parameter: DeadMansSwitch_timeout
   * Referenced by: '<S4>/Timeout in seconds'
   */
  0.2,

  /* Mask Parameter: LeadVehicle_vertOff
   * Referenced by: '<S293>/vertOff'
   */
  0.0,

  /* Mask Parameter: FollowerVehicle_vertOff
   * Referenced by: '<S97>/vertOff'
   */
  0.0,

  /* Mask Parameter: FollowerVehicle_w
   * Referenced by:
   *   '<S8>/vehicle model'
   *   '<S95>/w'
   *   '<S96>/w'
   *   '<S98>/w'
   *   '<S99>/w'
   */
  { 1.4, 1.4 },

  /* Mask Parameter: LeadVehicle_w
   * Referenced by:
   *   '<S11>/vehicle model'
   *   '<S291>/w'
   *   '<S292>/w'
   *   '<S294>/w'
   *   '<S295>/w'
   */
  { 1.4, 1.4 },

  /* Mask Parameter: FollowerVehicle_xdot_tol
   * Referenced by: '<S8>/vehicle model'
   */
  0.01,

  /* Mask Parameter: LeadVehicle_xdot_tol
   * Referenced by: '<S11>/vehicle model'
   */
  0.01,

  /* Mask Parameter: LeadVehicle_ydot_o
   * Referenced by: '<S378>/ydot_o'
   */
  0.0,

  /* Mask Parameter: FollowerVehicle_ydot_o
   * Referenced by: '<S182>/ydot_o'
   */
  0.0,

  /* Computed Parameter: Out1_Y0
   * Referenced by: '<S382>/Out1'
   */
  {
    {
      0.0,                             /* X */
      0.0,                             /* Y */
      0.0                              /* Z */
    },                                 /* Linear */

    {
      0.0,                             /* X */
      0.0,                             /* Y */
      0.0                              /* Z */
    }                                  /* Angular */
  },

  /* Computed Parameter: Constant_Value
   * Referenced by: '<S15>/Constant'
   */
  {
    {
      0.0,                             /* X */
      0.0,                             /* Y */
      0.0                              /* Z */
    },                                 /* Linear */

    {
      0.0,                             /* X */
      0.0,                             /* Y */
      0.0                              /* Z */
    }                                  /* Angular */
  },

  /* Computed Parameter: Constant_Value_b
   * Referenced by: '<S3>/Constant'
   */
  {
    {
      0.0,                             /* X */
      0.0,                             /* Y */
      0.0                              /* Z */
    },                                 /* Linear */

    {
      0.0,                             /* X */
      0.0,                             /* Y */
      0.0                              /* Z */
    }                                  /* Angular */
  },

  /* Computed Parameter: Constant_Value_k
   * Referenced by: '<S2>/Constant'
   */
  {
    {
      0.0,                             /* X */
      0.0,                             /* Y */
      0.0                              /* Z */
    },                                 /* Linear */

    {
      0.0,                             /* X */
      0.0,                             /* Y */
      0.0                              /* Z */
    }                                  /* Angular */
  },

  /* Expression: 0
   * Referenced by: '<S16>/Switch1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<Root>/Constant1'
   */
  0.0,

  /* Expression: Fxtire_sat
   * Referenced by: '<S8>/vehicle model'
   */
  70000.0,

  /* Expression: Fytire_sat
   * Referenced by: '<S8>/vehicle model'
   */
  70000.0,

  /* Expression: Fxtire_sat
   * Referenced by: '<S11>/vehicle model'
   */
  70000.0,

  /* Expression: Fytire_sat
   * Referenced by: '<S11>/vehicle model'
   */
  70000.0,

  /* Expression: 0
   * Referenced by: '<S5>/Constant'
   */
  0.0,

  /* Computed Parameter: Integrator_gainval
   * Referenced by: '<S22>/Integrator'
   */
  0.01,

  /* Expression: antiwindupUpperLimit
   * Referenced by: '<S22>/Integrator'
   */
  0.0,

  /* Expression: antiwindupLowerLimit
   * Referenced by: '<S22>/Integrator'
   */
  0.0,

  /* Expression: windupUpperLimit
   * Referenced by: '<S22>/Saturation'
   */
  0.0,

  /* Expression: windupLowerLimit
   * Referenced by: '<S22>/Saturation'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S6>/Constant'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S16>/Constant'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<Root>/Switch'
   */
  0.0,

  /* Computed Parameter: Integrator_gainval_p
   * Referenced by: '<S26>/Integrator'
   */
  0.01,

  /* Expression: antiwindupUpperLimit
   * Referenced by: '<S26>/Integrator'
   */
  0.0,

  /* Expression: antiwindupLowerLimit
   * Referenced by: '<S26>/Integrator'
   */
  0.0,

  /* Expression: windupUpperLimit
   * Referenced by: '<S26>/Saturation'
   */
  0.0,

  /* Expression: windupLowerLimit
   * Referenced by: '<S26>/Saturation'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<Root>/Constant12'
   */
  1.0,

  /* Expression: 3.748509500812156
   * Referenced by: '<S17>/reference_vel1'
   */
  3.7485095008121561,

  /* Expression: 1.5
   * Referenced by: '<S17>/max_accel1'
   */
  1.5,

  /* Expression: -1.5
   * Referenced by: '<S17>/max_decel1'
   */
  -1.5,

  /* Expression: 0.0
   * Referenced by: '<Root>/Steering angle'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S290>/Constant'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S290>/Constant7'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S290>/Constant2'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S94>/Constant'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S94>/Constant7'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S94>/Constant2'
   */
  0.0,

  /* Expression: [1.5 1.0 0.5]
   * Referenced by: '<S13>/decel'
   */
  { 1.5, 1.0, 0.5 },

  /* Expression: [0.7 0.7; 0.7 0.7]
   * Referenced by: '<Root>/Coefficient of Friction: Dry'
   */
  { 0.7, 0.7, 0.7, 0.7 },

  /* Expression: 0
   * Referenced by: '<Root>/Constant3'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<Root>/Constant4'
   */
  0.0,

  /* Expression: [0 0]
   * Referenced by: '<S45>/Crm'
   */
  { 0.0, 0.0 },

  /* Expression: [-1 1]
   * Referenced by: '<S45>/Crm'
   */
  { -1.0, 1.0 },

  /* Expression: ones(1,3)
   * Referenced by: '<S45>/Constant4'
   */
  { 1.0, 1.0, 1.0 },

  /* Expression: 0
   * Referenced by: '<S45>/Switch'
   */
  0.0,

  /* Expression: -1/2
   * Referenced by: '<S95>/Gain'
   */
  -0.5,

  /* Expression: 1/2
   * Referenced by: '<S96>/Gain'
   */
  0.5,

  /* Expression: -1/2
   * Referenced by: '<S98>/Gain'
   */
  -0.5,

  /* Expression: 1/2
   * Referenced by: '<S99>/Gain'
   */
  0.5,

  /* Expression: 0
   * Referenced by: '<S205>/Backlash'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S10>/index'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S10>/Switch'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S10>/Switch1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S174>/Constant'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S94>/Constant4'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S94>/Constant5'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S94>/Constant6'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S181>/lateral'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S186>/Backlash'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S9>/index'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S9>/Switch'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S9>/Switch1'
   */
  0.0,

  /* Expression: [0 0]
   * Referenced by: '<S241>/Crm'
   */
  { 0.0, 0.0 },

  /* Expression: [-1 1]
   * Referenced by: '<S241>/Crm'
   */
  { -1.0, 1.0 },

  /* Expression: ones(1,3)
   * Referenced by: '<S241>/Constant4'
   */
  { 1.0, 1.0, 1.0 },

  /* Expression: 0
   * Referenced by: '<S241>/Switch'
   */
  0.0,

  /* Expression: -1/2
   * Referenced by: '<S291>/Gain'
   */
  -0.5,

  /* Expression: 1/2
   * Referenced by: '<S292>/Gain'
   */
  0.5,

  /* Expression: -1/2
   * Referenced by: '<S294>/Gain'
   */
  -0.5,

  /* Expression: 1/2
   * Referenced by: '<S295>/Gain'
   */
  0.5,

  /* Expression: 0
   * Referenced by: '<S370>/Constant'
   */
  0.0,

  /* Expression: mu
   * Referenced by: '<S350>/Constant'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S290>/Constant4'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S290>/Constant5'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S290>/Constant6'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S377>/lateral'
   */
  0.0
};
