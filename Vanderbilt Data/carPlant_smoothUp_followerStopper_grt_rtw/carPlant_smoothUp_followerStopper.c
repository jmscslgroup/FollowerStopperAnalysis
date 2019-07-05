/*
 * carPlant_smoothUp_followerStopper.c
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

/* Block signals (default storage) */
B_carPlant_smoothUp_followerS_T carPlant_smoothUp_followerSto_B;

/* Continuous states */
X_carPlant_smoothUp_followerS_T carPlant_smoothUp_followerSto_X;

/* Block states (default storage) */
DW_carPlant_smoothUp_follower_T carPlant_smoothUp_followerSt_DW;

/* Real-time model */
RT_MODEL_carPlant_smoothUp_fo_T carPlant_smoothUp_followerSt_M_;
RT_MODEL_carPlant_smoothUp_fo_T *const carPlant_smoothUp_followerSt_M =
  &carPlant_smoothUp_followerSt_M_;

/* Forward declaration for local functions */
static void carPlant_s_automlvehdynftiresat(real_T Ftire_y, real_T b_Fxtire_sat,
  real_T b_Fytire_sat, real_T *Ftire_xs, real_T *Ftire_ys);
static void matlabCodegenHandle_matlabC_em4(robotics_slros_internal_bl_em_T *obj);
static void matlabCodegenHandle_matlabCodeg(robotics_slros_internal_block_T *obj);
static void matlabCodegenHandle_matlabCo_em(robotics_slros_internal_blo_e_T *obj);
real_T look1_binlcpw(real_T u0, const real_T bp0[], const real_T table[],
                     uint32_T maxIndex)
{
  real_T frac;
  uint32_T iRght;
  uint32_T iLeft;
  uint32_T bpIdx;

  /* Column-major Lookup 1-D
     Search method: 'binary'
     Use previous index: 'off'
     Interpolation method: 'Linear point-slope'
     Extrapolation method: 'Clip'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Clip'
     Use previous index: 'off'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u0 <= bp0[0U]) {
    iLeft = 0U;
    frac = 0.0;
  } else if (u0 < bp0[maxIndex]) {
    /* Binary Search */
    bpIdx = maxIndex >> 1U;
    iLeft = 0U;
    iRght = maxIndex;
    while (iRght - iLeft > 1U) {
      if (u0 < bp0[bpIdx]) {
        iRght = bpIdx;
      } else {
        iLeft = bpIdx;
      }

      bpIdx = (iRght + iLeft) >> 1U;
    }

    frac = (u0 - bp0[iLeft]) / (bp0[iLeft + 1U] - bp0[iLeft]);
  } else {
    iLeft = maxIndex - 1U;
    frac = 1.0;
  }

  /* Column-major Interpolation 1-D
     Interpolation method: 'Linear point-slope'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'portable wrapping'
   */
  return (table[iLeft + 1U] - table[iLeft]) * frac + table[iLeft];
}

real_T look1_binlxpw(real_T u0, const real_T bp0[], const real_T table[],
                     uint32_T maxIndex)
{
  real_T frac;
  uint32_T iRght;
  uint32_T iLeft;
  uint32_T bpIdx;

  /* Column-major Lookup 1-D
     Search method: 'binary'
     Use previous index: 'off'
     Interpolation method: 'Linear point-slope'
     Extrapolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Linear'
     Use previous index: 'off'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u0 <= bp0[0U]) {
    iLeft = 0U;
    frac = (u0 - bp0[0U]) / (bp0[1U] - bp0[0U]);
  } else if (u0 < bp0[maxIndex]) {
    /* Binary Search */
    bpIdx = maxIndex >> 1U;
    iLeft = 0U;
    iRght = maxIndex;
    while (iRght - iLeft > 1U) {
      if (u0 < bp0[bpIdx]) {
        iRght = bpIdx;
      } else {
        iLeft = bpIdx;
      }

      bpIdx = (iRght + iLeft) >> 1U;
    }

    frac = (u0 - bp0[iLeft]) / (bp0[iLeft + 1U] - bp0[iLeft]);
  } else {
    iLeft = maxIndex - 1U;
    frac = (u0 - bp0[maxIndex - 1U]) / (bp0[maxIndex] - bp0[maxIndex - 1U]);
  }

  /* Column-major Interpolation 1-D
     Interpolation method: 'Linear point-slope'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'portable wrapping'
   */
  return (table[iLeft + 1U] - table[iLeft]) * frac + table[iLeft];
}

/*
 * This function updates continuous states using the ODE3 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  /* Solver Matrices */
  static const real_T rt_ODE3_A[3] = {
    1.0/2.0, 3.0/4.0, 1.0
  };

  static const real_T rt_ODE3_B[3][3] = {
    { 1.0/2.0, 0.0, 0.0 },

    { 0.0, 3.0/4.0, 0.0 },

    { 2.0/9.0, 1.0/3.0, 4.0/9.0 }
  };

  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE3_IntgData *id = (ODE3_IntgData *)rtsiGetSolverData(si);
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T hB[3];
  int_T i;
  int_T nXc = 20;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(y, x,
                (uint_T)nXc*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  carPlant_smoothUp_followerStopper_derivatives();

  /* f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*)); */
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  carPlant_smoothUp_followerStopper_step();
  carPlant_smoothUp_followerStopper_derivatives();

  /* f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*)); */
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE3_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[1]);
  rtsiSetdX(si, f2);
  carPlant_smoothUp_followerStopper_step();
  carPlant_smoothUp_followerStopper_derivatives();

  /* tnew = t + hA(3);
     ynew = y + f*hB(:,3); */
  for (i = 0; i <= 2; i++) {
    hB[i] = h * rt_ODE3_B[2][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2]);
  }

  rtsiSetT(si, tnew);
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/*
 * Output and update for atomic system:
 *    '<S94>/COMB2I'
 *    '<S290>/COMB2I'
 */
void carPlant_smoothUp_follow_COMB2I(const real_T rtu_x[4],
  B_COMB2I_carPlant_smoothUp_fo_T *localB)
{
  real_T y_tmp;
  real_T y_tmp_0;
  y_tmp = sin(rtu_x[2]);
  y_tmp_0 = cos(rtu_x[2]);
  localB->y[0] = rtu_x[0] * y_tmp_0 - rtu_x[1] * y_tmp;
  localB->y[1] = rtu_x[0] * y_tmp + rtu_x[1] * y_tmp_0;
}

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T tmp;
  real_T tmp_0;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else {
    tmp = fabs(u0);
    tmp_0 = fabs(u1);
    if (rtIsInf(u1)) {
      if (tmp == 1.0) {
        y = 1.0;
      } else if (tmp > 1.0) {
        if (u1 > 0.0) {
          y = (rtInf);
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = (rtInf);
      }
    } else if (tmp_0 == 0.0) {
      y = 1.0;
    } else if (tmp_0 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = (rtNaN);
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

real_T rt_atan2d_snf(real_T u0, real_T u1)
{
  real_T y;
  int32_T u0_0;
  int32_T u1_0;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    if (u0 > 0.0) {
      u0_0 = 1;
    } else {
      u0_0 = -1;
    }

    if (u1 > 0.0) {
      u1_0 = 1;
    } else {
      u1_0 = -1;
    }

    y = atan2(u0_0, u1_0);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
}

/* Function for MATLAB Function: '<S8>/vehicle model' */
static void carPlant_s_automlvehdynftiresat(real_T Ftire_y, real_T b_Fxtire_sat,
  real_T b_Fytire_sat, real_T *Ftire_xs, real_T *Ftire_ys)
{
  real_T theta_Ftire;
  real_T b_a;
  real_T b_a_tmp;
  theta_Ftire = rt_atan2d_snf(0.0, Ftire_y);
  b_a_tmp = cos(theta_Ftire);
  b_a = b_Fxtire_sat * b_a_tmp;
  theta_Ftire = b_Fytire_sat * sin(theta_Ftire);
  theta_Ftire = b_Fxtire_sat * b_Fytire_sat / sqrt(b_a * b_a + theta_Ftire *
    theta_Ftire) * b_a_tmp;
  *Ftire_xs = 0.0;
  *Ftire_ys = Ftire_y;
  if (fabs(Ftire_y) > fabs(theta_Ftire)) {
    *Ftire_ys = theta_Ftire;
  }
}

static void matlabCodegenHandle_matlabC_em4(robotics_slros_internal_bl_em_T *obj)
{
  if (!obj->matlabCodegenIsDeleted) {
    obj->matlabCodegenIsDeleted = true;
  }
}

static void matlabCodegenHandle_matlabCodeg(robotics_slros_internal_block_T *obj)
{
  if (!obj->matlabCodegenIsDeleted) {
    obj->matlabCodegenIsDeleted = true;
  }
}

static void matlabCodegenHandle_matlabCo_em(robotics_slros_internal_blo_e_T *obj)
{
  if (!obj->matlabCodegenIsDeleted) {
    obj->matlabCodegenIsDeleted = true;
  }
}

/* Model step function */
void carPlant_smoothUp_followerStopper_step(void)
{
  /* local block i/o variables */
  real_T rtb_Product1_o;
  real_T dx1;
  real_T dx3;
  real_T delta_fl;
  real_T delta_fr;
  real_T delta_rl;
  real_T delta_rr;
  real_T alfa_fl;
  real_T alfa_fr;
  real_T alfa_rl;
  real_T alfa_rr;
  real_T Fz_fnom;
  real_T Fz_rnom;
  real_T yddot;
  real_T rdot;
  real_T b_Fy_rl;
  real_T b_Fy_rr;
  real_T z_data;
  real_T yddot_0;
  real_T rdot_0;
  real_T z_data_0;
  real_T z1_data;
  boolean_T b_varargout_1;
  int32_T iU;
  real_T rtb_Integrator_d;
  SL_Bus_carPlant_smoothUp_followe_Twist_domkw2 rtb_BusAssignment8;
  real_T rtb_VectorConcatenate_al[9];
  real_T rtb_uAPabsRT[6];
  int32_T i;
  real_T tmp[3];
  real_T rtb_VectorConcatenate1_k;
  real_T rtb_sincos_o1_idx_1;
  real_T rtb_sincos_o1_idx_2;
  real_T rtb_sincos_o2_idx_1;
  real_T rtb_sincos_o2_idx_2;
  real_T rtb_sincos_o1_idx_0;
  real_T rtb_UnaryMinus1_g_idx_2;
  real_T rtb_VectorConcatenate1_m_idx_0;
  real_T Fz_idx_0;
  real_T Fz_idx_1;
  real_T Fz_idx_2;
  real_T Fz_idx_3;
  real_T rtb_VectorConcatenate1_m_idx_0_;
  real_T b_Fy_fr_tmp;
  real_T yddot_tmp;
  real_T b_Fy_rl_tmp;
  real_T b_Fy_rl_tmp_0;
  real_T b_Fy_rr_tmp;
  real_T b_Fy_rr_tmp_0;
  real_T rtb_UnaryMinus1_g_idx_0_tmp;
  real_T rtb_sincos_o2_idx_0_tmp;
  if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
    /* set solver stop time */
    if (!(carPlant_smoothUp_followerSt_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&carPlant_smoothUp_followerSt_M->solverInfo,
                            ((carPlant_smoothUp_followerSt_M->Timing.clockTickH0
        + 1) * carPlant_smoothUp_followerSt_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&carPlant_smoothUp_followerSt_M->solverInfo,
                            ((carPlant_smoothUp_followerSt_M->Timing.clockTick0
        + 1) * carPlant_smoothUp_followerSt_M->Timing.stepSize0 +
        carPlant_smoothUp_followerSt_M->Timing.clockTickH0 *
        carPlant_smoothUp_followerSt_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(carPlant_smoothUp_followerSt_M)) {
    carPlant_smoothUp_followerSt_M->Timing.t[0] = rtsiGetT
      (&carPlant_smoothUp_followerSt_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
    /* MinMax: '<S19>/MinMax' incorporates:
     *  Constant: '<S19>/Time constant'
     *  Gain: '<S19>/Minimum sampling to time constant ratio'
     */
    rtb_Integrator_d = fmax
      (carPlant_smoothUp_followerSto_P.FilteredDerivativeDiscreteorC_g *
       carPlant_smoothUp_followerSto_B.Probe[0],
       carPlant_smoothUp_followerSto_P.FilteredDerivativeDiscreteor_hf);

    /* Fcn: '<S19>/Avoid Divide by Zero' */
    carPlant_smoothUp_followerSto_B.AvoidDividebyZero = (real_T)
      (rtb_Integrator_d == 0.0) * 2.2204460492503131e-16 + rtb_Integrator_d;
  }

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *)
      carPlant_smoothUp_followerSt_DW.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      carPlant_smoothUp_followerSt_DW.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex =
      carPlant_smoothUp_followerSt_DW.FromWorkspace_IWORK.PrevIndex;
    real_T t = carPlant_smoothUp_followerSt_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[9989]) {
      currTimeIndex = 9988;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    carPlant_smoothUp_followerSt_DW.FromWorkspace_IWORK.PrevIndex =
      currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          rtb_Product1_o = pDataValues[currTimeIndex];
        } else {
          rtb_Product1_o = pDataValues[currTimeIndex + 1];
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_Product1_o = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 9990;
      }
    }
  }

  /* SignalConversion: '<S378>/ConcatBufferAtVector Concatenate1In1' */
  /* Unit Conversion - from: m/s to: m/s
     Expression: output = (1*input) + (0) */
  carPlant_smoothUp_followerSto_B.VectorConcatenate1[0] = rtb_Product1_o;

  /* SignalConversion: '<S378>/ConcatBufferAtVector ConcatenateIn1' */
  carPlant_smoothUp_followerSto_B.VectorConcatenate[0] = rtb_Product1_o;
  if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
    /* Constant: '<S378>/ydot_o' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate[1] =
      carPlant_smoothUp_followerSto_P.LeadVehicle_ydot_o;

    /* Constant: '<S378>/psi_o' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate[2] =
      carPlant_smoothUp_followerSto_P.LeadVehicle_psi_o;

    /* Constant: '<S378>/r_o' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate[3] =
      carPlant_smoothUp_followerSto_P.LeadVehicle_r_o;
  }

  /* Integrator: '<S378>/Integrator' */
  if (carPlant_smoothUp_followerSt_DW.Integrator_IWORK != 0) {
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE[0] =
      carPlant_smoothUp_followerSto_B.VectorConcatenate[0];
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE[1] =
      carPlant_smoothUp_followerSto_B.VectorConcatenate[1];
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE[2] =
      carPlant_smoothUp_followerSto_B.VectorConcatenate[2];
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE[3] =
      carPlant_smoothUp_followerSto_B.VectorConcatenate[3];
  }

  /* SignalConversion: '<S378>/ConcatBufferAtVector Concatenate1In2' incorporates:
   *  Integrator: '<S378>/Integrator'
   */
  carPlant_smoothUp_followerSto_B.VectorConcatenate1[1] =
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE[1];
  carPlant_smoothUp_followerSto_B.VectorConcatenate1[2] =
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE[2];
  carPlant_smoothUp_followerSto_B.VectorConcatenate1[3] =
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE[3];
  if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
    /* DiscreteIntegrator: '<S22>/Integrator' incorporates:
     *  Constant: '<S5>/Constant'
     */
    if (carPlant_smoothUp_followerSt_DW.Integrator_IC_LOADING != 0) {
      carPlant_smoothUp_followerSt_DW.Integrator_DSTATE =
        carPlant_smoothUp_followerSto_B.VectorConcatenate1[0];
      if (carPlant_smoothUp_followerSt_DW.Integrator_DSTATE >=
          carPlant_smoothUp_followerSto_P.Integrator_UpperSat) {
        carPlant_smoothUp_followerSt_DW.Integrator_DSTATE =
          carPlant_smoothUp_followerSto_P.Integrator_UpperSat;
      } else {
        if (carPlant_smoothUp_followerSt_DW.Integrator_DSTATE <=
            carPlant_smoothUp_followerSto_P.Integrator_LowerSat) {
          carPlant_smoothUp_followerSt_DW.Integrator_DSTATE =
            carPlant_smoothUp_followerSto_P.Integrator_LowerSat;
        }
      }
    }

    if ((carPlant_smoothUp_followerSto_P.Constant_Value_f != 0.0) ||
        (carPlant_smoothUp_followerSt_DW.Integrator_PrevResetState != 0)) {
      carPlant_smoothUp_followerSt_DW.Integrator_DSTATE =
        carPlant_smoothUp_followerSto_B.VectorConcatenate1[0];
      if (carPlant_smoothUp_followerSt_DW.Integrator_DSTATE >=
          carPlant_smoothUp_followerSto_P.Integrator_UpperSat) {
        carPlant_smoothUp_followerSt_DW.Integrator_DSTATE =
          carPlant_smoothUp_followerSto_P.Integrator_UpperSat;
      } else {
        if (carPlant_smoothUp_followerSt_DW.Integrator_DSTATE <=
            carPlant_smoothUp_followerSto_P.Integrator_LowerSat) {
          carPlant_smoothUp_followerSt_DW.Integrator_DSTATE =
            carPlant_smoothUp_followerSto_P.Integrator_LowerSat;
        }
      }
    }

    if (carPlant_smoothUp_followerSt_DW.Integrator_DSTATE >=
        carPlant_smoothUp_followerSto_P.Integrator_UpperSat) {
      carPlant_smoothUp_followerSt_DW.Integrator_DSTATE =
        carPlant_smoothUp_followerSto_P.Integrator_UpperSat;
    } else {
      if (carPlant_smoothUp_followerSt_DW.Integrator_DSTATE <=
          carPlant_smoothUp_followerSto_P.Integrator_LowerSat) {
        carPlant_smoothUp_followerSt_DW.Integrator_DSTATE =
          carPlant_smoothUp_followerSto_P.Integrator_LowerSat;
      }
    }

    /* Saturate: '<S22>/Saturation' incorporates:
     *  DiscreteIntegrator: '<S22>/Integrator'
     */
    if (carPlant_smoothUp_followerSt_DW.Integrator_DSTATE >
        carPlant_smoothUp_followerSto_P.Saturation_UpperSat) {
      carPlant_smoothUp_followerSto_B.Saturation =
        carPlant_smoothUp_followerSto_P.Saturation_UpperSat;
    } else if (carPlant_smoothUp_followerSt_DW.Integrator_DSTATE <
               carPlant_smoothUp_followerSto_P.Saturation_LowerSat) {
      carPlant_smoothUp_followerSto_B.Saturation =
        carPlant_smoothUp_followerSto_P.Saturation_LowerSat;
    } else {
      carPlant_smoothUp_followerSto_B.Saturation =
        carPlant_smoothUp_followerSt_DW.Integrator_DSTATE;
    }

    /* End of Saturate: '<S22>/Saturation' */
  }

  /* Sum: '<S5>/Sum1' */
  rtb_Product1_o = carPlant_smoothUp_followerSto_B.VectorConcatenate1[0] -
    carPlant_smoothUp_followerSto_B.Saturation;

  /* Product: '<S5>/1//T' */
  carPlant_smoothUp_followerSto_B.uT = 1.0 /
    carPlant_smoothUp_followerSto_B.AvoidDividebyZero * rtb_Product1_o;

  /* Gain: '<S5>/Gain' */
  rtb_Product1_o =
    carPlant_smoothUp_followerSto_P.FilteredDerivativeDiscreteor_me *
    carPlant_smoothUp_followerSto_B.uT;

  /* Saturate: '<S5>/[A,B]' */
  if (rtb_Product1_o >
      carPlant_smoothUp_followerSto_P.FilteredDerivativeDiscreteorC_p) {
    carPlant_smoothUp_followerSto_B.AB =
      carPlant_smoothUp_followerSto_P.FilteredDerivativeDiscreteorC_p;
  } else if (rtb_Product1_o <
             carPlant_smoothUp_followerSto_P.FilteredDerivativeDiscreteorCon) {
    carPlant_smoothUp_followerSto_B.AB =
      carPlant_smoothUp_followerSto_P.FilteredDerivativeDiscreteorCon;
  } else {
    carPlant_smoothUp_followerSto_B.AB = rtb_Product1_o;
  }

  /* End of Saturate: '<S5>/[A,B]' */
  if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
    /* MinMax: '<S23>/MinMax' incorporates:
     *  Constant: '<S23>/Time constant'
     *  Gain: '<S23>/Minimum sampling to time constant ratio'
     */
    rtb_Integrator_d = fmax
      (carPlant_smoothUp_followerSto_P.FilteredDerivativeDiscreteorC_d *
       carPlant_smoothUp_followerSto_B.Probe_g[0],
       carPlant_smoothUp_followerSto_P.FilteredDerivativeDiscreteor_md);

    /* Fcn: '<S23>/Avoid Divide by Zero' */
    carPlant_smoothUp_followerSto_B.AvoidDividebyZero_l = (real_T)
      (rtb_Integrator_d == 0.0) * 2.2204460492503131e-16 + rtb_Integrator_d;

    /* Outputs for Atomic SubSystem: '<Root>/Subscribe' */
    /* MATLABSystem: '<S15>/SourceBlock' incorporates:
     *  Inport: '<S382>/In1'
     */
    b_varargout_1 = Sub_carPlant_smoothUp_followerStopper_335.getLatestMessage
      (&rtb_BusAssignment8);

    /* Outputs for Enabled SubSystem: '<S15>/Enabled Subsystem' incorporates:
     *  EnablePort: '<S382>/Enable'
     */
    if (b_varargout_1) {
      carPlant_smoothUp_followerSto_B.In1 = rtb_BusAssignment8;
    }

    /* End of MATLABSystem: '<S15>/SourceBlock' */
    /* End of Outputs for SubSystem: '<S15>/Enabled Subsystem' */
    /* End of Outputs for SubSystem: '<Root>/Subscribe' */

    /* Switch: '<Root>/Switch' incorporates:
     *  Constant: '<Root>/Constant1'
     *  Constant: '<S383>/Constant'
     *  Constant: '<S384>/Constant'
     *  Logic: '<S16>/Logical Operator'
     *  RelationalOperator: '<S383>/Compare'
     *  RelationalOperator: '<S384>/Compare'
     *  Switch: '<S16>/Switch'
     */
    if (carPlant_smoothUp_followerSto_B.In1.Linear.X >
        carPlant_smoothUp_followerSto_P.Switch_Threshold) {
      /* UnitConversion: '<S182>/Unit Conversion6' */
      carPlant_smoothUp_followerSto_B.UnitConversion6 =
        carPlant_smoothUp_followerSto_B.In1.Linear.X;
    } else if ((carPlant_smoothUp_followerSto_P.Constant1_Value >=
                -carPlant_smoothUp_followerSto_P.div0protectpoly_thresh) &&
               (carPlant_smoothUp_followerSto_P.Constant1_Value <=
                carPlant_smoothUp_followerSto_P.div0protectpoly_thresh)) {
      /* Switch: '<S16>/Switch1' incorporates:
       *  Constant: '<Root>/Constant1'
       *  Constant: '<S16>/Constant'
       *  Switch: '<S16>/Switch'
       *  UnaryMinus: '<S16>/Unary Minus'
       */
      if (carPlant_smoothUp_followerSto_P.Constant1_Value >=
          carPlant_smoothUp_followerSto_P.Switch1_Threshold) {
        rtb_Integrator_d = carPlant_smoothUp_followerSto_P.Constant_Value_g;
      } else {
        rtb_Integrator_d = -carPlant_smoothUp_followerSto_P.Constant_Value_g;
      }

      /* End of Switch: '<S16>/Switch1' */

      /* UnitConversion: '<S182>/Unit Conversion6' incorporates:
       *  Constant: '<Root>/Constant1'
       *  Fcn: '<S16>/Fcn'
       *  Product: '<S16>/Product'
       *  Switch: '<S16>/Switch'
       */
      carPlant_smoothUp_followerSto_B.UnitConversion6 = 0.002 / (3.0 -
        rt_powd_snf(carPlant_smoothUp_followerSto_P.Constant1_Value / 0.001, 2.0))
        * rtb_Integrator_d;
    } else {
      /* UnitConversion: '<S182>/Unit Conversion6' incorporates:
       *  Constant: '<Root>/Constant1'
       *  Switch: '<S16>/Switch'
       */
      carPlant_smoothUp_followerSto_B.UnitConversion6 =
        carPlant_smoothUp_followerSto_P.Constant1_Value;
    }

    /* End of Switch: '<Root>/Switch' */
    /* Unit Conversion - from: m/s to: m/s
       Expression: output = (1*input) + (0) */
  }

  /* SignalConversion: '<S182>/ConcatBufferAtVector Concatenate1In1' */
  carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[0] =
    carPlant_smoothUp_followerSto_B.UnitConversion6;
  if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
    /* SignalConversion: '<S182>/ConcatBufferAtVector ConcatenateIn1' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate_g[0] =
      carPlant_smoothUp_followerSto_B.UnitConversion6;

    /* Constant: '<S182>/ydot_o' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate_g[1] =
      carPlant_smoothUp_followerSto_P.FollowerVehicle_ydot_o;

    /* Constant: '<S182>/psi_o' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate_g[2] =
      carPlant_smoothUp_followerSto_P.FollowerVehicle_psi_o;

    /* Constant: '<S182>/r_o' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate_g[3] =
      carPlant_smoothUp_followerSto_P.FollowerVehicle_r_o;
  }

  /* Integrator: '<S182>/Integrator' */
  if (carPlant_smoothUp_followerSt_DW.Integrator_IWORK_h != 0) {
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE_j[0] =
      carPlant_smoothUp_followerSto_B.VectorConcatenate_g[0];
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE_j[1] =
      carPlant_smoothUp_followerSto_B.VectorConcatenate_g[1];
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE_j[2] =
      carPlant_smoothUp_followerSto_B.VectorConcatenate_g[2];
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE_j[3] =
      carPlant_smoothUp_followerSto_B.VectorConcatenate_g[3];
  }

  /* SignalConversion: '<S182>/ConcatBufferAtVector Concatenate1In2' incorporates:
   *  Integrator: '<S182>/Integrator'
   */
  carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[1] =
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE_j[1];
  carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[2] =
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE_j[2];
  carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[3] =
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE_j[3];
  if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
    /* DiscreteIntegrator: '<S26>/Integrator' incorporates:
     *  Constant: '<S6>/Constant'
     */
    if (carPlant_smoothUp_followerSt_DW.Integrator_IC_LOADING_g != 0) {
      carPlant_smoothUp_followerSt_DW.Integrator_DSTATE_m =
        carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[0];
      if (carPlant_smoothUp_followerSt_DW.Integrator_DSTATE_m >=
          carPlant_smoothUp_followerSto_P.Integrator_UpperSat_h) {
        carPlant_smoothUp_followerSt_DW.Integrator_DSTATE_m =
          carPlant_smoothUp_followerSto_P.Integrator_UpperSat_h;
      } else {
        if (carPlant_smoothUp_followerSt_DW.Integrator_DSTATE_m <=
            carPlant_smoothUp_followerSto_P.Integrator_LowerSat_j) {
          carPlant_smoothUp_followerSt_DW.Integrator_DSTATE_m =
            carPlant_smoothUp_followerSto_P.Integrator_LowerSat_j;
        }
      }
    }

    if ((carPlant_smoothUp_followerSto_P.Constant_Value_i != 0.0) ||
        (carPlant_smoothUp_followerSt_DW.Integrator_PrevResetState_c != 0)) {
      carPlant_smoothUp_followerSt_DW.Integrator_DSTATE_m =
        carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[0];
      if (carPlant_smoothUp_followerSt_DW.Integrator_DSTATE_m >=
          carPlant_smoothUp_followerSto_P.Integrator_UpperSat_h) {
        carPlant_smoothUp_followerSt_DW.Integrator_DSTATE_m =
          carPlant_smoothUp_followerSto_P.Integrator_UpperSat_h;
      } else {
        if (carPlant_smoothUp_followerSt_DW.Integrator_DSTATE_m <=
            carPlant_smoothUp_followerSto_P.Integrator_LowerSat_j) {
          carPlant_smoothUp_followerSt_DW.Integrator_DSTATE_m =
            carPlant_smoothUp_followerSto_P.Integrator_LowerSat_j;
        }
      }
    }

    if (carPlant_smoothUp_followerSt_DW.Integrator_DSTATE_m >=
        carPlant_smoothUp_followerSto_P.Integrator_UpperSat_h) {
      carPlant_smoothUp_followerSt_DW.Integrator_DSTATE_m =
        carPlant_smoothUp_followerSto_P.Integrator_UpperSat_h;
    } else {
      if (carPlant_smoothUp_followerSt_DW.Integrator_DSTATE_m <=
          carPlant_smoothUp_followerSto_P.Integrator_LowerSat_j) {
        carPlant_smoothUp_followerSt_DW.Integrator_DSTATE_m =
          carPlant_smoothUp_followerSto_P.Integrator_LowerSat_j;
      }
    }

    /* Saturate: '<S26>/Saturation' incorporates:
     *  DiscreteIntegrator: '<S26>/Integrator'
     */
    if (carPlant_smoothUp_followerSt_DW.Integrator_DSTATE_m >
        carPlant_smoothUp_followerSto_P.Saturation_UpperSat_h) {
      carPlant_smoothUp_followerSto_B.Saturation_l =
        carPlant_smoothUp_followerSto_P.Saturation_UpperSat_h;
    } else if (carPlant_smoothUp_followerSt_DW.Integrator_DSTATE_m <
               carPlant_smoothUp_followerSto_P.Saturation_LowerSat_i) {
      carPlant_smoothUp_followerSto_B.Saturation_l =
        carPlant_smoothUp_followerSto_P.Saturation_LowerSat_i;
    } else {
      carPlant_smoothUp_followerSto_B.Saturation_l =
        carPlant_smoothUp_followerSt_DW.Integrator_DSTATE_m;
    }

    /* End of Saturate: '<S26>/Saturation' */
  }

  /* Sum: '<S6>/Sum1' */
  rtb_Product1_o = carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[0] -
    carPlant_smoothUp_followerSto_B.Saturation_l;

  /* Product: '<S6>/1//T' */
  carPlant_smoothUp_followerSto_B.uT_k = 1.0 /
    carPlant_smoothUp_followerSto_B.AvoidDividebyZero_l * rtb_Product1_o;

  /* Gain: '<S6>/Gain' */
  rtb_Product1_o =
    carPlant_smoothUp_followerSto_P.FilteredDerivativeDiscreteorC_h *
    carPlant_smoothUp_followerSto_B.uT_k;

  /* Saturate: '<S6>/[A,B]' */
  if (rtb_Product1_o >
      carPlant_smoothUp_followerSto_P.FilteredDerivativeDiscreteorC_m) {
    carPlant_smoothUp_followerSto_B.AB_f =
      carPlant_smoothUp_followerSto_P.FilteredDerivativeDiscreteorC_m;
  } else if (rtb_Product1_o <
             carPlant_smoothUp_followerSto_P.FilteredDerivativeDiscreteorC_f) {
    carPlant_smoothUp_followerSto_B.AB_f =
      carPlant_smoothUp_followerSto_P.FilteredDerivativeDiscreteorC_f;
  } else {
    carPlant_smoothUp_followerSto_B.AB_f = rtb_Product1_o;
  }

  /* End of Saturate: '<S6>/[A,B]' */
  if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
    /* Scope: '<Root>/Acceleration' */
    if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
      StructLogVar *svar = (StructLogVar *)
        carPlant_smoothUp_followerSt_DW.Acceleration_PWORK.LoggedData[0];
      LogVar *var = svar->signals.values;

      /* time */
      {
        double locTime = (((carPlant_smoothUp_followerSt_M->Timing.clockTick1+
                            carPlant_smoothUp_followerSt_M->Timing.clockTickH1*
                            4294967296.0)) * 0.01);
        ;
        rt_UpdateLogVar((LogVar *)svar->time, &locTime, 0);
      }

      /* signals */
      {
        real_T up0[1];
        up0[0] = carPlant_smoothUp_followerSto_B.AB;
        rt_UpdateLogVar((LogVar *)var, up0, 0);
        var = var->next;
      }

      {
        real_T up1[1];
        up1[0] = carPlant_smoothUp_followerSto_B.AB_f;
        rt_UpdateLogVar((LogVar *)var, up1, 0);
      }
    }

    /* MATLAB Function: '<Root>/MATLAB Function' incorporates:
     *  Constant: '<S17>/max_accel1'
     *  Constant: '<S17>/max_decel1'
     *  Constant: '<S17>/reference_vel1'
     */
    if (!carPlant_smoothUp_followerSt_DW.y_not_empty) {
      carPlant_smoothUp_followerSt_DW.y =
        carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[0];
      carPlant_smoothUp_followerSt_DW.y_not_empty = true;
    }

    if (carPlant_smoothUp_followerSt_DW.y >
        carPlant_smoothUp_followerSto_P.max_accel1_Value * 0.01 +
        carPlant_smoothUp_followerSto_P.reference_vel1_Value) {
      carPlant_smoothUp_followerSt_DW.y = fmax
        (carPlant_smoothUp_followerSto_P.reference_vel1_Value,
         carPlant_smoothUp_followerSt_DW.y - fabs
         (carPlant_smoothUp_followerSto_P.max_decel1_Value) * 0.01);
    } else if (carPlant_smoothUp_followerSt_DW.y <
               carPlant_smoothUp_followerSto_P.reference_vel1_Value - fabs
               (carPlant_smoothUp_followerSto_P.max_decel1_Value) * 0.01) {
      carPlant_smoothUp_followerSt_DW.y = fmin
        (carPlant_smoothUp_followerSto_P.reference_vel1_Value,
         carPlant_smoothUp_followerSto_P.max_accel1_Value * 0.01 +
         carPlant_smoothUp_followerSt_DW.y);
    } else {
      carPlant_smoothUp_followerSt_DW.y =
        carPlant_smoothUp_followerSto_P.reference_vel1_Value;
    }

    /* BusAssignment: '<Root>/Bus Assignment7' incorporates:
     *  Constant: '<Root>/Steering angle'
     *  Constant: '<S2>/Constant'
     *  MATLAB Function: '<Root>/MATLAB Function'
     */
    carPlant_smoothUp_followerSto_B.BusAssignment7 =
      carPlant_smoothUp_followerSto_P.Constant_Value_k;
    carPlant_smoothUp_followerSto_B.BusAssignment7.Linear.X =
      carPlant_smoothUp_followerSt_DW.y;
    carPlant_smoothUp_followerSto_B.BusAssignment7.Angular.Z =
      carPlant_smoothUp_followerSto_P.Steeringangle_Value;

    /* Outputs for Atomic SubSystem: '<Root>/Dead Man's Switch' */
    /* MATLAB Function: '<S4>/timeout set to 0 output' incorporates:
     *  Constant: '<Root>/Constant12'
     */
    b_varargout_1 = true;
    if (!carPlant_smoothUp_followerSt_DW.sinceLastMsg_not_empty) {
      carPlant_smoothUp_followerSt_DW.sinceLastMsg =
        carPlant_smoothUp_followerSto_P.DeadMansSwitch_timeout /
        carPlant_smoothUp_followerSto_P.DeadMansSwitch_stepSize + 1.0;
      carPlant_smoothUp_followerSt_DW.sinceLastMsg_not_empty = true;
    }

    carPlant_smoothUp_followerSto_B.safeValue = 0.0;
    if (carPlant_smoothUp_followerSto_P.Constant12_Value == 1.0) {
      carPlant_smoothUp_followerSt_DW.sinceLastMsg = 0.0;
    } else {
      carPlant_smoothUp_followerSt_DW.sinceLastMsg++;
    }

    if (carPlant_smoothUp_followerSt_DW.sinceLastMsg <
        carPlant_smoothUp_followerSto_P.DeadMansSwitch_timeout /
        carPlant_smoothUp_followerSto_P.DeadMansSwitch_stepSize) {
      b_varargout_1 = false;
    }

    if (!b_varargout_1) {
      carPlant_smoothUp_followerSto_B.safeValue =
        carPlant_smoothUp_followerSto_B.BusAssignment7.Linear.X;
    }

    /* End of MATLAB Function: '<S4>/timeout set to 0 output' */
    /* End of Outputs for SubSystem: '<Root>/Dead Man's Switch' */

    /* Constant: '<S290>/X_o' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate3[0] =
      carPlant_smoothUp_followerSto_P.LeadVehicle_X_o;

    /* Constant: '<S290>/Y_o' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate3[1] =
      carPlant_smoothUp_followerSto_P.LeadVehicle_Y_o;
  }

  /* Integrator: '<S290>/Integrator' */
  if (carPlant_smoothUp_followerSt_DW.Integrator_IWORK_p != 0) {
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE_o[0] =
      carPlant_smoothUp_followerSto_B.VectorConcatenate3[0];
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE_o[1] =
      carPlant_smoothUp_followerSto_B.VectorConcatenate3[1];
  }

  /* Trigonometry: '<S314>/sincos' incorporates:
   *  Constant: '<S290>/Constant2'
   *  Constant: '<S290>/Constant7'
   *  SignalConversion: '<S297>/TmpSignal ConversionAtsincosInport1'
   *  Trigonometry: '<S242>/Trigonometric Function'
   */
  rtb_UnaryMinus1_g_idx_0_tmp = sin
    (carPlant_smoothUp_followerSto_B.VectorConcatenate1[2]);
  rtb_sincos_o2_idx_0_tmp = cos
    (carPlant_smoothUp_followerSto_B.VectorConcatenate1[2]);
  rtb_Integrator_d = sin(carPlant_smoothUp_followerSto_P.Constant7_Value);
  rtb_sincos_o2_idx_1 = cos(carPlant_smoothUp_followerSto_P.Constant7_Value);
  rtb_UnaryMinus1_g_idx_2 = sin(carPlant_smoothUp_followerSto_P.Constant2_Value);
  rtb_sincos_o2_idx_2 = cos(carPlant_smoothUp_followerSto_P.Constant2_Value);

  /* Fcn: '<S314>/Fcn11' incorporates:
   *  Trigonometry: '<S314>/sincos'
   */
  rtb_VectorConcatenate_al[0] = rtb_sincos_o2_idx_1 * rtb_sincos_o2_idx_0_tmp;

  /* Fcn: '<S314>/Fcn21' incorporates:
   *  Fcn: '<S314>/Fcn22'
   *  Trigonometry: '<S314>/sincos'
   */
  dx3 = rtb_UnaryMinus1_g_idx_2 * rtb_Integrator_d;
  rtb_VectorConcatenate_al[1] = dx3 * rtb_sincos_o2_idx_0_tmp -
    rtb_sincos_o2_idx_2 * rtb_UnaryMinus1_g_idx_0_tmp;

  /* Fcn: '<S314>/Fcn31' incorporates:
   *  Fcn: '<S314>/Fcn32'
   *  Trigonometry: '<S314>/sincos'
   */
  dx1 = rtb_sincos_o2_idx_2 * rtb_Integrator_d;
  rtb_VectorConcatenate_al[2] = dx1 * rtb_sincos_o2_idx_0_tmp +
    rtb_UnaryMinus1_g_idx_2 * rtb_UnaryMinus1_g_idx_0_tmp;

  /* Fcn: '<S314>/Fcn12' incorporates:
   *  Trigonometry: '<S314>/sincos'
   */
  rtb_VectorConcatenate_al[3] = rtb_sincos_o2_idx_1 *
    rtb_UnaryMinus1_g_idx_0_tmp;

  /* Fcn: '<S314>/Fcn22' incorporates:
   *  Trigonometry: '<S314>/sincos'
   */
  rtb_VectorConcatenate_al[4] = dx3 * rtb_UnaryMinus1_g_idx_0_tmp +
    rtb_sincos_o2_idx_2 * rtb_sincos_o2_idx_0_tmp;

  /* Fcn: '<S314>/Fcn32' incorporates:
   *  Trigonometry: '<S314>/sincos'
   */
  rtb_VectorConcatenate_al[5] = dx1 * rtb_UnaryMinus1_g_idx_0_tmp -
    rtb_UnaryMinus1_g_idx_2 * rtb_sincos_o2_idx_0_tmp;

  /* Fcn: '<S314>/Fcn13' */
  rtb_VectorConcatenate_al[6] = -rtb_Integrator_d;

  /* Fcn: '<S314>/Fcn23' */
  rtb_VectorConcatenate_al[7] = rtb_UnaryMinus1_g_idx_2 * rtb_sincos_o2_idx_1;

  /* Fcn: '<S314>/Fcn33' */
  rtb_VectorConcatenate_al[8] = rtb_sincos_o2_idx_2 * rtb_sincos_o2_idx_1;
  if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
    /* Constant: '<S293>/longOff' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate_l[0] =
      carPlant_smoothUp_followerSto_P.LeadVehicle_longOff;

    /* Constant: '<S293>/latOff' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate_l[1] =
      carPlant_smoothUp_followerSto_P.LeadVehicle_latOff;

    /* Constant: '<S293>/vertOff' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate_l[2] =
      carPlant_smoothUp_followerSto_P.LeadVehicle_vertOff;
  }

  /* Sum: '<S312>/Add' incorporates:
   *  Constant: '<S290>/Constant'
   *  Integrator: '<S290>/Integrator'
   */
  tmp[0] = carPlant_smoothUp_followerSto_X.Integrator_CSTATE_o[0];
  tmp[1] = carPlant_smoothUp_followerSto_X.Integrator_CSTATE_o[1];
  tmp[2] = carPlant_smoothUp_followerSto_P.Constant_Value_o;
  for (iU = 0; iU < 3; iU++) {
    /* Sum: '<S312>/Add' incorporates:
     *  Math: '<S312>/Transpose1'
     *  Product: '<S316>/Product'
     */
    carPlant_smoothUp_followerSto_B.Add[iU] = tmp[iU] +
      (rtb_VectorConcatenate_al[3 * iU + 2] *
       carPlant_smoothUp_followerSto_B.VectorConcatenate_l[2] +
       (rtb_VectorConcatenate_al[3 * iU + 1] *
        carPlant_smoothUp_followerSto_B.VectorConcatenate_l[1] +
        rtb_VectorConcatenate_al[3 * iU] *
        carPlant_smoothUp_followerSto_B.VectorConcatenate_l[0]));
  }

  if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
    /* Constant: '<S94>/X_o' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate3_a[0] =
      carPlant_smoothUp_followerSto_P.FollowerVehicle_X_o;

    /* Constant: '<S94>/Y_o' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate3_a[1] =
      carPlant_smoothUp_followerSto_P.FollowerVehicle_Y_o;
  }

  /* Integrator: '<S94>/Integrator' */
  if (carPlant_smoothUp_followerSt_DW.Integrator_IWORK_l != 0) {
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE_a[0] =
      carPlant_smoothUp_followerSto_B.VectorConcatenate3_a[0];
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE_a[1] =
      carPlant_smoothUp_followerSto_B.VectorConcatenate3_a[1];
  }

  /* Trigonometry: '<S118>/sincos' incorporates:
   *  Constant: '<S94>/Constant2'
   *  Constant: '<S94>/Constant7'
   *  SignalConversion: '<S101>/TmpSignal ConversionAtsincosInport1'
   *  Trigonometry: '<S46>/Trigonometric Function'
   */
  rtb_sincos_o1_idx_0 = sin
    (carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[2]);
  rtb_VectorConcatenate1_m_idx_0_ = cos
    (carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[2]);
  rtb_Integrator_d = sin(carPlant_smoothUp_followerSto_P.Constant7_Value_h);
  rtb_sincos_o2_idx_1 = cos(carPlant_smoothUp_followerSto_P.Constant7_Value_h);
  rtb_UnaryMinus1_g_idx_2 = sin
    (carPlant_smoothUp_followerSto_P.Constant2_Value_l);
  rtb_sincos_o2_idx_2 = cos(carPlant_smoothUp_followerSto_P.Constant2_Value_l);

  /* Fcn: '<S118>/Fcn11' incorporates:
   *  Trigonometry: '<S118>/sincos'
   */
  rtb_VectorConcatenate_al[0] = rtb_sincos_o2_idx_1 *
    rtb_VectorConcatenate1_m_idx_0_;

  /* Fcn: '<S118>/Fcn21' incorporates:
   *  Fcn: '<S118>/Fcn22'
   *  Trigonometry: '<S118>/sincos'
   */
  dx3 = rtb_UnaryMinus1_g_idx_2 * rtb_Integrator_d;
  rtb_VectorConcatenate_al[1] = dx3 * rtb_VectorConcatenate1_m_idx_0_ -
    rtb_sincos_o2_idx_2 * rtb_sincos_o1_idx_0;

  /* Fcn: '<S118>/Fcn31' incorporates:
   *  Fcn: '<S118>/Fcn32'
   *  Trigonometry: '<S118>/sincos'
   */
  dx1 = rtb_sincos_o2_idx_2 * rtb_Integrator_d;
  rtb_VectorConcatenate_al[2] = dx1 * rtb_VectorConcatenate1_m_idx_0_ +
    rtb_UnaryMinus1_g_idx_2 * rtb_sincos_o1_idx_0;

  /* Fcn: '<S118>/Fcn12' incorporates:
   *  Trigonometry: '<S118>/sincos'
   */
  rtb_VectorConcatenate_al[3] = rtb_sincos_o2_idx_1 * rtb_sincos_o1_idx_0;

  /* Fcn: '<S118>/Fcn22' incorporates:
   *  Trigonometry: '<S118>/sincos'
   */
  rtb_VectorConcatenate_al[4] = dx3 * rtb_sincos_o1_idx_0 + rtb_sincos_o2_idx_2 *
    rtb_VectorConcatenate1_m_idx_0_;

  /* Fcn: '<S118>/Fcn32' incorporates:
   *  Trigonometry: '<S118>/sincos'
   */
  rtb_VectorConcatenate_al[5] = dx1 * rtb_sincos_o1_idx_0 -
    rtb_UnaryMinus1_g_idx_2 * rtb_VectorConcatenate1_m_idx_0_;

  /* Fcn: '<S118>/Fcn13' */
  rtb_VectorConcatenate_al[6] = -rtb_Integrator_d;

  /* Fcn: '<S118>/Fcn23' */
  rtb_VectorConcatenate_al[7] = rtb_UnaryMinus1_g_idx_2 * rtb_sincos_o2_idx_1;

  /* Fcn: '<S118>/Fcn33' */
  rtb_VectorConcatenate_al[8] = rtb_sincos_o2_idx_2 * rtb_sincos_o2_idx_1;
  if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
    /* Constant: '<S97>/longOff' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate_n[0] =
      carPlant_smoothUp_followerSto_P.FollowerVehicle_longOff;

    /* Constant: '<S97>/latOff' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate_n[1] =
      carPlant_smoothUp_followerSto_P.FollowerVehicle_latOff;

    /* Constant: '<S97>/vertOff' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate_n[2] =
      carPlant_smoothUp_followerSto_P.FollowerVehicle_vertOff;
  }

  /* Sum: '<S116>/Add' incorporates:
   *  Constant: '<S94>/Constant'
   *  Integrator: '<S94>/Integrator'
   */
  tmp[0] = carPlant_smoothUp_followerSto_X.Integrator_CSTATE_a[0];
  tmp[1] = carPlant_smoothUp_followerSto_X.Integrator_CSTATE_a[1];
  tmp[2] = carPlant_smoothUp_followerSto_P.Constant_Value_e;
  for (iU = 0; iU < 3; iU++) {
    /* Sum: '<S116>/Add' incorporates:
     *  Math: '<S116>/Transpose1'
     *  Product: '<S120>/Product'
     */
    carPlant_smoothUp_followerSto_B.Add_b[iU] = tmp[iU] +
      (rtb_VectorConcatenate_al[3 * iU + 2] *
       carPlant_smoothUp_followerSto_B.VectorConcatenate_n[2] +
       (rtb_VectorConcatenate_al[3 * iU + 1] *
        carPlant_smoothUp_followerSto_B.VectorConcatenate_n[1] +
        rtb_VectorConcatenate_al[3 * iU] *
        carPlant_smoothUp_followerSto_B.VectorConcatenate_n[0]));
  }

  /* Sum: '<Root>/Subtract2' */
  carPlant_smoothUp_followerSto_B.relativeposition =
    carPlant_smoothUp_followerSto_B.Add[0] -
    carPlant_smoothUp_followerSto_B.Add_b[0];

  /* Sum: '<Root>/Subtract' */
  carPlant_smoothUp_followerSto_B.Subtract =
    carPlant_smoothUp_followerSto_B.VectorConcatenate1[0] -
    carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[0];

  /* MATLABSystem: '<S13>/Get Parameter1' */
  ParamGet_carPlant_smoothUp_followerStopper_99.get_parameter
    (&rtb_UnaryMinus1_g_idx_2);

  /* MATLABSystem: '<S13>/Get Parameter2' */
  ParamGet_carPlant_smoothUp_followerStopper_100.get_parameter(&dx3);

  /* MATLAB Function: '<S7>/MATLAB Function' incorporates:
   *  Constant: '<S13>/decel'
   *  MATLABSystem: '<S13>/Get Parameter1'
   *  MATLABSystem: '<S13>/Get Parameter2'
   */
  rtb_Integrator_d = fmin(carPlant_smoothUp_followerSto_B.safeValue, fmax
    (carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[0] +
     carPlant_smoothUp_followerSto_B.Subtract, 0.0));
  rtb_sincos_o2_idx_2 = fmin(carPlant_smoothUp_followerSto_B.Subtract, 0.0);
  rtb_sincos_o2_idx_1 = rtb_sincos_o2_idx_2 * rtb_sincos_o2_idx_2;
  dx1 = 1.0 / (2.0 * carPlant_smoothUp_followerSto_P.decel_Value[0]) *
    rtb_sincos_o2_idx_1 + rtb_UnaryMinus1_g_idx_2;
  rtb_UnaryMinus1_g_idx_2 = 1.0 / (2.0 *
    carPlant_smoothUp_followerSto_P.decel_Value[1]) * rtb_sincos_o2_idx_1 +
    (rtb_UnaryMinus1_g_idx_2 + dx3) / 2.0;
  dx3 += 1.0 / (2.0 * carPlant_smoothUp_followerSto_P.decel_Value[2]) *
    rtb_sincos_o2_idx_1;
  if (carPlant_smoothUp_followerSto_B.relativeposition <= dx1) {
    rtb_Integrator_d = 0.0;
  } else if (carPlant_smoothUp_followerSto_B.relativeposition <=
             rtb_UnaryMinus1_g_idx_2) {
    rtb_Integrator_d = fmin((carPlant_smoothUp_followerSto_B.relativeposition -
      dx1) * rtb_Integrator_d / (rtb_UnaryMinus1_g_idx_2 - dx1),
      carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[0] + 0.015);
  } else if (carPlant_smoothUp_followerSto_B.relativeposition <= dx3) {
    rtb_Integrator_d = fmin((carPlant_smoothUp_followerSto_B.safeValue -
      rtb_Integrator_d) * (carPlant_smoothUp_followerSto_B.relativeposition -
      rtb_UnaryMinus1_g_idx_2) / (dx3 - rtb_UnaryMinus1_g_idx_2) +
      rtb_Integrator_d, carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[0]
      + 0.015);
  } else {
    rtb_Integrator_d = fmin(carPlant_smoothUp_followerSto_B.safeValue,
      carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[0] + 0.015);
  }

  /* End of MATLAB Function: '<S7>/MATLAB Function' */

  /* BusAssignment: '<Root>/Bus Assignment8' incorporates:
   *  Constant: '<S3>/Constant'
   */
  rtb_BusAssignment8 = carPlant_smoothUp_followerSto_P.Constant_Value_b;
  rtb_BusAssignment8.Linear.X = rtb_Integrator_d;
  rtb_BusAssignment8.Angular.Z =
    carPlant_smoothUp_followerSto_B.BusAssignment7.Angular.Z;

  /* Outputs for Atomic SubSystem: '<Root>/Publish' */
  /* MATLABSystem: '<S14>/SinkBlock' */
  Pub_carPlant_smoothUp_followerStopper_334.publish(&rtb_BusAssignment8);

  /* End of Outputs for SubSystem: '<Root>/Publish' */
  if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
    /* Scope: '<Root>/Velocity' */
    if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
      StructLogVar *svar = (StructLogVar *)
        carPlant_smoothUp_followerSt_DW.Velocity_PWORK.LoggedData[0];
      LogVar *var = svar->signals.values;

      /* time */
      {
        double locTime = (((carPlant_smoothUp_followerSt_M->Timing.clockTick1+
                            carPlant_smoothUp_followerSt_M->Timing.clockTickH1*
                            4294967296.0)) * 0.01);
        ;
        rt_UpdateLogVar((LogVar *)svar->time, &locTime, 0);
      }

      /* signals */
      {
        real_T up0[1];
        up0[0] = carPlant_smoothUp_followerSto_B.VectorConcatenate1[0];
        rt_UpdateLogVar((LogVar *)var, up0, 0);
        var = var->next;
      }

      {
        real_T up1[1];
        up1[0] = carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[0];
        rt_UpdateLogVar((LogVar *)var, up1, 0);
      }
    }

    /* Scope: '<Root>/X-Position' */
    if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
      StructLogVar *svar = (StructLogVar *)
        carPlant_smoothUp_followerSt_DW.XPosition_PWORK.LoggedData[0];
      LogVar *var = svar->signals.values;

      /* time */
      {
        double locTime = (((carPlant_smoothUp_followerSt_M->Timing.clockTick1+
                            carPlant_smoothUp_followerSt_M->Timing.clockTickH1*
                            4294967296.0)) * 0.01);
        ;
        rt_UpdateLogVar((LogVar *)svar->time, &locTime, 0);
      }

      /* signals */
      {
        real_T up0[1];
        up0[0] = carPlant_smoothUp_followerSto_B.Add[0];
        rt_UpdateLogVar((LogVar *)var, up0, 0);
        var = var->next;
      }

      {
        real_T up1[1];
        up1[0] = carPlant_smoothUp_followerSto_B.Add_b[0];
        rt_UpdateLogVar((LogVar *)var, up1, 0);
      }
    }

    /* Scope: '<Root>/Y-Position' */
    if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
      StructLogVar *svar = (StructLogVar *)
        carPlant_smoothUp_followerSt_DW.YPosition_PWORK.LoggedData[0];
      LogVar *var = svar->signals.values;

      /* time */
      {
        double locTime = (((carPlant_smoothUp_followerSt_M->Timing.clockTick1+
                            carPlant_smoothUp_followerSt_M->Timing.clockTickH1*
                            4294967296.0)) * 0.01);
        ;
        rt_UpdateLogVar((LogVar *)svar->time, &locTime, 0);
      }

      /* signals */
      {
        real_T up0[1];
        up0[0] = carPlant_smoothUp_followerSto_B.Add[1];
        rt_UpdateLogVar((LogVar *)var, up0, 0);
        var = var->next;
      }

      {
        real_T up1[1];
        up1[0] = carPlant_smoothUp_followerSto_B.Add_b[1];
        rt_UpdateLogVar((LogVar *)var, up1, 0);
      }
    }

    /* SignalConversion: '<S185>/ConcatBufferAtVector Concatenate1In1' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate1_b[0] = 0.0;

    /* SignalConversion: '<S185>/ConcatBufferAtVector Concatenate1In2' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate1_b[1] = 0.0;

    /* SignalConversion: '<S185>/ConcatBufferAtVector Concatenate1In3' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate1_b[2] = 0.0;
  }

  /* Trigonometry: '<S46>/Trigonometric Function' */
  rtb_Product1_o = rtb_sincos_o1_idx_0;

  /* Sum: '<S46>/Add' incorporates:
   *  Product: '<S46>/Product1'
   *  Product: '<S46>/Product3'
   */
  rtb_sincos_o1_idx_0 = rtb_VectorConcatenate1_m_idx_0_ *
    carPlant_smoothUp_followerSto_B.VectorConcatenate1_b[0] + rtb_Product1_o *
    carPlant_smoothUp_followerSto_B.VectorConcatenate1_b[1];

  /* Product: '<S46>/Product' */
  rtb_Product1_o *= carPlant_smoothUp_followerSto_B.VectorConcatenate1_b[0];

  /* Sum: '<S45>/Add1' incorporates:
   *  SignalConversion: '<S29>/ConcatBufferAtVector Concatenate1In1'
   */
  rtb_VectorConcatenate1_k =
    carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[0] -
    rtb_sincos_o1_idx_0;
  rtb_VectorConcatenate1_m_idx_0 = rtb_VectorConcatenate1_k;

  /* Product: '<S45>/Product' */
  rtb_sincos_o1_idx_0 = rtb_VectorConcatenate1_k * rtb_VectorConcatenate1_k;

  /* Sum: '<S45>/Add1' incorporates:
   *  Product: '<S46>/Product2'
   *  SignalConversion: '<S29>/ConcatBufferAtVector Concatenate1In2'
   *  Sum: '<S46>/Add1'
   */
  rtb_VectorConcatenate1_k =
    carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[1] -
    (rtb_VectorConcatenate1_m_idx_0_ *
     carPlant_smoothUp_followerSto_B.VectorConcatenate1_b[1] - rtb_Product1_o);

  /* Sum: '<S45>/Sum of Elements' */
  iU = 0;

  /* Sqrt: '<S45>/Sqrt' incorporates:
   *  Product: '<S45>/Product'
   *  SignalConversion: '<S46>/ConcatBufferAtVector Concatenate2In3'
   *  Sum: '<S45>/Add1'
   *  Sum: '<S45>/Sum of Elements'
   */
  rtb_Integrator_d = sqrt((rtb_VectorConcatenate1_k * rtb_VectorConcatenate1_k +
    rtb_sincos_o1_idx_0) + (0.0 -
    carPlant_smoothUp_followerSto_B.VectorConcatenate1_b[2]) * (0.0 -
    carPlant_smoothUp_followerSto_B.VectorConcatenate1_b[2]));

  /* Product: '<S45>/Product2' */
  dx1 = rtb_Integrator_d * rtb_Integrator_d;
  if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
    /* Constant: '<S45>/Constant' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate_a[0] =
      carPlant_smoothUp_followerSto_P.FollowerVehicle_Cd;
  }

  /* Trigonometry: '<S45>/Trigonometric Function' */
  rtb_Integrator_d = rt_atan2d_snf(rtb_VectorConcatenate1_k,
    rtb_VectorConcatenate1_m_idx_0);

  /* Lookup_n-D: '<S45>/Cs' */
  carPlant_smoothUp_followerSto_B.VectorConcatenate_a[1] = look1_binlcpw
    (rtb_Integrator_d, carPlant_smoothUp_followerSto_P.FollowerVehicle_beta_w,
     carPlant_smoothUp_followerSto_P.FollowerVehicle_Cs, 30U);
  if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
    /* Constant: '<S45>/Constant1' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate_a[2] =
      carPlant_smoothUp_followerSto_P.FollowerVehicle_Cl;

    /* UnaryMinus: '<S45>/Unary Minus' incorporates:
     *  Constant: '<S45>/Constant4'
     */
    carPlant_smoothUp_followerSto_B.UnaryMinus[0] =
      -carPlant_smoothUp_followerSto_P.Constant4_Value_d[0];
    carPlant_smoothUp_followerSto_B.UnaryMinus[1] =
      -carPlant_smoothUp_followerSto_P.Constant4_Value_d[1];
    carPlant_smoothUp_followerSto_B.UnaryMinus[2] =
      -carPlant_smoothUp_followerSto_P.Constant4_Value_d[2];
  }

  /* Lookup_n-D: '<S45>/Crm' */
  carPlant_smoothUp_followerSto_B.VectorConcatenate_a[3] = look1_binlxpw
    (rtb_Integrator_d, carPlant_smoothUp_followerSto_P.Crm_bp01Data,
     carPlant_smoothUp_followerSto_P.Crm_tableData, 1U);

  /* Switch: '<S45>/Switch' incorporates:
   *  Constant: '<S45>/Constant4'
   */
  if (rtb_VectorConcatenate1_m_idx_0 >=
      carPlant_smoothUp_followerSto_P.Switch_Threshold_f) {
    rtb_VectorConcatenate1_m_idx_0 =
      carPlant_smoothUp_followerSto_P.Constant4_Value_d[0];
  } else {
    rtb_VectorConcatenate1_m_idx_0 = carPlant_smoothUp_followerSto_B.UnaryMinus
      [0];
  }

  /* Product: '<S45>/Product5' incorporates:
   *  Constant: '<S45>/Constant2'
   */
  carPlant_smoothUp_followerSto_B.VectorConcatenate_a[4] =
    rtb_VectorConcatenate1_m_idx_0 *
    carPlant_smoothUp_followerSto_P.FollowerVehicle_Cpm;

  /* Lookup_n-D: '<S45>/Cym' */
  carPlant_smoothUp_followerSto_B.VectorConcatenate_a[5] = look1_binlxpw
    (rtb_Integrator_d, carPlant_smoothUp_followerSto_P.FollowerVehicle_beta_w,
     carPlant_smoothUp_followerSto_P.FollowerVehicle_Cym, 30U);

  /* Gain: '<S45>/.5.*A.*Pabs.//R.//T' incorporates:
   *  Product: '<S45>/Product1'
   */
  rtb_Integrator_d = 0.5 * carPlant_smoothUp_followerSto_P.FollowerVehicle_Af *
    carPlant_smoothUp_followerSto_P.FollowerVehicle_Pabs /
    carPlant_smoothUp_followerSto_P.DragForce_R /
    carPlant_smoothUp_followerSto_P.FollowerVehicle_Tair;
  for (i = 0; i < 6; i++) {
    rtb_uAPabsRT[i] = dx1 *
      carPlant_smoothUp_followerSto_B.VectorConcatenate_a[i] * rtb_Integrator_d;
  }

  /* End of Gain: '<S45>/.5.*A.*Pabs.//R.//T' */

  /* Product: '<S45>/Product4' incorporates:
   *  Constant: '<S45>/Constant3'
   *  MATLAB Function: '<S8>/vehicle model'
   */
  rtb_VectorConcatenate1_m_idx_0_ =
    carPlant_smoothUp_followerSto_P.FollowerVehicle_a +
    carPlant_smoothUp_followerSto_P.FollowerVehicle_b;

  /* Sum: '<S8>/Add' incorporates:
   *  Constant: '<S45>/Constant3'
   *  Product: '<S45>/Product4'
   *  UnaryMinus: '<S29>/Unary Minus1'
   */
  rtb_sincos_o1_idx_0 = -(rtb_uAPabsRT[3] * rtb_VectorConcatenate1_m_idx_0_);

  /* Product: '<S45>/Product3' incorporates:
   *  UnaryMinus: '<S29>/Unary Minus'
   */
  rtb_VectorConcatenate1_m_idx_0 = -(rtb_VectorConcatenate1_m_idx_0 *
    rtb_uAPabsRT[0]);

  /* Sum: '<S8>/Add' incorporates:
   *  Constant: '<S45>/Constant3'
   *  Product: '<S45>/Product4'
   *  UnaryMinus: '<S29>/Unary Minus1'
   */
  rtb_sincos_o1_idx_1 = -(rtb_uAPabsRT[4] * rtb_VectorConcatenate1_m_idx_0_);

  /* Switch: '<S45>/Switch' incorporates:
   *  Constant: '<S45>/Constant4'
   */
  if (rtb_VectorConcatenate1_k >=
      carPlant_smoothUp_followerSto_P.Switch_Threshold_f) {
    rtb_VectorConcatenate1_k =
      carPlant_smoothUp_followerSto_P.Constant4_Value_d[1];
  } else {
    rtb_VectorConcatenate1_k = carPlant_smoothUp_followerSto_B.UnaryMinus[1];
  }

  /* Product: '<S45>/Product3' incorporates:
   *  UnaryMinus: '<S29>/Unary Minus'
   */
  rtb_sincos_o2_idx_1 = -(rtb_VectorConcatenate1_k * rtb_uAPabsRT[1]);

  /* Sum: '<S8>/Add' incorporates:
   *  Constant: '<S45>/Constant3'
   *  Product: '<S45>/Product4'
   *  UnaryMinus: '<S29>/Unary Minus1'
   */
  rtb_sincos_o1_idx_2 = -(rtb_uAPabsRT[5] * rtb_VectorConcatenate1_m_idx_0_);

  /* Switch: '<S45>/Switch' incorporates:
   *  Constant: '<S45>/Constant4'
   *  SignalConversion: '<S46>/ConcatBufferAtVector Concatenate2In3'
   *  Sum: '<S45>/Add1'
   */
  if (0.0 - carPlant_smoothUp_followerSto_B.VectorConcatenate1_b[2] >=
      carPlant_smoothUp_followerSto_P.Switch_Threshold_f) {
    rtb_Integrator_d = carPlant_smoothUp_followerSto_P.Constant4_Value_d[2];
  } else {
    rtb_Integrator_d = carPlant_smoothUp_followerSto_B.UnaryMinus[2];
  }

  /* UnaryMinus: '<S29>/Unary Minus' incorporates:
   *  Product: '<S45>/Product3'
   */
  rtb_VectorConcatenate1_k = -(rtb_Integrator_d * rtb_uAPabsRT[2]);
  if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
    /* SignalConversion: '<S41>/ConcatBufferAtVector Concatenate2In1' incorporates:
     *  Constant: '<S41>/Cyf'
     */
    carPlant_smoothUp_followerSto_B.VectorConcatenate2[0] =
      carPlant_smoothUp_followerSto_P.FollowerVehicle_Cy_f;

    /* SignalConversion: '<S41>/ConcatBufferAtVector Concatenate2In2' incorporates:
     *  Constant: '<S41>/Cyf'
     */
    carPlant_smoothUp_followerSto_B.VectorConcatenate2[1] =
      carPlant_smoothUp_followerSto_P.FollowerVehicle_Cy_f;

    /* SignalConversion: '<S41>/ConcatBufferAtVector Concatenate1In1' incorporates:
     *  Constant: '<S41>/Cyr'
     */
    carPlant_smoothUp_followerSto_B.VectorConcatenate1_bi[0] =
      carPlant_smoothUp_followerSto_P.FollowerVehicle_Cy_r;

    /* SignalConversion: '<S41>/ConcatBufferAtVector Concatenate1In2' incorporates:
     *  Constant: '<S41>/Cyr'
     */
    carPlant_smoothUp_followerSto_B.VectorConcatenate1_bi[1] =
      carPlant_smoothUp_followerSto_P.FollowerVehicle_Cy_r;

    /* Backlash: '<S205>/Backlash' incorporates:
     *  Constant: '<Root>/Constant4'
     */
    rtb_Integrator_d = carPlant_smoothUp_followerSto_P.KinematicSteering1_Db /
      2.0;
    if (carPlant_smoothUp_followerSto_P.Constant4_Value <
        carPlant_smoothUp_followerSt_DW.PrevY - rtb_Integrator_d) {
      carPlant_smoothUp_followerSto_B.Backlash =
        carPlant_smoothUp_followerSto_P.Constant4_Value + rtb_Integrator_d;
    } else if (carPlant_smoothUp_followerSto_P.Constant4_Value <=
               carPlant_smoothUp_followerSt_DW.PrevY + rtb_Integrator_d) {
      carPlant_smoothUp_followerSto_B.Backlash =
        carPlant_smoothUp_followerSt_DW.PrevY;
    } else {
      carPlant_smoothUp_followerSto_B.Backlash =
        carPlant_smoothUp_followerSto_P.Constant4_Value - rtb_Integrator_d;
    }

    /* End of Backlash: '<S205>/Backlash' */

    /* Saturate: '<S205>/Saturation' */
    if (carPlant_smoothUp_followerSto_B.Backlash >
        carPlant_smoothUp_followerSto_P.KinematicSteering1_StrgRng) {
      rtb_Integrator_d =
        carPlant_smoothUp_followerSto_P.KinematicSteering1_StrgRng;
    } else if (carPlant_smoothUp_followerSto_B.Backlash <
               -carPlant_smoothUp_followerSto_P.KinematicSteering1_StrgRng) {
      rtb_Integrator_d =
        -carPlant_smoothUp_followerSto_P.KinematicSteering1_StrgRng;
    } else {
      rtb_Integrator_d = carPlant_smoothUp_followerSto_B.Backlash;
    }

    /* End of Saturate: '<S205>/Saturation' */

    /* Trigonometry: '<S207>/Trigonometric Function' incorporates:
     *  Gain: '<S206>/Gain'
     */
    rtb_sincos_o2_idx_2 = tan(1.0 /
      carPlant_smoothUp_followerSto_P.KinematicSteering1_StrgRatio *
      rtb_Integrator_d);

    /* Product: '<S207>/Divide5' incorporates:
     *  Constant: '<S207>/Constant1'
     */
    rtb_Integrator_d =
      carPlant_smoothUp_followerSto_P.KinematicSteering1_TrckWdth / 2.0 *
      rtb_sincos_o2_idx_2;

    /* Sum: '<S207>/Add' incorporates:
     *  Constant: '<S207>/Constant'
     */
    dx1 = carPlant_smoothUp_followerSto_P.KinematicSteering1_WhlBase +
      rtb_Integrator_d;

    /* Switch: '<S212>/Switch' incorporates:
     *  Abs: '<S212>/Abs'
     *  Constant: '<S219>/Constant'
     *  Constant: '<S220>/Constant'
     *  Fcn: '<S212>/Fcn'
     *  Logic: '<S212>/Logical Operator'
     *  RelationalOperator: '<S219>/Compare'
     *  RelationalOperator: '<S220>/Compare'
     */
    if ((dx1 >= -carPlant_smoothUp_followerSto_P.div0protectabspoly3_thresh) &&
        (dx1 <= carPlant_smoothUp_followerSto_P.div0protectabspoly3_thresh)) {
      dx1 = 2.9802322387695312E-8 / (3.0 - rt_powd_snf(dx1 /
        1.4901161193847656e-8, 2.0));
    } else {
      dx1 = fabs(dx1);
    }

    /* End of Switch: '<S212>/Switch' */

    /* Product: '<S207>/Divide' incorporates:
     *  Constant: '<S207>/Constant'
     */
    rtb_sincos_o2_idx_2 *=
      carPlant_smoothUp_followerSto_P.KinematicSteering1_WhlBase;

    /* Trigonometry: '<S207>/Trigonometric Function1' incorporates:
     *  Product: '<S207>/Divide1'
     */
    rtb_UnaryMinus1_g_idx_2 = atan(1.0 / dx1 * rtb_sincos_o2_idx_2);

    /* Sum: '<S207>/Add1' incorporates:
     *  Constant: '<S207>/Constant'
     */
    rtb_Integrator_d =
      carPlant_smoothUp_followerSto_P.KinematicSteering1_WhlBase -
      rtb_Integrator_d;

    /* Switch: '<S209>/Switch' incorporates:
     *  Abs: '<S209>/Abs'
     *  Constant: '<S213>/Constant'
     *  Constant: '<S214>/Constant'
     *  Fcn: '<S209>/Fcn'
     *  Logic: '<S209>/Logical Operator'
     *  RelationalOperator: '<S213>/Compare'
     *  RelationalOperator: '<S214>/Compare'
     */
    if ((rtb_Integrator_d >=
         -carPlant_smoothUp_followerSto_P.div0protectabspoly_thresh) &&
        (rtb_Integrator_d <=
         carPlant_smoothUp_followerSto_P.div0protectabspoly_thresh)) {
      dx1 = 2.9802322387695312E-8 / (3.0 - rt_powd_snf(rtb_Integrator_d /
        1.4901161193847656e-8, 2.0));
    } else {
      dx1 = fabs(rtb_Integrator_d);
    }

    /* End of Switch: '<S209>/Switch' */

    /* Trigonometry: '<S207>/Trigonometric Function2' incorporates:
     *  Product: '<S207>/Divide2'
     */
    rtb_Integrator_d = atan(rtb_sincos_o2_idx_2 / dx1);

    /* Switch: '<S10>/Switch' incorporates:
     *  Constant: '<S10>/index'
     *  UnaryMinus: '<S10>/Unary Minus1'
     */
    if (carPlant_smoothUp_followerSto_P.index_Value >
        carPlant_smoothUp_followerSto_P.Switch_Threshold_j) {
      carPlant_smoothUp_followerSto_B.VectorConcatenate1_a[0] =
        rtb_UnaryMinus1_g_idx_2;
    } else {
      carPlant_smoothUp_followerSto_B.VectorConcatenate1_a[0] =
        -rtb_Integrator_d;
    }

    /* End of Switch: '<S10>/Switch' */

    /* Switch: '<S10>/Switch1' incorporates:
     *  Constant: '<S10>/index'
     *  UnaryMinus: '<S10>/Unary Minus'
     */
    if (carPlant_smoothUp_followerSto_P.index_Value >
        carPlant_smoothUp_followerSto_P.Switch1_Threshold_h) {
      carPlant_smoothUp_followerSto_B.VectorConcatenate1_a[1] = rtb_Integrator_d;
    } else {
      carPlant_smoothUp_followerSto_B.VectorConcatenate1_a[1] =
        -rtb_UnaryMinus1_g_idx_2;
    }

    /* End of Switch: '<S10>/Switch1' */

    /* SignalConversion: '<S174>/ConcatBufferAtVector Concatenate4In1' incorporates:
     *  Constant: '<S174>/Constant'
     */
    carPlant_smoothUp_followerSto_B.VectorConcatenate4[0] =
      carPlant_smoothUp_followerSto_P.Constant_Value_a;

    /* SignalConversion: '<S174>/ConcatBufferAtVector Concatenate4In2' incorporates:
     *  Constant: '<S174>/Constant'
     */
    carPlant_smoothUp_followerSto_B.VectorConcatenate4[1] =
      carPlant_smoothUp_followerSto_P.Constant_Value_a;

    /* SignalConversion: '<S160>/ConcatBufferAtVector Concatenate2In1' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate2_j[0] = 0.0;

    /* SignalConversion: '<S160>/ConcatBufferAtVector Concatenate2In2' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate2_j[1] = 0.0;

    /* SignalConversion: '<S160>/ConcatBufferAtVector Concatenate3In1' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate3_b[0] = 0.0;

    /* SignalConversion: '<S160>/ConcatBufferAtVector Concatenate3In2' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate3_b[1] = 0.0;

    /* SignalConversion: '<S170>/ConcatBufferAtVector Concatenate2In1' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate2_o[0] = 0.0;

    /* SignalConversion: '<S170>/ConcatBufferAtVector Concatenate2In2' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate2_o[1] = 0.0;

    /* SignalConversion: '<S170>/ConcatBufferAtVector Concatenate3In1' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate3_b5[0] = 0.0;

    /* SignalConversion: '<S170>/ConcatBufferAtVector Concatenate3In2' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate3_b5[1] = 0.0;
  }

  /* MATLAB Function: '<S8>/vehicle model' incorporates:
   *  Constant: '<Root>/Coefficient of Friction: Dry'
   *  Sum: '<S45>/Sum of Elements'
   */
  rtb_Integrator_d = carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[0];
  dx1 = carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[1];
  dx3 = carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[3];
  rtb_sincos_o2_idx_2 = carPlant_smoothUp_followerSto_P.FollowerVehicle_w[0] /
    2.0;
  rtb_UnaryMinus1_g_idx_2 = carPlant_smoothUp_followerSto_P.FollowerVehicle_w[1]
    / 2.0;
  rtb_sincos_o2_idx_2 = sqrt(carPlant_smoothUp_followerSto_P.FollowerVehicle_a *
    carPlant_smoothUp_followerSto_P.FollowerVehicle_a + rtb_sincos_o2_idx_2 *
    rtb_sincos_o2_idx_2) * carPlant_smoothUp_followerSto_B.VectorConcatenate1_i
    [3] + carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[1];
  delta_fl = carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[0] *
    carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[0];
  rtb_sincos_o2_idx_2 = sqrt(rtb_sincos_o2_idx_2 * rtb_sincos_o2_idx_2 +
    delta_fl);
  rtb_UnaryMinus1_g_idx_2 =
    carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[1] - sqrt
    (carPlant_smoothUp_followerSto_P.FollowerVehicle_b *
     carPlant_smoothUp_followerSto_P.FollowerVehicle_b + rtb_UnaryMinus1_g_idx_2
     * rtb_UnaryMinus1_g_idx_2) *
    carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[3];
  rtb_UnaryMinus1_g_idx_2 = sqrt(rtb_UnaryMinus1_g_idx_2 *
    rtb_UnaryMinus1_g_idx_2 + delta_fl);
  alfa_fl = fabs(carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[0]);
  if (alfa_fl < carPlant_smoothUp_followerSto_P.FollowerVehicle_xdot_tol) {
    iU = 1;
  }

  if (0 <= iU - 1) {
    z_data = alfa_fl / carPlant_smoothUp_followerSto_P.FollowerVehicle_xdot_tol;
    delta_fr = z_data;
  }

  if (0 <= iU - 1) {
    delta_fr = z_data * z_data;
  }

  i = iU - 1;
  delta_fl = delta_fr;
  for (iU = 0; iU <= i; iU++) {
    delta_fl = 2.0 * carPlant_smoothUp_followerSto_P.FollowerVehicle_xdot_tol /
      (3.0 - delta_fl);
    delta_fr = delta_fl;
  }

  alfa_rr = alfa_fl;
  if (alfa_fl < carPlant_smoothUp_followerSto_P.FollowerVehicle_xdot_tol) {
    alfa_rr = delta_fr;
  }

  iU = 0;
  if (carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[0] < 0.0) {
    iU = 1;
  }

  if (0 <= iU - 1) {
    z_data = -alfa_rr;
  }

  if (carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[0] < 0.0) {
    alfa_rr = z_data;
  }

  delta_fl = carPlant_smoothUp_followerSto_B.VectorConcatenate1_a[0];
  delta_fr = carPlant_smoothUp_followerSto_B.VectorConcatenate1_a[1];
  delta_rl = carPlant_smoothUp_followerSto_B.VectorConcatenate4[0];
  delta_rr = carPlant_smoothUp_followerSto_B.VectorConcatenate4[1];
  if (alfa_fl <= carPlant_smoothUp_followerSto_P.FollowerVehicle_xdot_tol) {
    Fz_rnom = tanh(4.0 * fabs
                   (carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[1]));
    alfa_fl = (atan((carPlant_smoothUp_followerSto_P.FollowerVehicle_a *
                     carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[3] +
                     carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[1]) /
                    (carPlant_smoothUp_followerSto_P.FollowerVehicle_w[0] / 2.0 *
                     carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[3] +
                     alfa_rr)) -
               carPlant_smoothUp_followerSto_B.VectorConcatenate1_a[0]) *
      Fz_rnom;
    alfa_fr = (atan((carPlant_smoothUp_followerSto_P.FollowerVehicle_a *
                     carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[3] +
                     carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[1]) /
                    (alfa_rr -
                     carPlant_smoothUp_followerSto_P.FollowerVehicle_w[0] / 2.0 *
                     carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[3])) -
               carPlant_smoothUp_followerSto_B.VectorConcatenate1_a[1]) *
      Fz_rnom;
    alfa_rl = (atan((carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[1] -
                     carPlant_smoothUp_followerSto_P.FollowerVehicle_b *
                     carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[3]) /
                    (carPlant_smoothUp_followerSto_P.FollowerVehicle_w[1] / 2.0 *
                     carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[3] +
                     alfa_rr)) -
               carPlant_smoothUp_followerSto_B.VectorConcatenate4[0]) * Fz_rnom;
    alfa_rr = (atan((carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[1] -
                     carPlant_smoothUp_followerSto_P.FollowerVehicle_b *
                     carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[3]) /
                    (alfa_rr -
                     carPlant_smoothUp_followerSto_P.FollowerVehicle_w[1] / 2.0 *
                     carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[3])) -
               carPlant_smoothUp_followerSto_B.VectorConcatenate4[1]) * Fz_rnom;
  } else {
    Fz_rnom = carPlant_smoothUp_followerSto_P.FollowerVehicle_w[0] / 2.0 *
      carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[3];
    alfa_fr = carPlant_smoothUp_followerSto_P.FollowerVehicle_a *
      carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[3] +
      carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[1];
    Fz_fnom = tanh(4.0 * carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[0]);
    alfa_fl = (atan(alfa_fr / (Fz_rnom + alfa_rr)) -
               carPlant_smoothUp_followerSto_B.VectorConcatenate1_a[0]) *
      Fz_fnom;
    alfa_fr = (atan(alfa_fr / (alfa_rr - Fz_rnom)) -
               carPlant_smoothUp_followerSto_B.VectorConcatenate1_a[1]) *
      Fz_fnom;
    Fz_rnom = carPlant_smoothUp_followerSto_P.FollowerVehicle_w[1] / 2.0 *
      carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[3];
    Fz_idx_0 = carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[1] -
      carPlant_smoothUp_followerSto_P.FollowerVehicle_b *
      carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[3];
    alfa_rl = (atan(Fz_idx_0 / (Fz_rnom + alfa_rr)) -
               carPlant_smoothUp_followerSto_B.VectorConcatenate4[0]) * Fz_fnom;
    alfa_rr = (atan(Fz_idx_0 / (alfa_rr - Fz_rnom)) -
               carPlant_smoothUp_followerSto_B.VectorConcatenate4[1]) * Fz_fnom;
  }

  Fz_idx_0 = 0.0;
  Fz_idx_1 = 0.0;
  Fz_idx_2 = 0.0;
  Fz_idx_3 = 0.0;
  for (iU = 0; iU < 11; iU++) {
    if (iU == 0) {
      Fz_rnom = -dx1 * dx3 * carPlant_smoothUp_followerSto_P.FollowerVehicle_h;
      yddot = rtb_VectorConcatenate1_m_idx_0 *
        carPlant_smoothUp_followerSto_P.FollowerVehicle_h;
      Fz_fnom = ((((carPlant_smoothUp_followerSto_P.FollowerVehicle_g *
                    carPlant_smoothUp_followerSto_P.FollowerVehicle_b - Fz_rnom)
                   * carPlant_smoothUp_followerSto_P.FollowerVehicle_m +
                   rtb_VectorConcatenate1_k *
                   carPlant_smoothUp_followerSto_P.FollowerVehicle_b) + yddot) -
                 rtb_sincos_o1_idx_1) / rtb_VectorConcatenate1_m_idx_0_;
      Fz_rnom = ((((Fz_rnom + carPlant_smoothUp_followerSto_P.FollowerVehicle_g *
                    carPlant_smoothUp_followerSto_P.FollowerVehicle_a) *
                   carPlant_smoothUp_followerSto_P.FollowerVehicle_m +
                   rtb_VectorConcatenate1_k *
                   carPlant_smoothUp_followerSto_P.FollowerVehicle_a) - yddot) +
                 rtb_sincos_o1_idx_1) / rtb_VectorConcatenate1_m_idx_0_;
      Fz_idx_1 = rtb_Integrator_d * dx3;
      Fz_idx_2 = rtb_sincos_o2_idx_1 *
        carPlant_smoothUp_followerSto_P.FollowerVehicle_h;
      yddot = (Fz_idx_1 * carPlant_smoothUp_followerSto_P.FollowerVehicle_m *
               carPlant_smoothUp_followerSto_P.FollowerVehicle_h - Fz_idx_2) +
        rtb_sincos_o1_idx_0;
      Fz_idx_0 = (yddot + Fz_fnom) *
        (carPlant_smoothUp_followerSto_P.FollowerVehicle_w[0] / 2.0 -
         carPlant_smoothUp_followerSto_P.FollowerVehicle_d) /
        (carPlant_smoothUp_followerSto_P.FollowerVehicle_w[0] / 2.0);
      Fz_idx_3 = (Fz_idx_1 * -carPlant_smoothUp_followerSto_P.FollowerVehicle_m *
                  carPlant_smoothUp_followerSto_P.FollowerVehicle_h + Fz_idx_2)
        + rtb_sincos_o1_idx_0;
      Fz_idx_1 = (Fz_idx_3 + Fz_fnom) *
        (carPlant_smoothUp_followerSto_P.FollowerVehicle_w[0] / 2.0 +
         carPlant_smoothUp_followerSto_P.FollowerVehicle_d) /
        (carPlant_smoothUp_followerSto_P.FollowerVehicle_w[0] / 2.0);
      Fz_idx_2 = (yddot + Fz_rnom) *
        (carPlant_smoothUp_followerSto_P.FollowerVehicle_w[1] / 2.0 -
         carPlant_smoothUp_followerSto_P.FollowerVehicle_d) /
        (carPlant_smoothUp_followerSto_P.FollowerVehicle_w[1] / 2.0);
      Fz_idx_3 = (Fz_idx_3 + Fz_rnom) *
        (carPlant_smoothUp_followerSto_P.FollowerVehicle_w[1] / 2.0 +
         carPlant_smoothUp_followerSto_P.FollowerVehicle_d) /
        (carPlant_smoothUp_followerSto_P.FollowerVehicle_w[1] / 2.0);
      if (Fz_idx_0 < 0.0) {
        Fz_idx_0 = 0.0;
      }

      if (Fz_idx_1 < 0.0) {
        Fz_idx_1 = 0.0;
      }

      if (Fz_idx_2 < 0.0) {
        Fz_idx_2 = 0.0;
      }

      if (Fz_idx_3 < 0.0) {
        Fz_idx_3 = 0.0;
      }
    }

    carPlant_s_automlvehdynftiresat
      (-carPlant_smoothUp_followerSto_B.VectorConcatenate2[0] / 2.0 * alfa_fl *
       carPlant_smoothUp_followerSto_P.CoefficientofFrictionDry_Value[0] *
       Fz_idx_0 / carPlant_smoothUp_followerSto_P.FollowerVehicle_Fznom,
       carPlant_smoothUp_followerSto_P.vehiclemodel_Fxtire_sat * Fz_idx_0 /
       carPlant_smoothUp_followerSto_P.FollowerVehicle_Fznom,
       carPlant_smoothUp_followerSto_P.vehiclemodel_Fytire_sat * Fz_idx_0 /
       carPlant_smoothUp_followerSto_P.FollowerVehicle_Fznom, &yddot, &rdot);
    carPlant_s_automlvehdynftiresat
      (-carPlant_smoothUp_followerSto_B.VectorConcatenate2[1] / 2.0 * alfa_fr *
       carPlant_smoothUp_followerSto_P.CoefficientofFrictionDry_Value[2] *
       Fz_idx_1 / carPlant_smoothUp_followerSto_P.FollowerVehicle_Fznom,
       carPlant_smoothUp_followerSto_P.vehiclemodel_Fxtire_sat * Fz_idx_1 /
       carPlant_smoothUp_followerSto_P.FollowerVehicle_Fznom,
       carPlant_smoothUp_followerSto_P.vehiclemodel_Fytire_sat * Fz_idx_1 /
       carPlant_smoothUp_followerSto_P.FollowerVehicle_Fznom, &yddot, &Fz_fnom);
    carPlant_s_automlvehdynftiresat
      (-carPlant_smoothUp_followerSto_B.VectorConcatenate1_bi[0] / 2.0 * alfa_rl
       * carPlant_smoothUp_followerSto_P.CoefficientofFrictionDry_Value[1] *
       Fz_idx_2 / carPlant_smoothUp_followerSto_P.FollowerVehicle_Fznom,
       carPlant_smoothUp_followerSto_P.vehiclemodel_Fxtire_sat * Fz_idx_2 /
       carPlant_smoothUp_followerSto_P.FollowerVehicle_Fznom,
       carPlant_smoothUp_followerSto_P.vehiclemodel_Fytire_sat * Fz_idx_2 /
       carPlant_smoothUp_followerSto_P.FollowerVehicle_Fznom, &yddot, &Fz_rnom);
    carPlant_s_automlvehdynftiresat
      (-carPlant_smoothUp_followerSto_B.VectorConcatenate1_bi[1] / 2.0 * alfa_rr
       * carPlant_smoothUp_followerSto_P.CoefficientofFrictionDry_Value[3] *
       Fz_idx_3 / carPlant_smoothUp_followerSto_P.FollowerVehicle_Fznom,
       carPlant_smoothUp_followerSto_P.vehiclemodel_Fxtire_sat * Fz_idx_3 /
       carPlant_smoothUp_followerSto_P.FollowerVehicle_Fznom,
       carPlant_smoothUp_followerSto_P.vehiclemodel_Fytire_sat * Fz_idx_3 /
       carPlant_smoothUp_followerSto_P.FollowerVehicle_Fznom, &yddot, &Fz_idx_0);
    Fz_idx_1 = cos(delta_fl);
    Fz_idx_2 = sin(delta_fl);
    Fz_idx_3 = cos(delta_fr);
    b_Fy_fr_tmp = sin(delta_fr);
    b_Fy_rl_tmp = cos(delta_rl);
    b_Fy_rl_tmp_0 = sin(delta_rl);
    b_Fy_rl = Fz_rnom * b_Fy_rl_tmp - 0.0 * b_Fy_rl_tmp_0;
    b_Fy_rr_tmp = cos(delta_rr);
    b_Fy_rr_tmp_0 = sin(delta_rr);
    b_Fy_rr = Fz_idx_0 * b_Fy_rr_tmp - 0.0 * b_Fy_rr_tmp_0;
    yddot_tmp = (rdot * Fz_idx_1 - 0.0 * Fz_idx_2) + (Fz_fnom * Fz_idx_3 - 0.0 *
      b_Fy_fr_tmp);
    yddot = (((yddot_tmp + b_Fy_rl) + b_Fy_rr) + rtb_sincos_o2_idx_1) /
      carPlant_smoothUp_followerSto_P.FollowerVehicle_m + -rtb_Integrator_d *
      dx3;
    rdot = (((((0.0 * Fz_idx_1 - rdot * Fz_idx_2) - (0.0 * Fz_idx_3 - Fz_fnom *
                b_Fy_fr_tmp)) *
              (carPlant_smoothUp_followerSto_P.FollowerVehicle_w[0] / 2.0) +
              (yddot_tmp * carPlant_smoothUp_followerSto_P.FollowerVehicle_a -
               (b_Fy_rl + b_Fy_rr) *
               carPlant_smoothUp_followerSto_P.FollowerVehicle_b)) + ((0.0 *
               b_Fy_rl_tmp - Fz_rnom * b_Fy_rl_tmp_0) - (0.0 * b_Fy_rr_tmp -
               Fz_idx_0 * b_Fy_rr_tmp_0)) *
             (carPlant_smoothUp_followerSto_P.FollowerVehicle_w[1] / 2.0)) +
            rtb_sincos_o1_idx_2) /
      carPlant_smoothUp_followerSto_P.FollowerVehicle_Izz;
    Fz_rnom = (0.0 - dx1 * dx3) *
      carPlant_smoothUp_followerSto_P.FollowerVehicle_h;
    Fz_fnom = ((((carPlant_smoothUp_followerSto_P.FollowerVehicle_g *
                  carPlant_smoothUp_followerSto_P.FollowerVehicle_b - Fz_rnom) *
                 carPlant_smoothUp_followerSto_P.FollowerVehicle_m +
                 rtb_VectorConcatenate1_k *
                 carPlant_smoothUp_followerSto_P.FollowerVehicle_b) +
                rtb_VectorConcatenate1_m_idx_0 *
                carPlant_smoothUp_followerSto_P.FollowerVehicle_h) -
               rtb_sincos_o1_idx_1) / rtb_VectorConcatenate1_m_idx_0_;
    Fz_rnom = ((((Fz_rnom + carPlant_smoothUp_followerSto_P.FollowerVehicle_g *
                  carPlant_smoothUp_followerSto_P.FollowerVehicle_a) *
                 carPlant_smoothUp_followerSto_P.FollowerVehicle_m +
                 rtb_VectorConcatenate1_k *
                 carPlant_smoothUp_followerSto_P.FollowerVehicle_a) -
                rtb_VectorConcatenate1_m_idx_0 *
                carPlant_smoothUp_followerSto_P.FollowerVehicle_h) +
               rtb_sincos_o1_idx_1) / rtb_VectorConcatenate1_m_idx_0_;
    Fz_idx_1 = rtb_Integrator_d * dx3 + yddot;
    Fz_idx_2 = (Fz_idx_1 * carPlant_smoothUp_followerSto_P.FollowerVehicle_m *
                carPlant_smoothUp_followerSto_P.FollowerVehicle_h -
                rtb_sincos_o2_idx_1 *
                carPlant_smoothUp_followerSto_P.FollowerVehicle_h) +
      rtb_sincos_o1_idx_0;
    Fz_idx_0 = (Fz_idx_2 + Fz_fnom) *
      (carPlant_smoothUp_followerSto_P.FollowerVehicle_w[0] / 2.0 -
       carPlant_smoothUp_followerSto_P.FollowerVehicle_d) /
      (carPlant_smoothUp_followerSto_P.FollowerVehicle_w[0] / 2.0);
    Fz_idx_3 = (Fz_idx_1 * -carPlant_smoothUp_followerSto_P.FollowerVehicle_m *
                carPlant_smoothUp_followerSto_P.FollowerVehicle_h +
                rtb_sincos_o2_idx_1 *
                carPlant_smoothUp_followerSto_P.FollowerVehicle_h) +
      rtb_sincos_o1_idx_0;
    Fz_idx_1 = (Fz_idx_3 + Fz_fnom) *
      (carPlant_smoothUp_followerSto_P.FollowerVehicle_w[0] / 2.0 +
       carPlant_smoothUp_followerSto_P.FollowerVehicle_d) /
      (carPlant_smoothUp_followerSto_P.FollowerVehicle_w[0] / 2.0);
    Fz_idx_2 = (Fz_idx_2 + Fz_rnom) *
      (carPlant_smoothUp_followerSto_P.FollowerVehicle_w[1] / 2.0 -
       carPlant_smoothUp_followerSto_P.FollowerVehicle_d) /
      (carPlant_smoothUp_followerSto_P.FollowerVehicle_w[1] / 2.0);
    Fz_idx_3 = (Fz_idx_3 + Fz_rnom) *
      (carPlant_smoothUp_followerSto_P.FollowerVehicle_w[1] / 2.0 +
       carPlant_smoothUp_followerSto_P.FollowerVehicle_d) /
      (carPlant_smoothUp_followerSto_P.FollowerVehicle_w[1] / 2.0);
    if (Fz_idx_0 < 0.0) {
      Fz_idx_0 = 0.0;
    }

    if (Fz_idx_1 < 0.0) {
      Fz_idx_1 = 0.0;
    }

    if (Fz_idx_2 < 0.0) {
      Fz_idx_2 = 0.0;
    }

    if (Fz_idx_3 < 0.0) {
      Fz_idx_3 = 0.0;
    }
  }

  carPlant_smoothUp_followerSto_B.stateDer_p[0] = 0.0;
  carPlant_smoothUp_followerSto_B.stateDer_p[1] = yddot;
  carPlant_smoothUp_followerSto_B.stateDer_p[2] =
    carPlant_smoothUp_followerSto_B.VectorConcatenate1_i[3];
  carPlant_smoothUp_followerSto_B.stateDer_p[3] = rdot;

  /* MATLAB Function: '<S94>/COMB2I' */
  carPlant_smoothUp_follow_COMB2I
    (carPlant_smoothUp_followerSto_B.VectorConcatenate1_i,
     &carPlant_smoothUp_followerSto_B.sf_COMB2I);
  if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
    /* SignalConversion: '<S178>/ConcatBufferAtVector ConcatenateIn1' incorporates:
     *  Constant: '<S178>/Constant1'
     */
    carPlant_smoothUp_followerSto_B.VectorConcatenate_o[0] =
      carPlant_smoothUp_followerSto_P.FollowerVehicle_sigma_f;

    /* SignalConversion: '<S178>/ConcatBufferAtVector ConcatenateIn2' incorporates:
     *  Constant: '<S178>/Constant1'
     */
    carPlant_smoothUp_followerSto_B.VectorConcatenate_o[1] =
      carPlant_smoothUp_followerSto_P.FollowerVehicle_sigma_f;

    /* SignalConversion: '<S178>/ConcatBufferAtVector ConcatenateIn3' incorporates:
     *  Constant: '<S178>/Constant2'
     */
    carPlant_smoothUp_followerSto_B.VectorConcatenate_o[2] =
      carPlant_smoothUp_followerSto_P.FollowerVehicle_sigma_r;

    /* SignalConversion: '<S178>/ConcatBufferAtVector ConcatenateIn4' incorporates:
     *  Constant: '<S178>/Constant2'
     */
    carPlant_smoothUp_followerSto_B.VectorConcatenate_o[3] =
      carPlant_smoothUp_followerSto_P.FollowerVehicle_sigma_r;

    /* Backlash: '<S186>/Backlash' incorporates:
     *  Constant: '<Root>/Constant3'
     */
    rtb_Integrator_d = carPlant_smoothUp_followerSto_P.KinematicSteering_Db /
      2.0;
    if (carPlant_smoothUp_followerSto_P.Constant3_Value <
        carPlant_smoothUp_followerSt_DW.PrevY_e - rtb_Integrator_d) {
      carPlant_smoothUp_followerSto_B.Backlash_a =
        carPlant_smoothUp_followerSto_P.Constant3_Value + rtb_Integrator_d;
    } else if (carPlant_smoothUp_followerSto_P.Constant3_Value <=
               carPlant_smoothUp_followerSt_DW.PrevY_e + rtb_Integrator_d) {
      carPlant_smoothUp_followerSto_B.Backlash_a =
        carPlant_smoothUp_followerSt_DW.PrevY_e;
    } else {
      carPlant_smoothUp_followerSto_B.Backlash_a =
        carPlant_smoothUp_followerSto_P.Constant3_Value - rtb_Integrator_d;
    }

    /* End of Backlash: '<S186>/Backlash' */

    /* Saturate: '<S186>/Saturation' */
    if (carPlant_smoothUp_followerSto_B.Backlash_a >
        carPlant_smoothUp_followerSto_P.KinematicSteering_StrgRng) {
      rtb_Integrator_d =
        carPlant_smoothUp_followerSto_P.KinematicSteering_StrgRng;
    } else if (carPlant_smoothUp_followerSto_B.Backlash_a <
               -carPlant_smoothUp_followerSto_P.KinematicSteering_StrgRng) {
      rtb_Integrator_d =
        -carPlant_smoothUp_followerSto_P.KinematicSteering_StrgRng;
    } else {
      rtb_Integrator_d = carPlant_smoothUp_followerSto_B.Backlash_a;
    }

    /* End of Saturate: '<S186>/Saturation' */

    /* Trigonometry: '<S188>/Trigonometric Function' incorporates:
     *  Gain: '<S187>/Gain'
     */
    dx1 = tan(1.0 / carPlant_smoothUp_followerSto_P.KinematicSteering_StrgRatio *
              rtb_Integrator_d);

    /* Product: '<S188>/Divide5' incorporates:
     *  Constant: '<S188>/Constant1'
     */
    rtb_Integrator_d =
      carPlant_smoothUp_followerSto_P.KinematicSteering_TrckWdth / 2.0 * dx1;

    /* Sum: '<S188>/Add' incorporates:
     *  Constant: '<S188>/Constant'
     */
    rdot = carPlant_smoothUp_followerSto_P.KinematicSteering_WhlBase +
      rtb_Integrator_d;

    /* Sum: '<S188>/Add1' incorporates:
     *  Constant: '<S188>/Constant'
     */
    yddot = carPlant_smoothUp_followerSto_P.KinematicSteering_WhlBase -
      rtb_Integrator_d;

    /* Product: '<S188>/Divide' incorporates:
     *  Constant: '<S188>/Constant'
     */
    dx1 *= carPlant_smoothUp_followerSto_P.KinematicSteering_WhlBase;

    /* Switch: '<S193>/Switch' incorporates:
     *  Abs: '<S193>/Abs'
     *  Constant: '<S200>/Constant'
     *  Constant: '<S201>/Constant'
     *  Fcn: '<S193>/Fcn'
     *  Logic: '<S193>/Logical Operator'
     *  RelationalOperator: '<S200>/Compare'
     *  RelationalOperator: '<S201>/Compare'
     */
    if ((rdot >= -carPlant_smoothUp_followerSto_P.div0protectabspoly3_thresh_g) &&
        (rdot <= carPlant_smoothUp_followerSto_P.div0protectabspoly3_thresh_g))
    {
      rtb_Integrator_d = 2.9802322387695312E-8 / (3.0 - rt_powd_snf(rdot /
        1.4901161193847656e-8, 2.0));
    } else {
      rtb_Integrator_d = fabs(rdot);
    }

    /* End of Switch: '<S193>/Switch' */

    /* Product: '<S188>/Divide1' */
    rdot = 1.0 / rtb_Integrator_d * dx1;

    /* Switch: '<S190>/Switch' incorporates:
     *  Abs: '<S190>/Abs'
     *  Constant: '<S194>/Constant'
     *  Constant: '<S195>/Constant'
     *  Fcn: '<S190>/Fcn'
     *  Logic: '<S190>/Logical Operator'
     *  RelationalOperator: '<S194>/Compare'
     *  RelationalOperator: '<S195>/Compare'
     */
    if ((yddot >= -carPlant_smoothUp_followerSto_P.div0protectabspoly_thresh_b) &&
        (yddot <= carPlant_smoothUp_followerSto_P.div0protectabspoly_thresh_b))
    {
      rtb_Integrator_d = 2.9802322387695312E-8 / (3.0 - rt_powd_snf(yddot /
        1.4901161193847656e-8, 2.0));
    } else {
      rtb_Integrator_d = fabs(yddot);
    }

    /* End of Switch: '<S190>/Switch' */

    /* Trigonometry: '<S188>/Trigonometric Function1' */
    yddot = atan(rdot);

    /* Trigonometry: '<S188>/Trigonometric Function2' incorporates:
     *  Product: '<S188>/Divide2'
     */
    rdot = atan(dx1 / rtb_Integrator_d);

    /* Switch: '<S9>/Switch' incorporates:
     *  Constant: '<S9>/index'
     *  UnaryMinus: '<S9>/Unary Minus1'
     */
    if (carPlant_smoothUp_followerSto_P.index_Value_d >
        carPlant_smoothUp_followerSto_P.Switch_Threshold_d) {
      carPlant_smoothUp_followerSto_B.VectorConcatenate_i[0] = yddot;
    } else {
      carPlant_smoothUp_followerSto_B.VectorConcatenate_i[0] = -rdot;
    }

    /* End of Switch: '<S9>/Switch' */

    /* Switch: '<S9>/Switch1' incorporates:
     *  Constant: '<S9>/index'
     *  UnaryMinus: '<S9>/Unary Minus'
     */
    if (carPlant_smoothUp_followerSto_P.index_Value_d >
        carPlant_smoothUp_followerSto_P.Switch1_Threshold_p) {
      carPlant_smoothUp_followerSto_B.VectorConcatenate_i[1] = rdot;
    } else {
      carPlant_smoothUp_followerSto_B.VectorConcatenate_i[1] = -yddot;
    }

    /* End of Switch: '<S9>/Switch1' */

    /* SignalConversion: '<S381>/ConcatBufferAtVector Concatenate1In1' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate1_iz[0] = 0.0;

    /* SignalConversion: '<S381>/ConcatBufferAtVector Concatenate1In2' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate1_iz[1] = 0.0;

    /* SignalConversion: '<S381>/ConcatBufferAtVector Concatenate1In3' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate1_iz[2] = 0.0;
  }

  /* Product: '<S181>/Product1' incorporates:
   *  Integrator: '<S181>/lateral'
   *  MATLAB Function: '<S8>/vehicle model'
   *  Sum: '<S181>/Sum'
   */
  carPlant_smoothUp_followerSto_B.Product1[0] = (alfa_fl -
    carPlant_smoothUp_followerSto_X.lateral_CSTATE[0]) * rtb_sincos_o2_idx_2 /
    carPlant_smoothUp_followerSto_B.VectorConcatenate_o[0];
  carPlant_smoothUp_followerSto_B.Product1[1] = (alfa_fr -
    carPlant_smoothUp_followerSto_X.lateral_CSTATE[1]) * rtb_sincos_o2_idx_2 /
    carPlant_smoothUp_followerSto_B.VectorConcatenate_o[1];
  carPlant_smoothUp_followerSto_B.Product1[2] = (alfa_rl -
    carPlant_smoothUp_followerSto_X.lateral_CSTATE[2]) * rtb_UnaryMinus1_g_idx_2
    / carPlant_smoothUp_followerSto_B.VectorConcatenate_o[2];
  carPlant_smoothUp_followerSto_B.Product1[3] = (alfa_rr -
    carPlant_smoothUp_followerSto_X.lateral_CSTATE[3]) * rtb_UnaryMinus1_g_idx_2
    / carPlant_smoothUp_followerSto_B.VectorConcatenate_o[3];

  /* Product: '<S242>/Product1' */
  rtb_Product1_o = rtb_UnaryMinus1_g_idx_0_tmp *
    carPlant_smoothUp_followerSto_B.VectorConcatenate1_iz[1];

  /* Sum: '<S241>/Add1' incorporates:
   *  Product: '<S242>/Product3'
   *  SignalConversion: '<S225>/ConcatBufferAtVector Concatenate1In1'
   *  Sum: '<S242>/Add'
   */
  rtb_VectorConcatenate1_k = carPlant_smoothUp_followerSto_B.VectorConcatenate1
    [0] - (rtb_sincos_o2_idx_0_tmp *
           carPlant_smoothUp_followerSto_B.VectorConcatenate1_iz[0] +
           rtb_Product1_o);
  rtb_VectorConcatenate1_m_idx_0 = rtb_VectorConcatenate1_k;

  /* Product: '<S241>/Product' */
  rtb_sincos_o1_idx_0 = rtb_VectorConcatenate1_k * rtb_VectorConcatenate1_k;

  /* Sum: '<S241>/Add1' incorporates:
   *  Product: '<S242>/Product'
   *  Product: '<S242>/Product2'
   *  SignalConversion: '<S225>/ConcatBufferAtVector Concatenate1In2'
   *  Sum: '<S242>/Add1'
   */
  rtb_VectorConcatenate1_k = carPlant_smoothUp_followerSto_B.VectorConcatenate1
    [1] - (rtb_sincos_o2_idx_0_tmp *
           carPlant_smoothUp_followerSto_B.VectorConcatenate1_iz[1] -
           rtb_UnaryMinus1_g_idx_0_tmp *
           carPlant_smoothUp_followerSto_B.VectorConcatenate1_iz[0]);

  /* Sum: '<S241>/Sum of Elements' */
  iU = 0;

  /* Sqrt: '<S241>/Sqrt' incorporates:
   *  Product: '<S241>/Product'
   *  SignalConversion: '<S242>/ConcatBufferAtVector Concatenate2In3'
   *  Sum: '<S241>/Add1'
   *  Sum: '<S241>/Sum of Elements'
   */
  rtb_Integrator_d = sqrt((rtb_VectorConcatenate1_k * rtb_VectorConcatenate1_k +
    rtb_sincos_o1_idx_0) + (0.0 -
    carPlant_smoothUp_followerSto_B.VectorConcatenate1_iz[2]) * (0.0 -
    carPlant_smoothUp_followerSto_B.VectorConcatenate1_iz[2]));

  /* Product: '<S241>/Product2' */
  yddot = rtb_Integrator_d * rtb_Integrator_d;
  if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
    /* Constant: '<S241>/Constant' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate_c[0] =
      carPlant_smoothUp_followerSto_P.LeadVehicle_Cd;
  }

  /* Trigonometry: '<S241>/Trigonometric Function' */
  rtb_Integrator_d = rt_atan2d_snf(rtb_VectorConcatenate1_k,
    rtb_VectorConcatenate1_m_idx_0);

  /* Lookup_n-D: '<S241>/Cs' */
  carPlant_smoothUp_followerSto_B.VectorConcatenate_c[1] = look1_binlcpw
    (rtb_Integrator_d, carPlant_smoothUp_followerSto_P.LeadVehicle_beta_w,
     carPlant_smoothUp_followerSto_P.LeadVehicle_Cs, 30U);
  if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
    /* Constant: '<S241>/Constant1' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate_c[2] =
      carPlant_smoothUp_followerSto_P.LeadVehicle_Cl;

    /* UnaryMinus: '<S241>/Unary Minus' incorporates:
     *  Constant: '<S241>/Constant4'
     */
    carPlant_smoothUp_followerSto_B.UnaryMinus_h[0] =
      -carPlant_smoothUp_followerSto_P.Constant4_Value_l[0];
    carPlant_smoothUp_followerSto_B.UnaryMinus_h[1] =
      -carPlant_smoothUp_followerSto_P.Constant4_Value_l[1];
    carPlant_smoothUp_followerSto_B.UnaryMinus_h[2] =
      -carPlant_smoothUp_followerSto_P.Constant4_Value_l[2];
  }

  /* Lookup_n-D: '<S241>/Crm' */
  carPlant_smoothUp_followerSto_B.VectorConcatenate_c[3] = look1_binlxpw
    (rtb_Integrator_d, carPlant_smoothUp_followerSto_P.Crm_bp01Data_p,
     carPlant_smoothUp_followerSto_P.Crm_tableData_n, 1U);

  /* Switch: '<S241>/Switch' incorporates:
   *  Constant: '<S241>/Constant4'
   */
  if (rtb_VectorConcatenate1_m_idx_0 >=
      carPlant_smoothUp_followerSto_P.Switch_Threshold_p) {
    rtb_VectorConcatenate1_m_idx_0 =
      carPlant_smoothUp_followerSto_P.Constant4_Value_l[0];
  } else {
    rtb_VectorConcatenate1_m_idx_0 =
      carPlant_smoothUp_followerSto_B.UnaryMinus_h[0];
  }

  /* Product: '<S241>/Product5' incorporates:
   *  Constant: '<S241>/Constant2'
   */
  carPlant_smoothUp_followerSto_B.VectorConcatenate_c[4] =
    rtb_VectorConcatenate1_m_idx_0 *
    carPlant_smoothUp_followerSto_P.LeadVehicle_Cpm;

  /* Lookup_n-D: '<S241>/Cym' */
  carPlant_smoothUp_followerSto_B.VectorConcatenate_c[5] = look1_binlxpw
    (rtb_Integrator_d, carPlant_smoothUp_followerSto_P.LeadVehicle_beta_w,
     carPlant_smoothUp_followerSto_P.LeadVehicle_Cym, 30U);

  /* Gain: '<S241>/.5.*A.*Pabs.//R.//T' incorporates:
   *  Product: '<S241>/Product1'
   */
  rtb_Integrator_d = 0.5 * carPlant_smoothUp_followerSto_P.LeadVehicle_Af *
    carPlant_smoothUp_followerSto_P.LeadVehicle_Pabs /
    carPlant_smoothUp_followerSto_P.DragForce_R_e /
    carPlant_smoothUp_followerSto_P.LeadVehicle_Tair;
  for (i = 0; i < 6; i++) {
    rtb_uAPabsRT[i] = yddot *
      carPlant_smoothUp_followerSto_B.VectorConcatenate_c[i] * rtb_Integrator_d;
  }

  /* End of Gain: '<S241>/.5.*A.*Pabs.//R.//T' */

  /* Product: '<S241>/Product4' incorporates:
   *  Constant: '<S241>/Constant3'
   *  MATLAB Function: '<S11>/vehicle model'
   */
  rtb_VectorConcatenate1_m_idx_0_ =
    carPlant_smoothUp_followerSto_P.LeadVehicle_a +
    carPlant_smoothUp_followerSto_P.LeadVehicle_b;

  /* Sum: '<S11>/Add' incorporates:
   *  Constant: '<S241>/Constant3'
   *  Product: '<S241>/Product4'
   *  UnaryMinus: '<S225>/Unary Minus1'
   */
  rtb_sincos_o1_idx_0 = -(rtb_uAPabsRT[3] * rtb_VectorConcatenate1_m_idx_0_);

  /* Product: '<S241>/Product3' incorporates:
   *  UnaryMinus: '<S225>/Unary Minus'
   */
  rtb_VectorConcatenate1_m_idx_0 = -(rtb_VectorConcatenate1_m_idx_0 *
    rtb_uAPabsRT[0]);

  /* Sum: '<S11>/Add' incorporates:
   *  Constant: '<S241>/Constant3'
   *  Product: '<S241>/Product4'
   *  UnaryMinus: '<S225>/Unary Minus1'
   */
  rtb_sincos_o1_idx_1 = -(rtb_uAPabsRT[4] * rtb_VectorConcatenate1_m_idx_0_);

  /* Switch: '<S241>/Switch' incorporates:
   *  Constant: '<S241>/Constant4'
   */
  if (rtb_VectorConcatenate1_k >=
      carPlant_smoothUp_followerSto_P.Switch_Threshold_p) {
    rtb_VectorConcatenate1_k =
      carPlant_smoothUp_followerSto_P.Constant4_Value_l[1];
  } else {
    rtb_VectorConcatenate1_k = carPlant_smoothUp_followerSto_B.UnaryMinus_h[1];
  }

  /* Product: '<S241>/Product3' incorporates:
   *  UnaryMinus: '<S225>/Unary Minus'
   */
  rtb_sincos_o2_idx_1 = -(rtb_VectorConcatenate1_k * rtb_uAPabsRT[1]);

  /* Sum: '<S11>/Add' incorporates:
   *  Constant: '<S241>/Constant3'
   *  Product: '<S241>/Product4'
   *  UnaryMinus: '<S225>/Unary Minus1'
   */
  rtb_sincos_o1_idx_2 = -(rtb_uAPabsRT[5] * rtb_VectorConcatenate1_m_idx_0_);

  /* Switch: '<S241>/Switch' incorporates:
   *  Constant: '<S241>/Constant4'
   *  SignalConversion: '<S242>/ConcatBufferAtVector Concatenate2In3'
   *  Sum: '<S241>/Add1'
   */
  if (0.0 - carPlant_smoothUp_followerSto_B.VectorConcatenate1_iz[2] >=
      carPlant_smoothUp_followerSto_P.Switch_Threshold_p) {
    rtb_Integrator_d = carPlant_smoothUp_followerSto_P.Constant4_Value_l[2];
  } else {
    rtb_Integrator_d = carPlant_smoothUp_followerSto_B.UnaryMinus_h[2];
  }

  /* UnaryMinus: '<S225>/Unary Minus' incorporates:
   *  Product: '<S241>/Product3'
   */
  rtb_VectorConcatenate1_k = -(rtb_Integrator_d * rtb_uAPabsRT[2]);
  if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
    /* SignalConversion: '<S237>/ConcatBufferAtVector Concatenate2In1' incorporates:
     *  Constant: '<S237>/Cyf'
     */
    carPlant_smoothUp_followerSto_B.VectorConcatenate2_i[0] =
      carPlant_smoothUp_followerSto_P.LeadVehicle_Cy_f;

    /* SignalConversion: '<S237>/ConcatBufferAtVector Concatenate2In2' incorporates:
     *  Constant: '<S237>/Cyf'
     */
    carPlant_smoothUp_followerSto_B.VectorConcatenate2_i[1] =
      carPlant_smoothUp_followerSto_P.LeadVehicle_Cy_f;

    /* SignalConversion: '<S237>/ConcatBufferAtVector Concatenate1In1' incorporates:
     *  Constant: '<S237>/Cyr'
     */
    carPlant_smoothUp_followerSto_B.VectorConcatenate1_m[0] =
      carPlant_smoothUp_followerSto_P.LeadVehicle_Cy_r;

    /* SignalConversion: '<S237>/ConcatBufferAtVector Concatenate1In2' incorporates:
     *  Constant: '<S237>/Cyr'
     */
    carPlant_smoothUp_followerSto_B.VectorConcatenate1_m[1] =
      carPlant_smoothUp_followerSto_P.LeadVehicle_Cy_r;

    /* SignalConversion: '<S370>/ConcatBufferAtVector Concatenate4In1' incorporates:
     *  Constant: '<S370>/Constant'
     */
    carPlant_smoothUp_followerSto_B.VectorConcatenate4_c[0] =
      carPlant_smoothUp_followerSto_P.Constant_Value_bd;

    /* SignalConversion: '<S370>/ConcatBufferAtVector Concatenate4In2' incorporates:
     *  Constant: '<S370>/Constant'
     */
    carPlant_smoothUp_followerSto_B.VectorConcatenate4_c[1] =
      carPlant_smoothUp_followerSto_P.Constant_Value_bd;

    /* Concatenate: '<S350>/Vector Concatenate4' incorporates:
     *  Constant: '<S350>/Constant'
     */
    carPlant_smoothUp_followerSto_B.VectorConcatenate4_b[0] =
      carPlant_smoothUp_followerSto_P.Constant_Value_j;
    carPlant_smoothUp_followerSto_B.VectorConcatenate4_b[1] =
      carPlant_smoothUp_followerSto_P.Constant_Value_j;
    carPlant_smoothUp_followerSto_B.VectorConcatenate4_b[2] =
      carPlant_smoothUp_followerSto_P.Constant_Value_j;
    carPlant_smoothUp_followerSto_B.VectorConcatenate4_b[3] =
      carPlant_smoothUp_followerSto_P.Constant_Value_j;

    /* SignalConversion: '<S356>/ConcatBufferAtVector Concatenate2In1' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate2_h[0] = 0.0;

    /* SignalConversion: '<S356>/ConcatBufferAtVector Concatenate2In2' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate2_h[1] = 0.0;

    /* SignalConversion: '<S356>/ConcatBufferAtVector Concatenate3In1' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate3_m[0] = 0.0;

    /* SignalConversion: '<S356>/ConcatBufferAtVector Concatenate3In2' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate3_m[1] = 0.0;

    /* SignalConversion: '<S366>/ConcatBufferAtVector Concatenate2In1' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate2_c[0] = 0.0;

    /* SignalConversion: '<S366>/ConcatBufferAtVector Concatenate2In2' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate2_c[1] = 0.0;

    /* SignalConversion: '<S366>/ConcatBufferAtVector Concatenate3In1' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate3_f[0] = 0.0;

    /* SignalConversion: '<S366>/ConcatBufferAtVector Concatenate3In2' */
    carPlant_smoothUp_followerSto_B.VectorConcatenate3_f[1] = 0.0;
  }

  /* MATLAB Function: '<S11>/vehicle model' incorporates:
   *  Sum: '<S241>/Sum of Elements'
   */
  rtb_Integrator_d = carPlant_smoothUp_followerSto_B.VectorConcatenate1[0];
  dx1 = carPlant_smoothUp_followerSto_B.VectorConcatenate1[1];
  dx3 = carPlant_smoothUp_followerSto_B.VectorConcatenate1[3];
  rtb_sincos_o2_idx_2 = carPlant_smoothUp_followerSto_P.LeadVehicle_w[0] / 2.0;
  rtb_UnaryMinus1_g_idx_2 = carPlant_smoothUp_followerSto_P.LeadVehicle_w[1] /
    2.0;
  rtb_sincos_o2_idx_2 = sqrt(carPlant_smoothUp_followerSto_P.LeadVehicle_a *
    carPlant_smoothUp_followerSto_P.LeadVehicle_a + rtb_sincos_o2_idx_2 *
    rtb_sincos_o2_idx_2) * carPlant_smoothUp_followerSto_B.VectorConcatenate1[3]
    + carPlant_smoothUp_followerSto_B.VectorConcatenate1[1];
  delta_fl = carPlant_smoothUp_followerSto_B.VectorConcatenate1[0] *
    carPlant_smoothUp_followerSto_B.VectorConcatenate1[0];
  rtb_sincos_o2_idx_2 = sqrt(rtb_sincos_o2_idx_2 * rtb_sincos_o2_idx_2 +
    delta_fl);
  rtb_UnaryMinus1_g_idx_2 = carPlant_smoothUp_followerSto_B.VectorConcatenate1[1]
    - sqrt(carPlant_smoothUp_followerSto_P.LeadVehicle_b *
           carPlant_smoothUp_followerSto_P.LeadVehicle_b +
           rtb_UnaryMinus1_g_idx_2 * rtb_UnaryMinus1_g_idx_2) *
    carPlant_smoothUp_followerSto_B.VectorConcatenate1[3];
  rtb_UnaryMinus1_g_idx_2 = sqrt(rtb_UnaryMinus1_g_idx_2 *
    rtb_UnaryMinus1_g_idx_2 + delta_fl);
  alfa_fl = fabs(carPlant_smoothUp_followerSto_B.VectorConcatenate1[0]);
  if (alfa_fl < carPlant_smoothUp_followerSto_P.LeadVehicle_xdot_tol) {
    iU = 1;
  }

  if (0 <= iU - 1) {
    z_data_0 = alfa_fl / carPlant_smoothUp_followerSto_P.LeadVehicle_xdot_tol;
    z1_data = z_data;
  }

  if (0 <= iU - 1) {
    z1_data = z_data_0 * z_data_0;
  }

  i = iU - 1;
  delta_fl = z1_data;
  for (iU = 0; iU <= i; iU++) {
    delta_fl = 2.0 * carPlant_smoothUp_followerSto_P.LeadVehicle_xdot_tol / (3.0
      - delta_fl);
    z1_data = delta_fl;
  }

  alfa_rr = alfa_fl;
  if (alfa_fl < carPlant_smoothUp_followerSto_P.LeadVehicle_xdot_tol) {
    alfa_rr = z1_data;
  }

  iU = 0;
  if (carPlant_smoothUp_followerSto_B.VectorConcatenate1[0] < 0.0) {
    iU = 1;
  }

  if (0 <= iU - 1) {
    z_data_0 = -alfa_rr;
  }

  if (carPlant_smoothUp_followerSto_B.VectorConcatenate1[0] < 0.0) {
    alfa_rr = z_data_0;
  }

  delta_fl = carPlant_smoothUp_followerSto_B.VectorConcatenate_i[0];
  delta_fr = carPlant_smoothUp_followerSto_B.VectorConcatenate_i[1];
  delta_rl = carPlant_smoothUp_followerSto_B.VectorConcatenate4_c[0];
  delta_rr = carPlant_smoothUp_followerSto_B.VectorConcatenate4_c[1];
  if (alfa_fl <= carPlant_smoothUp_followerSto_P.LeadVehicle_xdot_tol) {
    Fz_rnom = tanh(4.0 * fabs
                   (carPlant_smoothUp_followerSto_B.VectorConcatenate1[1]));
    alfa_fl = (atan((carPlant_smoothUp_followerSto_P.LeadVehicle_a *
                     carPlant_smoothUp_followerSto_B.VectorConcatenate1[3] +
                     carPlant_smoothUp_followerSto_B.VectorConcatenate1[1]) /
                    (carPlant_smoothUp_followerSto_P.LeadVehicle_w[0] / 2.0 *
                     carPlant_smoothUp_followerSto_B.VectorConcatenate1[3] +
                     alfa_rr)) -
               carPlant_smoothUp_followerSto_B.VectorConcatenate_i[0]) * Fz_rnom;
    alfa_fr = (atan((carPlant_smoothUp_followerSto_P.LeadVehicle_a *
                     carPlant_smoothUp_followerSto_B.VectorConcatenate1[3] +
                     carPlant_smoothUp_followerSto_B.VectorConcatenate1[1]) /
                    (alfa_rr - carPlant_smoothUp_followerSto_P.LeadVehicle_w[0] /
                     2.0 * carPlant_smoothUp_followerSto_B.VectorConcatenate1[3]))
               - carPlant_smoothUp_followerSto_B.VectorConcatenate_i[1]) *
      Fz_rnom;
    alfa_rl = (atan((carPlant_smoothUp_followerSto_B.VectorConcatenate1[1] -
                     carPlant_smoothUp_followerSto_P.LeadVehicle_b *
                     carPlant_smoothUp_followerSto_B.VectorConcatenate1[3]) /
                    (carPlant_smoothUp_followerSto_P.LeadVehicle_w[1] / 2.0 *
                     carPlant_smoothUp_followerSto_B.VectorConcatenate1[3] +
                     alfa_rr)) -
               carPlant_smoothUp_followerSto_B.VectorConcatenate4_c[0]) *
      Fz_rnom;
    alfa_rr = (atan((carPlant_smoothUp_followerSto_B.VectorConcatenate1[1] -
                     carPlant_smoothUp_followerSto_P.LeadVehicle_b *
                     carPlant_smoothUp_followerSto_B.VectorConcatenate1[3]) /
                    (alfa_rr - carPlant_smoothUp_followerSto_P.LeadVehicle_w[1] /
                     2.0 * carPlant_smoothUp_followerSto_B.VectorConcatenate1[3]))
               - carPlant_smoothUp_followerSto_B.VectorConcatenate4_c[1]) *
      Fz_rnom;
  } else {
    Fz_rnom = carPlant_smoothUp_followerSto_P.LeadVehicle_w[0] / 2.0 *
      carPlant_smoothUp_followerSto_B.VectorConcatenate1[3];
    alfa_fr = carPlant_smoothUp_followerSto_P.LeadVehicle_a *
      carPlant_smoothUp_followerSto_B.VectorConcatenate1[3] +
      carPlant_smoothUp_followerSto_B.VectorConcatenate1[1];
    Fz_fnom = tanh(4.0 * carPlant_smoothUp_followerSto_B.VectorConcatenate1[0]);
    alfa_fl = (atan(alfa_fr / (Fz_rnom + alfa_rr)) -
               carPlant_smoothUp_followerSto_B.VectorConcatenate_i[0]) * Fz_fnom;
    alfa_fr = (atan(alfa_fr / (alfa_rr - Fz_rnom)) -
               carPlant_smoothUp_followerSto_B.VectorConcatenate_i[1]) * Fz_fnom;
    Fz_rnom = carPlant_smoothUp_followerSto_P.LeadVehicle_w[1] / 2.0 *
      carPlant_smoothUp_followerSto_B.VectorConcatenate1[3];
    Fz_idx_0 = carPlant_smoothUp_followerSto_B.VectorConcatenate1[1] -
      carPlant_smoothUp_followerSto_P.LeadVehicle_b *
      carPlant_smoothUp_followerSto_B.VectorConcatenate1[3];
    alfa_rl = (atan(Fz_idx_0 / (Fz_rnom + alfa_rr)) -
               carPlant_smoothUp_followerSto_B.VectorConcatenate4_c[0]) *
      Fz_fnom;
    alfa_rr = (atan(Fz_idx_0 / (alfa_rr - Fz_rnom)) -
               carPlant_smoothUp_followerSto_B.VectorConcatenate4_c[1]) *
      Fz_fnom;
  }

  Fz_idx_0 = 0.0;
  Fz_idx_1 = 0.0;
  Fz_idx_2 = 0.0;
  Fz_idx_3 = 0.0;
  for (iU = 0; iU < 11; iU++) {
    if (iU == 0) {
      Fz_rnom = -dx1 * dx3 * carPlant_smoothUp_followerSto_P.LeadVehicle_h;
      yddot = rtb_VectorConcatenate1_m_idx_0 *
        carPlant_smoothUp_followerSto_P.LeadVehicle_h;
      Fz_fnom = ((((carPlant_smoothUp_followerSto_P.LeadVehicle_g *
                    carPlant_smoothUp_followerSto_P.LeadVehicle_b - Fz_rnom) *
                   carPlant_smoothUp_followerSto_P.LeadVehicle_m +
                   rtb_VectorConcatenate1_k *
                   carPlant_smoothUp_followerSto_P.LeadVehicle_b) + yddot) -
                 rtb_sincos_o1_idx_1) / rtb_VectorConcatenate1_m_idx_0_;
      Fz_rnom = ((((Fz_rnom + carPlant_smoothUp_followerSto_P.LeadVehicle_g *
                    carPlant_smoothUp_followerSto_P.LeadVehicle_a) *
                   carPlant_smoothUp_followerSto_P.LeadVehicle_m +
                   rtb_VectorConcatenate1_k *
                   carPlant_smoothUp_followerSto_P.LeadVehicle_a) - yddot) +
                 rtb_sincos_o1_idx_1) / rtb_VectorConcatenate1_m_idx_0_;
      Fz_idx_1 = rtb_Integrator_d * dx3;
      Fz_idx_2 = rtb_sincos_o2_idx_1 *
        carPlant_smoothUp_followerSto_P.LeadVehicle_h;
      yddot = (Fz_idx_1 * carPlant_smoothUp_followerSto_P.LeadVehicle_m *
               carPlant_smoothUp_followerSto_P.LeadVehicle_h - Fz_idx_2) +
        rtb_sincos_o1_idx_0;
      Fz_idx_0 = (yddot + Fz_fnom) *
        (carPlant_smoothUp_followerSto_P.LeadVehicle_w[0] / 2.0 -
         carPlant_smoothUp_followerSto_P.LeadVehicle_d) /
        (carPlant_smoothUp_followerSto_P.LeadVehicle_w[0] / 2.0);
      Fz_idx_3 = (Fz_idx_1 * -carPlant_smoothUp_followerSto_P.LeadVehicle_m *
                  carPlant_smoothUp_followerSto_P.LeadVehicle_h + Fz_idx_2) +
        rtb_sincos_o1_idx_0;
      Fz_idx_1 = (Fz_idx_3 + Fz_fnom) *
        (carPlant_smoothUp_followerSto_P.LeadVehicle_w[0] / 2.0 +
         carPlant_smoothUp_followerSto_P.LeadVehicle_d) /
        (carPlant_smoothUp_followerSto_P.LeadVehicle_w[0] / 2.0);
      Fz_idx_2 = (yddot + Fz_rnom) *
        (carPlant_smoothUp_followerSto_P.LeadVehicle_w[1] / 2.0 -
         carPlant_smoothUp_followerSto_P.LeadVehicle_d) /
        (carPlant_smoothUp_followerSto_P.LeadVehicle_w[1] / 2.0);
      Fz_idx_3 = (Fz_idx_3 + Fz_rnom) *
        (carPlant_smoothUp_followerSto_P.LeadVehicle_w[1] / 2.0 +
         carPlant_smoothUp_followerSto_P.LeadVehicle_d) /
        (carPlant_smoothUp_followerSto_P.LeadVehicle_w[1] / 2.0);
      if (Fz_idx_0 < 0.0) {
        Fz_idx_0 = 0.0;
      }

      if (Fz_idx_1 < 0.0) {
        Fz_idx_1 = 0.0;
      }

      if (Fz_idx_2 < 0.0) {
        Fz_idx_2 = 0.0;
      }

      if (Fz_idx_3 < 0.0) {
        Fz_idx_3 = 0.0;
      }
    }

    carPlant_s_automlvehdynftiresat
      (-carPlant_smoothUp_followerSto_B.VectorConcatenate2_i[0] / 2.0 * alfa_fl *
       carPlant_smoothUp_followerSto_B.VectorConcatenate4_b[0] * Fz_idx_0 /
       carPlant_smoothUp_followerSto_P.LeadVehicle_Fznom,
       carPlant_smoothUp_followerSto_P.vehiclemodel_Fxtire_sat_d * Fz_idx_0 /
       carPlant_smoothUp_followerSto_P.LeadVehicle_Fznom,
       carPlant_smoothUp_followerSto_P.vehiclemodel_Fytire_sat_e * Fz_idx_0 /
       carPlant_smoothUp_followerSto_P.LeadVehicle_Fznom, &yddot, &rdot);
    carPlant_s_automlvehdynftiresat
      (-carPlant_smoothUp_followerSto_B.VectorConcatenate2_i[1] / 2.0 * alfa_fr *
       carPlant_smoothUp_followerSto_B.VectorConcatenate4_b[2] * Fz_idx_1 /
       carPlant_smoothUp_followerSto_P.LeadVehicle_Fznom,
       carPlant_smoothUp_followerSto_P.vehiclemodel_Fxtire_sat_d * Fz_idx_1 /
       carPlant_smoothUp_followerSto_P.LeadVehicle_Fznom,
       carPlant_smoothUp_followerSto_P.vehiclemodel_Fytire_sat_e * Fz_idx_1 /
       carPlant_smoothUp_followerSto_P.LeadVehicle_Fznom, &yddot, &Fz_fnom);
    carPlant_s_automlvehdynftiresat
      (-carPlant_smoothUp_followerSto_B.VectorConcatenate1_m[0] / 2.0 * alfa_rl *
       carPlant_smoothUp_followerSto_B.VectorConcatenate4_b[1] * Fz_idx_2 /
       carPlant_smoothUp_followerSto_P.LeadVehicle_Fznom,
       carPlant_smoothUp_followerSto_P.vehiclemodel_Fxtire_sat_d * Fz_idx_2 /
       carPlant_smoothUp_followerSto_P.LeadVehicle_Fznom,
       carPlant_smoothUp_followerSto_P.vehiclemodel_Fytire_sat_e * Fz_idx_2 /
       carPlant_smoothUp_followerSto_P.LeadVehicle_Fznom, &yddot, &Fz_rnom);
    carPlant_s_automlvehdynftiresat
      (-carPlant_smoothUp_followerSto_B.VectorConcatenate1_m[1] / 2.0 * alfa_rr *
       carPlant_smoothUp_followerSto_B.VectorConcatenate4_b[3] * Fz_idx_3 /
       carPlant_smoothUp_followerSto_P.LeadVehicle_Fznom,
       carPlant_smoothUp_followerSto_P.vehiclemodel_Fxtire_sat_d * Fz_idx_3 /
       carPlant_smoothUp_followerSto_P.LeadVehicle_Fznom,
       carPlant_smoothUp_followerSto_P.vehiclemodel_Fytire_sat_e * Fz_idx_3 /
       carPlant_smoothUp_followerSto_P.LeadVehicle_Fznom, &yddot, &Fz_idx_0);
    Fz_idx_1 = cos(delta_fl);
    Fz_idx_2 = sin(delta_fl);
    Fz_idx_3 = cos(delta_fr);
    b_Fy_fr_tmp = sin(delta_fr);
    b_Fy_rl_tmp = cos(delta_rl);
    b_Fy_rl_tmp_0 = sin(delta_rl);
    b_Fy_rl = Fz_rnom * b_Fy_rl_tmp - 0.0 * b_Fy_rl_tmp_0;
    b_Fy_rr_tmp = cos(delta_rr);
    b_Fy_rr_tmp_0 = sin(delta_rr);
    b_Fy_rr = Fz_idx_0 * b_Fy_rr_tmp - 0.0 * b_Fy_rr_tmp_0;
    yddot_tmp = (rdot * Fz_idx_1 - 0.0 * Fz_idx_2) + (Fz_fnom * Fz_idx_3 - 0.0 *
      b_Fy_fr_tmp);
    yddot_0 = (((yddot_tmp + b_Fy_rl) + b_Fy_rr) + rtb_sincos_o2_idx_1) /
      carPlant_smoothUp_followerSto_P.LeadVehicle_m + -rtb_Integrator_d * dx3;
    rdot_0 = (((((0.0 * Fz_idx_1 - rdot * Fz_idx_2) - (0.0 * Fz_idx_3 - Fz_fnom *
      b_Fy_fr_tmp)) * (carPlant_smoothUp_followerSto_P.LeadVehicle_w[0] / 2.0) +
                (yddot_tmp * carPlant_smoothUp_followerSto_P.LeadVehicle_a -
                 (b_Fy_rl + b_Fy_rr) *
                 carPlant_smoothUp_followerSto_P.LeadVehicle_b)) + ((0.0 *
      b_Fy_rl_tmp - Fz_rnom * b_Fy_rl_tmp_0) - (0.0 * b_Fy_rr_tmp - Fz_idx_0 *
      b_Fy_rr_tmp_0)) * (carPlant_smoothUp_followerSto_P.LeadVehicle_w[1] / 2.0))
              + rtb_sincos_o1_idx_2) /
      carPlant_smoothUp_followerSto_P.LeadVehicle_Izz;
    Fz_rnom = (0.0 - dx1 * dx3) * carPlant_smoothUp_followerSto_P.LeadVehicle_h;
    Fz_fnom = ((((carPlant_smoothUp_followerSto_P.LeadVehicle_g *
                  carPlant_smoothUp_followerSto_P.LeadVehicle_b - Fz_rnom) *
                 carPlant_smoothUp_followerSto_P.LeadVehicle_m +
                 rtb_VectorConcatenate1_k *
                 carPlant_smoothUp_followerSto_P.LeadVehicle_b) +
                rtb_VectorConcatenate1_m_idx_0 *
                carPlant_smoothUp_followerSto_P.LeadVehicle_h) -
               rtb_sincos_o1_idx_1) / rtb_VectorConcatenate1_m_idx_0_;
    Fz_rnom = ((((Fz_rnom + carPlant_smoothUp_followerSto_P.LeadVehicle_g *
                  carPlant_smoothUp_followerSto_P.LeadVehicle_a) *
                 carPlant_smoothUp_followerSto_P.LeadVehicle_m +
                 rtb_VectorConcatenate1_k *
                 carPlant_smoothUp_followerSto_P.LeadVehicle_a) -
                rtb_VectorConcatenate1_m_idx_0 *
                carPlant_smoothUp_followerSto_P.LeadVehicle_h) +
               rtb_sincos_o1_idx_1) / rtb_VectorConcatenate1_m_idx_0_;
    Fz_idx_1 = rtb_Integrator_d * dx3 + yddot_0;
    Fz_idx_2 = (Fz_idx_1 * carPlant_smoothUp_followerSto_P.LeadVehicle_m *
                carPlant_smoothUp_followerSto_P.LeadVehicle_h -
                rtb_sincos_o2_idx_1 *
                carPlant_smoothUp_followerSto_P.LeadVehicle_h) +
      rtb_sincos_o1_idx_0;
    Fz_idx_0 = (Fz_idx_2 + Fz_fnom) *
      (carPlant_smoothUp_followerSto_P.LeadVehicle_w[0] / 2.0 -
       carPlant_smoothUp_followerSto_P.LeadVehicle_d) /
      (carPlant_smoothUp_followerSto_P.LeadVehicle_w[0] / 2.0);
    Fz_idx_3 = (Fz_idx_1 * -carPlant_smoothUp_followerSto_P.LeadVehicle_m *
                carPlant_smoothUp_followerSto_P.LeadVehicle_h +
                rtb_sincos_o2_idx_1 *
                carPlant_smoothUp_followerSto_P.LeadVehicle_h) +
      rtb_sincos_o1_idx_0;
    Fz_idx_1 = (Fz_idx_3 + Fz_fnom) *
      (carPlant_smoothUp_followerSto_P.LeadVehicle_w[0] / 2.0 +
       carPlant_smoothUp_followerSto_P.LeadVehicle_d) /
      (carPlant_smoothUp_followerSto_P.LeadVehicle_w[0] / 2.0);
    Fz_idx_2 = (Fz_idx_2 + Fz_rnom) *
      (carPlant_smoothUp_followerSto_P.LeadVehicle_w[1] / 2.0 -
       carPlant_smoothUp_followerSto_P.LeadVehicle_d) /
      (carPlant_smoothUp_followerSto_P.LeadVehicle_w[1] / 2.0);
    Fz_idx_3 = (Fz_idx_3 + Fz_rnom) *
      (carPlant_smoothUp_followerSto_P.LeadVehicle_w[1] / 2.0 +
       carPlant_smoothUp_followerSto_P.LeadVehicle_d) /
      (carPlant_smoothUp_followerSto_P.LeadVehicle_w[1] / 2.0);
    if (Fz_idx_0 < 0.0) {
      Fz_idx_0 = 0.0;
    }

    if (Fz_idx_1 < 0.0) {
      Fz_idx_1 = 0.0;
    }

    if (Fz_idx_2 < 0.0) {
      Fz_idx_2 = 0.0;
    }

    if (Fz_idx_3 < 0.0) {
      Fz_idx_3 = 0.0;
    }
  }

  carPlant_smoothUp_followerSto_B.stateDer[0] = 0.0;
  carPlant_smoothUp_followerSto_B.stateDer[1] = yddot_0;
  carPlant_smoothUp_followerSto_B.stateDer[2] =
    carPlant_smoothUp_followerSto_B.VectorConcatenate1[3];
  carPlant_smoothUp_followerSto_B.stateDer[3] = rdot_0;

  /* MATLAB Function: '<S290>/COMB2I' */
  carPlant_smoothUp_follow_COMB2I
    (carPlant_smoothUp_followerSto_B.VectorConcatenate1,
     &carPlant_smoothUp_followerSto_B.sf_COMB2I_f);
  if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
    /* SignalConversion: '<S374>/ConcatBufferAtVector ConcatenateIn1' incorporates:
     *  Constant: '<S374>/Constant1'
     */
    carPlant_smoothUp_followerSto_B.VectorConcatenate_e[0] =
      carPlant_smoothUp_followerSto_P.LeadVehicle_sigma_f;

    /* SignalConversion: '<S374>/ConcatBufferAtVector ConcatenateIn2' incorporates:
     *  Constant: '<S374>/Constant1'
     */
    carPlant_smoothUp_followerSto_B.VectorConcatenate_e[1] =
      carPlant_smoothUp_followerSto_P.LeadVehicle_sigma_f;

    /* SignalConversion: '<S374>/ConcatBufferAtVector ConcatenateIn3' incorporates:
     *  Constant: '<S374>/Constant2'
     */
    carPlant_smoothUp_followerSto_B.VectorConcatenate_e[2] =
      carPlant_smoothUp_followerSto_P.LeadVehicle_sigma_r;

    /* SignalConversion: '<S374>/ConcatBufferAtVector ConcatenateIn4' incorporates:
     *  Constant: '<S374>/Constant2'
     */
    carPlant_smoothUp_followerSto_B.VectorConcatenate_e[3] =
      carPlant_smoothUp_followerSto_P.LeadVehicle_sigma_r;
  }

  /* Product: '<S377>/Product1' incorporates:
   *  Integrator: '<S377>/lateral'
   *  MATLAB Function: '<S11>/vehicle model'
   *  Sum: '<S377>/Sum'
   */
  carPlant_smoothUp_followerSto_B.Product1_l[0] = (alfa_fl -
    carPlant_smoothUp_followerSto_X.lateral_CSTATE_f[0]) * rtb_sincos_o2_idx_2 /
    carPlant_smoothUp_followerSto_B.VectorConcatenate_e[0];
  carPlant_smoothUp_followerSto_B.Product1_l[1] = (alfa_fr -
    carPlant_smoothUp_followerSto_X.lateral_CSTATE_f[1]) * rtb_sincos_o2_idx_2 /
    carPlant_smoothUp_followerSto_B.VectorConcatenate_e[1];
  carPlant_smoothUp_followerSto_B.Product1_l[2] = (alfa_rl -
    carPlant_smoothUp_followerSto_X.lateral_CSTATE_f[2]) *
    rtb_UnaryMinus1_g_idx_2 /
    carPlant_smoothUp_followerSto_B.VectorConcatenate_e[2];
  carPlant_smoothUp_followerSto_B.Product1_l[3] = (alfa_rr -
    carPlant_smoothUp_followerSto_X.lateral_CSTATE_f[3]) *
    rtb_UnaryMinus1_g_idx_2 /
    carPlant_smoothUp_followerSto_B.VectorConcatenate_e[3];
  if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
    /* Matfile logging */
    rt_UpdateTXYLogVars(carPlant_smoothUp_followerSt_M->rtwLogInfo,
                        (carPlant_smoothUp_followerSt_M->Timing.t));
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
    /* Update for Integrator: '<S378>/Integrator' */
    carPlant_smoothUp_followerSt_DW.Integrator_IWORK = 0;
    if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
      /* Update for DiscreteIntegrator: '<S22>/Integrator' incorporates:
       *  Constant: '<S5>/Constant'
       */
      carPlant_smoothUp_followerSt_DW.Integrator_IC_LOADING = 0U;
      carPlant_smoothUp_followerSt_DW.Integrator_DSTATE +=
        carPlant_smoothUp_followerSto_P.Integrator_gainval *
        carPlant_smoothUp_followerSto_B.uT;
      if (carPlant_smoothUp_followerSt_DW.Integrator_DSTATE >=
          carPlant_smoothUp_followerSto_P.Integrator_UpperSat) {
        carPlant_smoothUp_followerSt_DW.Integrator_DSTATE =
          carPlant_smoothUp_followerSto_P.Integrator_UpperSat;
      } else {
        if (carPlant_smoothUp_followerSt_DW.Integrator_DSTATE <=
            carPlant_smoothUp_followerSto_P.Integrator_LowerSat) {
          carPlant_smoothUp_followerSt_DW.Integrator_DSTATE =
            carPlant_smoothUp_followerSto_P.Integrator_LowerSat;
        }
      }

      if (carPlant_smoothUp_followerSto_P.Constant_Value_f > 0.0) {
        carPlant_smoothUp_followerSt_DW.Integrator_PrevResetState = 1;
      } else if (carPlant_smoothUp_followerSto_P.Constant_Value_f < 0.0) {
        carPlant_smoothUp_followerSt_DW.Integrator_PrevResetState = -1;
      } else if (carPlant_smoothUp_followerSto_P.Constant_Value_f == 0.0) {
        carPlant_smoothUp_followerSt_DW.Integrator_PrevResetState = 0;
      } else {
        carPlant_smoothUp_followerSt_DW.Integrator_PrevResetState = 2;
      }

      /* End of Update for DiscreteIntegrator: '<S22>/Integrator' */

      /* Update for DiscreteIntegrator: '<S26>/Integrator' incorporates:
       *  Constant: '<S6>/Constant'
       */
      carPlant_smoothUp_followerSt_DW.Integrator_IC_LOADING_g = 0U;
      carPlant_smoothUp_followerSt_DW.Integrator_DSTATE_m +=
        carPlant_smoothUp_followerSto_P.Integrator_gainval_p *
        carPlant_smoothUp_followerSto_B.uT_k;
      if (carPlant_smoothUp_followerSt_DW.Integrator_DSTATE_m >=
          carPlant_smoothUp_followerSto_P.Integrator_UpperSat_h) {
        carPlant_smoothUp_followerSt_DW.Integrator_DSTATE_m =
          carPlant_smoothUp_followerSto_P.Integrator_UpperSat_h;
      } else {
        if (carPlant_smoothUp_followerSt_DW.Integrator_DSTATE_m <=
            carPlant_smoothUp_followerSto_P.Integrator_LowerSat_j) {
          carPlant_smoothUp_followerSt_DW.Integrator_DSTATE_m =
            carPlant_smoothUp_followerSto_P.Integrator_LowerSat_j;
        }
      }

      if (carPlant_smoothUp_followerSto_P.Constant_Value_i > 0.0) {
        carPlant_smoothUp_followerSt_DW.Integrator_PrevResetState_c = 1;
      } else if (carPlant_smoothUp_followerSto_P.Constant_Value_i < 0.0) {
        carPlant_smoothUp_followerSt_DW.Integrator_PrevResetState_c = -1;
      } else if (carPlant_smoothUp_followerSto_P.Constant_Value_i == 0.0) {
        carPlant_smoothUp_followerSt_DW.Integrator_PrevResetState_c = 0;
      } else {
        carPlant_smoothUp_followerSt_DW.Integrator_PrevResetState_c = 2;
      }

      /* End of Update for DiscreteIntegrator: '<S26>/Integrator' */
    }

    /* Update for Integrator: '<S182>/Integrator' */
    carPlant_smoothUp_followerSt_DW.Integrator_IWORK_h = 0;

    /* Update for Integrator: '<S290>/Integrator' */
    carPlant_smoothUp_followerSt_DW.Integrator_IWORK_p = 0;

    /* Update for Integrator: '<S94>/Integrator' */
    carPlant_smoothUp_followerSt_DW.Integrator_IWORK_l = 0;
    if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
      /* Update for Backlash: '<S205>/Backlash' */
      carPlant_smoothUp_followerSt_DW.PrevY =
        carPlant_smoothUp_followerSto_B.Backlash;

      /* Update for Backlash: '<S186>/Backlash' */
      carPlant_smoothUp_followerSt_DW.PrevY_e =
        carPlant_smoothUp_followerSto_B.Backlash_a;
    }
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(carPlant_smoothUp_followerSt_M)) {
    /* signal main to stop simulation */
    {                                  /* Sample time: [0.0s, 0.0s] */
      if ((rtmGetTFinal(carPlant_smoothUp_followerSt_M)!=-1) &&
          !((rtmGetTFinal(carPlant_smoothUp_followerSt_M)-
             (((carPlant_smoothUp_followerSt_M->Timing.clockTick1+
                carPlant_smoothUp_followerSt_M->Timing.clockTickH1* 4294967296.0))
              * 0.01)) > (((carPlant_smoothUp_followerSt_M->Timing.clockTick1+
                            carPlant_smoothUp_followerSt_M->Timing.clockTickH1*
                            4294967296.0)) * 0.01) * (DBL_EPSILON))) {
        rtmSetErrorStatus(carPlant_smoothUp_followerSt_M, "Simulation finished");
      }
    }

    rt_ertODEUpdateContinuousStates(&carPlant_smoothUp_followerSt_M->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++carPlant_smoothUp_followerSt_M->Timing.clockTick0)) {
      ++carPlant_smoothUp_followerSt_M->Timing.clockTickH0;
    }

    carPlant_smoothUp_followerSt_M->Timing.t[0] = rtsiGetSolverStopTime
      (&carPlant_smoothUp_followerSt_M->solverInfo);

    {
      /* Update absolute timer for sample time: [0.01s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.01, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       * Timer of this task consists of two 32 bit unsigned integers.
       * The two integers represent the low bits Timing.clockTick1 and the high bits
       * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
       */
      carPlant_smoothUp_followerSt_M->Timing.clockTick1++;
      if (!carPlant_smoothUp_followerSt_M->Timing.clockTick1) {
        carPlant_smoothUp_followerSt_M->Timing.clockTickH1++;
      }
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void carPlant_smoothUp_followerStopper_derivatives(void)
{
  XDot_carPlant_smoothUp_follow_T *_rtXdot;
  _rtXdot = ((XDot_carPlant_smoothUp_follow_T *)
             carPlant_smoothUp_followerSt_M->derivs);

  /* Derivatives for Integrator: '<S378>/Integrator' */
  _rtXdot->Integrator_CSTATE[0] = carPlant_smoothUp_followerSto_B.stateDer[0];

  /* Derivatives for Integrator: '<S182>/Integrator' */
  _rtXdot->Integrator_CSTATE_j[0] = carPlant_smoothUp_followerSto_B.stateDer_p[0];

  /* Derivatives for Integrator: '<S378>/Integrator' */
  _rtXdot->Integrator_CSTATE[1] = carPlant_smoothUp_followerSto_B.stateDer[1];

  /* Derivatives for Integrator: '<S182>/Integrator' */
  _rtXdot->Integrator_CSTATE_j[1] = carPlant_smoothUp_followerSto_B.stateDer_p[1];

  /* Derivatives for Integrator: '<S378>/Integrator' */
  _rtXdot->Integrator_CSTATE[2] = carPlant_smoothUp_followerSto_B.stateDer[2];

  /* Derivatives for Integrator: '<S182>/Integrator' */
  _rtXdot->Integrator_CSTATE_j[2] = carPlant_smoothUp_followerSto_B.stateDer_p[2];

  /* Derivatives for Integrator: '<S378>/Integrator' */
  _rtXdot->Integrator_CSTATE[3] = carPlant_smoothUp_followerSto_B.stateDer[3];

  /* Derivatives for Integrator: '<S182>/Integrator' */
  _rtXdot->Integrator_CSTATE_j[3] = carPlant_smoothUp_followerSto_B.stateDer_p[3];

  /* Derivatives for Integrator: '<S290>/Integrator' */
  _rtXdot->Integrator_CSTATE_o[0] =
    carPlant_smoothUp_followerSto_B.sf_COMB2I_f.y[0];

  /* Derivatives for Integrator: '<S94>/Integrator' */
  _rtXdot->Integrator_CSTATE_a[0] = carPlant_smoothUp_followerSto_B.sf_COMB2I.y
    [0];

  /* Derivatives for Integrator: '<S290>/Integrator' */
  _rtXdot->Integrator_CSTATE_o[1] =
    carPlant_smoothUp_followerSto_B.sf_COMB2I_f.y[1];

  /* Derivatives for Integrator: '<S94>/Integrator' */
  _rtXdot->Integrator_CSTATE_a[1] = carPlant_smoothUp_followerSto_B.sf_COMB2I.y
    [1];

  /* Derivatives for Integrator: '<S181>/lateral' */
  _rtXdot->lateral_CSTATE[0] = carPlant_smoothUp_followerSto_B.Product1[0];

  /* Derivatives for Integrator: '<S377>/lateral' */
  _rtXdot->lateral_CSTATE_f[0] = carPlant_smoothUp_followerSto_B.Product1_l[0];

  /* Derivatives for Integrator: '<S181>/lateral' */
  _rtXdot->lateral_CSTATE[1] = carPlant_smoothUp_followerSto_B.Product1[1];

  /* Derivatives for Integrator: '<S377>/lateral' */
  _rtXdot->lateral_CSTATE_f[1] = carPlant_smoothUp_followerSto_B.Product1_l[1];

  /* Derivatives for Integrator: '<S181>/lateral' */
  _rtXdot->lateral_CSTATE[2] = carPlant_smoothUp_followerSto_B.Product1[2];

  /* Derivatives for Integrator: '<S377>/lateral' */
  _rtXdot->lateral_CSTATE_f[2] = carPlant_smoothUp_followerSto_B.Product1_l[2];

  /* Derivatives for Integrator: '<S181>/lateral' */
  _rtXdot->lateral_CSTATE[3] = carPlant_smoothUp_followerSto_B.Product1[3];

  /* Derivatives for Integrator: '<S377>/lateral' */
  _rtXdot->lateral_CSTATE_f[3] = carPlant_smoothUp_followerSto_B.Product1_l[3];
}

/* Model initialize function */
void carPlant_smoothUp_followerStopper_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  carPlant_smoothUp_followerSto_P.Integrator_UpperSat = rtInf;
  carPlant_smoothUp_followerSto_P.Integrator_LowerSat = rtMinusInf;
  carPlant_smoothUp_followerSto_P.Saturation_UpperSat = rtInf;
  carPlant_smoothUp_followerSto_P.Saturation_LowerSat = rtMinusInf;
  carPlant_smoothUp_followerSto_P.Integrator_UpperSat_h = rtInf;
  carPlant_smoothUp_followerSto_P.Integrator_LowerSat_j = rtMinusInf;
  carPlant_smoothUp_followerSto_P.Saturation_UpperSat_h = rtInf;
  carPlant_smoothUp_followerSto_P.Saturation_LowerSat_i = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)carPlant_smoothUp_followerSt_M, 0,
                sizeof(RT_MODEL_carPlant_smoothUp_fo_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&carPlant_smoothUp_followerSt_M->solverInfo,
                          &carPlant_smoothUp_followerSt_M->Timing.simTimeStep);
    rtsiSetTPtr(&carPlant_smoothUp_followerSt_M->solverInfo, &rtmGetTPtr
                (carPlant_smoothUp_followerSt_M));
    rtsiSetStepSizePtr(&carPlant_smoothUp_followerSt_M->solverInfo,
                       &carPlant_smoothUp_followerSt_M->Timing.stepSize0);
    rtsiSetdXPtr(&carPlant_smoothUp_followerSt_M->solverInfo,
                 &carPlant_smoothUp_followerSt_M->derivs);
    rtsiSetContStatesPtr(&carPlant_smoothUp_followerSt_M->solverInfo, (real_T **)
                         &carPlant_smoothUp_followerSt_M->contStates);
    rtsiSetNumContStatesPtr(&carPlant_smoothUp_followerSt_M->solverInfo,
      &carPlant_smoothUp_followerSt_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&carPlant_smoothUp_followerSt_M->solverInfo,
      &carPlant_smoothUp_followerSt_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr
      (&carPlant_smoothUp_followerSt_M->solverInfo,
       &carPlant_smoothUp_followerSt_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr
      (&carPlant_smoothUp_followerSt_M->solverInfo,
       &carPlant_smoothUp_followerSt_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&carPlant_smoothUp_followerSt_M->solverInfo,
                          (&rtmGetErrorStatus(carPlant_smoothUp_followerSt_M)));
    rtsiSetRTModelPtr(&carPlant_smoothUp_followerSt_M->solverInfo,
                      carPlant_smoothUp_followerSt_M);
  }

  rtsiSetSimTimeStep(&carPlant_smoothUp_followerSt_M->solverInfo,
                     MAJOR_TIME_STEP);
  carPlant_smoothUp_followerSt_M->intgData.y =
    carPlant_smoothUp_followerSt_M->odeY;
  carPlant_smoothUp_followerSt_M->intgData.f[0] =
    carPlant_smoothUp_followerSt_M->odeF[0];
  carPlant_smoothUp_followerSt_M->intgData.f[1] =
    carPlant_smoothUp_followerSt_M->odeF[1];
  carPlant_smoothUp_followerSt_M->intgData.f[2] =
    carPlant_smoothUp_followerSt_M->odeF[2];
  carPlant_smoothUp_followerSt_M->contStates = ((X_carPlant_smoothUp_followerS_T
    *) &carPlant_smoothUp_followerSto_X);
  rtsiSetSolverData(&carPlant_smoothUp_followerSt_M->solverInfo, (void *)
                    &carPlant_smoothUp_followerSt_M->intgData);
  rtsiSetSolverName(&carPlant_smoothUp_followerSt_M->solverInfo,"ode3");
  rtmSetTPtr(carPlant_smoothUp_followerSt_M,
             &carPlant_smoothUp_followerSt_M->Timing.tArray[0]);
  rtmSetTFinal(carPlant_smoothUp_followerSt_M, -1);
  carPlant_smoothUp_followerSt_M->Timing.stepSize0 = 0.01;
  rtmSetFirstInitCond(carPlant_smoothUp_followerSt_M, 1);

  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;
    rt_DataLoggingInfo.loggingInterval = NULL;
    carPlant_smoothUp_followerSt_M->rtwLogInfo = &rt_DataLoggingInfo;
  }

  /* Setup for data logging */
  {
    rtliSetLogXSignalInfo(carPlant_smoothUp_followerSt_M->rtwLogInfo, (NULL));
    rtliSetLogXSignalPtrs(carPlant_smoothUp_followerSt_M->rtwLogInfo, (NULL));
    rtliSetLogT(carPlant_smoothUp_followerSt_M->rtwLogInfo, "tout");
    rtliSetLogX(carPlant_smoothUp_followerSt_M->rtwLogInfo, "");
    rtliSetLogXFinal(carPlant_smoothUp_followerSt_M->rtwLogInfo, "");
    rtliSetLogVarNameModifier(carPlant_smoothUp_followerSt_M->rtwLogInfo, "rt_");
    rtliSetLogFormat(carPlant_smoothUp_followerSt_M->rtwLogInfo, 4);
    rtliSetLogMaxRows(carPlant_smoothUp_followerSt_M->rtwLogInfo, 0);
    rtliSetLogDecimation(carPlant_smoothUp_followerSt_M->rtwLogInfo, 1);
    rtliSetLogY(carPlant_smoothUp_followerSt_M->rtwLogInfo, "");
    rtliSetLogYSignalInfo(carPlant_smoothUp_followerSt_M->rtwLogInfo, (NULL));
    rtliSetLogYSignalPtrs(carPlant_smoothUp_followerSt_M->rtwLogInfo, (NULL));
  }

  /* block I/O */
  (void) memset(((void *) &carPlant_smoothUp_followerSto_B), 0,
                sizeof(B_carPlant_smoothUp_followerS_T));

  /* states (continuous) */
  {
    (void) memset((void *)&carPlant_smoothUp_followerSto_X, 0,
                  sizeof(X_carPlant_smoothUp_followerS_T));
  }

  /* states (dwork) */
  (void) memset((void *)&carPlant_smoothUp_followerSt_DW, 0,
                sizeof(DW_carPlant_smoothUp_follower_T));

  /* Matfile logging */
  rt_StartDataLoggingWithStartTime(carPlant_smoothUp_followerSt_M->rtwLogInfo,
    0.0, rtmGetTFinal(carPlant_smoothUp_followerSt_M),
    carPlant_smoothUp_followerSt_M->Timing.stepSize0, (&rtmGetErrorStatus
    (carPlant_smoothUp_followerSt_M)));

  {
    static const char_T tmp[7] = { 'c', 'm', 'd', '_', 'v', 'e', 'l' };

    static const char_T tmp_0[12] = { '/', 'd', 'x', '_', 'a', 'c', 't', 'i',
      'v', 'a', 't', 'e' };

    static const char_T tmp_1[7] = { '/', 'd', 'x', '_', 'm', 'i', 'n' };

    char_T tmp_2[13];
    char_T tmp_3[8];
    int32_T i;

    /* Start for Probe: '<S19>/Probe' */
    carPlant_smoothUp_followerSto_B.Probe[0] = 0.01;
    carPlant_smoothUp_followerSto_B.Probe[1] = 0.0;

    /* Start for FromWorkspace: '<Root>/From Workspace' */
    {
      static real_T pTimeValues0[] = { 0.0, 0.033333333333333333,
        0.066666666666666666, 0.1, 0.13333333333333333, 0.16666666666666666, 0.2,
        0.23333333333333334, 0.26666666666666666, 0.3, 0.33333333333333331,
        0.36666666666666664, 0.4, 0.43333333333333335, 0.46666666666666667, 0.5,
        0.53333333333333333, 0.56666666666666665, 0.6, 0.6333333333333333,
        0.66666666666666663, 0.7, 0.73333333333333328, 0.76666666666666661, 0.8,
        0.83333333333333337, 0.8666666666666667, 0.9, 0.93333333333333335,
        0.96666666666666667, 1.0, 1.0333333333333332, 1.0666666666666667, 1.1,
        1.1333333333333333, 1.1666666666666667, 1.2, 1.2333333333333334,
        1.2666666666666666, 1.3, 1.3333333333333333, 1.3666666666666667, 1.4,
        1.4333333333333333, 1.4666666666666666, 1.5, 1.5333333333333332,
        1.5666666666666667, 1.6, 1.6333333333333333, 1.6666666666666667, 1.7,
        1.7333333333333334, 1.7666666666666666, 1.8, 1.8333333333333333,
        1.8666666666666667, 1.9, 1.9333333333333333, 1.9666666666666666, 2.0,
        2.0333333333333332, 2.0666666666666664, 2.1, 2.1333333333333333,
        2.1666666666666665, 2.2, 2.2333333333333334, 2.2666666666666666, 2.3,
        2.3333333333333335, 2.3666666666666667, 2.4, 2.4333333333333331,
        2.4666666666666668, 2.5, 2.5333333333333332, 2.5666666666666664, 2.6,
        2.6333333333333333, 2.6666666666666665, 2.7, 2.7333333333333334,
        2.7666666666666666, 2.8, 2.8333333333333335, 2.8666666666666667, 2.9,
        2.9333333333333331, 2.9666666666666668, 3.0, 3.0333333333333332,
        3.0666666666666664, 3.1, 3.1333333333333333, 3.1666666666666665, 3.2,
        3.2333333333333334, 3.2666666666666666, 3.3, 3.3333333333333335,
        3.3666666666666667, 3.4, 3.4333333333333331, 3.4666666666666668, 3.5,
        3.5333333333333332, 3.5666666666666664, 3.6, 3.6333333333333333,
        3.6666666666666665, 3.6999999999999997, 3.7333333333333334,
        3.7666666666666666, 3.8, 3.8333333333333335, 3.8666666666666667, 3.9,
        3.9333333333333331, 3.9666666666666668, 4.0, 4.0333333333333332,
        4.0666666666666664, 4.1, 4.1333333333333329, 4.166666666666667, 4.2,
        4.2333333333333334, 4.2666666666666666, 4.3, 4.333333333333333,
        4.3666666666666663, 4.4, 4.4333333333333336, 4.4666666666666668, 4.5,
        4.5333333333333332, 4.5666666666666664, 4.6, 4.6333333333333329,
        4.666666666666667, 4.7, 4.7333333333333334, 4.7666666666666666, 4.8,
        4.833333333333333, 4.8666666666666663, 4.9, 4.9333333333333336,
        4.9666666666666668, 5.0, 5.0333333333333332, 5.0666666666666664, 5.1,
        5.1333333333333329, 5.166666666666667, 5.2, 5.2333333333333334,
        5.2666666666666666, 5.3, 5.333333333333333, 5.3666666666666663, 5.4,
        5.4333333333333336, 5.4666666666666668, 5.5, 5.5333333333333332,
        5.5666666666666664, 5.6, 5.6333333333333329, 5.666666666666667, 5.7,
        5.7333333333333334, 5.7666666666666666, 5.8, 5.833333333333333,
        5.8666666666666663, 5.9, 5.9333333333333336, 5.9666666666666668, 6.0,
        6.0333333333333332, 6.0666666666666664, 6.1, 6.1333333333333329,
        6.166666666666667, 6.2, 6.2333333333333334, 6.2666666666666666, 6.3,
        6.333333333333333, 6.3666666666666663, 6.4, 6.4333333333333336,
        6.4666666666666668, 6.5, 6.5333333333333332, 6.5666666666666664, 6.6,
        6.6333333333333329, 6.666666666666667, 6.7, 6.7333333333333334,
        6.7666666666666666, 6.8, 6.833333333333333, 6.8666666666666663,
        6.8999999999999995, 6.9333333333333336, 6.9666666666666668, 7.0,
        7.0333333333333332, 7.0666666666666664, 7.1, 7.1333333333333329,
        7.166666666666667, 7.2, 7.2333333333333334, 7.2666666666666666, 7.3,
        7.333333333333333, 7.3666666666666663, 7.3999999999999995,
        7.4333333333333336, 7.4666666666666668, 7.5, 7.5333333333333332,
        7.5666666666666664, 7.6, 7.6333333333333329, 7.666666666666667, 7.7,
        7.7333333333333334, 7.7666666666666666, 7.8, 7.833333333333333,
        7.8666666666666663, 7.8999999999999995, 7.9333333333333336,
        7.9666666666666668, 8.0, 8.0333333333333332, 8.0666666666666664, 8.1,
        8.1333333333333329, 8.1666666666666661, 8.2, 8.2333333333333325,
        8.2666666666666657, 8.3, 8.3333333333333339, 8.3666666666666671, 8.4,
        8.4333333333333336, 8.4666666666666668, 8.5, 8.5333333333333332,
        8.5666666666666664, 8.6, 8.6333333333333329, 8.6666666666666661, 8.7,
        8.7333333333333325, 8.7666666666666657, 8.8, 8.8333333333333339,
        8.8666666666666671, 8.9, 8.9333333333333336, 8.9666666666666668, 9.0,
        9.0333333333333332, 9.0666666666666664, 9.1, 9.1333333333333329,
        9.1666666666666661, 9.2, 9.2333333333333325, 9.2666666666666657, 9.3,
        9.3333333333333339, 9.3666666666666671, 9.4, 9.4333333333333336,
        9.4666666666666668, 9.5, 9.5333333333333332, 9.5666666666666664, 9.6,
        9.6333333333333329, 9.6666666666666661, 9.7, 9.7333333333333325,
        9.7666666666666657, 9.8, 9.8333333333333339, 9.8666666666666671, 9.9,
        9.9333333333333336, 9.9666666666666668, 10.0, 10.033333333333333,
        10.066666666666666, 10.1, 10.133333333333333, 10.166666666666666, 10.2,
        10.233333333333333, 10.266666666666666, 10.3, 10.333333333333334,
        10.366666666666667, 10.4, 10.433333333333334, 10.466666666666667, 10.5,
        10.533333333333333, 10.566666666666666, 10.6, 10.633333333333333,
        10.666666666666666, 10.7, 10.733333333333333, 10.766666666666666, 10.8,
        10.833333333333334, 10.866666666666667, 10.9, 10.933333333333334,
        10.966666666666667, 11.0, 11.033333333333333, 11.066666666666666, 11.1,
        11.133333333333333, 11.166666666666666, 11.2, 11.233333333333333,
        11.266666666666666, 11.3, 11.333333333333334, 11.366666666666667, 11.4,
        11.433333333333334, 11.466666666666667, 11.5, 11.533333333333333,
        11.566666666666666, 11.6, 11.633333333333333, 11.666666666666666, 11.7,
        11.733333333333333, 11.766666666666666, 11.8, 11.833333333333334,
        11.866666666666667, 11.9, 11.933333333333334, 11.966666666666667, 12.0,
        12.033333333333333, 12.066666666666666, 12.1, 12.133333333333333,
        12.166666666666666, 12.2, 12.233333333333333, 12.266666666666666, 12.3,
        12.333333333333334, 12.366666666666667, 12.4, 12.433333333333334,
        12.466666666666667, 12.5, 12.533333333333333, 12.566666666666666, 12.6,
        12.633333333333333, 12.666666666666666, 12.7, 12.733333333333333,
        12.766666666666666, 12.8, 12.833333333333334, 12.866666666666667, 12.9,
        12.933333333333334, 12.966666666666667, 13.0, 13.033333333333333,
        13.066666666666666, 13.1, 13.133333333333333, 13.166666666666666, 13.2,
        13.233333333333333, 13.266666666666666, 13.299999999999999,
        13.333333333333334, 13.366666666666667, 13.4, 13.433333333333334,
        13.466666666666667, 13.5, 13.533333333333333, 13.566666666666666, 13.6,
        13.633333333333333, 13.666666666666666, 13.7, 13.733333333333333,
        13.766666666666666, 13.799999999999999, 13.833333333333334,
        13.866666666666667, 13.9, 13.933333333333334, 13.966666666666667, 14.0,
        14.033333333333333, 14.066666666666666, 14.1, 14.133333333333333,
        14.166666666666666, 14.2, 14.233333333333333, 14.266666666666666,
        14.299999999999999, 14.333333333333334, 14.366666666666667, 14.4,
        14.433333333333334, 14.466666666666667, 14.5, 14.533333333333333,
        14.566666666666666, 14.6, 14.633333333333333, 14.666666666666666, 14.7,
        14.733333333333333, 14.766666666666666, 14.799999999999999,
        14.833333333333334, 14.866666666666667, 14.9, 14.933333333333334,
        14.966666666666667, 15.0, 15.033333333333333, 15.066666666666666, 15.1,
        15.133333333333333, 15.166666666666666, 15.2, 15.233333333333333,
        15.266666666666666, 15.299999999999999, 15.333333333333334,
        15.366666666666667, 15.4, 15.433333333333334, 15.466666666666667, 15.5,
        15.533333333333333, 15.566666666666666, 15.6, 15.633333333333333,
        15.666666666666666, 15.7, 15.733333333333333, 15.766666666666666,
        15.799999999999999, 15.833333333333334, 15.866666666666667, 15.9,
        15.933333333333334, 15.966666666666667, 16.0, 16.033333333333331,
        16.066666666666666, 16.1, 16.133333333333333, 16.166666666666668, 16.2,
        16.233333333333334, 16.266666666666666, 16.3, 16.333333333333332,
        16.366666666666667, 16.4, 16.433333333333334, 16.466666666666665, 16.5,
        16.533333333333331, 16.566666666666666, 16.6, 16.633333333333333,
        16.666666666666668, 16.7, 16.733333333333334, 16.766666666666666, 16.8,
        16.833333333333332, 16.866666666666667, 16.9, 16.933333333333334,
        16.966666666666665, 17.0, 17.033333333333331, 17.066666666666666, 17.1,
        17.133333333333333, 17.166666666666668, 17.2, 17.233333333333334,
        17.266666666666666, 17.3, 17.333333333333332, 17.366666666666667, 17.4,
        17.433333333333334, 17.466666666666665, 17.5, 17.533333333333331,
        17.566666666666666, 17.6, 17.633333333333333, 17.666666666666668, 17.7,
        17.733333333333334, 17.766666666666666, 17.8, 17.833333333333332,
        17.866666666666667, 17.9, 17.933333333333334, 17.966666666666665, 18.0,
        18.033333333333331, 18.066666666666666, 18.1, 18.133333333333333,
        18.166666666666668, 18.2, 18.233333333333334, 18.266666666666666, 18.3,
        18.333333333333332, 18.366666666666667, 18.4, 18.433333333333334,
        18.466666666666665, 18.5, 18.533333333333331, 18.566666666666666, 18.6,
        18.633333333333333, 18.666666666666668, 18.7, 18.733333333333334,
        18.766666666666666, 18.8, 18.833333333333332, 18.866666666666667, 18.9,
        18.933333333333334, 18.966666666666665, 19.0, 19.033333333333331,
        19.066666666666666, 19.1, 19.133333333333333, 19.166666666666668, 19.2,
        19.233333333333334, 19.266666666666666, 19.3, 19.333333333333332,
        19.366666666666667, 19.4, 19.433333333333334, 19.466666666666665, 19.5,
        19.533333333333331, 19.566666666666666, 19.6, 19.633333333333333,
        19.666666666666668, 19.7, 19.733333333333334, 19.766666666666666, 19.8,
        19.833333333333332, 19.866666666666667, 19.9, 19.933333333333334,
        19.966666666666665, 20.0, 20.033333333333331, 20.066666666666666, 20.1,
        20.133333333333333, 20.166666666666668, 20.2, 20.233333333333334,
        20.266666666666666, 20.3, 20.333333333333332, 20.366666666666667, 20.4,
        20.433333333333334, 20.466666666666665, 20.5, 20.533333333333331,
        20.566666666666666, 20.6, 20.633333333333333, 20.666666666666668, 20.7,
        20.733333333333334, 20.766666666666666, 20.8, 20.833333333333332,
        20.866666666666667, 20.9, 20.933333333333334, 20.966666666666665, 21.0,
        21.033333333333331, 21.066666666666666, 21.1, 21.133333333333333,
        21.166666666666668, 21.2, 21.233333333333334, 21.266666666666666, 21.3,
        21.333333333333332, 21.366666666666667, 21.4, 21.433333333333334,
        21.466666666666665, 21.5, 21.533333333333331, 21.566666666666666, 21.6,
        21.633333333333333, 21.666666666666668, 21.7, 21.733333333333334,
        21.766666666666666, 21.8, 21.833333333333332, 21.866666666666667, 21.9,
        21.933333333333334, 21.966666666666665, 22.0, 22.033333333333331,
        22.066666666666666, 22.1, 22.133333333333333, 22.166666666666668, 22.2,
        22.233333333333334, 22.266666666666666, 22.3, 22.333333333333332,
        22.366666666666667, 22.4, 22.433333333333334, 22.466666666666665, 22.5,
        22.533333333333331, 22.566666666666666, 22.6, 22.633333333333333,
        22.666666666666668, 22.7, 22.733333333333334, 22.766666666666666, 22.8,
        22.833333333333332, 22.866666666666667, 22.9, 22.933333333333334,
        22.966666666666665, 23.0, 23.033333333333331, 23.066666666666666, 23.1,
        23.133333333333333, 23.166666666666668, 23.2, 23.233333333333334,
        23.266666666666666, 23.3, 23.333333333333332, 23.366666666666667, 23.4,
        23.433333333333334, 23.466666666666665, 23.5, 23.533333333333331,
        23.566666666666666, 23.6, 23.633333333333333, 23.666666666666668, 23.7,
        23.733333333333334, 23.766666666666666, 23.8, 23.833333333333332,
        23.866666666666667, 23.9, 23.933333333333334, 23.966666666666665, 24.0,
        24.033333333333331, 24.066666666666666, 24.1, 24.133333333333333,
        24.166666666666668, 24.2, 24.233333333333334, 24.266666666666666, 24.3,
        24.333333333333332, 24.366666666666667, 24.4, 24.433333333333334,
        24.466666666666665, 24.5, 24.533333333333331, 24.566666666666666, 24.6,
        24.633333333333333, 24.666666666666668, 24.7, 24.733333333333334,
        24.766666666666666, 24.8, 24.833333333333332, 24.866666666666667, 24.9,
        24.933333333333334, 24.966666666666665, 25.0, 25.033333333333331,
        25.066666666666666, 25.1, 25.133333333333333, 25.166666666666668, 25.2,
        25.233333333333334, 25.266666666666666, 25.3, 25.333333333333332,
        25.366666666666667, 25.4, 25.433333333333334, 25.466666666666665, 25.5,
        25.533333333333331, 25.566666666666666, 25.6, 25.633333333333333,
        25.666666666666668, 25.7, 25.733333333333334, 25.766666666666666, 25.8,
        25.833333333333332, 25.866666666666667, 25.9, 25.933333333333334,
        25.966666666666665, 26.0, 26.033333333333331, 26.066666666666666,
        26.099999999999998, 26.133333333333333, 26.166666666666668, 26.2,
        26.233333333333334, 26.266666666666666, 26.3, 26.333333333333332,
        26.366666666666667, 26.4, 26.433333333333334, 26.466666666666665, 26.5,
        26.533333333333331, 26.566666666666666, 26.599999999999998,
        26.633333333333333, 26.666666666666668, 26.7, 26.733333333333334,
        26.766666666666666, 26.8, 26.833333333333332, 26.866666666666667, 26.9,
        26.933333333333334, 26.966666666666665, 27.0, 27.033333333333331,
        27.066666666666666, 27.099999999999998, 27.133333333333333,
        27.166666666666668, 27.2, 27.233333333333334, 27.266666666666666, 27.3,
        27.333333333333332, 27.366666666666667, 27.4, 27.433333333333334,
        27.466666666666665, 27.5, 27.533333333333331, 27.566666666666666,
        27.599999999999998, 27.633333333333333, 27.666666666666668, 27.7,
        27.733333333333334, 27.766666666666666, 27.8, 27.833333333333332,
        27.866666666666667, 27.9, 27.933333333333334, 27.966666666666665, 28.0,
        28.033333333333331, 28.066666666666666, 28.099999999999998,
        28.133333333333333, 28.166666666666668, 28.2, 28.233333333333334,
        28.266666666666666, 28.3, 28.333333333333332, 28.366666666666667, 28.4,
        28.433333333333334, 28.466666666666665, 28.5, 28.533333333333331,
        28.566666666666666, 28.599999999999998, 28.633333333333333,
        28.666666666666668, 28.7, 28.733333333333334, 28.766666666666666, 28.8,
        28.833333333333332, 28.866666666666667, 28.9, 28.933333333333334,
        28.966666666666665, 29.0, 29.033333333333331, 29.066666666666666,
        29.099999999999998, 29.133333333333333, 29.166666666666668, 29.2,
        29.233333333333334, 29.266666666666666, 29.3, 29.333333333333332,
        29.366666666666667, 29.4, 29.433333333333334, 29.466666666666665, 29.5,
        29.533333333333331, 29.566666666666666, 29.599999999999998,
        29.633333333333333, 29.666666666666668, 29.7, 29.733333333333334,
        29.766666666666666, 29.8, 29.833333333333332, 29.866666666666667, 29.9,
        29.933333333333334, 29.966666666666665, 30.0, 30.033333333333331,
        30.066666666666666, 30.099999999999998, 30.133333333333333,
        30.166666666666668, 30.2, 30.233333333333334, 30.266666666666666, 30.3,
        30.333333333333332, 30.366666666666667, 30.4, 30.433333333333334,
        30.466666666666665, 30.5, 30.533333333333331, 30.566666666666666,
        30.599999999999998, 30.633333333333333, 30.666666666666668, 30.7,
        30.733333333333334, 30.766666666666666, 30.8, 30.833333333333332,
        30.866666666666667, 30.9, 30.933333333333334, 30.966666666666665, 31.0,
        31.033333333333331, 31.066666666666666, 31.099999999999998,
        31.133333333333333, 31.166666666666668, 31.2, 31.233333333333334,
        31.266666666666666, 31.3, 31.333333333333332, 31.366666666666667, 31.4,
        31.433333333333334, 31.466666666666665, 31.5, 31.533333333333331,
        31.566666666666666, 31.599999999999998, 31.633333333333333,
        31.666666666666668, 31.7, 31.733333333333334, 31.766666666666666, 31.8,
        31.833333333333332, 31.866666666666667, 31.9, 31.933333333333334,
        31.966666666666665, 32.0, 32.033333333333331, 32.066666666666663, 32.1,
        32.133333333333333, 32.166666666666664, 32.2, 32.233333333333334,
        32.266666666666666, 32.3, 32.333333333333336, 32.366666666666667, 32.4,
        32.43333333333333, 32.466666666666669, 32.5, 32.533333333333331,
        32.566666666666663, 32.6, 32.633333333333333, 32.666666666666664, 32.7,
        32.733333333333334, 32.766666666666666, 32.8, 32.833333333333336,
        32.866666666666667, 32.9, 32.93333333333333, 32.966666666666669, 33.0,
        33.033333333333331, 33.066666666666663, 33.1, 33.133333333333333,
        33.166666666666664, 33.2, 33.233333333333334, 33.266666666666666, 33.3,
        33.333333333333336, 33.366666666666667, 33.4, 33.43333333333333,
        33.466666666666669, 33.5, 33.533333333333331, 33.566666666666663, 33.6,
        33.633333333333333, 33.666666666666664, 33.7, 33.733333333333334,
        33.766666666666666, 33.8, 33.833333333333336, 33.866666666666667, 33.9,
        33.93333333333333, 33.966666666666669, 34.0, 34.033333333333331,
        34.066666666666663, 34.1, 34.133333333333333, 34.166666666666664, 34.2,
        34.233333333333334, 34.266666666666666, 34.3, 34.333333333333336,
        34.366666666666667, 34.4, 34.43333333333333, 34.466666666666669, 34.5,
        34.533333333333331, 34.566666666666663, 34.6, 34.633333333333333,
        34.666666666666664, 34.7, 34.733333333333334, 34.766666666666666, 34.8,
        34.833333333333336, 34.866666666666667, 34.9, 34.93333333333333,
        34.966666666666669, 35.0, 35.033333333333331, 35.066666666666663, 35.1,
        35.133333333333333, 35.166666666666664, 35.2, 35.233333333333334,
        35.266666666666666, 35.3, 35.333333333333336, 35.366666666666667, 35.4,
        35.43333333333333, 35.466666666666669, 35.5, 35.533333333333331,
        35.566666666666663, 35.6, 35.633333333333333, 35.666666666666664, 35.7,
        35.733333333333334, 35.766666666666666, 35.8, 35.833333333333336,
        35.866666666666667, 35.9, 35.93333333333333, 35.966666666666669, 36.0,
        36.033333333333331, 36.066666666666663, 36.1, 36.133333333333333,
        36.166666666666664, 36.2, 36.233333333333334, 36.266666666666666, 36.3,
        36.333333333333336, 36.366666666666667, 36.4, 36.43333333333333,
        36.466666666666669, 36.5, 36.533333333333331, 36.566666666666663, 36.6,
        36.633333333333333, 36.666666666666664, 36.7, 36.733333333333334,
        36.766666666666666, 36.8, 36.833333333333336, 36.866666666666667, 36.9,
        36.93333333333333, 36.966666666666669, 37.0, 37.033333333333331,
        37.066666666666663, 37.1, 37.133333333333333, 37.166666666666664, 37.2,
        37.233333333333334, 37.266666666666666, 37.3, 37.333333333333336,
        37.366666666666667, 37.4, 37.43333333333333, 37.466666666666669, 37.5,
        37.533333333333331, 37.566666666666663, 37.6, 37.633333333333333,
        37.666666666666664, 37.7, 37.733333333333334, 37.766666666666666, 37.8,
        37.833333333333336, 37.866666666666667, 37.9, 37.93333333333333,
        37.966666666666669, 38.0, 38.033333333333331, 38.066666666666663, 38.1,
        38.133333333333333, 38.166666666666664, 38.2, 38.233333333333334,
        38.266666666666666, 38.3, 38.333333333333336, 38.366666666666667, 38.4,
        38.43333333333333, 38.466666666666669, 38.5, 38.533333333333331,
        38.566666666666663, 38.6, 38.633333333333333, 38.666666666666664, 38.7,
        38.733333333333334, 38.766666666666666, 38.8, 38.833333333333336,
        38.866666666666667, 38.9, 38.93333333333333, 38.966666666666669, 39.0,
        39.033333333333331, 39.066666666666663, 39.1, 39.133333333333333,
        39.166666666666664, 39.2, 39.233333333333334, 39.266666666666666, 39.3,
        39.333333333333336, 39.366666666666667, 39.4, 39.43333333333333,
        39.466666666666669, 39.5, 39.533333333333331, 39.566666666666663, 39.6,
        39.633333333333333, 39.666666666666664, 39.7, 39.733333333333334,
        39.766666666666666, 39.8, 39.833333333333336, 39.866666666666667, 39.9,
        39.93333333333333, 39.966666666666669, 40.0, 40.033333333333331,
        40.066666666666663, 40.1, 40.133333333333333, 40.166666666666664, 40.2,
        40.233333333333334, 40.266666666666666, 40.3, 40.333333333333336,
        40.366666666666667, 40.4, 40.43333333333333, 40.466666666666669, 40.5,
        40.533333333333331, 40.566666666666663, 40.6, 40.633333333333333,
        40.666666666666664, 40.7, 40.733333333333334, 40.766666666666666, 40.8,
        40.833333333333336, 40.866666666666667, 40.9, 40.93333333333333,
        40.966666666666669, 41.0, 41.033333333333331, 41.066666666666663, 41.1,
        41.133333333333333, 41.166666666666664, 41.2, 41.233333333333334,
        41.266666666666666, 41.3, 41.333333333333336, 41.366666666666667, 41.4,
        41.43333333333333, 41.466666666666669, 41.5, 41.533333333333331,
        41.566666666666663, 41.6, 41.633333333333333, 41.666666666666664, 41.7,
        41.733333333333334, 41.766666666666666, 41.8, 41.833333333333336,
        41.866666666666667, 41.9, 41.93333333333333, 41.966666666666669, 42.0,
        42.033333333333331, 42.066666666666663, 42.1, 42.133333333333333,
        42.166666666666664, 42.2, 42.233333333333334, 42.266666666666666, 42.3,
        42.333333333333336, 42.366666666666667, 42.4, 42.43333333333333,
        42.466666666666669, 42.5, 42.533333333333331, 42.566666666666663, 42.6,
        42.633333333333333, 42.666666666666664, 42.7, 42.733333333333334,
        42.766666666666666, 42.8, 42.833333333333336, 42.866666666666667, 42.9,
        42.93333333333333, 42.966666666666669, 43.0, 43.033333333333331,
        43.066666666666663, 43.1, 43.133333333333333, 43.166666666666664, 43.2,
        43.233333333333334, 43.266666666666666, 43.3, 43.333333333333336,
        43.366666666666667, 43.4, 43.43333333333333, 43.466666666666669, 43.5,
        43.533333333333331, 43.566666666666663, 43.6, 43.633333333333333,
        43.666666666666664, 43.7, 43.733333333333334, 43.766666666666666, 43.8,
        43.833333333333336, 43.866666666666667, 43.9, 43.93333333333333,
        43.966666666666669, 44.0, 44.033333333333331, 44.066666666666663, 44.1,
        44.133333333333333, 44.166666666666664, 44.2, 44.233333333333334,
        44.266666666666666, 44.3, 44.333333333333336, 44.366666666666667, 44.4,
        44.43333333333333, 44.466666666666669, 44.5, 44.533333333333331,
        44.566666666666663, 44.6, 44.633333333333333, 44.666666666666664, 44.7,
        44.733333333333334, 44.766666666666666, 44.8, 44.833333333333336,
        44.866666666666667, 44.9, 44.93333333333333, 44.966666666666669, 45.0,
        45.033333333333331, 45.066666666666663, 45.1, 45.133333333333333,
        45.166666666666664, 45.2, 45.233333333333334, 45.266666666666666, 45.3,
        45.333333333333336, 45.366666666666667, 45.4, 45.43333333333333,
        45.466666666666669, 45.5, 45.533333333333331, 45.566666666666663, 45.6,
        45.633333333333333, 45.666666666666664, 45.7, 45.733333333333334,
        45.766666666666666, 45.8, 45.833333333333336, 45.866666666666667, 45.9,
        45.93333333333333, 45.966666666666669, 46.0, 46.033333333333331,
        46.066666666666663, 46.1, 46.133333333333333, 46.166666666666664, 46.2,
        46.233333333333334, 46.266666666666666, 46.3, 46.333333333333336,
        46.366666666666667, 46.4, 46.43333333333333, 46.466666666666669, 46.5,
        46.533333333333331, 46.566666666666663, 46.6, 46.633333333333333,
        46.666666666666664, 46.7, 46.733333333333334, 46.766666666666666, 46.8,
        46.833333333333336, 46.866666666666667, 46.9, 46.93333333333333,
        46.966666666666669, 47.0, 47.033333333333331, 47.066666666666663, 47.1,
        47.133333333333333, 47.166666666666664, 47.2, 47.233333333333334,
        47.266666666666666, 47.3, 47.333333333333336, 47.366666666666667, 47.4,
        47.43333333333333, 47.466666666666669, 47.5, 47.533333333333331,
        47.566666666666663, 47.6, 47.633333333333333, 47.666666666666664, 47.7,
        47.733333333333334, 47.766666666666666, 47.8, 47.833333333333336,
        47.866666666666667, 47.9, 47.93333333333333, 47.966666666666669, 48.0,
        48.033333333333331, 48.066666666666663, 48.1, 48.133333333333333,
        48.166666666666664, 48.2, 48.233333333333334, 48.266666666666666, 48.3,
        48.333333333333336, 48.366666666666667, 48.4, 48.43333333333333,
        48.466666666666669, 48.5, 48.533333333333331, 48.566666666666663, 48.6,
        48.633333333333333, 48.666666666666664, 48.7, 48.733333333333334,
        48.766666666666666, 48.8, 48.833333333333336, 48.866666666666667, 48.9,
        48.93333333333333, 48.966666666666669, 49.0, 49.033333333333331,
        49.066666666666663, 49.1, 49.133333333333333, 49.166666666666664, 49.2,
        49.233333333333334, 49.266666666666666, 49.3, 49.333333333333336,
        49.366666666666667, 49.4, 49.43333333333333, 49.466666666666669, 49.5,
        49.533333333333331, 49.566666666666663, 49.6, 49.633333333333333,
        49.666666666666664, 49.7, 49.733333333333334, 49.766666666666666, 49.8,
        49.833333333333336, 49.866666666666667, 49.9, 49.93333333333333,
        49.966666666666669, 50.0, 50.033333333333331, 50.066666666666663, 50.1,
        50.133333333333333, 50.166666666666664, 50.2, 50.233333333333334,
        50.266666666666666, 50.3, 50.333333333333336, 50.366666666666667, 50.4,
        50.43333333333333, 50.466666666666669, 50.5, 50.533333333333331,
        50.566666666666663, 50.6, 50.633333333333333, 50.666666666666664, 50.7,
        50.733333333333334, 50.766666666666666, 50.8, 50.833333333333336,
        50.866666666666667, 50.9, 50.93333333333333, 50.966666666666669, 51.0,
        51.033333333333331, 51.066666666666663, 51.1, 51.133333333333333,
        51.166666666666664, 51.2, 51.233333333333334, 51.266666666666666, 51.3,
        51.333333333333336, 51.366666666666667, 51.4, 51.43333333333333,
        51.466666666666669, 51.5, 51.533333333333331, 51.566666666666663, 51.6,
        51.633333333333333, 51.666666666666664, 51.699999999999996,
        51.733333333333334, 51.766666666666666, 51.8, 51.833333333333336,
        51.866666666666667, 51.9, 51.93333333333333, 51.966666666666669, 52.0,
        52.033333333333331, 52.066666666666663, 52.1, 52.133333333333333,
        52.166666666666664, 52.199999999999996, 52.233333333333334,
        52.266666666666666, 52.3, 52.333333333333336, 52.366666666666667, 52.4,
        52.43333333333333, 52.466666666666669, 52.5, 52.533333333333331,
        52.566666666666663, 52.6, 52.633333333333333, 52.666666666666664,
        52.699999999999996, 52.733333333333334, 52.766666666666666, 52.8,
        52.833333333333336, 52.866666666666667, 52.9, 52.93333333333333,
        52.966666666666669, 53.0, 53.033333333333331, 53.066666666666663, 53.1,
        53.133333333333333, 53.166666666666664, 53.199999999999996,
        53.233333333333334, 53.266666666666666, 53.3, 53.333333333333336,
        53.366666666666667, 53.4, 53.43333333333333, 53.466666666666669, 53.5,
        53.533333333333331, 53.566666666666663, 53.6, 53.633333333333333,
        53.666666666666664, 53.699999999999996, 53.733333333333334,
        53.766666666666666, 53.8, 53.833333333333336, 53.866666666666667, 53.9,
        53.93333333333333, 53.966666666666669, 54.0, 54.033333333333331,
        54.066666666666663, 54.1, 54.133333333333333, 54.166666666666664,
        54.199999999999996, 54.233333333333334, 54.266666666666666, 54.3,
        54.333333333333336, 54.366666666666667, 54.4, 54.43333333333333,
        54.466666666666669, 54.5, 54.533333333333331, 54.566666666666663, 54.6,
        54.633333333333333, 54.666666666666664, 54.699999999999996,
        54.733333333333334, 54.766666666666666, 54.8, 54.833333333333336,
        54.866666666666667, 54.9, 54.93333333333333, 54.966666666666669, 55.0,
        55.033333333333331, 55.066666666666663, 55.1, 55.133333333333333,
        55.166666666666664, 55.199999999999996, 55.233333333333334,
        55.266666666666666, 55.3, 55.333333333333336, 55.366666666666667, 55.4,
        55.43333333333333, 55.466666666666669, 55.5, 55.533333333333331,
        55.566666666666663, 55.6, 55.633333333333333, 55.666666666666664,
        55.699999999999996, 55.733333333333334, 55.766666666666666, 55.8,
        55.833333333333336, 55.866666666666667, 55.9, 55.93333333333333,
        55.966666666666669, 56.0, 56.033333333333331, 56.066666666666663, 56.1,
        56.133333333333333, 56.166666666666664, 56.199999999999996,
        56.233333333333334, 56.266666666666666, 56.3, 56.333333333333336,
        56.366666666666667, 56.4, 56.43333333333333, 56.466666666666669, 56.5,
        56.533333333333331, 56.566666666666663, 56.6, 56.633333333333333,
        56.666666666666664, 56.699999999999996, 56.733333333333334,
        56.766666666666666, 56.8, 56.833333333333336, 56.866666666666667, 56.9,
        56.93333333333333, 56.966666666666669, 57.0, 57.033333333333331,
        57.066666666666663, 57.1, 57.133333333333333, 57.166666666666664,
        57.199999999999996, 57.233333333333334, 57.266666666666666, 57.3,
        57.333333333333336, 57.366666666666667, 57.4, 57.43333333333333,
        57.466666666666669, 57.5, 57.533333333333331, 57.566666666666663, 57.6,
        57.633333333333333, 57.666666666666664, 57.699999999999996,
        57.733333333333334, 57.766666666666666, 57.8, 57.833333333333336,
        57.866666666666667, 57.9, 57.93333333333333, 57.966666666666669, 58.0,
        58.033333333333331, 58.066666666666663, 58.1, 58.133333333333333,
        58.166666666666664, 58.199999999999996, 58.233333333333334,
        58.266666666666666, 58.3, 58.333333333333336, 58.366666666666667, 58.4,
        58.43333333333333, 58.466666666666669, 58.5, 58.533333333333331,
        58.566666666666663, 58.6, 58.633333333333333, 58.666666666666664,
        58.699999999999996, 58.733333333333334, 58.766666666666666, 58.8,
        58.833333333333336, 58.866666666666667, 58.9, 58.93333333333333,
        58.966666666666669, 59.0, 59.033333333333331, 59.066666666666663, 59.1,
        59.133333333333333, 59.166666666666664, 59.199999999999996,
        59.233333333333334, 59.266666666666666, 59.3, 59.333333333333336,
        59.366666666666667, 59.4, 59.43333333333333, 59.466666666666669, 59.5,
        59.533333333333331, 59.566666666666663, 59.6, 59.633333333333333,
        59.666666666666664, 59.699999999999996, 59.733333333333334,
        59.766666666666666, 59.8, 59.833333333333336, 59.866666666666667, 59.9,
        59.93333333333333, 59.966666666666669, 60.0, 60.033333333333331,
        60.066666666666663, 60.1, 60.133333333333333, 60.166666666666664,
        60.199999999999996, 60.233333333333334, 60.266666666666666, 60.3,
        60.333333333333336, 60.366666666666667, 60.4, 60.43333333333333,
        60.466666666666669, 60.5, 60.533333333333331, 60.566666666666663, 60.6,
        60.633333333333333, 60.666666666666664, 60.699999999999996,
        60.733333333333334, 60.766666666666666, 60.8, 60.833333333333336,
        60.866666666666667, 60.9, 60.93333333333333, 60.966666666666669, 61.0,
        61.033333333333331, 61.066666666666663, 61.1, 61.133333333333333,
        61.166666666666664, 61.199999999999996, 61.233333333333334,
        61.266666666666666, 61.3, 61.333333333333336, 61.366666666666667, 61.4,
        61.43333333333333, 61.466666666666669, 61.5, 61.533333333333331,
        61.566666666666663, 61.6, 61.633333333333333, 61.666666666666664,
        61.699999999999996, 61.733333333333334, 61.766666666666666, 61.8,
        61.833333333333336, 61.866666666666667, 61.9, 61.93333333333333,
        61.966666666666669, 62.0, 62.033333333333331, 62.066666666666663, 62.1,
        62.133333333333333, 62.166666666666664, 62.199999999999996,
        62.233333333333334, 62.266666666666666, 62.3, 62.333333333333336,
        62.366666666666667, 62.4, 62.43333333333333, 62.466666666666669, 62.5,
        62.533333333333331, 62.566666666666663, 62.6, 62.633333333333333,
        62.666666666666664, 62.699999999999996, 62.733333333333334,
        62.766666666666666, 62.8, 62.833333333333336, 62.866666666666667, 62.9,
        62.93333333333333, 62.966666666666669, 63.0, 63.033333333333331,
        63.066666666666663, 63.1, 63.133333333333333, 63.166666666666664,
        63.199999999999996, 63.233333333333334, 63.266666666666666, 63.3,
        63.333333333333336, 63.366666666666667, 63.4, 63.43333333333333,
        63.466666666666669, 63.5, 63.533333333333331, 63.566666666666663, 63.6,
        63.633333333333333, 63.666666666666664, 63.699999999999996,
        63.733333333333334, 63.766666666666666, 63.8, 63.833333333333336,
        63.866666666666667, 63.9, 63.93333333333333, 63.966666666666669, 64.0,
        64.033333333333331, 64.066666666666663, 64.1, 64.133333333333326,
        64.166666666666671, 64.2, 64.233333333333334, 64.266666666666666, 64.3,
        64.333333333333329, 64.36666666666666, 64.4, 64.433333333333337,
        64.466666666666669, 64.5, 64.533333333333331, 64.566666666666663, 64.6,
        64.633333333333326, 64.666666666666671, 64.7, 64.733333333333334,
        64.766666666666666, 64.8, 64.833333333333329, 64.86666666666666, 64.9,
        64.933333333333337, 64.966666666666669, 65.0, 65.033333333333331,
        65.066666666666663, 65.1, 65.133333333333326, 65.166666666666671, 65.2,
        65.233333333333334, 65.266666666666666, 65.3, 65.333333333333329,
        65.36666666666666, 65.4, 65.433333333333337, 65.466666666666669, 65.5,
        65.533333333333331, 65.566666666666663, 65.6, 65.633333333333326,
        65.666666666666671, 65.7, 65.733333333333334, 65.766666666666666, 65.8,
        65.833333333333329, 65.86666666666666, 65.9, 65.933333333333337,
        65.966666666666669, 66.0, 66.033333333333331, 66.066666666666663, 66.1,
        66.133333333333326, 66.166666666666671, 66.2, 66.233333333333334,
        66.266666666666666, 66.3, 66.333333333333329, 66.36666666666666, 66.4,
        66.433333333333337, 66.466666666666669, 66.5, 66.533333333333331,
        66.566666666666663, 66.6, 66.633333333333326, 66.666666666666671, 66.7,
        66.733333333333334, 66.766666666666666, 66.8, 66.833333333333329,
        66.86666666666666, 66.9, 66.933333333333337, 66.966666666666669, 67.0,
        67.033333333333331, 67.066666666666663, 67.1, 67.133333333333326,
        67.166666666666671, 67.2, 67.233333333333334, 67.266666666666666, 67.3,
        67.333333333333329, 67.36666666666666, 67.4, 67.433333333333337,
        67.466666666666669, 67.5, 67.533333333333331, 67.566666666666663, 67.6,
        67.633333333333326, 67.666666666666671, 67.7, 67.733333333333334,
        67.766666666666666, 67.8, 67.833333333333329, 67.86666666666666, 67.9,
        67.933333333333337, 67.966666666666669, 68.0, 68.033333333333331,
        68.066666666666663, 68.1, 68.133333333333326, 68.166666666666671, 68.2,
        68.233333333333334, 68.266666666666666, 68.3, 68.333333333333329,
        68.36666666666666, 68.4, 68.433333333333337, 68.466666666666669, 68.5,
        68.533333333333331, 68.566666666666663, 68.6, 68.633333333333326,
        68.666666666666671, 68.7, 68.733333333333334, 68.766666666666666, 68.8,
        68.833333333333329, 68.86666666666666, 68.9, 68.933333333333337,
        68.966666666666669, 69.0, 69.033333333333331, 69.066666666666663, 69.1,
        69.133333333333326, 69.166666666666671, 69.2, 69.233333333333334,
        69.266666666666666, 69.3, 69.333333333333329, 69.36666666666666, 69.4,
        69.433333333333337, 69.466666666666669, 69.5, 69.533333333333331,
        69.566666666666663, 69.6, 69.633333333333326, 69.666666666666671, 69.7,
        69.733333333333334, 69.766666666666666, 69.8, 69.833333333333329,
        69.86666666666666, 69.9, 69.933333333333337, 69.966666666666669, 70.0,
        70.033333333333331, 70.066666666666663, 70.1, 70.133333333333326,
        70.166666666666671, 70.2, 70.233333333333334, 70.266666666666666, 70.3,
        70.333333333333329, 70.36666666666666, 70.4, 70.433333333333337,
        70.466666666666669, 70.5, 70.533333333333331, 70.566666666666663, 70.6,
        70.633333333333326, 70.666666666666671, 70.7, 70.733333333333334,
        70.766666666666666, 70.8, 70.833333333333329, 70.86666666666666, 70.9,
        70.933333333333337, 70.966666666666669, 71.0, 71.033333333333331,
        71.066666666666663, 71.1, 71.133333333333326, 71.166666666666671, 71.2,
        71.233333333333334, 71.266666666666666, 71.3, 71.333333333333329,
        71.36666666666666, 71.4, 71.433333333333337, 71.466666666666669, 71.5,
        71.533333333333331, 71.566666666666663, 71.6, 71.633333333333326,
        71.666666666666671, 71.7, 71.733333333333334, 71.766666666666666, 71.8,
        71.833333333333329, 71.86666666666666, 71.9, 71.933333333333337,
        71.966666666666669, 72.0, 72.033333333333331, 72.066666666666663, 72.1,
        72.133333333333326, 72.166666666666671, 72.2, 72.233333333333334,
        72.266666666666666, 72.3, 72.333333333333329, 72.36666666666666, 72.4,
        72.433333333333337, 72.466666666666669, 72.5, 72.533333333333331,
        72.566666666666663, 72.6, 72.633333333333326, 72.666666666666671, 72.7,
        72.733333333333334, 72.766666666666666, 72.8, 72.833333333333329,
        72.86666666666666, 72.9, 72.933333333333337, 72.966666666666669, 73.0,
        73.033333333333331, 73.066666666666663, 73.1, 73.133333333333326,
        73.166666666666671, 73.2, 73.233333333333334, 73.266666666666666, 73.3,
        73.333333333333329, 73.36666666666666, 73.4, 73.433333333333337,
        73.466666666666669, 73.5, 73.533333333333331, 73.566666666666663, 73.6,
        73.633333333333326, 73.666666666666671, 73.7, 73.733333333333334,
        73.766666666666666, 73.8, 73.833333333333329, 73.86666666666666, 73.9,
        73.933333333333337, 73.966666666666669, 74.0, 74.033333333333331,
        74.066666666666663, 74.1, 74.133333333333326, 74.166666666666671, 74.2,
        74.233333333333334, 74.266666666666666, 74.3, 74.333333333333329,
        74.36666666666666, 74.4, 74.433333333333337, 74.466666666666669, 74.5,
        74.533333333333331, 74.566666666666663, 74.6, 74.633333333333326,
        74.666666666666671, 74.7, 74.733333333333334, 74.766666666666666, 74.8,
        74.833333333333329, 74.86666666666666, 74.9, 74.933333333333337,
        74.966666666666669, 75.0, 75.033333333333331, 75.066666666666663, 75.1,
        75.133333333333326, 75.166666666666671, 75.2, 75.233333333333334,
        75.266666666666666, 75.3, 75.333333333333329, 75.36666666666666, 75.4,
        75.433333333333337, 75.466666666666669, 75.5, 75.533333333333331,
        75.566666666666663, 75.6, 75.633333333333326, 75.666666666666671, 75.7,
        75.733333333333334, 75.766666666666666, 75.8, 75.833333333333329,
        75.86666666666666, 75.9, 75.933333333333337, 75.966666666666669, 76.0,
        76.033333333333331, 76.066666666666663, 76.1, 76.133333333333326,
        76.166666666666671, 76.2, 76.233333333333334, 76.266666666666666, 76.3,
        76.333333333333329, 76.36666666666666, 76.4, 76.433333333333337,
        76.466666666666669, 76.5, 76.533333333333331, 76.566666666666663, 76.6,
        76.633333333333326, 76.666666666666671, 76.7, 76.733333333333334,
        76.766666666666666, 76.8, 76.833333333333329, 76.86666666666666, 76.9,
        76.933333333333337, 76.966666666666669, 77.0, 77.033333333333331,
        77.066666666666663, 77.1, 77.133333333333326, 77.166666666666671, 77.2,
        77.233333333333334, 77.266666666666666, 77.3, 77.333333333333329,
        77.36666666666666, 77.4, 77.433333333333337, 77.466666666666669, 77.5,
        77.533333333333331, 77.566666666666663, 77.6, 77.633333333333326,
        77.666666666666671, 77.7, 77.733333333333334, 77.766666666666666, 77.8,
        77.833333333333329, 77.86666666666666, 77.9, 77.933333333333337,
        77.966666666666669, 78.0, 78.033333333333331, 78.066666666666663, 78.1,
        78.133333333333326, 78.166666666666671, 78.2, 78.233333333333334,
        78.266666666666666, 78.3, 78.333333333333329, 78.36666666666666, 78.4,
        78.433333333333337, 78.466666666666669, 78.5, 78.533333333333331,
        78.566666666666663, 78.6, 78.633333333333326, 78.666666666666671, 78.7,
        78.733333333333334, 78.766666666666666, 78.8, 78.833333333333329,
        78.86666666666666, 78.9, 78.933333333333337, 78.966666666666669, 79.0,
        79.033333333333331, 79.066666666666663, 79.1, 79.133333333333326,
        79.166666666666671, 79.2, 79.233333333333334, 79.266666666666666, 79.3,
        79.333333333333329, 79.36666666666666, 79.4, 79.433333333333337,
        79.466666666666669, 79.5, 79.533333333333331, 79.566666666666663, 79.6,
        79.633333333333326, 79.666666666666671, 79.7, 79.733333333333334,
        79.766666666666666, 79.8, 79.833333333333329, 79.86666666666666, 79.9,
        79.933333333333337, 79.966666666666669, 80.0, 80.033333333333331,
        80.066666666666663, 80.1, 80.133333333333326, 80.166666666666671, 80.2,
        80.233333333333334, 80.266666666666666, 80.3, 80.333333333333329,
        80.36666666666666, 80.4, 80.433333333333337, 80.466666666666669, 80.5,
        80.533333333333331, 80.566666666666663, 80.6, 80.633333333333326,
        80.666666666666671, 80.7, 80.733333333333334, 80.766666666666666, 80.8,
        80.833333333333329, 80.86666666666666, 80.9, 80.933333333333337,
        80.966666666666669, 81.0, 81.033333333333331, 81.066666666666663, 81.1,
        81.133333333333326, 81.166666666666671, 81.2, 81.233333333333334,
        81.266666666666666, 81.3, 81.333333333333329, 81.36666666666666, 81.4,
        81.433333333333337, 81.466666666666669, 81.5, 81.533333333333331,
        81.566666666666663, 81.6, 81.633333333333326, 81.666666666666671, 81.7,
        81.733333333333334, 81.766666666666666, 81.8, 81.833333333333329,
        81.86666666666666, 81.9, 81.933333333333337, 81.966666666666669, 82.0,
        82.033333333333331, 82.066666666666663, 82.1, 82.133333333333326,
        82.166666666666671, 82.2, 82.233333333333334, 82.266666666666666, 82.3,
        82.333333333333329, 82.36666666666666, 82.4, 82.433333333333337,
        82.466666666666669, 82.5, 82.533333333333331, 82.566666666666663, 82.6,
        82.633333333333326, 82.666666666666671, 82.7, 82.733333333333334,
        82.766666666666666, 82.8, 82.833333333333329, 82.86666666666666, 82.9,
        82.933333333333337, 82.966666666666669, 83.0, 83.033333333333331,
        83.066666666666663, 83.1, 83.133333333333326, 83.166666666666671, 83.2,
        83.233333333333334, 83.266666666666666, 83.3, 83.333333333333329,
        83.36666666666666, 83.4, 83.433333333333337, 83.466666666666669, 83.5,
        83.533333333333331, 83.566666666666663, 83.6, 83.633333333333326,
        83.666666666666671, 83.7, 83.733333333333334, 83.766666666666666, 83.8,
        83.833333333333329, 83.86666666666666, 83.9, 83.933333333333337,
        83.966666666666669, 84.0, 84.033333333333331, 84.066666666666663, 84.1,
        84.133333333333326, 84.166666666666671, 84.2, 84.233333333333334,
        84.266666666666666, 84.3, 84.333333333333329, 84.36666666666666, 84.4,
        84.433333333333337, 84.466666666666669, 84.5, 84.533333333333331,
        84.566666666666663, 84.6, 84.633333333333326, 84.666666666666671, 84.7,
        84.733333333333334, 84.766666666666666, 84.8, 84.833333333333329,
        84.86666666666666, 84.9, 84.933333333333337, 84.966666666666669, 85.0,
        85.033333333333331, 85.066666666666663, 85.1, 85.133333333333326,
        85.166666666666671, 85.2, 85.233333333333334, 85.266666666666666, 85.3,
        85.333333333333329, 85.36666666666666, 85.4, 85.433333333333337,
        85.466666666666669, 85.5, 85.533333333333331, 85.566666666666663, 85.6,
        85.633333333333326, 85.666666666666671, 85.7, 85.733333333333334,
        85.766666666666666, 85.8, 85.833333333333329, 85.86666666666666, 85.9,
        85.933333333333337, 85.966666666666669, 86.0, 86.033333333333331,
        86.066666666666663, 86.1, 86.133333333333326, 86.166666666666671, 86.2,
        86.233333333333334, 86.266666666666666, 86.3, 86.333333333333329,
        86.36666666666666, 86.4, 86.433333333333337, 86.466666666666669, 86.5,
        86.533333333333331, 86.566666666666663, 86.6, 86.633333333333326,
        86.666666666666671, 86.7, 86.733333333333334, 86.766666666666666, 86.8,
        86.833333333333329, 86.86666666666666, 86.9, 86.933333333333337,
        86.966666666666669, 87.0, 87.033333333333331, 87.066666666666663, 87.1,
        87.133333333333326, 87.166666666666671, 87.2, 87.233333333333334,
        87.266666666666666, 87.3, 87.333333333333329, 87.36666666666666, 87.4,
        87.433333333333337, 87.466666666666669, 87.5, 87.533333333333331,
        87.566666666666663, 87.6, 87.633333333333326, 87.666666666666671, 87.7,
        87.733333333333334, 87.766666666666666, 87.8, 87.833333333333329,
        87.86666666666666, 87.9, 87.933333333333337, 87.966666666666669, 88.0,
        88.033333333333331, 88.066666666666663, 88.1, 88.133333333333326,
        88.166666666666671, 88.2, 88.233333333333334, 88.266666666666666, 88.3,
        88.333333333333329, 88.36666666666666, 88.4, 88.433333333333337,
        88.466666666666669, 88.5, 88.533333333333331, 88.566666666666663, 88.6,
        88.633333333333326, 88.666666666666671, 88.7, 88.733333333333334,
        88.766666666666666, 88.8, 88.833333333333329, 88.86666666666666, 88.9,
        88.933333333333337, 88.966666666666669, 89.0, 89.033333333333331,
        89.066666666666663, 89.1, 89.133333333333326, 89.166666666666671, 89.2,
        89.233333333333334, 89.266666666666666, 89.3, 89.333333333333329,
        89.36666666666666, 89.4, 89.433333333333337, 89.466666666666669, 89.5,
        89.533333333333331, 89.566666666666663, 89.6, 89.633333333333326,
        89.666666666666671, 89.7, 89.733333333333334, 89.766666666666666, 89.8,
        89.833333333333329, 89.86666666666666, 89.9, 89.933333333333337,
        89.966666666666669, 90.0, 90.033333333333331, 90.066666666666663, 90.1,
        90.133333333333326, 90.166666666666671, 90.2, 90.233333333333334,
        90.266666666666666, 90.3, 90.333333333333329, 90.36666666666666, 90.4,
        90.433333333333337, 90.466666666666669, 90.5, 90.533333333333331,
        90.566666666666663, 90.6, 90.633333333333326, 90.666666666666671, 90.7,
        90.733333333333334, 90.766666666666666, 90.8, 90.833333333333329,
        90.86666666666666, 90.9, 90.933333333333337, 90.966666666666669, 91.0,
        91.033333333333331, 91.066666666666663, 91.1, 91.133333333333326,
        91.166666666666671, 91.2, 91.233333333333334, 91.266666666666666, 91.3,
        91.333333333333329, 91.36666666666666, 91.4, 91.433333333333337,
        91.466666666666669, 91.5, 91.533333333333331, 91.566666666666663, 91.6,
        91.633333333333326, 91.666666666666671, 91.7, 91.733333333333334,
        91.766666666666666, 91.8, 91.833333333333329, 91.86666666666666, 91.9,
        91.933333333333337, 91.966666666666669, 92.0, 92.033333333333331,
        92.066666666666663, 92.1, 92.133333333333326, 92.166666666666671, 92.2,
        92.233333333333334, 92.266666666666666, 92.3, 92.333333333333329,
        92.36666666666666, 92.4, 92.433333333333337, 92.466666666666669, 92.5,
        92.533333333333331, 92.566666666666663, 92.6, 92.633333333333326,
        92.666666666666671, 92.7, 92.733333333333334, 92.766666666666666, 92.8,
        92.833333333333329, 92.86666666666666, 92.9, 92.933333333333337,
        92.966666666666669, 93.0, 93.033333333333331, 93.066666666666663, 93.1,
        93.133333333333326, 93.166666666666671, 93.2, 93.233333333333334,
        93.266666666666666, 93.3, 93.333333333333329, 93.36666666666666, 93.4,
        93.433333333333337, 93.466666666666669, 93.5, 93.533333333333331,
        93.566666666666663, 93.6, 93.633333333333326, 93.666666666666671, 93.7,
        93.733333333333334, 93.766666666666666, 93.8, 93.833333333333329,
        93.86666666666666, 93.9, 93.933333333333337, 93.966666666666669, 94.0,
        94.033333333333331, 94.066666666666663, 94.1, 94.133333333333326,
        94.166666666666671, 94.2, 94.233333333333334, 94.266666666666666, 94.3,
        94.333333333333329, 94.36666666666666, 94.4, 94.433333333333337,
        94.466666666666669, 94.5, 94.533333333333331, 94.566666666666663, 94.6,
        94.633333333333326, 94.666666666666671, 94.7, 94.733333333333334,
        94.766666666666666, 94.8, 94.833333333333329, 94.86666666666666, 94.9,
        94.933333333333337, 94.966666666666669, 95.0, 95.033333333333331,
        95.066666666666663, 95.1, 95.133333333333326, 95.166666666666671, 95.2,
        95.233333333333334, 95.266666666666666, 95.3, 95.333333333333329,
        95.36666666666666, 95.4, 95.433333333333337, 95.466666666666669, 95.5,
        95.533333333333331, 95.566666666666663, 95.6, 95.633333333333326,
        95.666666666666671, 95.7, 95.733333333333334, 95.766666666666666, 95.8,
        95.833333333333329, 95.86666666666666, 95.9, 95.933333333333337,
        95.966666666666669, 96.0, 96.033333333333331, 96.066666666666663, 96.1,
        96.133333333333326, 96.166666666666671, 96.2, 96.233333333333334,
        96.266666666666666, 96.3, 96.333333333333329, 96.36666666666666, 96.4,
        96.433333333333337, 96.466666666666669, 96.5, 96.533333333333331,
        96.566666666666663, 96.6, 96.633333333333326, 96.666666666666671, 96.7,
        96.733333333333334, 96.766666666666666, 96.8, 96.833333333333329,
        96.86666666666666, 96.9, 96.933333333333337, 96.966666666666669, 97.0,
        97.033333333333331, 97.066666666666663, 97.1, 97.133333333333326,
        97.166666666666671, 97.2, 97.233333333333334, 97.266666666666666, 97.3,
        97.333333333333329, 97.36666666666666, 97.4, 97.433333333333337,
        97.466666666666669, 97.5, 97.533333333333331, 97.566666666666663, 97.6,
        97.633333333333326, 97.666666666666671, 97.7, 97.733333333333334,
        97.766666666666666, 97.8, 97.833333333333329, 97.86666666666666, 97.9,
        97.933333333333337, 97.966666666666669, 98.0, 98.033333333333331,
        98.066666666666663, 98.1, 98.133333333333326, 98.166666666666671, 98.2,
        98.233333333333334, 98.266666666666666, 98.3, 98.333333333333329,
        98.36666666666666, 98.4, 98.433333333333337, 98.466666666666669, 98.5,
        98.533333333333331, 98.566666666666663, 98.6, 98.633333333333326,
        98.666666666666671, 98.7, 98.733333333333334, 98.766666666666666, 98.8,
        98.833333333333329, 98.86666666666666, 98.9, 98.933333333333337,
        98.966666666666669, 99.0, 99.033333333333331, 99.066666666666663, 99.1,
        99.133333333333326, 99.166666666666671, 99.2, 99.233333333333334,
        99.266666666666666, 99.3, 99.333333333333329, 99.36666666666666, 99.4,
        99.433333333333337, 99.466666666666669, 99.5, 99.533333333333331,
        99.566666666666663, 99.6, 99.633333333333326, 99.666666666666671, 99.7,
        99.733333333333334, 99.766666666666666, 99.8, 99.833333333333329,
        99.86666666666666, 99.9, 99.933333333333337, 99.966666666666669, 100.0,
        100.03333333333333, 100.06666666666666, 100.1, 100.13333333333333,
        100.16666666666667, 100.2, 100.23333333333333, 100.26666666666667, 100.3,
        100.33333333333333, 100.36666666666666, 100.4, 100.43333333333334,
        100.46666666666667, 100.5, 100.53333333333333, 100.56666666666666, 100.6,
        100.63333333333333, 100.66666666666667, 100.7, 100.73333333333333,
        100.76666666666667, 100.8, 100.83333333333333, 100.86666666666666, 100.9,
        100.93333333333334, 100.96666666666667, 101.0, 101.03333333333333,
        101.06666666666666, 101.1, 101.13333333333333, 101.16666666666667, 101.2,
        101.23333333333333, 101.26666666666667, 101.3, 101.33333333333333,
        101.36666666666666, 101.4, 101.43333333333334, 101.46666666666667, 101.5,
        101.53333333333333, 101.56666666666666, 101.6, 101.63333333333333,
        101.66666666666667, 101.7, 101.73333333333333, 101.76666666666667, 101.8,
        101.83333333333333, 101.86666666666666, 101.9, 101.93333333333334,
        101.96666666666667, 102.0, 102.03333333333333, 102.06666666666666, 102.1,
        102.13333333333333, 102.16666666666667, 102.2, 102.23333333333333,
        102.26666666666667, 102.3, 102.33333333333333, 102.36666666666666, 102.4,
        102.43333333333334, 102.46666666666667, 102.5, 102.53333333333333,
        102.56666666666666, 102.6, 102.63333333333333, 102.66666666666667, 102.7,
        102.73333333333333, 102.76666666666667, 102.8, 102.83333333333333,
        102.86666666666666, 102.89999999999999, 102.93333333333334,
        102.96666666666667, 103.0, 103.03333333333333, 103.06666666666666, 103.1,
        103.13333333333333, 103.16666666666667, 103.2, 103.23333333333333,
        103.26666666666667, 103.3, 103.33333333333333, 103.36666666666666,
        103.39999999999999, 103.43333333333334, 103.46666666666667, 103.5,
        103.53333333333333, 103.56666666666666, 103.6, 103.63333333333333,
        103.66666666666667, 103.7, 103.73333333333333, 103.76666666666667, 103.8,
        103.83333333333333, 103.86666666666666, 103.89999999999999,
        103.93333333333334, 103.96666666666667, 104.0, 104.03333333333333,
        104.06666666666666, 104.1, 104.13333333333333, 104.16666666666667, 104.2,
        104.23333333333333, 104.26666666666667, 104.3, 104.33333333333333,
        104.36666666666666, 104.39999999999999, 104.43333333333334,
        104.46666666666667, 104.5, 104.53333333333333, 104.56666666666666, 104.6,
        104.63333333333333, 104.66666666666667, 104.7, 104.73333333333333,
        104.76666666666667, 104.8, 104.83333333333333, 104.86666666666666,
        104.89999999999999, 104.93333333333334, 104.96666666666667, 105.0,
        105.03333333333333, 105.06666666666666, 105.1, 105.13333333333333,
        105.16666666666667, 105.2, 105.23333333333333, 105.26666666666667, 105.3,
        105.33333333333333, 105.36666666666666, 105.39999999999999,
        105.43333333333334, 105.46666666666667, 105.5, 105.53333333333333,
        105.56666666666666, 105.6, 105.63333333333333, 105.66666666666667, 105.7,
        105.73333333333333, 105.76666666666667, 105.8, 105.83333333333333,
        105.86666666666666, 105.89999999999999, 105.93333333333334,
        105.96666666666667, 106.0, 106.03333333333333, 106.06666666666666, 106.1,
        106.13333333333333, 106.16666666666667, 106.2, 106.23333333333333,
        106.26666666666667, 106.3, 106.33333333333333, 106.36666666666666,
        106.39999999999999, 106.43333333333334, 106.46666666666667, 106.5,
        106.53333333333333, 106.56666666666666, 106.6, 106.63333333333333,
        106.66666666666667, 106.7, 106.73333333333333, 106.76666666666667, 106.8,
        106.83333333333333, 106.86666666666666, 106.89999999999999,
        106.93333333333334, 106.96666666666667, 107.0, 107.03333333333333,
        107.06666666666666, 107.1, 107.13333333333333, 107.16666666666667, 107.2,
        107.23333333333333, 107.26666666666667, 107.3, 107.33333333333333,
        107.36666666666666, 107.39999999999999, 107.43333333333334,
        107.46666666666667, 107.5, 107.53333333333333, 107.56666666666666, 107.6,
        107.63333333333333, 107.66666666666667, 107.7, 107.73333333333333,
        107.76666666666667, 107.8, 107.83333333333333, 107.86666666666666,
        107.89999999999999, 107.93333333333334, 107.96666666666667, 108.0,
        108.03333333333333, 108.06666666666666, 108.1, 108.13333333333333,
        108.16666666666667, 108.2, 108.23333333333333, 108.26666666666667, 108.3,
        108.33333333333333, 108.36666666666666, 108.39999999999999,
        108.43333333333334, 108.46666666666667, 108.5, 108.53333333333333,
        108.56666666666666, 108.6, 108.63333333333333, 108.66666666666667, 108.7,
        108.73333333333333, 108.76666666666667, 108.8, 108.83333333333333,
        108.86666666666666, 108.89999999999999, 108.93333333333334,
        108.96666666666667, 109.0, 109.03333333333333, 109.06666666666666, 109.1,
        109.13333333333333, 109.16666666666667, 109.2, 109.23333333333333,
        109.26666666666667, 109.3, 109.33333333333333, 109.36666666666666,
        109.39999999999999, 109.43333333333334, 109.46666666666667, 109.5,
        109.53333333333333, 109.56666666666666, 109.6, 109.63333333333333,
        109.66666666666667, 109.7, 109.73333333333333, 109.76666666666667, 109.8,
        109.83333333333333, 109.86666666666666, 109.89999999999999,
        109.93333333333334, 109.96666666666667, 110.0, 110.03333333333333,
        110.06666666666666, 110.1, 110.13333333333333, 110.16666666666667, 110.2,
        110.23333333333333, 110.26666666666667, 110.3, 110.33333333333333,
        110.36666666666666, 110.39999999999999, 110.43333333333334,
        110.46666666666667, 110.5, 110.53333333333333, 110.56666666666666, 110.6,
        110.63333333333333, 110.66666666666667, 110.7, 110.73333333333333,
        110.76666666666667, 110.8, 110.83333333333333, 110.86666666666666,
        110.89999999999999, 110.93333333333334, 110.96666666666667, 111.0,
        111.03333333333333, 111.06666666666666, 111.1, 111.13333333333333,
        111.16666666666667, 111.2, 111.23333333333333, 111.26666666666667, 111.3,
        111.33333333333333, 111.36666666666666, 111.39999999999999,
        111.43333333333334, 111.46666666666667, 111.5, 111.53333333333333,
        111.56666666666666, 111.6, 111.63333333333333, 111.66666666666667, 111.7,
        111.73333333333333, 111.76666666666667, 111.8, 111.83333333333333,
        111.86666666666666, 111.89999999999999, 111.93333333333334,
        111.96666666666667, 112.0, 112.03333333333333, 112.06666666666666, 112.1,
        112.13333333333333, 112.16666666666667, 112.2, 112.23333333333333,
        112.26666666666667, 112.3, 112.33333333333333, 112.36666666666666,
        112.39999999999999, 112.43333333333334, 112.46666666666667, 112.5,
        112.53333333333333, 112.56666666666666, 112.6, 112.63333333333333,
        112.66666666666667, 112.7, 112.73333333333333, 112.76666666666667, 112.8,
        112.83333333333333, 112.86666666666666, 112.89999999999999,
        112.93333333333334, 112.96666666666667, 113.0, 113.03333333333333,
        113.06666666666666, 113.1, 113.13333333333333, 113.16666666666667, 113.2,
        113.23333333333333, 113.26666666666667, 113.3, 113.33333333333333,
        113.36666666666666, 113.39999999999999, 113.43333333333334,
        113.46666666666667, 113.5, 113.53333333333333, 113.56666666666666, 113.6,
        113.63333333333333, 113.66666666666667, 113.7, 113.73333333333333,
        113.76666666666667, 113.8, 113.83333333333333, 113.86666666666666,
        113.89999999999999, 113.93333333333334, 113.96666666666667, 114.0,
        114.03333333333333, 114.06666666666666, 114.1, 114.13333333333333,
        114.16666666666667, 114.2, 114.23333333333333, 114.26666666666667, 114.3,
        114.33333333333333, 114.36666666666666, 114.39999999999999,
        114.43333333333334, 114.46666666666667, 114.5, 114.53333333333333,
        114.56666666666666, 114.6, 114.63333333333333, 114.66666666666667, 114.7,
        114.73333333333333, 114.76666666666667, 114.8, 114.83333333333333,
        114.86666666666666, 114.89999999999999, 114.93333333333334,
        114.96666666666667, 115.0, 115.03333333333333, 115.06666666666666, 115.1,
        115.13333333333333, 115.16666666666667, 115.2, 115.23333333333333,
        115.26666666666667, 115.3, 115.33333333333333, 115.36666666666666,
        115.39999999999999, 115.43333333333334, 115.46666666666667, 115.5,
        115.53333333333333, 115.56666666666666, 115.6, 115.63333333333333,
        115.66666666666667, 115.7, 115.73333333333333, 115.76666666666667, 115.8,
        115.83333333333333, 115.86666666666666, 115.89999999999999,
        115.93333333333334, 115.96666666666667, 116.0, 116.03333333333333,
        116.06666666666666, 116.1, 116.13333333333333, 116.16666666666667, 116.2,
        116.23333333333333, 116.26666666666667, 116.3, 116.33333333333333,
        116.36666666666666, 116.39999999999999, 116.43333333333334,
        116.46666666666667, 116.5, 116.53333333333333, 116.56666666666666, 116.6,
        116.63333333333333, 116.66666666666667, 116.7, 116.73333333333333,
        116.76666666666667, 116.8, 116.83333333333333, 116.86666666666666,
        116.89999999999999, 116.93333333333334, 116.96666666666667, 117.0,
        117.03333333333333, 117.06666666666666, 117.1, 117.13333333333333,
        117.16666666666667, 117.2, 117.23333333333333, 117.26666666666667, 117.3,
        117.33333333333333, 117.36666666666666, 117.39999999999999,
        117.43333333333334, 117.46666666666667, 117.5, 117.53333333333333,
        117.56666666666666, 117.6, 117.63333333333333, 117.66666666666667, 117.7,
        117.73333333333333, 117.76666666666667, 117.8, 117.83333333333333,
        117.86666666666666, 117.89999999999999, 117.93333333333334,
        117.96666666666667, 118.0, 118.03333333333333, 118.06666666666666, 118.1,
        118.13333333333333, 118.16666666666667, 118.2, 118.23333333333333,
        118.26666666666667, 118.3, 118.33333333333333, 118.36666666666666,
        118.39999999999999, 118.43333333333334, 118.46666666666667, 118.5,
        118.53333333333333, 118.56666666666666, 118.6, 118.63333333333333,
        118.66666666666667, 118.7, 118.73333333333333, 118.76666666666667, 118.8,
        118.83333333333333, 118.86666666666666, 118.89999999999999,
        118.93333333333334, 118.96666666666667, 119.0, 119.03333333333333,
        119.06666666666666, 119.1, 119.13333333333333, 119.16666666666667, 119.2,
        119.23333333333333, 119.26666666666667, 119.3, 119.33333333333333,
        119.36666666666666, 119.39999999999999, 119.43333333333334,
        119.46666666666667, 119.5, 119.53333333333333, 119.56666666666666, 119.6,
        119.63333333333333, 119.66666666666667, 119.7, 119.73333333333333,
        119.76666666666667, 119.8, 119.83333333333333, 119.86666666666666,
        119.89999999999999, 119.93333333333334, 119.96666666666667, 120.0,
        120.03333333333333, 120.06666666666666, 120.1, 120.13333333333333,
        120.16666666666667, 120.2, 120.23333333333333, 120.26666666666667, 120.3,
        120.33333333333333, 120.36666666666666, 120.39999999999999,
        120.43333333333334, 120.46666666666667, 120.5, 120.53333333333333,
        120.56666666666666, 120.6, 120.63333333333333, 120.66666666666667, 120.7,
        120.73333333333333, 120.76666666666667, 120.8, 120.83333333333333,
        120.86666666666666, 120.89999999999999, 120.93333333333334,
        120.96666666666667, 121.0, 121.03333333333333, 121.06666666666666, 121.1,
        121.13333333333333, 121.16666666666667, 121.2, 121.23333333333333,
        121.26666666666667, 121.3, 121.33333333333333, 121.36666666666666,
        121.39999999999999, 121.43333333333334, 121.46666666666667, 121.5,
        121.53333333333333, 121.56666666666666, 121.6, 121.63333333333333,
        121.66666666666667, 121.7, 121.73333333333333, 121.76666666666667, 121.8,
        121.83333333333333, 121.86666666666666, 121.89999999999999,
        121.93333333333334, 121.96666666666667, 122.0, 122.03333333333333,
        122.06666666666666, 122.1, 122.13333333333333, 122.16666666666667, 122.2,
        122.23333333333333, 122.26666666666667, 122.3, 122.33333333333333,
        122.36666666666666, 122.39999999999999, 122.43333333333334,
        122.46666666666667, 122.5, 122.53333333333333, 122.56666666666666, 122.6,
        122.63333333333333, 122.66666666666667, 122.7, 122.73333333333333,
        122.76666666666667, 122.8, 122.83333333333333, 122.86666666666666,
        122.89999999999999, 122.93333333333334, 122.96666666666667, 123.0,
        123.03333333333333, 123.06666666666666, 123.1, 123.13333333333333,
        123.16666666666667, 123.2, 123.23333333333333, 123.26666666666667, 123.3,
        123.33333333333333, 123.36666666666666, 123.39999999999999,
        123.43333333333334, 123.46666666666667, 123.5, 123.53333333333333,
        123.56666666666666, 123.6, 123.63333333333333, 123.66666666666667, 123.7,
        123.73333333333333, 123.76666666666667, 123.8, 123.83333333333333,
        123.86666666666666, 123.89999999999999, 123.93333333333334,
        123.96666666666667, 124.0, 124.03333333333333, 124.06666666666666, 124.1,
        124.13333333333333, 124.16666666666667, 124.2, 124.23333333333333,
        124.26666666666667, 124.3, 124.33333333333333, 124.36666666666666,
        124.39999999999999, 124.43333333333334, 124.46666666666667, 124.5,
        124.53333333333333, 124.56666666666666, 124.6, 124.63333333333333,
        124.66666666666667, 124.7, 124.73333333333333, 124.76666666666667, 124.8,
        124.83333333333333, 124.86666666666666, 124.89999999999999,
        124.93333333333334, 124.96666666666667, 125.0, 125.03333333333333,
        125.06666666666666, 125.1, 125.13333333333333, 125.16666666666667, 125.2,
        125.23333333333333, 125.26666666666667, 125.3, 125.33333333333333,
        125.36666666666666, 125.39999999999999, 125.43333333333334,
        125.46666666666667, 125.5, 125.53333333333333, 125.56666666666666, 125.6,
        125.63333333333333, 125.66666666666667, 125.7, 125.73333333333333,
        125.76666666666667, 125.8, 125.83333333333333, 125.86666666666666,
        125.89999999999999, 125.93333333333334, 125.96666666666667, 126.0,
        126.03333333333333, 126.06666666666666, 126.1, 126.13333333333333,
        126.16666666666667, 126.2, 126.23333333333333, 126.26666666666667, 126.3,
        126.33333333333333, 126.36666666666666, 126.39999999999999,
        126.43333333333334, 126.46666666666667, 126.5, 126.53333333333333,
        126.56666666666666, 126.6, 126.63333333333333, 126.66666666666667, 126.7,
        126.73333333333333, 126.76666666666667, 126.8, 126.83333333333333,
        126.86666666666666, 126.89999999999999, 126.93333333333334,
        126.96666666666667, 127.0, 127.03333333333333, 127.06666666666666, 127.1,
        127.13333333333333, 127.16666666666667, 127.2, 127.23333333333333,
        127.26666666666667, 127.3, 127.33333333333333, 127.36666666666666,
        127.39999999999999, 127.43333333333334, 127.46666666666667, 127.5,
        127.53333333333333, 127.56666666666666, 127.6, 127.63333333333333,
        127.66666666666667, 127.7, 127.73333333333333, 127.76666666666667, 127.8,
        127.83333333333333, 127.86666666666666, 127.89999999999999,
        127.93333333333334, 127.96666666666667, 128.0, 128.03333333333333,
        128.06666666666666, 128.1, 128.13333333333333, 128.16666666666666, 128.2,
        128.23333333333332, 128.26666666666665, 128.3, 128.33333333333334,
        128.36666666666667, 128.4, 128.43333333333334, 128.46666666666667, 128.5,
        128.53333333333333, 128.56666666666666, 128.6, 128.63333333333333,
        128.66666666666666, 128.7, 128.73333333333332, 128.76666666666665, 128.8,
        128.83333333333334, 128.86666666666667, 128.9, 128.93333333333334,
        128.96666666666667, 129.0, 129.03333333333333, 129.06666666666666, 129.1,
        129.13333333333333, 129.16666666666666, 129.2, 129.23333333333332,
        129.26666666666665, 129.3, 129.33333333333334, 129.36666666666667, 129.4,
        129.43333333333334, 129.46666666666667, 129.5, 129.53333333333333,
        129.56666666666666, 129.6, 129.63333333333333, 129.66666666666666, 129.7,
        129.73333333333332, 129.76666666666665, 129.8, 129.83333333333334,
        129.86666666666667, 129.9, 129.93333333333334, 129.96666666666667, 130.0,
        130.03333333333333, 130.06666666666666, 130.1, 130.13333333333333,
        130.16666666666666, 130.2, 130.23333333333332, 130.26666666666665, 130.3,
        130.33333333333334, 130.36666666666667, 130.4, 130.43333333333334,
        130.46666666666667, 130.5, 130.53333333333333, 130.56666666666666, 130.6,
        130.63333333333333, 130.66666666666666, 130.7, 130.73333333333332,
        130.76666666666665, 130.8, 130.83333333333334, 130.86666666666667, 130.9,
        130.93333333333334, 130.96666666666667, 131.0, 131.03333333333333,
        131.06666666666666, 131.1, 131.13333333333333, 131.16666666666666, 131.2,
        131.23333333333332, 131.26666666666665, 131.3, 131.33333333333334,
        131.36666666666667, 131.4, 131.43333333333334, 131.46666666666667, 131.5,
        131.53333333333333, 131.56666666666666, 131.6, 131.63333333333333,
        131.66666666666666, 131.7, 131.73333333333332, 131.76666666666665, 131.8,
        131.83333333333334, 131.86666666666667, 131.9, 131.93333333333334,
        131.96666666666667, 132.0, 132.03333333333333, 132.06666666666666, 132.1,
        132.13333333333333, 132.16666666666666, 132.2, 132.23333333333332,
        132.26666666666665, 132.3, 132.33333333333334, 132.36666666666667, 132.4,
        132.43333333333334, 132.46666666666667, 132.5, 132.53333333333333,
        132.56666666666666, 132.6, 132.63333333333333, 132.66666666666666, 132.7,
        132.73333333333332, 132.76666666666665, 132.8, 132.83333333333334,
        132.86666666666667, 132.9, 132.93333333333334, 132.96666666666667, 133.0,
        133.03333333333333, 133.06666666666666, 133.1, 133.13333333333333,
        133.16666666666666, 133.2, 133.23333333333332, 133.26666666666665, 133.3,
        133.33333333333334, 133.36666666666667, 133.4, 133.43333333333334,
        133.46666666666667, 133.5, 133.53333333333333, 133.56666666666666, 133.6,
        133.63333333333333, 133.66666666666666, 133.7, 133.73333333333332,
        133.76666666666665, 133.8, 133.83333333333334, 133.86666666666667, 133.9,
        133.93333333333334, 133.96666666666667, 134.0, 134.03333333333333,
        134.06666666666666, 134.1, 134.13333333333333, 134.16666666666666, 134.2,
        134.23333333333332, 134.26666666666665, 134.3, 134.33333333333334,
        134.36666666666667, 134.4, 134.43333333333334, 134.46666666666667, 134.5,
        134.53333333333333, 134.56666666666666, 134.6, 134.63333333333333,
        134.66666666666666, 134.7, 134.73333333333332, 134.76666666666665, 134.8,
        134.83333333333334, 134.86666666666667, 134.9, 134.93333333333334,
        134.96666666666667, 135.0, 135.03333333333333, 135.06666666666666, 135.1,
        135.13333333333333, 135.16666666666666, 135.2, 135.23333333333332,
        135.26666666666665, 135.3, 135.33333333333334, 135.36666666666667, 135.4,
        135.43333333333334, 135.46666666666667, 135.5, 135.53333333333333,
        135.56666666666666, 135.6, 135.63333333333333, 135.66666666666666, 135.7,
        135.73333333333332, 135.76666666666665, 135.8, 135.83333333333334,
        135.86666666666667, 135.9, 135.93333333333334, 135.96666666666667, 136.0,
        136.03333333333333, 136.06666666666666, 136.1, 136.13333333333333,
        136.16666666666666, 136.2, 136.23333333333332, 136.26666666666665, 136.3,
        136.33333333333334, 136.36666666666667, 136.4, 136.43333333333334,
        136.46666666666667, 136.5, 136.53333333333333, 136.56666666666666, 136.6,
        136.63333333333333, 136.66666666666666, 136.7, 136.73333333333332,
        136.76666666666665, 136.8, 136.83333333333334, 136.86666666666667, 136.9,
        136.93333333333334, 136.96666666666667, 137.0, 137.03333333333333,
        137.06666666666666, 137.1, 137.13333333333333, 137.16666666666666, 137.2,
        137.23333333333332, 137.26666666666665, 137.3, 137.33333333333334,
        137.36666666666667, 137.4, 137.43333333333334, 137.46666666666667, 137.5,
        137.53333333333333, 137.56666666666666, 137.6, 137.63333333333333,
        137.66666666666666, 137.7, 137.73333333333332, 137.76666666666665, 137.8,
        137.83333333333334, 137.86666666666667, 137.9, 137.93333333333334,
        137.96666666666667, 138.0, 138.03333333333333, 138.06666666666666, 138.1,
        138.13333333333333, 138.16666666666666, 138.2, 138.23333333333332,
        138.26666666666665, 138.3, 138.33333333333334, 138.36666666666667, 138.4,
        138.43333333333334, 138.46666666666667, 138.5, 138.53333333333333,
        138.56666666666666, 138.6, 138.63333333333333, 138.66666666666666, 138.7,
        138.73333333333332, 138.76666666666665, 138.8, 138.83333333333334,
        138.86666666666667, 138.9, 138.93333333333334, 138.96666666666667, 139.0,
        139.03333333333333, 139.06666666666666, 139.1, 139.13333333333333,
        139.16666666666666, 139.2, 139.23333333333332, 139.26666666666665, 139.3,
        139.33333333333334, 139.36666666666667, 139.4, 139.43333333333334,
        139.46666666666667, 139.5, 139.53333333333333, 139.56666666666666, 139.6,
        139.63333333333333, 139.66666666666666, 139.7, 139.73333333333332,
        139.76666666666665, 139.8, 139.83333333333334, 139.86666666666667, 139.9,
        139.93333333333334, 139.96666666666667, 140.0, 140.03333333333333,
        140.06666666666666, 140.1, 140.13333333333333, 140.16666666666666, 140.2,
        140.23333333333332, 140.26666666666665, 140.3, 140.33333333333334,
        140.36666666666667, 140.4, 140.43333333333334, 140.46666666666667, 140.5,
        140.53333333333333, 140.56666666666666, 140.6, 140.63333333333333,
        140.66666666666666, 140.7, 140.73333333333332, 140.76666666666665, 140.8,
        140.83333333333334, 140.86666666666667, 140.9, 140.93333333333334,
        140.96666666666667, 141.0, 141.03333333333333, 141.06666666666666, 141.1,
        141.13333333333333, 141.16666666666666, 141.2, 141.23333333333332,
        141.26666666666665, 141.3, 141.33333333333334, 141.36666666666667, 141.4,
        141.43333333333334, 141.46666666666667, 141.5, 141.53333333333333,
        141.56666666666666, 141.6, 141.63333333333333, 141.66666666666666, 141.7,
        141.73333333333332, 141.76666666666665, 141.8, 141.83333333333334,
        141.86666666666667, 141.9, 141.93333333333334, 141.96666666666667, 142.0,
        142.03333333333333, 142.06666666666666, 142.1, 142.13333333333333,
        142.16666666666666, 142.2, 142.23333333333332, 142.26666666666665, 142.3,
        142.33333333333334, 142.36666666666667, 142.4, 142.43333333333334,
        142.46666666666667, 142.5, 142.53333333333333, 142.56666666666666, 142.6,
        142.63333333333333, 142.66666666666666, 142.7, 142.73333333333332,
        142.76666666666665, 142.8, 142.83333333333334, 142.86666666666667, 142.9,
        142.93333333333334, 142.96666666666667, 143.0, 143.03333333333333,
        143.06666666666666, 143.1, 143.13333333333333, 143.16666666666666, 143.2,
        143.23333333333332, 143.26666666666665, 143.3, 143.33333333333334,
        143.36666666666667, 143.4, 143.43333333333334, 143.46666666666667, 143.5,
        143.53333333333333, 143.56666666666666, 143.6, 143.63333333333333,
        143.66666666666666, 143.7, 143.73333333333332, 143.76666666666665, 143.8,
        143.83333333333334, 143.86666666666667, 143.9, 143.93333333333334,
        143.96666666666667, 144.0, 144.03333333333333, 144.06666666666666, 144.1,
        144.13333333333333, 144.16666666666666, 144.2, 144.23333333333332,
        144.26666666666665, 144.3, 144.33333333333334, 144.36666666666667, 144.4,
        144.43333333333334, 144.46666666666667, 144.5, 144.53333333333333,
        144.56666666666666, 144.6, 144.63333333333333, 144.66666666666666, 144.7,
        144.73333333333332, 144.76666666666665, 144.8, 144.83333333333334,
        144.86666666666667, 144.9, 144.93333333333334, 144.96666666666667, 145.0,
        145.03333333333333, 145.06666666666666, 145.1, 145.13333333333333,
        145.16666666666666, 145.2, 145.23333333333332, 145.26666666666665, 145.3,
        145.33333333333334, 145.36666666666667, 145.4, 145.43333333333334,
        145.46666666666667, 145.5, 145.53333333333333, 145.56666666666666, 145.6,
        145.63333333333333, 145.66666666666666, 145.7, 145.73333333333332,
        145.76666666666665, 145.8, 145.83333333333334, 145.86666666666667, 145.9,
        145.93333333333334, 145.96666666666667, 146.0, 146.03333333333333,
        146.06666666666666, 146.1, 146.13333333333333, 146.16666666666666, 146.2,
        146.23333333333332, 146.26666666666665, 146.3, 146.33333333333334,
        146.36666666666667, 146.4, 146.43333333333334, 146.46666666666667, 146.5,
        146.53333333333333, 146.56666666666666, 146.6, 146.63333333333333,
        146.66666666666666, 146.7, 146.73333333333332, 146.76666666666665, 146.8,
        146.83333333333334, 146.86666666666667, 146.9, 146.93333333333334,
        146.96666666666667, 147.0, 147.03333333333333, 147.06666666666666, 147.1,
        147.13333333333333, 147.16666666666666, 147.2, 147.23333333333332,
        147.26666666666665, 147.3, 147.33333333333334, 147.36666666666667, 147.4,
        147.43333333333334, 147.46666666666667, 147.5, 147.53333333333333,
        147.56666666666666, 147.6, 147.63333333333333, 147.66666666666666, 147.7,
        147.73333333333332, 147.76666666666665, 147.8, 147.83333333333334,
        147.86666666666667, 147.9, 147.93333333333334, 147.96666666666667, 148.0,
        148.03333333333333, 148.06666666666666, 148.1, 148.13333333333333,
        148.16666666666666, 148.2, 148.23333333333332, 148.26666666666665, 148.3,
        148.33333333333334, 148.36666666666667, 148.4, 148.43333333333334,
        148.46666666666667, 148.5, 148.53333333333333, 148.56666666666666, 148.6,
        148.63333333333333, 148.66666666666666, 148.7, 148.73333333333332,
        148.76666666666665, 148.8, 148.83333333333334, 148.86666666666667, 148.9,
        148.93333333333334, 148.96666666666667, 149.0, 149.03333333333333,
        149.06666666666666, 149.1, 149.13333333333333, 149.16666666666666, 149.2,
        149.23333333333332, 149.26666666666665, 149.3, 149.33333333333334,
        149.36666666666667, 149.4, 149.43333333333334, 149.46666666666667, 149.5,
        149.53333333333333, 149.56666666666666, 149.6, 149.63333333333333,
        149.66666666666666, 149.7, 149.73333333333332, 149.76666666666665, 149.8,
        149.83333333333334, 149.86666666666667, 149.9, 149.93333333333334,
        149.96666666666667, 150.0, 150.03333333333333, 150.06666666666666, 150.1,
        150.13333333333333, 150.16666666666666, 150.2, 150.23333333333332,
        150.26666666666665, 150.3, 150.33333333333334, 150.36666666666667, 150.4,
        150.43333333333334, 150.46666666666667, 150.5, 150.53333333333333,
        150.56666666666666, 150.6, 150.63333333333333, 150.66666666666666, 150.7,
        150.73333333333332, 150.76666666666665, 150.8, 150.83333333333334,
        150.86666666666667, 150.9, 150.93333333333334, 150.96666666666667, 151.0,
        151.03333333333333, 151.06666666666666, 151.1, 151.13333333333333,
        151.16666666666666, 151.2, 151.23333333333332, 151.26666666666665, 151.3,
        151.33333333333334, 151.36666666666667, 151.4, 151.43333333333334,
        151.46666666666667, 151.5, 151.53333333333333, 151.56666666666666, 151.6,
        151.63333333333333, 151.66666666666666, 151.7, 151.73333333333332,
        151.76666666666665, 151.8, 151.83333333333334, 151.86666666666667, 151.9,
        151.93333333333334, 151.96666666666667, 152.0, 152.03333333333333,
        152.06666666666666, 152.1, 152.13333333333333, 152.16666666666666, 152.2,
        152.23333333333332, 152.26666666666665, 152.3, 152.33333333333334,
        152.36666666666667, 152.4, 152.43333333333334, 152.46666666666667, 152.5,
        152.53333333333333, 152.56666666666666, 152.6, 152.63333333333333,
        152.66666666666666, 152.7, 152.73333333333332, 152.76666666666665, 152.8,
        152.83333333333334, 152.86666666666667, 152.9, 152.93333333333334,
        152.96666666666667, 153.0, 153.03333333333333, 153.06666666666666, 153.1,
        153.13333333333333, 153.16666666666666, 153.2, 153.23333333333332,
        153.26666666666665, 153.3, 153.33333333333334, 153.36666666666667, 153.4,
        153.43333333333334, 153.46666666666667, 153.5, 153.53333333333333,
        153.56666666666666, 153.6, 153.63333333333333, 153.66666666666666, 153.7,
        153.73333333333332, 153.76666666666665, 153.8, 153.83333333333334,
        153.86666666666667, 153.9, 153.93333333333334, 153.96666666666667, 154.0,
        154.03333333333333, 154.06666666666666, 154.1, 154.13333333333333,
        154.16666666666666, 154.2, 154.23333333333332, 154.26666666666665, 154.3,
        154.33333333333334, 154.36666666666667, 154.4, 154.43333333333334,
        154.46666666666667, 154.5, 154.53333333333333, 154.56666666666666, 154.6,
        154.63333333333333, 154.66666666666666, 154.7, 154.73333333333332,
        154.76666666666665, 154.8, 154.83333333333334, 154.86666666666667, 154.9,
        154.93333333333334, 154.96666666666667, 155.0, 155.03333333333333,
        155.06666666666666, 155.1, 155.13333333333333, 155.16666666666666, 155.2,
        155.23333333333332, 155.26666666666665, 155.3, 155.33333333333334,
        155.36666666666667, 155.4, 155.43333333333334, 155.46666666666667, 155.5,
        155.53333333333333, 155.56666666666666, 155.6, 155.63333333333333,
        155.66666666666666, 155.7, 155.73333333333332, 155.76666666666665, 155.8,
        155.83333333333334, 155.86666666666667, 155.9, 155.93333333333334,
        155.96666666666667, 156.0, 156.03333333333333, 156.06666666666666, 156.1,
        156.13333333333333, 156.16666666666666, 156.2, 156.23333333333332,
        156.26666666666665, 156.3, 156.33333333333334, 156.36666666666667, 156.4,
        156.43333333333334, 156.46666666666667, 156.5, 156.53333333333333,
        156.56666666666666, 156.6, 156.63333333333333, 156.66666666666666, 156.7,
        156.73333333333332, 156.76666666666665, 156.8, 156.83333333333334,
        156.86666666666667, 156.9, 156.93333333333334, 156.96666666666667, 157.0,
        157.03333333333333, 157.06666666666666, 157.1, 157.13333333333333,
        157.16666666666666, 157.2, 157.23333333333332, 157.26666666666665, 157.3,
        157.33333333333334, 157.36666666666667, 157.4, 157.43333333333334,
        157.46666666666667, 157.5, 157.53333333333333, 157.56666666666666, 157.6,
        157.63333333333333, 157.66666666666666, 157.7, 157.73333333333332,
        157.76666666666665, 157.8, 157.83333333333334, 157.86666666666667, 157.9,
        157.93333333333334, 157.96666666666667, 158.0, 158.03333333333333,
        158.06666666666666, 158.1, 158.13333333333333, 158.16666666666666, 158.2,
        158.23333333333332, 158.26666666666665, 158.3, 158.33333333333334,
        158.36666666666667, 158.4, 158.43333333333334, 158.46666666666667, 158.5,
        158.53333333333333, 158.56666666666666, 158.6, 158.63333333333333,
        158.66666666666666, 158.7, 158.73333333333332, 158.76666666666665, 158.8,
        158.83333333333334, 158.86666666666667, 158.9, 158.93333333333334,
        158.96666666666667, 159.0, 159.03333333333333, 159.06666666666666, 159.1,
        159.13333333333333, 159.16666666666666, 159.2, 159.23333333333332,
        159.26666666666665, 159.3, 159.33333333333334, 159.36666666666667, 159.4,
        159.43333333333334, 159.46666666666667, 159.5, 159.53333333333333,
        159.56666666666666, 159.6, 159.63333333333333, 159.66666666666666, 159.7,
        159.73333333333332, 159.76666666666665, 159.8, 159.83333333333334,
        159.86666666666667, 159.9, 159.93333333333334, 159.96666666666667, 160.0,
        160.03333333333333, 160.06666666666666, 160.1, 160.13333333333333,
        160.16666666666666, 160.2, 160.23333333333332, 160.26666666666665, 160.3,
        160.33333333333334, 160.36666666666667, 160.4, 160.43333333333334,
        160.46666666666667, 160.5, 160.53333333333333, 160.56666666666666, 160.6,
        160.63333333333333, 160.66666666666666, 160.7, 160.73333333333332,
        160.76666666666665, 160.8, 160.83333333333334, 160.86666666666667, 160.9,
        160.93333333333334, 160.96666666666667, 161.0, 161.03333333333333,
        161.06666666666666, 161.1, 161.13333333333333, 161.16666666666666, 161.2,
        161.23333333333332, 161.26666666666665, 161.3, 161.33333333333334,
        161.36666666666667, 161.4, 161.43333333333334, 161.46666666666667, 161.5,
        161.53333333333333, 161.56666666666666, 161.6, 161.63333333333333,
        161.66666666666666, 161.7, 161.73333333333332, 161.76666666666665, 161.8,
        161.83333333333334, 161.86666666666667, 161.9, 161.93333333333334,
        161.96666666666667, 162.0, 162.03333333333333, 162.06666666666666, 162.1,
        162.13333333333333, 162.16666666666666, 162.2, 162.23333333333332,
        162.26666666666665, 162.3, 162.33333333333334, 162.36666666666667, 162.4,
        162.43333333333334, 162.46666666666667, 162.5, 162.53333333333333,
        162.56666666666666, 162.6, 162.63333333333333, 162.66666666666666, 162.7,
        162.73333333333332, 162.76666666666665, 162.8, 162.83333333333334,
        162.86666666666667, 162.9, 162.93333333333334, 162.96666666666667, 163.0,
        163.03333333333333, 163.06666666666666, 163.1, 163.13333333333333,
        163.16666666666666, 163.2, 163.23333333333332, 163.26666666666665, 163.3,
        163.33333333333334, 163.36666666666667, 163.4, 163.43333333333334,
        163.46666666666667, 163.5, 163.53333333333333, 163.56666666666666, 163.6,
        163.63333333333333, 163.66666666666666, 163.7, 163.73333333333332,
        163.76666666666665, 163.8, 163.83333333333334, 163.86666666666667, 163.9,
        163.93333333333334, 163.96666666666667, 164.0, 164.03333333333333,
        164.06666666666666, 164.1, 164.13333333333333, 164.16666666666666, 164.2,
        164.23333333333332, 164.26666666666665, 164.3, 164.33333333333334,
        164.36666666666667, 164.4, 164.43333333333334, 164.46666666666667, 164.5,
        164.53333333333333, 164.56666666666666, 164.6, 164.63333333333333,
        164.66666666666666, 164.7, 164.73333333333332, 164.76666666666665, 164.8,
        164.83333333333334, 164.86666666666667, 164.9, 164.93333333333334,
        164.96666666666667, 165.0, 165.03333333333333, 165.06666666666666, 165.1,
        165.13333333333333, 165.16666666666666, 165.2, 165.23333333333332,
        165.26666666666665, 165.3, 165.33333333333334, 165.36666666666667, 165.4,
        165.43333333333334, 165.46666666666667, 165.5, 165.53333333333333,
        165.56666666666666, 165.6, 165.63333333333333, 165.66666666666666, 165.7,
        165.73333333333332, 165.76666666666665, 165.8, 165.83333333333334,
        165.86666666666667, 165.9, 165.93333333333334, 165.96666666666667, 166.0,
        166.03333333333333, 166.06666666666666, 166.1, 166.13333333333333,
        166.16666666666666, 166.2, 166.23333333333332, 166.26666666666665, 166.3,
        166.33333333333334, 166.36666666666667, 166.4, 166.43333333333334,
        166.46666666666667, 166.5, 166.53333333333333, 166.56666666666666, 166.6,
        166.63333333333333, 166.66666666666666, 166.7, 166.73333333333332,
        166.76666666666665, 166.8, 166.83333333333334, 166.86666666666667, 166.9,
        166.93333333333334, 166.96666666666667, 167.0, 167.03333333333333,
        167.06666666666666, 167.1, 167.13333333333333, 167.16666666666666, 167.2,
        167.23333333333332, 167.26666666666665, 167.3, 167.33333333333334,
        167.36666666666667, 167.4, 167.43333333333334, 167.46666666666667, 167.5,
        167.53333333333333, 167.56666666666666, 167.6, 167.63333333333333,
        167.66666666666666, 167.7, 167.73333333333332, 167.76666666666665, 167.8,
        167.83333333333334, 167.86666666666667, 167.9, 167.93333333333334,
        167.96666666666667, 168.0, 168.03333333333333, 168.06666666666666, 168.1,
        168.13333333333333, 168.16666666666666, 168.2, 168.23333333333332,
        168.26666666666665, 168.3, 168.33333333333334, 168.36666666666667, 168.4,
        168.43333333333334, 168.46666666666667, 168.5, 168.53333333333333,
        168.56666666666666, 168.6, 168.63333333333333, 168.66666666666666, 168.7,
        168.73333333333332, 168.76666666666665, 168.8, 168.83333333333334,
        168.86666666666667, 168.9, 168.93333333333334, 168.96666666666667, 169.0,
        169.03333333333333, 169.06666666666666, 169.1, 169.13333333333333,
        169.16666666666666, 169.2, 169.23333333333332, 169.26666666666665, 169.3,
        169.33333333333334, 169.36666666666667, 169.4, 169.43333333333334,
        169.46666666666667, 169.5, 169.53333333333333, 169.56666666666666, 169.6,
        169.63333333333333, 169.66666666666666, 169.7, 169.73333333333332,
        169.76666666666665, 169.8, 169.83333333333334, 169.86666666666667, 169.9,
        169.93333333333334, 169.96666666666667, 170.0, 170.03333333333333,
        170.06666666666666, 170.1, 170.13333333333333, 170.16666666666666, 170.2,
        170.23333333333332, 170.26666666666665, 170.3, 170.33333333333334,
        170.36666666666667, 170.4, 170.43333333333334, 170.46666666666667, 170.5,
        170.53333333333333, 170.56666666666666, 170.6, 170.63333333333333,
        170.66666666666666, 170.7, 170.73333333333332, 170.76666666666665, 170.8,
        170.83333333333334, 170.86666666666667, 170.9, 170.93333333333334,
        170.96666666666667, 171.0, 171.03333333333333, 171.06666666666666, 171.1,
        171.13333333333333, 171.16666666666666, 171.2, 171.23333333333332,
        171.26666666666665, 171.3, 171.33333333333334, 171.36666666666667, 171.4,
        171.43333333333334, 171.46666666666667, 171.5, 171.53333333333333,
        171.56666666666666, 171.6, 171.63333333333333, 171.66666666666666, 171.7,
        171.73333333333332, 171.76666666666665, 171.8, 171.83333333333334,
        171.86666666666667, 171.9, 171.93333333333334, 171.96666666666667, 172.0,
        172.03333333333333, 172.06666666666666, 172.1, 172.13333333333333,
        172.16666666666666, 172.2, 172.23333333333332, 172.26666666666665, 172.3,
        172.33333333333334, 172.36666666666667, 172.4, 172.43333333333334,
        172.46666666666667, 172.5, 172.53333333333333, 172.56666666666666, 172.6,
        172.63333333333333, 172.66666666666666, 172.7, 172.73333333333332,
        172.76666666666665, 172.8, 172.83333333333334, 172.86666666666667, 172.9,
        172.93333333333334, 172.96666666666667, 173.0, 173.03333333333333,
        173.06666666666666, 173.1, 173.13333333333333, 173.16666666666666, 173.2,
        173.23333333333332, 173.26666666666665, 173.3, 173.33333333333334,
        173.36666666666667, 173.4, 173.43333333333334, 173.46666666666667, 173.5,
        173.53333333333333, 173.56666666666666, 173.6, 173.63333333333333,
        173.66666666666666, 173.7, 173.73333333333332, 173.76666666666665, 173.8,
        173.83333333333334, 173.86666666666667, 173.9, 173.93333333333334,
        173.96666666666667, 174.0, 174.03333333333333, 174.06666666666666, 174.1,
        174.13333333333333, 174.16666666666666, 174.2, 174.23333333333332,
        174.26666666666665, 174.3, 174.33333333333334, 174.36666666666667, 174.4,
        174.43333333333334, 174.46666666666667, 174.5, 174.53333333333333,
        174.56666666666666, 174.6, 174.63333333333333, 174.66666666666666, 174.7,
        174.73333333333332, 174.76666666666665, 174.8, 174.83333333333334,
        174.86666666666667, 174.9, 174.93333333333334, 174.96666666666667, 175.0,
        175.03333333333333, 175.06666666666666, 175.1, 175.13333333333333,
        175.16666666666666, 175.2, 175.23333333333332, 175.26666666666665, 175.3,
        175.33333333333334, 175.36666666666667, 175.4, 175.43333333333334,
        175.46666666666667, 175.5, 175.53333333333333, 175.56666666666666, 175.6,
        175.63333333333333, 175.66666666666666, 175.7, 175.73333333333332,
        175.76666666666665, 175.8, 175.83333333333334, 175.86666666666667, 175.9,
        175.93333333333334, 175.96666666666667, 176.0, 176.03333333333333,
        176.06666666666666, 176.1, 176.13333333333333, 176.16666666666666, 176.2,
        176.23333333333332, 176.26666666666665, 176.3, 176.33333333333334,
        176.36666666666667, 176.4, 176.43333333333334, 176.46666666666667, 176.5,
        176.53333333333333, 176.56666666666666, 176.6, 176.63333333333333,
        176.66666666666666, 176.7, 176.73333333333332, 176.76666666666665, 176.8,
        176.83333333333334, 176.86666666666667, 176.9, 176.93333333333334,
        176.96666666666667, 177.0, 177.03333333333333, 177.06666666666666, 177.1,
        177.13333333333333, 177.16666666666666, 177.2, 177.23333333333332,
        177.26666666666665, 177.3, 177.33333333333334, 177.36666666666667, 177.4,
        177.43333333333334, 177.46666666666667, 177.5, 177.53333333333333,
        177.56666666666666, 177.6, 177.63333333333333, 177.66666666666666, 177.7,
        177.73333333333332, 177.76666666666665, 177.8, 177.83333333333334,
        177.86666666666667, 177.9, 177.93333333333334, 177.96666666666667, 178.0,
        178.03333333333333, 178.06666666666666, 178.1, 178.13333333333333,
        178.16666666666666, 178.2, 178.23333333333332, 178.26666666666665, 178.3,
        178.33333333333334, 178.36666666666667, 178.4, 178.43333333333334,
        178.46666666666667, 178.5, 178.53333333333333, 178.56666666666666, 178.6,
        178.63333333333333, 178.66666666666666, 178.7, 178.73333333333332,
        178.76666666666665, 178.8, 178.83333333333334, 178.86666666666667, 178.9,
        178.93333333333334, 178.96666666666667, 179.0, 179.03333333333333,
        179.06666666666666, 179.1, 179.13333333333333, 179.16666666666666, 179.2,
        179.23333333333332, 179.26666666666665, 179.3, 179.33333333333334,
        179.36666666666667, 179.4, 179.43333333333334, 179.46666666666667, 179.5,
        179.53333333333333, 179.56666666666666, 179.6, 179.63333333333333,
        179.66666666666666, 179.7, 179.73333333333332, 179.76666666666665, 179.8,
        179.83333333333334, 179.86666666666667, 179.9, 179.93333333333334,
        179.96666666666667, 180.0, 180.03333333333333, 180.06666666666666, 180.1,
        180.13333333333333, 180.16666666666666, 180.2, 180.23333333333332,
        180.26666666666665, 180.3, 180.33333333333334, 180.36666666666667, 180.4,
        180.43333333333334, 180.46666666666667, 180.5, 180.53333333333333,
        180.56666666666666, 180.6, 180.63333333333333, 180.66666666666666, 180.7,
        180.73333333333332, 180.76666666666665, 180.8, 180.83333333333334,
        180.86666666666667, 180.9, 180.93333333333334, 180.96666666666667, 181.0,
        181.03333333333333, 181.06666666666666, 181.1, 181.13333333333333,
        181.16666666666666, 181.2, 181.23333333333332, 181.26666666666665, 181.3,
        181.33333333333334, 181.36666666666667, 181.4, 181.43333333333334,
        181.46666666666667, 181.5, 181.53333333333333, 181.56666666666666, 181.6,
        181.63333333333333, 181.66666666666666, 181.7, 181.73333333333332,
        181.76666666666665, 181.8, 181.83333333333334, 181.86666666666667, 181.9,
        181.93333333333334, 181.96666666666667, 182.0, 182.03333333333333,
        182.06666666666666, 182.1, 182.13333333333333, 182.16666666666666, 182.2,
        182.23333333333332, 182.26666666666665, 182.3, 182.33333333333334,
        182.36666666666667, 182.4, 182.43333333333334, 182.46666666666667, 182.5,
        182.53333333333333, 182.56666666666666, 182.6, 182.63333333333333,
        182.66666666666666, 182.7, 182.73333333333332, 182.76666666666665, 182.8,
        182.83333333333334, 182.86666666666667, 182.9, 182.93333333333334,
        182.96666666666667, 183.0, 183.03333333333333, 183.06666666666666, 183.1,
        183.13333333333333, 183.16666666666666, 183.2, 183.23333333333332,
        183.26666666666665, 183.3, 183.33333333333334, 183.36666666666667, 183.4,
        183.43333333333334, 183.46666666666667, 183.5, 183.53333333333333,
        183.56666666666666, 183.6, 183.63333333333333, 183.66666666666666, 183.7,
        183.73333333333332, 183.76666666666665, 183.8, 183.83333333333334,
        183.86666666666667, 183.9, 183.93333333333334, 183.96666666666667, 184.0,
        184.03333333333333, 184.06666666666666, 184.1, 184.13333333333333,
        184.16666666666666, 184.2, 184.23333333333332, 184.26666666666665, 184.3,
        184.33333333333334, 184.36666666666667, 184.4, 184.43333333333334,
        184.46666666666667, 184.5, 184.53333333333333, 184.56666666666666, 184.6,
        184.63333333333333, 184.66666666666666, 184.7, 184.73333333333332,
        184.76666666666665, 184.8, 184.83333333333334, 184.86666666666667, 184.9,
        184.93333333333334, 184.96666666666667, 185.0, 185.03333333333333,
        185.06666666666666, 185.1, 185.13333333333333, 185.16666666666666, 185.2,
        185.23333333333332, 185.26666666666665, 185.3, 185.33333333333334,
        185.36666666666667, 185.4, 185.43333333333334, 185.46666666666667, 185.5,
        185.53333333333333, 185.56666666666666, 185.6, 185.63333333333333,
        185.66666666666666, 185.7, 185.73333333333332, 185.76666666666665, 185.8,
        185.83333333333334, 185.86666666666667, 185.9, 185.93333333333334,
        185.96666666666667, 186.0, 186.03333333333333, 186.06666666666666, 186.1,
        186.13333333333333, 186.16666666666666, 186.2, 186.23333333333332,
        186.26666666666665, 186.3, 186.33333333333334, 186.36666666666667, 186.4,
        186.43333333333334, 186.46666666666667, 186.5, 186.53333333333333,
        186.56666666666666, 186.6, 186.63333333333333, 186.66666666666666, 186.7,
        186.73333333333332, 186.76666666666665, 186.8, 186.83333333333334,
        186.86666666666667, 186.9, 186.93333333333334, 186.96666666666667, 187.0,
        187.03333333333333, 187.06666666666666, 187.1, 187.13333333333333,
        187.16666666666666, 187.2, 187.23333333333332, 187.26666666666665, 187.3,
        187.33333333333334, 187.36666666666667, 187.4, 187.43333333333334,
        187.46666666666667, 187.5, 187.53333333333333, 187.56666666666666, 187.6,
        187.63333333333333, 187.66666666666666, 187.7, 187.73333333333332,
        187.76666666666665, 187.8, 187.83333333333334, 187.86666666666667, 187.9,
        187.93333333333334, 187.96666666666667, 188.0, 188.03333333333333,
        188.06666666666666, 188.1, 188.13333333333333, 188.16666666666666, 188.2,
        188.23333333333332, 188.26666666666665, 188.3, 188.33333333333334,
        188.36666666666667, 188.4, 188.43333333333334, 188.46666666666667, 188.5,
        188.53333333333333, 188.56666666666666, 188.6, 188.63333333333333,
        188.66666666666666, 188.7, 188.73333333333332, 188.76666666666665, 188.8,
        188.83333333333334, 188.86666666666667, 188.9, 188.93333333333334,
        188.96666666666667, 189.0, 189.03333333333333, 189.06666666666666, 189.1,
        189.13333333333333, 189.16666666666666, 189.2, 189.23333333333332,
        189.26666666666665, 189.3, 189.33333333333334, 189.36666666666667, 189.4,
        189.43333333333334, 189.46666666666667, 189.5, 189.53333333333333,
        189.56666666666666, 189.6, 189.63333333333333, 189.66666666666666, 189.7,
        189.73333333333332, 189.76666666666665, 189.8, 189.83333333333334,
        189.86666666666667, 189.9, 189.93333333333334, 189.96666666666667, 190.0,
        190.03333333333333, 190.06666666666666, 190.1, 190.13333333333333,
        190.16666666666666, 190.2, 190.23333333333332, 190.26666666666665, 190.3,
        190.33333333333334, 190.36666666666667, 190.4, 190.43333333333334,
        190.46666666666667, 190.5, 190.53333333333333, 190.56666666666666, 190.6,
        190.63333333333333, 190.66666666666666, 190.7, 190.73333333333332,
        190.76666666666665, 190.8, 190.83333333333334, 190.86666666666667, 190.9,
        190.93333333333334, 190.96666666666667, 191.0, 191.03333333333333,
        191.06666666666666, 191.1, 191.13333333333333, 191.16666666666666, 191.2,
        191.23333333333332, 191.26666666666665, 191.3, 191.33333333333334,
        191.36666666666667, 191.4, 191.43333333333334, 191.46666666666667, 191.5,
        191.53333333333333, 191.56666666666666, 191.6, 191.63333333333333,
        191.66666666666666, 191.7, 191.73333333333332, 191.76666666666665, 191.8,
        191.83333333333334, 191.86666666666667, 191.9, 191.93333333333334,
        191.96666666666667, 192.0, 192.03333333333333, 192.06666666666666, 192.1,
        192.13333333333333, 192.16666666666666, 192.2, 192.23333333333332,
        192.26666666666665, 192.3, 192.33333333333334, 192.36666666666667, 192.4,
        192.43333333333334, 192.46666666666667, 192.5, 192.53333333333333,
        192.56666666666666, 192.6, 192.63333333333333, 192.66666666666666, 192.7,
        192.73333333333332, 192.76666666666665, 192.8, 192.83333333333334,
        192.86666666666667, 192.9, 192.93333333333334, 192.96666666666667, 193.0,
        193.03333333333333, 193.06666666666666, 193.1, 193.13333333333333,
        193.16666666666666, 193.2, 193.23333333333332, 193.26666666666665, 193.3,
        193.33333333333334, 193.36666666666667, 193.4, 193.43333333333334,
        193.46666666666667, 193.5, 193.53333333333333, 193.56666666666666, 193.6,
        193.63333333333333, 193.66666666666666, 193.7, 193.73333333333332,
        193.76666666666665, 193.8, 193.83333333333334, 193.86666666666667, 193.9,
        193.93333333333334, 193.96666666666667, 194.0, 194.03333333333333,
        194.06666666666666, 194.1, 194.13333333333333, 194.16666666666666, 194.2,
        194.23333333333332, 194.26666666666665, 194.3, 194.33333333333334,
        194.36666666666667, 194.4, 194.43333333333334, 194.46666666666667, 194.5,
        194.53333333333333, 194.56666666666666, 194.6, 194.63333333333333,
        194.66666666666666, 194.7, 194.73333333333332, 194.76666666666665, 194.8,
        194.83333333333334, 194.86666666666667, 194.9, 194.93333333333334,
        194.96666666666667, 195.0, 195.03333333333333, 195.06666666666666, 195.1,
        195.13333333333333, 195.16666666666666, 195.2, 195.23333333333332,
        195.26666666666665, 195.3, 195.33333333333334, 195.36666666666667, 195.4,
        195.43333333333334, 195.46666666666667, 195.5, 195.53333333333333,
        195.56666666666666, 195.6, 195.63333333333333, 195.66666666666666, 195.7,
        195.73333333333332, 195.76666666666665, 195.8, 195.83333333333334,
        195.86666666666667, 195.9, 195.93333333333334, 195.96666666666667, 196.0,
        196.03333333333333, 196.06666666666666, 196.1, 196.13333333333333,
        196.16666666666666, 196.2, 196.23333333333332, 196.26666666666665, 196.3,
        196.33333333333334, 196.36666666666667, 196.4, 196.43333333333334,
        196.46666666666667, 196.5, 196.53333333333333, 196.56666666666666, 196.6,
        196.63333333333333, 196.66666666666666, 196.7, 196.73333333333332,
        196.76666666666665, 196.8, 196.83333333333334, 196.86666666666667, 196.9,
        196.93333333333334, 196.96666666666667, 197.0, 197.03333333333333,
        197.06666666666666, 197.1, 197.13333333333333, 197.16666666666666, 197.2,
        197.23333333333332, 197.26666666666665, 197.3, 197.33333333333334,
        197.36666666666667, 197.4, 197.43333333333334, 197.46666666666667, 197.5,
        197.53333333333333, 197.56666666666666, 197.6, 197.63333333333333,
        197.66666666666666, 197.7, 197.73333333333332, 197.76666666666665, 197.8,
        197.83333333333334, 197.86666666666667, 197.9, 197.93333333333334,
        197.96666666666667, 198.0, 198.03333333333333, 198.06666666666666, 198.1,
        198.13333333333333, 198.16666666666666, 198.2, 198.23333333333332,
        198.26666666666665, 198.3, 198.33333333333334, 198.36666666666667, 198.4,
        198.43333333333334, 198.46666666666667, 198.5, 198.53333333333333,
        198.56666666666666, 198.6, 198.63333333333333, 198.66666666666666, 198.7,
        198.73333333333332, 198.76666666666665, 198.8, 198.83333333333334,
        198.86666666666667, 198.9, 198.93333333333334, 198.96666666666667, 199.0,
        199.03333333333333, 199.06666666666666, 199.1, 199.13333333333333,
        199.16666666666666, 199.2, 199.23333333333332, 199.26666666666665, 199.3,
        199.33333333333334, 199.36666666666667, 199.4, 199.43333333333334,
        199.46666666666667, 199.5, 199.53333333333333, 199.56666666666666, 199.6,
        199.63333333333333, 199.66666666666666, 199.7, 199.73333333333332,
        199.76666666666665, 199.8, 199.83333333333334, 199.86666666666667, 199.9,
        199.93333333333334, 199.96666666666667, 200.0, 200.03333333333333,
        200.06666666666666, 200.1, 200.13333333333333, 200.16666666666666, 200.2,
        200.23333333333332, 200.26666666666665, 200.3, 200.33333333333334,
        200.36666666666667, 200.4, 200.43333333333334, 200.46666666666667, 200.5,
        200.53333333333333, 200.56666666666666, 200.6, 200.63333333333333,
        200.66666666666666, 200.7, 200.73333333333332, 200.76666666666665, 200.8,
        200.83333333333334, 200.86666666666667, 200.9, 200.93333333333334,
        200.96666666666667, 201.0, 201.03333333333333, 201.06666666666666, 201.1,
        201.13333333333333, 201.16666666666666, 201.2, 201.23333333333332,
        201.26666666666665, 201.3, 201.33333333333334, 201.36666666666667, 201.4,
        201.43333333333334, 201.46666666666667, 201.5, 201.53333333333333,
        201.56666666666666, 201.6, 201.63333333333333, 201.66666666666666, 201.7,
        201.73333333333332, 201.76666666666665, 201.8, 201.83333333333334,
        201.86666666666667, 201.9, 201.93333333333334, 201.96666666666667, 202.0,
        202.03333333333333, 202.06666666666666, 202.1, 202.13333333333333,
        202.16666666666666, 202.2, 202.23333333333332, 202.26666666666665, 202.3,
        202.33333333333334, 202.36666666666667, 202.4, 202.43333333333334,
        202.46666666666667, 202.5, 202.53333333333333, 202.56666666666666, 202.6,
        202.63333333333333, 202.66666666666666, 202.7, 202.73333333333332,
        202.76666666666665, 202.8, 202.83333333333334, 202.86666666666667, 202.9,
        202.93333333333334, 202.96666666666667, 203.0, 203.03333333333333,
        203.06666666666666, 203.1, 203.13333333333333, 203.16666666666666, 203.2,
        203.23333333333332, 203.26666666666665, 203.3, 203.33333333333334,
        203.36666666666667, 203.4, 203.43333333333334, 203.46666666666667, 203.5,
        203.53333333333333, 203.56666666666666, 203.6, 203.63333333333333,
        203.66666666666666, 203.7, 203.73333333333332, 203.76666666666665, 203.8,
        203.83333333333334, 203.86666666666667, 203.9, 203.93333333333334,
        203.96666666666667, 204.0, 204.03333333333333, 204.06666666666666, 204.1,
        204.13333333333333, 204.16666666666666, 204.2, 204.23333333333332,
        204.26666666666665, 204.3, 204.33333333333334, 204.36666666666667, 204.4,
        204.43333333333334, 204.46666666666667, 204.5, 204.53333333333333,
        204.56666666666666, 204.6, 204.63333333333333, 204.66666666666666, 204.7,
        204.73333333333332, 204.76666666666665, 204.8, 204.83333333333334,
        204.86666666666667, 204.9, 204.93333333333334, 204.96666666666667, 205.0,
        205.03333333333333, 205.06666666666666, 205.1, 205.13333333333333,
        205.16666666666666, 205.2, 205.23333333333332, 205.26666666666665,
        205.29999999999998, 205.33333333333334, 205.36666666666667, 205.4,
        205.43333333333334, 205.46666666666667, 205.5, 205.53333333333333,
        205.56666666666666, 205.6, 205.63333333333333, 205.66666666666666, 205.7,
        205.73333333333332, 205.76666666666665, 205.79999999999998,
        205.83333333333334, 205.86666666666667, 205.9, 205.93333333333334,
        205.96666666666667, 206.0, 206.03333333333333, 206.06666666666666, 206.1,
        206.13333333333333, 206.16666666666666, 206.2, 206.23333333333332,
        206.26666666666665, 206.29999999999998, 206.33333333333334,
        206.36666666666667, 206.4, 206.43333333333334, 206.46666666666667, 206.5,
        206.53333333333333, 206.56666666666666, 206.6, 206.63333333333333,
        206.66666666666666, 206.7, 206.73333333333332, 206.76666666666665,
        206.79999999999998, 206.83333333333334, 206.86666666666667, 206.9,
        206.93333333333334, 206.96666666666667, 207.0, 207.03333333333333,
        207.06666666666666, 207.1, 207.13333333333333, 207.16666666666666, 207.2,
        207.23333333333332, 207.26666666666665, 207.29999999999998,
        207.33333333333334, 207.36666666666667, 207.4, 207.43333333333334,
        207.46666666666667, 207.5, 207.53333333333333, 207.56666666666666, 207.6,
        207.63333333333333, 207.66666666666666, 207.7, 207.73333333333332,
        207.76666666666665, 207.79999999999998, 207.83333333333334,
        207.86666666666667, 207.9, 207.93333333333334, 207.96666666666667, 208.0,
        208.03333333333333, 208.06666666666666, 208.1, 208.13333333333333,
        208.16666666666666, 208.2, 208.23333333333332, 208.26666666666665,
        208.29999999999998, 208.33333333333334, 208.36666666666667, 208.4,
        208.43333333333334, 208.46666666666667, 208.5, 208.53333333333333,
        208.56666666666666, 208.6, 208.63333333333333, 208.66666666666666, 208.7,
        208.73333333333332, 208.76666666666665, 208.79999999999998,
        208.83333333333334, 208.86666666666667, 208.9, 208.93333333333334,
        208.96666666666667, 209.0, 209.03333333333333, 209.06666666666666, 209.1,
        209.13333333333333, 209.16666666666666, 209.2, 209.23333333333332,
        209.26666666666665, 209.29999999999998, 209.33333333333334,
        209.36666666666667, 209.4, 209.43333333333334, 209.46666666666667, 209.5,
        209.53333333333333, 209.56666666666666, 209.6, 209.63333333333333,
        209.66666666666666, 209.7, 209.73333333333332, 209.76666666666665,
        209.79999999999998, 209.83333333333334, 209.86666666666667, 209.9,
        209.93333333333334, 209.96666666666667, 210.0, 210.03333333333333,
        210.06666666666666, 210.1, 210.13333333333333, 210.16666666666666, 210.2,
        210.23333333333332, 210.26666666666665, 210.29999999999998,
        210.33333333333334, 210.36666666666667, 210.4, 210.43333333333334,
        210.46666666666667, 210.5, 210.53333333333333, 210.56666666666666, 210.6,
        210.63333333333333, 210.66666666666666, 210.7, 210.73333333333332,
        210.76666666666665, 210.79999999999998, 210.83333333333334,
        210.86666666666667, 210.9, 210.93333333333334, 210.96666666666667, 211.0,
        211.03333333333333, 211.06666666666666, 211.1, 211.13333333333333,
        211.16666666666666, 211.2, 211.23333333333332, 211.26666666666665,
        211.29999999999998, 211.33333333333334, 211.36666666666667, 211.4,
        211.43333333333334, 211.46666666666667, 211.5, 211.53333333333333,
        211.56666666666666, 211.6, 211.63333333333333, 211.66666666666666, 211.7,
        211.73333333333332, 211.76666666666665, 211.79999999999998,
        211.83333333333334, 211.86666666666667, 211.9, 211.93333333333334,
        211.96666666666667, 212.0, 212.03333333333333, 212.06666666666666, 212.1,
        212.13333333333333, 212.16666666666666, 212.2, 212.23333333333332,
        212.26666666666665, 212.29999999999998, 212.33333333333334,
        212.36666666666667, 212.4, 212.43333333333334, 212.46666666666667, 212.5,
        212.53333333333333, 212.56666666666666, 212.6, 212.63333333333333,
        212.66666666666666, 212.7, 212.73333333333332, 212.76666666666665,
        212.79999999999998, 212.83333333333334, 212.86666666666667, 212.9,
        212.93333333333334, 212.96666666666667, 213.0, 213.03333333333333,
        213.06666666666666, 213.1, 213.13333333333333, 213.16666666666666, 213.2,
        213.23333333333332, 213.26666666666665, 213.29999999999998,
        213.33333333333334, 213.36666666666667, 213.4, 213.43333333333334,
        213.46666666666667, 213.5, 213.53333333333333, 213.56666666666666, 213.6,
        213.63333333333333, 213.66666666666666, 213.7, 213.73333333333332,
        213.76666666666665, 213.79999999999998, 213.83333333333334,
        213.86666666666667, 213.9, 213.93333333333334, 213.96666666666667, 214.0,
        214.03333333333333, 214.06666666666666, 214.1, 214.13333333333333,
        214.16666666666666, 214.2, 214.23333333333332, 214.26666666666665,
        214.29999999999998, 214.33333333333334, 214.36666666666667, 214.4,
        214.43333333333334, 214.46666666666667, 214.5, 214.53333333333333,
        214.56666666666666, 214.6, 214.63333333333333, 214.66666666666666, 214.7,
        214.73333333333332, 214.76666666666665, 214.79999999999998,
        214.83333333333334, 214.86666666666667, 214.9, 214.93333333333334,
        214.96666666666667, 215.0, 215.03333333333333, 215.06666666666666, 215.1,
        215.13333333333333, 215.16666666666666, 215.2, 215.23333333333332,
        215.26666666666665, 215.29999999999998, 215.33333333333334,
        215.36666666666667, 215.4, 215.43333333333334, 215.46666666666667, 215.5,
        215.53333333333333, 215.56666666666666, 215.6, 215.63333333333333,
        215.66666666666666, 215.7, 215.73333333333332, 215.76666666666665,
        215.79999999999998, 215.83333333333334, 215.86666666666667, 215.9,
        215.93333333333334, 215.96666666666667, 216.0, 216.03333333333333,
        216.06666666666666, 216.1, 216.13333333333333, 216.16666666666666, 216.2,
        216.23333333333332, 216.26666666666665, 216.29999999999998,
        216.33333333333334, 216.36666666666667, 216.4, 216.43333333333334,
        216.46666666666667, 216.5, 216.53333333333333, 216.56666666666666, 216.6,
        216.63333333333333, 216.66666666666666, 216.7, 216.73333333333332,
        216.76666666666665, 216.79999999999998, 216.83333333333334,
        216.86666666666667, 216.9, 216.93333333333334, 216.96666666666667, 217.0,
        217.03333333333333, 217.06666666666666, 217.1, 217.13333333333333,
        217.16666666666666, 217.2, 217.23333333333332, 217.26666666666665,
        217.29999999999998, 217.33333333333334, 217.36666666666667, 217.4,
        217.43333333333334, 217.46666666666667, 217.5, 217.53333333333333,
        217.56666666666666, 217.6, 217.63333333333333, 217.66666666666666, 217.7,
        217.73333333333332, 217.76666666666665, 217.79999999999998,
        217.83333333333334, 217.86666666666667, 217.9, 217.93333333333334,
        217.96666666666667, 218.0, 218.03333333333333, 218.06666666666666, 218.1,
        218.13333333333333, 218.16666666666666, 218.2, 218.23333333333332,
        218.26666666666665, 218.29999999999998, 218.33333333333334,
        218.36666666666667, 218.4, 218.43333333333334, 218.46666666666667, 218.5,
        218.53333333333333, 218.56666666666666, 218.6, 218.63333333333333,
        218.66666666666666, 218.7, 218.73333333333332, 218.76666666666665,
        218.79999999999998, 218.83333333333334, 218.86666666666667, 218.9,
        218.93333333333334, 218.96666666666667, 219.0, 219.03333333333333,
        219.06666666666666, 219.1, 219.13333333333333, 219.16666666666666, 219.2,
        219.23333333333332, 219.26666666666665, 219.29999999999998,
        219.33333333333334, 219.36666666666667, 219.4, 219.43333333333334,
        219.46666666666667, 219.5, 219.53333333333333, 219.56666666666666, 219.6,
        219.63333333333333, 219.66666666666666, 219.7, 219.73333333333332,
        219.76666666666665, 219.79999999999998, 219.83333333333334,
        219.86666666666667, 219.9, 219.93333333333334, 219.96666666666667, 220.0,
        220.03333333333333, 220.06666666666666, 220.1, 220.13333333333333,
        220.16666666666666, 220.2, 220.23333333333332, 220.26666666666665,
        220.29999999999998, 220.33333333333334, 220.36666666666667, 220.4,
        220.43333333333334, 220.46666666666667, 220.5, 220.53333333333333,
        220.56666666666666, 220.6, 220.63333333333333, 220.66666666666666, 220.7,
        220.73333333333332, 220.76666666666665, 220.79999999999998,
        220.83333333333334, 220.86666666666667, 220.9, 220.93333333333334,
        220.96666666666667, 221.0, 221.03333333333333, 221.06666666666666, 221.1,
        221.13333333333333, 221.16666666666666, 221.2, 221.23333333333332,
        221.26666666666665, 221.29999999999998, 221.33333333333334,
        221.36666666666667, 221.4, 221.43333333333334, 221.46666666666667, 221.5,
        221.53333333333333, 221.56666666666666, 221.6, 221.63333333333333,
        221.66666666666666, 221.7, 221.73333333333332, 221.76666666666665,
        221.79999999999998, 221.83333333333334, 221.86666666666667, 221.9,
        221.93333333333334, 221.96666666666667, 222.0, 222.03333333333333,
        222.06666666666666, 222.1, 222.13333333333333, 222.16666666666666, 222.2,
        222.23333333333332, 222.26666666666665, 222.29999999999998,
        222.33333333333334, 222.36666666666667, 222.4, 222.43333333333334,
        222.46666666666667, 222.5, 222.53333333333333, 222.56666666666666, 222.6,
        222.63333333333333, 222.66666666666666, 222.7, 222.73333333333332,
        222.76666666666665, 222.79999999999998, 222.83333333333334,
        222.86666666666667, 222.9, 222.93333333333334, 222.96666666666667, 223.0,
        223.03333333333333, 223.06666666666666, 223.1, 223.13333333333333,
        223.16666666666666, 223.2, 223.23333333333332, 223.26666666666665,
        223.29999999999998, 223.33333333333334, 223.36666666666667, 223.4,
        223.43333333333334, 223.46666666666667, 223.5, 223.53333333333333,
        223.56666666666666, 223.6, 223.63333333333333, 223.66666666666666, 223.7,
        223.73333333333332, 223.76666666666665, 223.79999999999998,
        223.83333333333334, 223.86666666666667, 223.9, 223.93333333333334,
        223.96666666666667, 224.0, 224.03333333333333, 224.06666666666666, 224.1,
        224.13333333333333, 224.16666666666666, 224.2, 224.23333333333332,
        224.26666666666665, 224.29999999999998, 224.33333333333334,
        224.36666666666667, 224.4, 224.43333333333334, 224.46666666666667, 224.5,
        224.53333333333333, 224.56666666666666, 224.6, 224.63333333333333,
        224.66666666666666, 224.7, 224.73333333333332, 224.76666666666665,
        224.79999999999998, 224.83333333333334, 224.86666666666667, 224.9,
        224.93333333333334, 224.96666666666667, 225.0, 225.03333333333333,
        225.06666666666666, 225.1, 225.13333333333333, 225.16666666666666, 225.2,
        225.23333333333332, 225.26666666666665, 225.29999999999998,
        225.33333333333334, 225.36666666666667, 225.4, 225.43333333333334,
        225.46666666666667, 225.5, 225.53333333333333, 225.56666666666666, 225.6,
        225.63333333333333, 225.66666666666666, 225.7, 225.73333333333332,
        225.76666666666665, 225.79999999999998, 225.83333333333334,
        225.86666666666667, 225.9, 225.93333333333334, 225.96666666666667, 226.0,
        226.03333333333333, 226.06666666666666, 226.1, 226.13333333333333,
        226.16666666666666, 226.2, 226.23333333333332, 226.26666666666665,
        226.29999999999998, 226.33333333333334, 226.36666666666667, 226.4,
        226.43333333333334, 226.46666666666667, 226.5, 226.53333333333333,
        226.56666666666666, 226.6, 226.63333333333333, 226.66666666666666, 226.7,
        226.73333333333332, 226.76666666666665, 226.79999999999998,
        226.83333333333334, 226.86666666666667, 226.9, 226.93333333333334,
        226.96666666666667, 227.0, 227.03333333333333, 227.06666666666666, 227.1,
        227.13333333333333, 227.16666666666666, 227.2, 227.23333333333332,
        227.26666666666665, 227.29999999999998, 227.33333333333334,
        227.36666666666667, 227.4, 227.43333333333334, 227.46666666666667, 227.5,
        227.53333333333333, 227.56666666666666, 227.6, 227.63333333333333,
        227.66666666666666, 227.7, 227.73333333333332, 227.76666666666665,
        227.79999999999998, 227.83333333333334, 227.86666666666667, 227.9,
        227.93333333333334, 227.96666666666667, 228.0, 228.03333333333333,
        228.06666666666666, 228.1, 228.13333333333333, 228.16666666666666, 228.2,
        228.23333333333332, 228.26666666666665, 228.29999999999998,
        228.33333333333334, 228.36666666666667, 228.4, 228.43333333333334,
        228.46666666666667, 228.5, 228.53333333333333, 228.56666666666666, 228.6,
        228.63333333333333, 228.66666666666666, 228.7, 228.73333333333332,
        228.76666666666665, 228.79999999999998, 228.83333333333334,
        228.86666666666667, 228.9, 228.93333333333334, 228.96666666666667, 229.0,
        229.03333333333333, 229.06666666666666, 229.1, 229.13333333333333,
        229.16666666666666, 229.2, 229.23333333333332, 229.26666666666665,
        229.29999999999998, 229.33333333333334, 229.36666666666667, 229.4,
        229.43333333333334, 229.46666666666667, 229.5, 229.53333333333333,
        229.56666666666666, 229.6, 229.63333333333333, 229.66666666666666, 229.7,
        229.73333333333332, 229.76666666666665, 229.79999999999998,
        229.83333333333334, 229.86666666666667, 229.9, 229.93333333333334,
        229.96666666666667, 230.0, 230.03333333333333, 230.06666666666666, 230.1,
        230.13333333333333, 230.16666666666666, 230.2, 230.23333333333332,
        230.26666666666665, 230.29999999999998, 230.33333333333334,
        230.36666666666667, 230.4, 230.43333333333334, 230.46666666666667, 230.5,
        230.53333333333333, 230.56666666666666, 230.6, 230.63333333333333,
        230.66666666666666, 230.7, 230.73333333333332, 230.76666666666665,
        230.79999999999998, 230.83333333333334, 230.86666666666667, 230.9,
        230.93333333333334, 230.96666666666667, 231.0, 231.03333333333333,
        231.06666666666666, 231.1, 231.13333333333333, 231.16666666666666, 231.2,
        231.23333333333332, 231.26666666666665, 231.29999999999998,
        231.33333333333334, 231.36666666666667, 231.4, 231.43333333333334,
        231.46666666666667, 231.5, 231.53333333333333, 231.56666666666666, 231.6,
        231.63333333333333, 231.66666666666666, 231.7, 231.73333333333332,
        231.76666666666665, 231.79999999999998, 231.83333333333334,
        231.86666666666667, 231.9, 231.93333333333334, 231.96666666666667, 232.0,
        232.03333333333333, 232.06666666666666, 232.1, 232.13333333333333,
        232.16666666666666, 232.2, 232.23333333333332, 232.26666666666665,
        232.29999999999998, 232.33333333333334, 232.36666666666667, 232.4,
        232.43333333333334, 232.46666666666667, 232.5, 232.53333333333333,
        232.56666666666666, 232.6, 232.63333333333333, 232.66666666666666, 232.7,
        232.73333333333332, 232.76666666666665, 232.79999999999998,
        232.83333333333334, 232.86666666666667, 232.9, 232.93333333333334,
        232.96666666666667, 233.0, 233.03333333333333, 233.06666666666666, 233.1,
        233.13333333333333, 233.16666666666666, 233.2, 233.23333333333332,
        233.26666666666665, 233.29999999999998, 233.33333333333334,
        233.36666666666667, 233.4, 233.43333333333334, 233.46666666666667, 233.5,
        233.53333333333333, 233.56666666666666, 233.6, 233.63333333333333,
        233.66666666666666, 233.7, 233.73333333333332, 233.76666666666665,
        233.79999999999998, 233.83333333333334, 233.86666666666667, 233.9,
        233.93333333333334, 233.96666666666667, 234.0, 234.03333333333333,
        234.06666666666666, 234.1, 234.13333333333333, 234.16666666666666, 234.2,
        234.23333333333332, 234.26666666666665, 234.29999999999998,
        234.33333333333334, 234.36666666666667, 234.4, 234.43333333333334,
        234.46666666666667, 234.5, 234.53333333333333, 234.56666666666666, 234.6,
        234.63333333333333, 234.66666666666666, 234.7, 234.73333333333332,
        234.76666666666665, 234.79999999999998, 234.83333333333334,
        234.86666666666667, 234.9, 234.93333333333334, 234.96666666666667, 235.0,
        235.03333333333333, 235.06666666666666, 235.1, 235.13333333333333,
        235.16666666666666, 235.2, 235.23333333333332, 235.26666666666665,
        235.29999999999998, 235.33333333333334, 235.36666666666667, 235.4,
        235.43333333333334, 235.46666666666667, 235.5, 235.53333333333333,
        235.56666666666666, 235.6, 235.63333333333333, 235.66666666666666, 235.7,
        235.73333333333332, 235.76666666666665, 235.79999999999998,
        235.83333333333334, 235.86666666666667, 235.9, 235.93333333333334,
        235.96666666666667, 236.0, 236.03333333333333, 236.06666666666666, 236.1,
        236.13333333333333, 236.16666666666666, 236.2, 236.23333333333332,
        236.26666666666665, 236.29999999999998, 236.33333333333334,
        236.36666666666667, 236.4, 236.43333333333334, 236.46666666666667, 236.5,
        236.53333333333333, 236.56666666666666, 236.6, 236.63333333333333,
        236.66666666666666, 236.7, 236.73333333333332, 236.76666666666665,
        236.79999999999998, 236.83333333333334, 236.86666666666667, 236.9,
        236.93333333333334, 236.96666666666667, 237.0, 237.03333333333333,
        237.06666666666666, 237.1, 237.13333333333333, 237.16666666666666, 237.2,
        237.23333333333332, 237.26666666666665, 237.29999999999998,
        237.33333333333334, 237.36666666666667, 237.4, 237.43333333333334,
        237.46666666666667, 237.5, 237.53333333333333, 237.56666666666666, 237.6,
        237.63333333333333, 237.66666666666666, 237.7, 237.73333333333332,
        237.76666666666665, 237.79999999999998, 237.83333333333334,
        237.86666666666667, 237.9, 237.93333333333334, 237.96666666666667, 238.0,
        238.03333333333333, 238.06666666666666, 238.1, 238.13333333333333,
        238.16666666666666, 238.2, 238.23333333333332, 238.26666666666665,
        238.29999999999998, 238.33333333333334, 238.36666666666667, 238.4,
        238.43333333333334, 238.46666666666667, 238.5, 238.53333333333333,
        238.56666666666666, 238.6, 238.63333333333333, 238.66666666666666, 238.7,
        238.73333333333332, 238.76666666666665, 238.79999999999998,
        238.83333333333334, 238.86666666666667, 238.9, 238.93333333333334,
        238.96666666666667, 239.0, 239.03333333333333, 239.06666666666666, 239.1,
        239.13333333333333, 239.16666666666666, 239.2, 239.23333333333332,
        239.26666666666665, 239.29999999999998, 239.33333333333334,
        239.36666666666667, 239.4, 239.43333333333334, 239.46666666666667, 239.5,
        239.53333333333333, 239.56666666666666, 239.6, 239.63333333333333,
        239.66666666666666, 239.7, 239.73333333333332, 239.76666666666665,
        239.79999999999998, 239.83333333333334, 239.86666666666667, 239.9,
        239.93333333333334, 239.96666666666667, 240.0, 240.03333333333333,
        240.06666666666666, 240.1, 240.13333333333333, 240.16666666666666, 240.2,
        240.23333333333332, 240.26666666666665, 240.29999999999998,
        240.33333333333334, 240.36666666666667, 240.4, 240.43333333333334,
        240.46666666666667, 240.5, 240.53333333333333, 240.56666666666666, 240.6,
        240.63333333333333, 240.66666666666666, 240.7, 240.73333333333332,
        240.76666666666665, 240.79999999999998, 240.83333333333334,
        240.86666666666667, 240.9, 240.93333333333334, 240.96666666666667, 241.0,
        241.03333333333333, 241.06666666666666, 241.1, 241.13333333333333,
        241.16666666666666, 241.2, 241.23333333333332, 241.26666666666665,
        241.29999999999998, 241.33333333333334, 241.36666666666667, 241.4,
        241.43333333333334, 241.46666666666667, 241.5, 241.53333333333333,
        241.56666666666666, 241.6, 241.63333333333333, 241.66666666666666, 241.7,
        241.73333333333332, 241.76666666666665, 241.79999999999998,
        241.83333333333334, 241.86666666666667, 241.9, 241.93333333333334,
        241.96666666666667, 242.0, 242.03333333333333, 242.06666666666666, 242.1,
        242.13333333333333, 242.16666666666666, 242.2, 242.23333333333332,
        242.26666666666665, 242.29999999999998, 242.33333333333334,
        242.36666666666667, 242.4, 242.43333333333334, 242.46666666666667, 242.5,
        242.53333333333333, 242.56666666666666, 242.6, 242.63333333333333,
        242.66666666666666, 242.7, 242.73333333333332, 242.76666666666665,
        242.79999999999998, 242.83333333333334, 242.86666666666667, 242.9,
        242.93333333333334, 242.96666666666667, 243.0, 243.03333333333333,
        243.06666666666666, 243.1, 243.13333333333333, 243.16666666666666, 243.2,
        243.23333333333332, 243.26666666666665, 243.29999999999998,
        243.33333333333334, 243.36666666666667, 243.4, 243.43333333333334,
        243.46666666666667, 243.5, 243.53333333333333, 243.56666666666666, 243.6,
        243.63333333333333, 243.66666666666666, 243.7, 243.73333333333332,
        243.76666666666665, 243.79999999999998, 243.83333333333334,
        243.86666666666667, 243.9, 243.93333333333334, 243.96666666666667, 244.0,
        244.03333333333333, 244.06666666666666, 244.1, 244.13333333333333,
        244.16666666666666, 244.2, 244.23333333333332, 244.26666666666665,
        244.29999999999998, 244.33333333333334, 244.36666666666667, 244.4,
        244.43333333333334, 244.46666666666667, 244.5, 244.53333333333333,
        244.56666666666666, 244.6, 244.63333333333333, 244.66666666666666, 244.7,
        244.73333333333332, 244.76666666666665, 244.79999999999998,
        244.83333333333334, 244.86666666666667, 244.9, 244.93333333333334,
        244.96666666666667, 245.0, 245.03333333333333, 245.06666666666666, 245.1,
        245.13333333333333, 245.16666666666666, 245.2, 245.23333333333332,
        245.26666666666665, 245.29999999999998, 245.33333333333334,
        245.36666666666667, 245.4, 245.43333333333334, 245.46666666666667, 245.5,
        245.53333333333333, 245.56666666666666, 245.6, 245.63333333333333,
        245.66666666666666, 245.7, 245.73333333333332, 245.76666666666665,
        245.79999999999998, 245.83333333333334, 245.86666666666667, 245.9,
        245.93333333333334, 245.96666666666667, 246.0, 246.03333333333333,
        246.06666666666666, 246.1, 246.13333333333333, 246.16666666666666, 246.2,
        246.23333333333332, 246.26666666666665, 246.29999999999998,
        246.33333333333334, 246.36666666666667, 246.4, 246.43333333333334,
        246.46666666666667, 246.5, 246.53333333333333, 246.56666666666666, 246.6,
        246.63333333333333, 246.66666666666666, 246.7, 246.73333333333332,
        246.76666666666665, 246.79999999999998, 246.83333333333334,
        246.86666666666667, 246.9, 246.93333333333334, 246.96666666666667, 247.0,
        247.03333333333333, 247.06666666666666, 247.1, 247.13333333333333,
        247.16666666666666, 247.2, 247.23333333333332, 247.26666666666665,
        247.29999999999998, 247.33333333333334, 247.36666666666667, 247.4,
        247.43333333333334, 247.46666666666667, 247.5, 247.53333333333333,
        247.56666666666666, 247.6, 247.63333333333333, 247.66666666666666, 247.7,
        247.73333333333332, 247.76666666666665, 247.79999999999998,
        247.83333333333334, 247.86666666666667, 247.9, 247.93333333333334,
        247.96666666666667, 248.0, 248.03333333333333, 248.06666666666666, 248.1,
        248.13333333333333, 248.16666666666666, 248.2, 248.23333333333332,
        248.26666666666665, 248.29999999999998, 248.33333333333334,
        248.36666666666667, 248.4, 248.43333333333334, 248.46666666666667, 248.5,
        248.53333333333333, 248.56666666666666, 248.6, 248.63333333333333,
        248.66666666666666, 248.7, 248.73333333333332, 248.76666666666665,
        248.79999999999998, 248.83333333333334, 248.86666666666667, 248.9,
        248.93333333333334, 248.96666666666667, 249.0, 249.03333333333333,
        249.06666666666666, 249.1, 249.13333333333333, 249.16666666666666, 249.2,
        249.23333333333332, 249.26666666666665, 249.29999999999998,
        249.33333333333334, 249.36666666666667, 249.4, 249.43333333333334,
        249.46666666666667, 249.5, 249.53333333333333, 249.56666666666666, 249.6,
        249.63333333333333, 249.66666666666666, 249.7, 249.73333333333332,
        249.76666666666665, 249.79999999999998, 249.83333333333334,
        249.86666666666667, 249.9, 249.93333333333334, 249.96666666666667, 250.0,
        250.03333333333333, 250.06666666666666, 250.1, 250.13333333333333,
        250.16666666666666, 250.2, 250.23333333333332, 250.26666666666665,
        250.29999999999998, 250.33333333333334, 250.36666666666667, 250.4,
        250.43333333333334, 250.46666666666667, 250.5, 250.53333333333333,
        250.56666666666666, 250.6, 250.63333333333333, 250.66666666666666, 250.7,
        250.73333333333332, 250.76666666666665, 250.79999999999998,
        250.83333333333334, 250.86666666666667, 250.9, 250.93333333333334,
        250.96666666666667, 251.0, 251.03333333333333, 251.06666666666666, 251.1,
        251.13333333333333, 251.16666666666666, 251.2, 251.23333333333332,
        251.26666666666665, 251.29999999999998, 251.33333333333334,
        251.36666666666667, 251.4, 251.43333333333334, 251.46666666666667, 251.5,
        251.53333333333333, 251.56666666666666, 251.6, 251.63333333333333,
        251.66666666666666, 251.7, 251.73333333333332, 251.76666666666665,
        251.79999999999998, 251.83333333333334, 251.86666666666667, 251.9,
        251.93333333333334, 251.96666666666667, 252.0, 252.03333333333333,
        252.06666666666666, 252.1, 252.13333333333333, 252.16666666666666, 252.2,
        252.23333333333332, 252.26666666666665, 252.29999999999998,
        252.33333333333334, 252.36666666666667, 252.4, 252.43333333333334,
        252.46666666666667, 252.5, 252.53333333333333, 252.56666666666666, 252.6,
        252.63333333333333, 252.66666666666666, 252.7, 252.73333333333332,
        252.76666666666665, 252.79999999999998, 252.83333333333334,
        252.86666666666667, 252.9, 252.93333333333334, 252.96666666666667, 253.0,
        253.03333333333333, 253.06666666666666, 253.1, 253.13333333333333,
        253.16666666666666, 253.2, 253.23333333333332, 253.26666666666665,
        253.29999999999998, 253.33333333333334, 253.36666666666667, 253.4,
        253.43333333333334, 253.46666666666667, 253.5, 253.53333333333333,
        253.56666666666666, 253.6, 253.63333333333333, 253.66666666666666, 253.7,
        253.73333333333332, 253.76666666666665, 253.79999999999998,
        253.83333333333334, 253.86666666666667, 253.9, 253.93333333333334,
        253.96666666666667, 254.0, 254.03333333333333, 254.06666666666666, 254.1,
        254.13333333333333, 254.16666666666666, 254.2, 254.23333333333332,
        254.26666666666665, 254.29999999999998, 254.33333333333334,
        254.36666666666667, 254.4, 254.43333333333334, 254.46666666666667, 254.5,
        254.53333333333333, 254.56666666666666, 254.6, 254.63333333333333,
        254.66666666666666, 254.7, 254.73333333333332, 254.76666666666665,
        254.79999999999998, 254.83333333333334, 254.86666666666667, 254.9,
        254.93333333333334, 254.96666666666667, 255.0, 255.03333333333333,
        255.06666666666666, 255.1, 255.13333333333333, 255.16666666666666, 255.2,
        255.23333333333332, 255.26666666666665, 255.29999999999998,
        255.33333333333334, 255.36666666666667, 255.4, 255.43333333333334,
        255.46666666666667, 255.5, 255.53333333333333, 255.56666666666666, 255.6,
        255.63333333333333, 255.66666666666666, 255.7, 255.73333333333332,
        255.76666666666665, 255.79999999999998, 255.83333333333334,
        255.86666666666667, 255.9, 255.93333333333334, 255.96666666666667, 256.0,
        256.0333333333333, 256.06666666666666, 256.1, 256.13333333333333,
        256.16666666666669, 256.2, 256.23333333333335, 256.26666666666665, 256.3,
        256.33333333333331, 256.36666666666667, 256.4, 256.43333333333334,
        256.46666666666664, 256.5, 256.5333333333333, 256.56666666666666, 256.6,
        256.63333333333333, 256.66666666666669, 256.7, 256.73333333333335,
        256.76666666666665, 256.8, 256.83333333333331, 256.86666666666667, 256.9,
        256.93333333333334, 256.96666666666664, 257.0, 257.0333333333333,
        257.06666666666666, 257.1, 257.13333333333333, 257.16666666666669, 257.2,
        257.23333333333335, 257.26666666666665, 257.3, 257.33333333333331,
        257.36666666666667, 257.4, 257.43333333333334, 257.46666666666664, 257.5,
        257.5333333333333, 257.56666666666666, 257.6, 257.63333333333333,
        257.66666666666669, 257.7, 257.73333333333335, 257.76666666666665, 257.8,
        257.83333333333331, 257.86666666666667, 257.9, 257.93333333333334,
        257.96666666666664, 258.0, 258.0333333333333, 258.06666666666666, 258.1,
        258.13333333333333, 258.16666666666669, 258.2, 258.23333333333335,
        258.26666666666665, 258.3, 258.33333333333331, 258.36666666666667, 258.4,
        258.43333333333334, 258.46666666666664, 258.5, 258.5333333333333,
        258.56666666666666, 258.6, 258.63333333333333, 258.66666666666669, 258.7,
        258.73333333333335, 258.76666666666665, 258.8, 258.83333333333331,
        258.86666666666667, 258.9, 258.93333333333334, 258.96666666666664, 259.0,
        259.0333333333333, 259.06666666666666, 259.1, 259.13333333333333,
        259.16666666666669, 259.2, 259.23333333333335, 259.26666666666665, 259.3,
        259.33333333333331, 259.36666666666667, 259.4, 259.43333333333334,
        259.46666666666664, 259.5, 259.5333333333333, 259.56666666666666, 259.6,
        259.63333333333333, 259.66666666666669, 259.7, 259.73333333333335,
        259.76666666666665, 259.8, 259.83333333333331, 259.86666666666667, 259.9,
        259.93333333333334, 259.96666666666664, 260.0, 260.0333333333333,
        260.06666666666666, 260.1, 260.13333333333333, 260.16666666666669, 260.2,
        260.23333333333335, 260.26666666666665, 260.3, 260.33333333333331,
        260.36666666666667, 260.4, 260.43333333333334, 260.46666666666664, 260.5,
        260.5333333333333, 260.56666666666666, 260.6, 260.63333333333333,
        260.66666666666669, 260.7, 260.73333333333335, 260.76666666666665, 260.8,
        260.83333333333331, 260.86666666666667, 260.9, 260.93333333333334,
        260.96666666666664, 261.0, 261.0333333333333, 261.06666666666666, 261.1,
        261.13333333333333, 261.16666666666669, 261.2, 261.23333333333335,
        261.26666666666665, 261.3, 261.33333333333331, 261.36666666666667, 261.4,
        261.43333333333334, 261.46666666666664, 261.5, 261.5333333333333,
        261.56666666666666, 261.6, 261.63333333333333, 261.66666666666669, 261.7,
        261.73333333333335, 261.76666666666665, 261.8, 261.83333333333331,
        261.86666666666667, 261.9, 261.93333333333334, 261.96666666666664, 262.0,
        262.0333333333333, 262.06666666666666, 262.1, 262.13333333333333,
        262.16666666666669, 262.2, 262.23333333333335, 262.26666666666665, 262.3,
        262.33333333333331, 262.36666666666667, 262.4, 262.43333333333334,
        262.46666666666664, 262.5, 262.5333333333333, 262.56666666666666, 262.6,
        262.63333333333333, 262.66666666666669, 262.7, 262.73333333333335,
        262.76666666666665, 262.8, 262.83333333333331, 262.86666666666667, 262.9,
        262.93333333333334, 262.96666666666664, 263.0, 263.0333333333333,
        263.06666666666666, 263.1, 263.13333333333333, 263.16666666666669, 263.2,
        263.23333333333335, 263.26666666666665, 263.3, 263.33333333333331,
        263.36666666666667, 263.4, 263.43333333333334, 263.46666666666664, 263.5,
        263.5333333333333, 263.56666666666666, 263.6, 263.63333333333333,
        263.66666666666669, 263.7, 263.73333333333335, 263.76666666666665, 263.8,
        263.83333333333331, 263.86666666666667, 263.9, 263.93333333333334,
        263.96666666666664, 264.0, 264.0333333333333, 264.06666666666666, 264.1,
        264.13333333333333, 264.16666666666669, 264.2, 264.23333333333335,
        264.26666666666665, 264.3, 264.33333333333331, 264.36666666666667, 264.4,
        264.43333333333334, 264.46666666666664, 264.5, 264.5333333333333,
        264.56666666666666, 264.6, 264.63333333333333, 264.66666666666669, 264.7,
        264.73333333333335, 264.76666666666665, 264.8, 264.83333333333331,
        264.86666666666667, 264.9, 264.93333333333334, 264.96666666666664, 265.0,
        265.0333333333333, 265.06666666666666, 265.1, 265.13333333333333,
        265.16666666666669, 265.2, 265.23333333333335, 265.26666666666665, 265.3,
        265.33333333333331, 265.36666666666667, 265.4, 265.43333333333334,
        265.46666666666664, 265.5, 265.5333333333333, 265.56666666666666, 265.6,
        265.63333333333333, 265.66666666666669, 265.7, 265.73333333333335,
        265.76666666666665, 265.8, 265.83333333333331, 265.86666666666667, 265.9,
        265.93333333333334, 265.96666666666664, 266.0, 266.0333333333333,
        266.06666666666666, 266.1, 266.13333333333333, 266.16666666666669, 266.2,
        266.23333333333335, 266.26666666666665, 266.3, 266.33333333333331,
        266.36666666666667, 266.4, 266.43333333333334, 266.46666666666664, 266.5,
        266.5333333333333, 266.56666666666666, 266.6, 266.63333333333333,
        266.66666666666669, 266.7, 266.73333333333335, 266.76666666666665, 266.8,
        266.83333333333331, 266.86666666666667, 266.9, 266.93333333333334,
        266.96666666666664, 267.0, 267.0333333333333, 267.06666666666666, 267.1,
        267.13333333333333, 267.16666666666669, 267.2, 267.23333333333335,
        267.26666666666665, 267.3, 267.33333333333331, 267.36666666666667, 267.4,
        267.43333333333334, 267.46666666666664, 267.5, 267.5333333333333,
        267.56666666666666, 267.6, 267.63333333333333, 267.66666666666669, 267.7,
        267.73333333333335, 267.76666666666665, 267.8, 267.83333333333331,
        267.86666666666667, 267.9, 267.93333333333334, 267.96666666666664, 268.0,
        268.0333333333333, 268.06666666666666, 268.1, 268.13333333333333,
        268.16666666666669, 268.2, 268.23333333333335, 268.26666666666665, 268.3,
        268.33333333333331, 268.36666666666667, 268.4, 268.43333333333334,
        268.46666666666664, 268.5, 268.5333333333333, 268.56666666666666, 268.6,
        268.63333333333333, 268.66666666666669, 268.7, 268.73333333333335,
        268.76666666666665, 268.8, 268.83333333333331, 268.86666666666667, 268.9,
        268.93333333333334, 268.96666666666664, 269.0, 269.0333333333333,
        269.06666666666666, 269.1, 269.13333333333333, 269.16666666666669, 269.2,
        269.23333333333335, 269.26666666666665, 269.3, 269.33333333333331,
        269.36666666666667, 269.4, 269.43333333333334, 269.46666666666664, 269.5,
        269.5333333333333, 269.56666666666666, 269.6, 269.63333333333333,
        269.66666666666669, 269.7, 269.73333333333335, 269.76666666666665, 269.8,
        269.83333333333331, 269.86666666666667, 269.9, 269.93333333333334,
        269.96666666666664, 270.0, 270.0333333333333, 270.06666666666666, 270.1,
        270.13333333333333, 270.16666666666669, 270.2, 270.23333333333335,
        270.26666666666665, 270.3, 270.33333333333331, 270.36666666666667, 270.4,
        270.43333333333334, 270.46666666666664, 270.5, 270.5333333333333,
        270.56666666666666, 270.6, 270.63333333333333, 270.66666666666669, 270.7,
        270.73333333333335, 270.76666666666665, 270.8, 270.83333333333331,
        270.86666666666667, 270.9, 270.93333333333334, 270.96666666666664, 271.0,
        271.0333333333333, 271.06666666666666, 271.1, 271.13333333333333,
        271.16666666666669, 271.2, 271.23333333333335, 271.26666666666665, 271.3,
        271.33333333333331, 271.36666666666667, 271.4, 271.43333333333334,
        271.46666666666664, 271.5, 271.5333333333333, 271.56666666666666, 271.6,
        271.63333333333333, 271.66666666666669, 271.7, 271.73333333333335,
        271.76666666666665, 271.8, 271.83333333333331, 271.86666666666667, 271.9,
        271.93333333333334, 271.96666666666664, 272.0, 272.0333333333333,
        272.06666666666666, 272.1, 272.13333333333333, 272.16666666666669, 272.2,
        272.23333333333335, 272.26666666666665, 272.3, 272.33333333333331,
        272.36666666666667, 272.4, 272.43333333333334, 272.46666666666664, 272.5,
        272.5333333333333, 272.56666666666666, 272.6, 272.63333333333333,
        272.66666666666669, 272.7, 272.73333333333335, 272.76666666666665, 272.8,
        272.83333333333331, 272.86666666666667, 272.9, 272.93333333333334,
        272.96666666666664, 273.0, 273.0333333333333, 273.06666666666666, 273.1,
        273.13333333333333, 273.16666666666669, 273.2, 273.23333333333335,
        273.26666666666665, 273.3, 273.33333333333331, 273.36666666666667, 273.4,
        273.43333333333334, 273.46666666666664, 273.5, 273.5333333333333,
        273.56666666666666, 273.6, 273.63333333333333, 273.66666666666669, 273.7,
        273.73333333333335, 273.76666666666665, 273.8, 273.83333333333331,
        273.86666666666667, 273.9, 273.93333333333334, 273.96666666666664, 274.0,
        274.0333333333333, 274.06666666666666, 274.1, 274.13333333333333,
        274.16666666666669, 274.2, 274.23333333333335, 274.26666666666665, 274.3,
        274.33333333333331, 274.36666666666667, 274.4, 274.43333333333334,
        274.46666666666664, 274.5, 274.5333333333333, 274.56666666666666, 274.6,
        274.63333333333333, 274.66666666666669, 274.7, 274.73333333333335,
        274.76666666666665, 274.8, 274.83333333333331, 274.86666666666667, 274.9,
        274.93333333333334, 274.96666666666664, 275.0, 275.0333333333333,
        275.06666666666666, 275.1, 275.13333333333333, 275.16666666666669, 275.2,
        275.23333333333335, 275.26666666666665, 275.3, 275.33333333333331,
        275.36666666666667, 275.4, 275.43333333333334, 275.46666666666664, 275.5,
        275.5333333333333, 275.56666666666666, 275.6, 275.63333333333333,
        275.66666666666669, 275.7, 275.73333333333335, 275.76666666666665, 275.8,
        275.83333333333331, 275.86666666666667, 275.9, 275.93333333333334,
        275.96666666666664, 276.0, 276.0333333333333, 276.06666666666666, 276.1,
        276.13333333333333, 276.16666666666669, 276.2, 276.23333333333335,
        276.26666666666665, 276.3, 276.33333333333331, 276.36666666666667, 276.4,
        276.43333333333334, 276.46666666666664, 276.5, 276.5333333333333,
        276.56666666666666, 276.6, 276.63333333333333, 276.66666666666669, 276.7,
        276.73333333333335, 276.76666666666665, 276.8, 276.83333333333331,
        276.86666666666667, 276.9, 276.93333333333334, 276.96666666666664, 277.0,
        277.0333333333333, 277.06666666666666, 277.1, 277.13333333333333,
        277.16666666666669, 277.2, 277.23333333333335, 277.26666666666665, 277.3,
        277.33333333333331, 277.36666666666667, 277.4, 277.43333333333334,
        277.46666666666664, 277.5, 277.5333333333333, 277.56666666666666, 277.6,
        277.63333333333333, 277.66666666666669, 277.7, 277.73333333333335,
        277.76666666666665, 277.8, 277.83333333333331, 277.86666666666667, 277.9,
        277.93333333333334, 277.96666666666664, 278.0, 278.0333333333333,
        278.06666666666666, 278.1, 278.13333333333333, 278.16666666666669, 278.2,
        278.23333333333335, 278.26666666666665, 278.3, 278.33333333333331,
        278.36666666666667, 278.4, 278.43333333333334, 278.46666666666664, 278.5,
        278.5333333333333, 278.56666666666666, 278.6, 278.63333333333333,
        278.66666666666669, 278.7, 278.73333333333335, 278.76666666666665, 278.8,
        278.83333333333331, 278.86666666666667, 278.9, 278.93333333333334,
        278.96666666666664, 279.0, 279.0333333333333, 279.06666666666666, 279.1,
        279.13333333333333, 279.16666666666669, 279.2, 279.23333333333335,
        279.26666666666665, 279.3, 279.33333333333331, 279.36666666666667, 279.4,
        279.43333333333334, 279.46666666666664, 279.5, 279.5333333333333,
        279.56666666666666, 279.6, 279.63333333333333, 279.66666666666669, 279.7,
        279.73333333333335, 279.76666666666665, 279.8, 279.83333333333331,
        279.86666666666667, 279.9, 279.93333333333334, 279.96666666666664, 280.0,
        280.0333333333333, 280.06666666666666, 280.1, 280.13333333333333,
        280.16666666666669, 280.2, 280.23333333333335, 280.26666666666665, 280.3,
        280.33333333333331, 280.36666666666667, 280.4, 280.43333333333334,
        280.46666666666664, 280.5, 280.5333333333333, 280.56666666666666, 280.6,
        280.63333333333333, 280.66666666666669, 280.7, 280.73333333333335,
        280.76666666666665, 280.8, 280.83333333333331, 280.86666666666667, 280.9,
        280.93333333333334, 280.96666666666664, 281.0, 281.0333333333333,
        281.06666666666666, 281.1, 281.13333333333333, 281.16666666666669, 281.2,
        281.23333333333335, 281.26666666666665, 281.3, 281.33333333333331,
        281.36666666666667, 281.4, 281.43333333333334, 281.46666666666664, 281.5,
        281.5333333333333, 281.56666666666666, 281.6, 281.63333333333333,
        281.66666666666669, 281.7, 281.73333333333335, 281.76666666666665, 281.8,
        281.83333333333331, 281.86666666666667, 281.9, 281.93333333333334,
        281.96666666666664, 282.0, 282.0333333333333, 282.06666666666666, 282.1,
        282.13333333333333, 282.16666666666669, 282.2, 282.23333333333335,
        282.26666666666665, 282.3, 282.33333333333331, 282.36666666666667, 282.4,
        282.43333333333334, 282.46666666666664, 282.5, 282.5333333333333,
        282.56666666666666, 282.6, 282.63333333333333, 282.66666666666669, 282.7,
        282.73333333333335, 282.76666666666665, 282.8, 282.83333333333331,
        282.86666666666667, 282.9, 282.93333333333334, 282.96666666666664, 283.0,
        283.0333333333333, 283.06666666666666, 283.1, 283.13333333333333,
        283.16666666666669, 283.2, 283.23333333333335, 283.26666666666665, 283.3,
        283.33333333333331, 283.36666666666667, 283.4, 283.43333333333334,
        283.46666666666664, 283.5, 283.5333333333333, 283.56666666666666, 283.6,
        283.63333333333333, 283.66666666666669, 283.7, 283.73333333333335,
        283.76666666666665, 283.8, 283.83333333333331, 283.86666666666667, 283.9,
        283.93333333333334, 283.96666666666664, 284.0, 284.0333333333333,
        284.06666666666666, 284.1, 284.13333333333333, 284.16666666666669, 284.2,
        284.23333333333335, 284.26666666666665, 284.3, 284.33333333333331,
        284.36666666666667, 284.4, 284.43333333333334, 284.46666666666664, 284.5,
        284.5333333333333, 284.56666666666666, 284.6, 284.63333333333333,
        284.66666666666669, 284.7, 284.73333333333335, 284.76666666666665, 284.8,
        284.83333333333331, 284.86666666666667, 284.9, 284.93333333333334,
        284.96666666666664, 285.0, 285.0333333333333, 285.06666666666666, 285.1,
        285.13333333333333, 285.16666666666669, 285.2, 285.23333333333335,
        285.26666666666665, 285.3, 285.33333333333331, 285.36666666666667, 285.4,
        285.43333333333334, 285.46666666666664, 285.5, 285.5333333333333,
        285.56666666666666, 285.6, 285.63333333333333, 285.66666666666669, 285.7,
        285.73333333333335, 285.76666666666665, 285.8, 285.83333333333331,
        285.86666666666667, 285.9, 285.93333333333334, 285.96666666666664, 286.0,
        286.0333333333333, 286.06666666666666, 286.1, 286.13333333333333,
        286.16666666666669, 286.2, 286.23333333333335, 286.26666666666665, 286.3,
        286.33333333333331, 286.36666666666667, 286.4, 286.43333333333334,
        286.46666666666664, 286.5, 286.5333333333333, 286.56666666666666, 286.6,
        286.63333333333333, 286.66666666666669, 286.7, 286.73333333333335,
        286.76666666666665, 286.8, 286.83333333333331, 286.86666666666667, 286.9,
        286.93333333333334, 286.96666666666664, 287.0, 287.0333333333333,
        287.06666666666666, 287.1, 287.13333333333333, 287.16666666666669, 287.2,
        287.23333333333335, 287.26666666666665, 287.3, 287.33333333333331,
        287.36666666666667, 287.4, 287.43333333333334, 287.46666666666664, 287.5,
        287.5333333333333, 287.56666666666666, 287.6, 287.63333333333333,
        287.66666666666669, 287.7, 287.73333333333335, 287.76666666666665, 287.8,
        287.83333333333331, 287.86666666666667, 287.9, 287.93333333333334,
        287.96666666666664, 288.0, 288.0333333333333, 288.06666666666666, 288.1,
        288.13333333333333, 288.16666666666669, 288.2, 288.23333333333335,
        288.26666666666665, 288.3, 288.33333333333331, 288.36666666666667, 288.4,
        288.43333333333334, 288.46666666666664, 288.5, 288.5333333333333,
        288.56666666666666, 288.6, 288.63333333333333, 288.66666666666669, 288.7,
        288.73333333333335, 288.76666666666665, 288.8, 288.83333333333331,
        288.86666666666667, 288.9, 288.93333333333334, 288.96666666666664, 289.0,
        289.0333333333333, 289.06666666666666, 289.1, 289.13333333333333,
        289.16666666666669, 289.2, 289.23333333333335, 289.26666666666665, 289.3,
        289.33333333333331, 289.36666666666667, 289.4, 289.43333333333334,
        289.46666666666664, 289.5, 289.5333333333333, 289.56666666666666, 289.6,
        289.63333333333333, 289.66666666666669, 289.7, 289.73333333333335,
        289.76666666666665, 289.8, 289.83333333333331, 289.86666666666667, 289.9,
        289.93333333333334, 289.96666666666664, 290.0, 290.0333333333333,
        290.06666666666666, 290.1, 290.13333333333333, 290.16666666666669, 290.2,
        290.23333333333335, 290.26666666666665, 290.3, 290.33333333333331,
        290.36666666666667, 290.4, 290.43333333333334, 290.46666666666664, 290.5,
        290.5333333333333, 290.56666666666666, 290.6, 290.63333333333333,
        290.66666666666669, 290.7, 290.73333333333335, 290.76666666666665, 290.8,
        290.83333333333331, 290.86666666666667, 290.9, 290.93333333333334,
        290.96666666666664, 291.0, 291.0333333333333, 291.06666666666666, 291.1,
        291.13333333333333, 291.16666666666669, 291.2, 291.23333333333335,
        291.26666666666665, 291.3, 291.33333333333331, 291.36666666666667, 291.4,
        291.43333333333334, 291.46666666666664, 291.5, 291.5333333333333,
        291.56666666666666, 291.6, 291.63333333333333, 291.66666666666669, 291.7,
        291.73333333333335, 291.76666666666665, 291.8, 291.83333333333331,
        291.86666666666667, 291.9, 291.93333333333334, 291.96666666666664, 292.0,
        292.0333333333333, 292.06666666666666, 292.1, 292.13333333333333,
        292.16666666666669, 292.2, 292.23333333333335, 292.26666666666665, 292.3,
        292.33333333333331, 292.36666666666667, 292.4, 292.43333333333334,
        292.46666666666664, 292.5, 292.5333333333333, 292.56666666666666, 292.6,
        292.63333333333333, 292.66666666666669, 292.7, 292.73333333333335,
        292.76666666666665, 292.8, 292.83333333333331, 292.86666666666667, 292.9,
        292.93333333333334, 292.96666666666664, 293.0, 293.0333333333333,
        293.06666666666666, 293.1, 293.13333333333333, 293.16666666666669, 293.2,
        293.23333333333335, 293.26666666666665, 293.3, 293.33333333333331,
        293.36666666666667, 293.4, 293.43333333333334, 293.46666666666664, 293.5,
        293.5333333333333, 293.56666666666666, 293.6, 293.63333333333333,
        293.66666666666669, 293.7, 293.73333333333335, 293.76666666666665, 293.8,
        293.83333333333331, 293.86666666666667, 293.9, 293.93333333333334,
        293.96666666666664, 294.0, 294.0333333333333, 294.06666666666666, 294.1,
        294.13333333333333, 294.16666666666669, 294.2, 294.23333333333335,
        294.26666666666665, 294.3, 294.33333333333331, 294.36666666666667, 294.4,
        294.43333333333334, 294.46666666666664, 294.5, 294.5333333333333,
        294.56666666666666, 294.6, 294.63333333333333, 294.66666666666669, 294.7,
        294.73333333333335, 294.76666666666665, 294.8, 294.83333333333331,
        294.86666666666667, 294.9, 294.93333333333334, 294.96666666666664, 295.0,
        295.0333333333333, 295.06666666666666, 295.1, 295.13333333333333,
        295.16666666666669, 295.2, 295.23333333333335, 295.26666666666665, 295.3,
        295.33333333333331, 295.36666666666667, 295.4, 295.43333333333334,
        295.46666666666664, 295.5, 295.5333333333333, 295.56666666666666, 295.6,
        295.63333333333333, 295.66666666666669, 295.7, 295.73333333333335,
        295.76666666666665, 295.8, 295.83333333333331, 295.86666666666667, 295.9,
        295.93333333333334, 295.96666666666664, 296.0, 296.0333333333333,
        296.06666666666666, 296.1, 296.13333333333333, 296.16666666666669, 296.2,
        296.23333333333335, 296.26666666666665, 296.3, 296.33333333333331,
        296.36666666666667, 296.4, 296.43333333333334, 296.46666666666664, 296.5,
        296.5333333333333, 296.56666666666666, 296.6, 296.63333333333333,
        296.66666666666669, 296.7, 296.73333333333335, 296.76666666666665, 296.8,
        296.83333333333331, 296.86666666666667, 296.9, 296.93333333333334,
        296.96666666666664, 297.0, 297.0333333333333, 297.06666666666666, 297.1,
        297.13333333333333, 297.16666666666669, 297.2, 297.23333333333335,
        297.26666666666665, 297.3, 297.33333333333331, 297.36666666666667, 297.4,
        297.43333333333334, 297.46666666666664, 297.5, 297.5333333333333,
        297.56666666666666, 297.6, 297.63333333333333, 297.66666666666669, 297.7,
        297.73333333333335, 297.76666666666665, 297.8, 297.83333333333331,
        297.86666666666667, 297.9, 297.93333333333334, 297.96666666666664, 298.0,
        298.0333333333333, 298.06666666666666, 298.1, 298.13333333333333,
        298.16666666666669, 298.2, 298.23333333333335, 298.26666666666665, 298.3,
        298.33333333333331, 298.36666666666667, 298.4, 298.43333333333334,
        298.46666666666664, 298.5, 298.5333333333333, 298.56666666666666, 298.6,
        298.63333333333333, 298.66666666666669, 298.7, 298.73333333333335,
        298.76666666666665, 298.8, 298.83333333333331, 298.86666666666667, 298.9,
        298.93333333333334, 298.96666666666664, 299.0, 299.0333333333333,
        299.06666666666666, 299.1, 299.13333333333333, 299.16666666666669, 299.2,
        299.23333333333335, 299.26666666666665, 299.3, 299.33333333333331,
        299.36666666666667, 299.4, 299.43333333333334, 299.46666666666664, 299.5,
        299.5333333333333, 299.56666666666666, 299.6, 299.63333333333333,
        299.66666666666669, 299.7, 299.73333333333335, 299.76666666666665, 299.8,
        299.83333333333331, 299.86666666666667, 299.9, 299.93333333333334,
        299.96666666666664, 300.0, 300.0333333333333, 300.06666666666666, 300.1,
        300.13333333333333, 300.16666666666669, 300.2, 300.23333333333335,
        300.26666666666665, 300.3, 300.33333333333331, 300.36666666666667, 300.4,
        300.43333333333334, 300.46666666666664, 300.5, 300.5333333333333,
        300.56666666666666, 300.6, 300.63333333333333, 300.66666666666669, 300.7,
        300.73333333333335, 300.76666666666665, 300.8, 300.83333333333331,
        300.86666666666667, 300.9, 300.93333333333334, 300.96666666666664, 301.0,
        301.0333333333333, 301.06666666666666, 301.1, 301.13333333333333,
        301.16666666666669, 301.2, 301.23333333333335, 301.26666666666665, 301.3,
        301.33333333333331, 301.36666666666667, 301.4, 301.43333333333334,
        301.46666666666664, 301.5, 301.5333333333333, 301.56666666666666, 301.6,
        301.63333333333333, 301.66666666666669, 301.7, 301.73333333333335,
        301.76666666666665, 301.8, 301.83333333333331, 301.86666666666667, 301.9,
        301.93333333333334, 301.96666666666664, 302.0, 302.0333333333333,
        302.06666666666666, 302.1, 302.13333333333333, 302.16666666666669, 302.2,
        302.23333333333335, 302.26666666666665, 302.3, 302.33333333333331,
        302.36666666666667, 302.4, 302.43333333333334, 302.46666666666664, 302.5,
        302.5333333333333, 302.56666666666666, 302.6, 302.63333333333333,
        302.66666666666669, 302.7, 302.73333333333335, 302.76666666666665, 302.8,
        302.83333333333331, 302.86666666666667, 302.9, 302.93333333333334,
        302.96666666666664, 303.0, 303.0333333333333, 303.06666666666666, 303.1,
        303.13333333333333, 303.16666666666669, 303.2, 303.23333333333335,
        303.26666666666665, 303.3, 303.33333333333331, 303.36666666666667, 303.4,
        303.43333333333334, 303.46666666666664, 303.5, 303.5333333333333,
        303.56666666666666, 303.6, 303.63333333333333, 303.66666666666669, 303.7,
        303.73333333333335, 303.76666666666665, 303.8, 303.83333333333331,
        303.86666666666667, 303.9, 303.93333333333334, 303.96666666666664, 304.0,
        304.0333333333333, 304.06666666666666, 304.1, 304.13333333333333,
        304.16666666666669, 304.2, 304.23333333333335, 304.26666666666665, 304.3,
        304.33333333333331, 304.36666666666667, 304.4, 304.43333333333334,
        304.46666666666664, 304.5, 304.5333333333333, 304.56666666666666, 304.6,
        304.63333333333333, 304.66666666666669, 304.7, 304.73333333333335,
        304.76666666666665, 304.8, 304.83333333333331, 304.86666666666667, 304.9,
        304.93333333333334, 304.96666666666664, 305.0, 305.0333333333333,
        305.06666666666666, 305.1, 305.13333333333333, 305.16666666666669, 305.2,
        305.23333333333335, 305.26666666666665, 305.3, 305.33333333333331,
        305.36666666666667, 305.4, 305.43333333333334, 305.46666666666664, 305.5,
        305.5333333333333, 305.56666666666666, 305.6, 305.63333333333333,
        305.66666666666669, 305.7, 305.73333333333335, 305.76666666666665, 305.8,
        305.83333333333331, 305.86666666666667, 305.9, 305.93333333333334,
        305.96666666666664, 306.0, 306.0333333333333, 306.06666666666666, 306.1,
        306.13333333333333, 306.16666666666669, 306.2, 306.23333333333335,
        306.26666666666665, 306.3, 306.33333333333331, 306.36666666666667, 306.4,
        306.43333333333334, 306.46666666666664, 306.5, 306.5333333333333,
        306.56666666666666, 306.6, 306.63333333333333, 306.66666666666669, 306.7,
        306.73333333333335, 306.76666666666665, 306.8, 306.83333333333331,
        306.86666666666667, 306.9, 306.93333333333334, 306.96666666666664, 307.0,
        307.0333333333333, 307.06666666666666, 307.1, 307.13333333333333,
        307.16666666666669, 307.2, 307.23333333333335, 307.26666666666665, 307.3,
        307.33333333333331, 307.36666666666667, 307.4, 307.43333333333334,
        307.46666666666664, 307.5, 307.5333333333333, 307.56666666666666, 307.6,
        307.63333333333333, 307.66666666666669, 307.7, 307.73333333333335,
        307.76666666666665, 307.8, 307.83333333333331, 307.86666666666667, 307.9,
        307.93333333333334, 307.96666666666664, 308.0, 308.0333333333333,
        308.06666666666666, 308.1, 308.13333333333333, 308.16666666666669, 308.2,
        308.23333333333335, 308.26666666666665, 308.3, 308.33333333333331,
        308.36666666666667, 308.4, 308.43333333333334, 308.46666666666664, 308.5,
        308.5333333333333, 308.56666666666666, 308.6, 308.63333333333333,
        308.66666666666669, 308.7, 308.73333333333335, 308.76666666666665, 308.8,
        308.83333333333331, 308.86666666666667, 308.9, 308.93333333333334,
        308.96666666666664, 309.0, 309.0333333333333, 309.06666666666666, 309.1,
        309.13333333333333, 309.16666666666669, 309.2, 309.23333333333335,
        309.26666666666665, 309.3, 309.33333333333331, 309.36666666666667, 309.4,
        309.43333333333334, 309.46666666666664, 309.5, 309.5333333333333,
        309.56666666666666, 309.6, 309.63333333333333, 309.66666666666669, 309.7,
        309.73333333333335, 309.76666666666665, 309.8, 309.83333333333331,
        309.86666666666667, 309.9, 309.93333333333334, 309.96666666666664, 310.0,
        310.0333333333333, 310.06666666666666, 310.1, 310.13333333333333,
        310.16666666666669, 310.2, 310.23333333333335, 310.26666666666665, 310.3,
        310.33333333333331, 310.36666666666667, 310.4, 310.43333333333334,
        310.46666666666664, 310.5, 310.5333333333333, 310.56666666666666, 310.6,
        310.63333333333333, 310.66666666666669, 310.7, 310.73333333333335,
        310.76666666666665, 310.8, 310.83333333333331, 310.86666666666667, 310.9,
        310.93333333333334, 310.96666666666664, 311.0, 311.0333333333333,
        311.06666666666666, 311.1, 311.13333333333333, 311.16666666666669, 311.2,
        311.23333333333335, 311.26666666666665, 311.3, 311.33333333333331,
        311.36666666666667, 311.4, 311.43333333333334, 311.46666666666664, 311.5,
        311.5333333333333, 311.56666666666666, 311.6, 311.63333333333333,
        311.66666666666669, 311.7, 311.73333333333335, 311.76666666666665, 311.8,
        311.83333333333331, 311.86666666666667, 311.9, 311.93333333333334,
        311.96666666666664, 312.0, 312.0333333333333, 312.06666666666666, 312.1,
        312.13333333333333, 312.16666666666669, 312.2, 312.23333333333335,
        312.26666666666665, 312.3, 312.33333333333331, 312.36666666666667, 312.4,
        312.43333333333334, 312.46666666666664, 312.5, 312.5333333333333,
        312.56666666666666, 312.6, 312.63333333333333, 312.66666666666669, 312.7,
        312.73333333333335, 312.76666666666665, 312.8, 312.83333333333331,
        312.86666666666667, 312.9, 312.93333333333334, 312.96666666666664, 313.0,
        313.0333333333333, 313.06666666666666, 313.1, 313.13333333333333,
        313.16666666666669, 313.2, 313.23333333333335, 313.26666666666665, 313.3,
        313.33333333333331, 313.36666666666667, 313.4, 313.43333333333334,
        313.46666666666664, 313.5, 313.5333333333333, 313.56666666666666, 313.6,
        313.63333333333333, 313.66666666666669, 313.7, 313.73333333333335,
        313.76666666666665, 313.8, 313.83333333333331, 313.86666666666667, 313.9,
        313.93333333333334, 313.96666666666664, 314.0, 314.0333333333333,
        314.06666666666666, 314.1, 314.13333333333333, 314.16666666666669, 314.2,
        314.23333333333335, 314.26666666666665, 314.3, 314.33333333333331,
        314.36666666666667, 314.4, 314.43333333333334, 314.46666666666664, 314.5,
        314.5333333333333, 314.56666666666666, 314.6, 314.63333333333333,
        314.66666666666669, 314.7, 314.73333333333335, 314.76666666666665, 314.8,
        314.83333333333331, 314.86666666666667, 314.9, 314.93333333333334,
        314.96666666666664, 315.0, 315.0333333333333, 315.06666666666666, 315.1,
        315.13333333333333, 315.16666666666669, 315.2, 315.23333333333335,
        315.26666666666665, 315.3, 315.33333333333331, 315.36666666666667, 315.4,
        315.43333333333334, 315.46666666666664, 315.5, 315.5333333333333,
        315.56666666666666, 315.6, 315.63333333333333, 315.66666666666669, 315.7,
        315.73333333333335, 315.76666666666665, 315.8, 315.83333333333331,
        315.86666666666667, 315.9, 315.93333333333334, 315.96666666666664, 316.0,
        316.0333333333333, 316.06666666666666, 316.1, 316.13333333333333,
        316.16666666666669, 316.2, 316.23333333333335, 316.26666666666665, 316.3,
        316.33333333333331, 316.36666666666667, 316.4, 316.43333333333334,
        316.46666666666664, 316.5, 316.5333333333333, 316.56666666666666, 316.6,
        316.63333333333333, 316.66666666666669, 316.7, 316.73333333333335,
        316.76666666666665, 316.8, 316.83333333333331, 316.86666666666667, 316.9,
        316.93333333333334, 316.96666666666664, 317.0, 317.0333333333333,
        317.06666666666666, 317.1, 317.13333333333333, 317.16666666666669, 317.2,
        317.23333333333335, 317.26666666666665, 317.3, 317.33333333333331,
        317.36666666666667, 317.4, 317.43333333333334, 317.46666666666664, 317.5,
        317.5333333333333, 317.56666666666666, 317.6, 317.63333333333333,
        317.66666666666669, 317.7, 317.73333333333335, 317.76666666666665, 317.8,
        317.83333333333331, 317.86666666666667, 317.9, 317.93333333333334,
        317.96666666666664, 318.0, 318.0333333333333, 318.06666666666666, 318.1,
        318.13333333333333, 318.16666666666669, 318.2, 318.23333333333335,
        318.26666666666665, 318.3, 318.33333333333331, 318.36666666666667, 318.4,
        318.43333333333334, 318.46666666666664, 318.5, 318.5333333333333,
        318.56666666666666, 318.6, 318.63333333333333, 318.66666666666669, 318.7,
        318.73333333333335, 318.76666666666665, 318.8, 318.83333333333331,
        318.86666666666667, 318.9, 318.93333333333334, 318.96666666666664, 319.0,
        319.0333333333333, 319.06666666666666, 319.1, 319.13333333333333,
        319.16666666666669, 319.2, 319.23333333333335, 319.26666666666665, 319.3,
        319.33333333333331, 319.36666666666667, 319.4, 319.43333333333334,
        319.46666666666664, 319.5, 319.5333333333333, 319.56666666666666, 319.6,
        319.63333333333333, 319.66666666666669, 319.7, 319.73333333333335,
        319.76666666666665, 319.8, 319.83333333333331, 319.86666666666667, 319.9,
        319.93333333333334, 319.96666666666664, 320.0, 320.0333333333333,
        320.06666666666666, 320.1, 320.13333333333333, 320.16666666666669, 320.2,
        320.23333333333335, 320.26666666666665, 320.3, 320.33333333333331,
        320.36666666666667, 320.4, 320.43333333333334, 320.46666666666664, 320.5,
        320.5333333333333, 320.56666666666666, 320.6, 320.63333333333333,
        320.66666666666669, 320.7, 320.73333333333335, 320.76666666666665, 320.8,
        320.83333333333331, 320.86666666666667, 320.9, 320.93333333333334,
        320.96666666666664, 321.0, 321.0333333333333, 321.06666666666666, 321.1,
        321.13333333333333, 321.16666666666669, 321.2, 321.23333333333335,
        321.26666666666665, 321.3, 321.33333333333331, 321.36666666666667, 321.4,
        321.43333333333334, 321.46666666666664, 321.5, 321.5333333333333,
        321.56666666666666, 321.6, 321.63333333333333, 321.66666666666669, 321.7,
        321.73333333333335, 321.76666666666665, 321.8, 321.83333333333331,
        321.86666666666667, 321.9, 321.93333333333334, 321.96666666666664, 322.0,
        322.0333333333333, 322.06666666666666, 322.1, 322.13333333333333,
        322.16666666666669, 322.2, 322.23333333333335, 322.26666666666665, 322.3,
        322.33333333333331, 322.36666666666667, 322.4, 322.43333333333334,
        322.46666666666664, 322.5, 322.5333333333333, 322.56666666666666, 322.6,
        322.63333333333333, 322.66666666666669, 322.7, 322.73333333333335,
        322.76666666666665, 322.8, 322.83333333333331, 322.86666666666667, 322.9,
        322.93333333333334, 322.96666666666664, 323.0, 323.0333333333333,
        323.06666666666666, 323.1, 323.13333333333333, 323.16666666666669, 323.2,
        323.23333333333335, 323.26666666666665, 323.3, 323.33333333333331,
        323.36666666666667, 323.4, 323.43333333333334, 323.46666666666664, 323.5,
        323.5333333333333, 323.56666666666666, 323.6, 323.63333333333333,
        323.66666666666669, 323.7, 323.73333333333335, 323.76666666666665, 323.8,
        323.83333333333331, 323.86666666666667, 323.9, 323.93333333333334,
        323.96666666666664, 324.0, 324.0333333333333, 324.06666666666666, 324.1,
        324.13333333333333, 324.16666666666669, 324.2, 324.23333333333335,
        324.26666666666665, 324.3, 324.33333333333331, 324.36666666666667, 324.4,
        324.43333333333334, 324.46666666666664, 324.5, 324.5333333333333,
        324.56666666666666, 324.6, 324.63333333333333, 324.66666666666669, 324.7,
        324.73333333333335, 324.76666666666665, 324.8, 324.83333333333331,
        324.86666666666667, 324.9, 324.93333333333334, 324.96666666666664, 325.0,
        325.0333333333333, 325.06666666666666, 325.1, 325.13333333333333,
        325.16666666666669, 325.2, 325.23333333333335, 325.26666666666665, 325.3,
        325.33333333333331, 325.36666666666667, 325.4, 325.43333333333334,
        325.46666666666664, 325.5, 325.5333333333333, 325.56666666666666, 325.6,
        325.63333333333333, 325.66666666666669, 325.7, 325.73333333333335,
        325.76666666666665, 325.8, 325.83333333333331, 325.86666666666667, 325.9,
        325.93333333333334, 325.96666666666664, 326.0, 326.0333333333333,
        326.06666666666666, 326.1, 326.13333333333333, 326.16666666666669, 326.2,
        326.23333333333335, 326.26666666666665, 326.3, 326.33333333333331,
        326.36666666666667, 326.4, 326.43333333333334, 326.46666666666664, 326.5,
        326.5333333333333, 326.56666666666666, 326.6, 326.63333333333333,
        326.66666666666669, 326.7, 326.73333333333335, 326.76666666666665, 326.8,
        326.83333333333331, 326.86666666666667, 326.9, 326.93333333333334,
        326.96666666666664, 327.0, 327.0333333333333, 327.06666666666666, 327.1,
        327.13333333333333, 327.16666666666669, 327.2, 327.23333333333335,
        327.26666666666665, 327.3, 327.33333333333331, 327.36666666666667, 327.4,
        327.43333333333334, 327.46666666666664, 327.5, 327.5333333333333,
        327.56666666666666, 327.6, 327.63333333333333, 327.66666666666669, 327.7,
        327.73333333333335, 327.76666666666665, 327.8, 327.83333333333331,
        327.86666666666667, 327.9, 327.93333333333334, 327.96666666666664, 328.0,
        328.0333333333333, 328.06666666666666, 328.1, 328.13333333333333,
        328.16666666666669, 328.2, 328.23333333333335, 328.26666666666665, 328.3,
        328.33333333333331, 328.36666666666667, 328.4, 328.43333333333334,
        328.46666666666664, 328.5, 328.5333333333333, 328.56666666666666, 328.6,
        328.63333333333333, 328.66666666666669, 328.7, 328.73333333333335,
        328.76666666666665, 328.8, 328.83333333333331, 328.86666666666667, 328.9,
        328.93333333333334, 328.96666666666664, 329.0, 329.0333333333333,
        329.06666666666666, 329.1, 329.13333333333333, 329.16666666666669, 329.2,
        329.23333333333335, 329.26666666666665, 329.3, 329.33333333333331,
        329.36666666666667, 329.4, 329.43333333333334, 329.46666666666664, 329.5,
        329.5333333333333, 329.56666666666666, 329.6, 329.63333333333333,
        329.66666666666669, 329.7, 329.73333333333335, 329.76666666666665, 329.8,
        329.83333333333331, 329.86666666666667, 329.9, 329.93333333333334,
        329.96666666666664, 330.0, 330.0333333333333, 330.06666666666666, 330.1,
        330.13333333333333, 330.16666666666669, 330.2, 330.23333333333335,
        330.26666666666665, 330.3, 330.33333333333331, 330.36666666666667, 330.4,
        330.43333333333334, 330.46666666666664, 330.5, 330.5333333333333,
        330.56666666666666, 330.6, 330.63333333333333, 330.66666666666669, 330.7,
        330.73333333333335, 330.76666666666665, 330.8, 330.83333333333331,
        330.86666666666667, 330.9, 330.93333333333334, 330.96666666666664, 331.0,
        331.0333333333333, 331.06666666666666, 331.1, 331.13333333333333,
        331.16666666666669, 331.2, 331.23333333333335, 331.26666666666665, 331.3,
        331.33333333333331, 331.36666666666667, 331.4, 331.43333333333334,
        331.46666666666664, 331.5, 331.5333333333333, 331.56666666666666, 331.6,
        331.63333333333333, 331.66666666666669, 331.7, 331.73333333333335,
        331.76666666666665, 331.8, 331.83333333333331, 331.86666666666667, 331.9,
        331.93333333333334, 331.96666666666664, 332.0, 332.0333333333333,
        332.06666666666666, 332.1, 332.13333333333333, 332.16666666666669, 332.2,
        332.23333333333335, 332.26666666666665, 332.3, 332.33333333333331,
        332.36666666666667, 332.4, 332.43333333333334, 332.46666666666664, 332.5,
        332.5333333333333, 332.56666666666666, 332.6, 332.63333333333333,
        332.66666666666669, 332.7, 332.73333333333335, 332.76666666666665, 332.8,
        332.83333333333331, 332.86666666666667, 332.9, 332.93333333333334,
        332.96666666666664 } ;

      static real_T pDataValues0[] = { -0.012216242412705467,
        -0.012085013502739082, -0.011950937288450139, -0.011814075134701065,
        -0.011674488392667524, -0.011532238431811625, -0.011387386598062932,
        -0.011239994265395127, -0.011090122787678071, -0.010937833517754571,
        -0.010783187831933105, -0.01062624708364268, -0.010467072630060411,
        -0.010305725824473201, -0.010142268040187474, -0.00997676063290598,
        -0.0098092649618131168, -0.0096398423869281658, -0.00946855426539273,
        -0.0092954619591978371, -0.0091206268360544147, -0.0089441102405807219,
        -0.0087659735502931754, -0.0085862781118706231, -0.00840508529481809,
        -0.0082224564476440776, -0.0080384529376154052, -0.0078531361296363435,
        -0.0076665673772247124, -0.0074788080428511722, -0.0072899194764453014,
        -0.0070999630648493871, -0.0069090001325733423, -0.0067170920662640316,
        -0.00652430021975903, -0.0063306859477663142, -0.0061363106090084364,
        -0.0059412355708587947, -0.0057455221979907396, -0.0055492318387528883,
        -0.005352425846361088, -0.0051551656021686654, -0.0049575124659284615,
        -0.004759527767905781, -0.0045612728958488421, -0.0043628092038616667,
        -0.0041641980516438065, -0.0039655007958216979, -0.0037667787954553953,
        -0.0035680934153781165, -0.0033695060233363, -0.0031710779538407525,
        -0.00297287059355611, -0.0027749452940817325, -0.0025773634048320559,
        -0.0023801862981009951, -0.0021834753287386377, -0.0019872918709040737,
        -0.0017916972564435698, -0.0015967528666571706, -0.001402520058313427,
        -0.0012090601880565541, -0.0010164346191032804, -0.00082470470491813189,
        -0.00063393181667592492, -0.00044417730604706792,
        -0.00025550253186068484, -6.79688664284583E-5, 0.00011836235261242227,
        0.00030342974309149862, 0.00048717195205938549, 0.00066952762209027063,
        0.00085043539433726648, 0.0010298339093495193, 0.0012076617893086456,
        0.0013838577043578939, 0.0015583602771229723, 0.0017311081410653619,
        0.0019020399441060909, 0.0020710943377721916, 0.0022382099383477785,
        0.0024033254045185974, 0.0025663793701191694, 0.0027273104696590219,
        0.0028860573528888271, 0.0030425586523819006, 0.0031967530099308271,
        0.0033485790657650175, 0.0034979754659580911, 0.0036448808346634123,
        0.0037892338290035359, 0.0039309730829906195, 0.0040700372283986042,
        0.0042063649147649988, 0.0043398947808626024, 0.004470565460064079,
        0.0045983156055484659, 0.0047230838414513794, 0.0048448088175057176,
        0.0049634291762323821, 0.0050788835393333906, 0.0051911105734525551,
        0.0053000489020504846, 0.0054056371631894251, 0.00550781400189493,
        0.00560651806235768, 0.0057016879810056782, 0.0057932623884227041,
        0.005881179938906903, 0.0059653792663105689, 0.0060457990101703283,
        0.0061223778108576718, 0.0061950543038769125, 0.0062637671368470762,
        0.0063284549536390972, 0.006389056373610205, 0.00644551005411386,
        0.0064977546385769254, 0.0065457287514902891, 0.00658937104370594,
        0.0066286201484189291, 0.0066634147110100645, 0.0066936933725081161,
        0.0067193947681509561, 0.0067404575391094368, 0.0067568203257728516,
        0.0067684217702357564, 0.0067752005029042793, 0.0067770951769219566,
        0.0067740444268338485, 0.0067659868913949612, 0.0067528612104616817,
        0.0067346060255068854, 0.0067111599706848124, 0.0066824616967409173,
        0.0066484498410801724, 0.0066090630314974436, 0.0065642399203191624,
        0.0065139191437957123, 0.0064580393444124522, 0.0063965391568387995,
        0.0063293572247502755, 0.0062564321844505114, 0.006177702680343366,
        0.0060931073539549829, 0.0060025848359401852, 0.0059060737771154774,
        0.0058035128115107937, 0.0056948405701539966, 0.0055799957189606193,
        0.005458916873131177, 0.0053315426786859649, 0.0051978117833328813,
        0.0050576628142451211, 0.0049110344251879108, 0.0047578652505464447,
        0.00459809392511447, 0.0044316590956228353, 0.0042584993964212151,
        0.00407855346901799, 0.0038917599617605, 0.0036980575036338241,
        0.003497384736412798, 0.0032896803028492523, 0.0030748828384652427,
        0.0028529310165304007, 0.0026237633536396751, 0.0023873188304968035,
        0.0021435351712639737, 0.0018923535233739133, 0.0016337056089277267,
        0.0013675491459430634, 0.0010937749944339753, 0.00081271962383056289,
        0.00052582599003869941, 0.00023498268678423593, -5.7988585488997841E-5,
        -0.00035124012091714289, -0.00064293362479261819,
        -0.00093122738677564685, -0.001214280955323947, -0.0014902534211991276,
        -0.0017573040173068865, -0.00201359194672789, -0.0022572764202699485,
        -0.0024865166421150788, -0.0026994718159656463, -0.0028943011582072235,
        -0.0030691638653834842, -0.0032222191470055008, -0.0033516262101329635,
        -0.0034555442625894164, -0.0035321325144366152, -0.0035795501645985331,
        -0.003595956421715841, -0.0035795104966673936, -0.0035283715964418158,
        -0.0034406989260382568, -0.003314651685588635, -0.0031483890894889958,
        -0.0029400703450654935, -0.0026878546530007077, -0.0023899012288275488,
        -0.0020443692705463053, -0.00164941798637131, -0.0012032065875366995,
        -0.00070389428033834688, -0.00014964026364695129, 0.00046139624359398281,
        0.0011310560512001208, 0.0018611799304757487, 0.0026536086911654868,
        0.0035101831287143023, 0.0044327440192759313, 0.0054231321745845272,
        0.0064831883778993, 0.0076147534253756793, 0.0088196681027774872,
        0.010099773217433531, 0.011456909556937272, 0.012892917906981436,
        0.01440963906654597, 0.016008913835783217, 0.017692582994293059,
        0.019462487347556848, 0.021320467683345685, 0.023268364792876639,
        0.025308019475840292, 0.027441272517574222, 0.029669964713558841,
        0.031995936865438786, 0.034421029762104052, 0.036947084181928734,
        0.039575940940324156, 0.042309440821579759, 0.045149424615917873,
        0.048097733122211946, 0.05115620712610161, 0.054326687429799875,
        0.057611014820402345, 0.061011030094615629, 0.064528574045344081,
        0.068165487465403671, 0.071923611142104243, 0.07580478587808441,
        0.079810852466924112, 0.083943651686019555, 0.0882050243505352,
        0.092596811248258928, 0.097120853153071188, 0.1017789908782166,
        0.10657306522343685, 0.11150491695367788, 0.11657638688692672,
        0.12178931580960438, 0.12714554451259066, 0.13264691379035345,
        0.13829526443363044, 0.14409243724141921, 0.1500402730020772,
        0.15614061251696443, 0.16239529656118479, 0.16880616595061509,
        0.1753750614662033, 0.18210382389334095, 0.18899429404747595,
        0.19604831270628581, 0.20326772066070548, 0.21065435871467103,
        0.21821006765282724, 0.22593668827609018, 0.23383606136992169,
        0.2419100277308148, 0.25016042815659473, 0.25858910342966812,
        0.26719789435887348, 0.27598864171835391, 0.28496318632091611,
        0.29412336894967372, 0.30347103038927248, 0.3130080114537579,
        0.322736152924713, 0.33265729558648383, 0.34277328024826426,
        0.35308594769695517, 0.36359713872642052, 0.37430869412759327,
        0.38522245469593625, 0.39634026122661004, 0.40766395450461496,
        0.419195375338416, 0.43093636450002337, 0.44288876280697215,
        0.45505441103738276, 0.46743514998055541, 0.48003282043911488,
        0.49284926321222672, 0.50588631907963766, 0.5191458288348213,
        0.5326296332831878, 0.54633957320533211, 0.560277489405989,
        0.57444522267285891, 0.58884461379604, 0.60347750356856122,
        0.618345732792582, 0.63345114225462262, 0.648795572748029,
        0.66438086506878935, 0.68020886000467662, 0.69628139835670166,
        0.7126003209167413, 0.729167468468968, 0.74598468181648336,
        0.76305380174997151, 0.780376669066441, 0.79795512454540274,
        0.81579100899363433, 0.83388616320928166, 0.85224242796288463,
        0.8708616440684861, 0.889745652314114, 0.90889629348825263,
        0.928315408396977, 0.9480048377947361, 0.96796642258042032,
        0.98820200326569818, 1.0087134214148499, 1.0295025156857067,
        1.0505711327899228, 1.0719210971433528, 1.0935542948504244,
        1.1154724417525093, 1.1376776913626194, 1.1601692811219917,
        1.1829392073009473, 1.2059765500929196, 1.2292708273789772,
        1.2528113867639457, 1.2765876375205449, 1.3005889666715691,
        1.3248047692764908, 1.3492244374530187, 1.3738373644246464,
        1.3986329430000892, 1.4236005661304016, 1.4487296267241998,
        1.4740095176974195, 1.4994296319787868, 1.5249793624705765,
        1.5506481021046572, 1.5764252437897348, 1.6023001804506458,
        1.6282623050044807, 1.6543010103546891, 1.6804056894371708,
        1.7065657351691024, 1.7327705404533471, 1.7590094982212061,
        1.7852720013847077, 1.8115474428646356, 1.837825215573214,
        1.8640947124378007, 1.8903453263685568, 1.9165664502844377,
        1.9427474771093736, 1.9688777997484406, 1.9949468111323456,
        2.02094390417393, 2.046858471787746, 2.0726799068936566,
        2.0983976024208215, 2.1240009512619626, 2.1494793463605326,
        2.17482218062364, 2.2000188469632218, 2.2250587383072231,
        2.2499312475695712, 2.274625767665527, 2.2991316915122169,
        2.3234384120388785, 2.3475353221493616, 2.3714118147679146,
        2.395057282812119, 2.4184611191963294, 2.4416127168472355,
        2.4645014686743583, 2.4871167675934291, 2.5094480065342935,
        2.5314845783985529, 2.5532158761192663, 2.5746312926060626,
        2.5957202207770456, 2.6164720535554857, 2.6368761838498593,
        2.6569220045883033, 2.6765989086769717, 2.6958962890482945,
        2.7148035386090981, 2.7333100502772378, 2.7514052169741534,
        2.7690784316257853, 2.7863190871298031, 2.80311657642149,
        2.81946029240543, 2.8353396280409022, 2.8507439760983924,
        2.86566272993171, 2.8800852812461133, 2.8940010263077651,
        2.9073993487672123, 2.9202696672132946, 2.9326013034818281,
        2.9443838465470114, 2.9556061986405959, 2.9662618376889061,
        2.9763556040263524, 2.985896913723824, 2.9948944960548585,
        3.0033573474900983, 3.011294367683341, 3.0187144912778514,
        3.0256266402776539, 3.0320397412550832, 3.0379627191191112,
        3.0434044993936293, 3.0483740073725478, 3.0528801684425546,
        3.0569319079534965, 3.0605381512433349, 3.0637078236935,
        3.0664498506569982, 3.0687731574726778, 3.07068666951966,
        3.0721993121439373, 3.0733200107012886, 3.07405769055751,
        3.0744212770646522, 3.0744196955720442, 3.0740618714509371,
        3.0733567300513891, 3.0723131967168524, 3.0709401968368564,
        3.0692466557405371, 3.0672414987888956, 3.0649336513514411,
        3.0623320387780884, 3.0594455864206012, 3.05628321964517,
        3.052853863800113, 3.0491664442458983, 3.0452298863457123,
        3.0410531154499871, 3.0366450569171128, 3.0320146360977347,
        3.0271707783670103, 3.0221224090586571, 3.01687845354597,
        3.011447837187553, 3.0058394853236066, 3.0000623233243306,
        2.9941252765509074, 2.9880372703446341, 2.9818072300784348,
        2.9754440810971623, 2.9689567487635244, 2.9623541584418795,
        2.9556452354722862, 2.9488389052235604, 2.9419440930462062,
        2.934969724310426, 2.9279247243605773, 2.920818018552545,
        2.9136585322536392, 2.9064551908080372, 2.8992169195864204,
        2.8919526439401824, 2.8846712892175148, 2.8773817807902016,
        2.8700930440050505, 2.8628140042257844, 2.8555535868045028,
        2.8483207170995417, 2.8411243204667356, 2.8339733222686125,
        2.8268766478410714, 2.8198432226219441, 2.8128819717280176,
        2.8060018211600655, 2.7992116945314209, 2.7925205220196068,
        2.7859372156066713, 2.7794707376816343, 2.773129911504463,
        2.7669239179283251, 2.7608595552138762, 2.7549377050729169,
        2.7491568665763655, 2.7435158964261239, 2.7380135121843021,
        2.7326484818477597, 2.727419555166942, 2.7223254884736034,
        2.7173650357533239, 2.7125369518230178, 2.7078399911960025,
        2.7032729084972393, 2.6988344583140842, 2.6945233952541061,
        2.6903384739078637, 2.6862784488646989, 2.6823420747291942,
        2.6785281060985509, 2.6748352975741971, 2.6712624037391475,
        2.667808179201439, 2.6644713785563661, 2.6612507563935051,
        2.6581450673184031, 2.6551530659249871, 2.6522735068093715,
        2.6495051445680233, 2.6468467337945882, 2.6442970290965455,
        2.6418547850631229, 2.639518756281265, 2.6372876973716433,
        2.6351603629111229, 2.6331355075007949, 2.6312118857448135,
        2.6293882522283867, 2.6276633615626381, 2.6260359683322059,
        2.6245048271371827, 2.6230686925769846, 2.6217263192454348,
        2.6204764617427321, 2.6193178746612054, 2.6182493126015496,
        2.6172695301596529, 2.6163772819292279, 2.6155713225111752,
        2.6148504065028968, 2.6142132884974161, 2.6136587230891717,
        2.6131854648807278, 2.6127922684721403, 2.6124778884520605,
        2.6122410794179718, 2.6120805959677282, 2.6119951927058076,
        2.6119836242200827, 2.6120446451039485, 2.6121770099684563,
        2.6123794734040748, 2.6126507899953171, 2.6129897143521776,
        2.6133950010753089, 2.613865404751321, 2.614399679976307,
        2.6149965813529712, 2.6156548634822583, 2.6163732809472067,
        2.6171505883576449, 2.6179855403054928, 2.6188768913815523,
        2.6198233961956157, 2.6208238093314318, 2.6218768853963228,
        2.6229813789785896, 2.6241360446778783, 2.625339637096352,
        2.6265909108270411, 2.6278886204604128, 2.6292315205990406,
        2.6306183658462339, 2.6320479107854404, 2.6335189100245509,
        2.6350301181556457, 2.6365802897744488, 2.6381681794792389,
        2.6397925418613171, 2.641452131534328, 2.6431457030803678,
        2.6448720110930157, 2.6466298101774779, 2.6484178549320263,
        2.6502348999527476, 2.65207969982483, 2.6539510091595191,
        2.6558475825511172, 2.6577681745896786, 2.659711539875476,
        2.6616764330057214, 2.6636616085770464, 2.6656658211853661,
        2.66768782542882, 2.6697263759111949, 2.6717802272137883,
        2.6738481339397695, 2.6759288506913967, 2.6780211320611351,
        2.6801237326493319, 2.6822354070460457, 2.684354909853762,
        2.6864809956666234, 2.6886124190802105, 2.6907479346960903,
        2.6928862971088132, 2.6950262609150424, 2.6971665807058995,
        2.6993060110850835, 2.7014433066549963, 2.7035772219954963,
        2.7057065117216172, 2.7078299304159015, 2.70994623267927,
        2.7120541731206678, 2.7141525063135954, 2.7162399868733997,
        2.7183153693940474, 2.7203774084625207, 2.7224248586885933,
        2.7244564746593012, 2.7264710109772947, 2.7284672222359871,
        2.7304438630317938, 2.7323996879678449, 2.7343334516299662,
        2.7362439086267449, 2.7381298135469052, 2.7399899209882239,
        2.7418229855533389, 2.7436277618335323, 2.7454030044223923,
        2.7471474679228356, 2.7488599069348139, 2.7505390760445274,
        2.7521837298556937, 2.7537926229680303, 2.7553645099700317,
        2.756898145464489, 2.7583922840459882, 2.7598456803108009,
        2.7612570888580286, 2.7626252642800182, 2.7639489611857195,
        2.7652269341479476, 2.7664579378048613, 2.7676407266538403,
        2.768774055580546, 2.7698566783568266, 2.7708873518744173,
        2.7718648243863981, 2.7727878679763687, 2.7736551934686657,
        2.7744659199491926, 2.775220180338112, 2.7759185158603734,
        2.7765614064297757, 2.7771493558145286, 2.7776828591558722,
        2.7781624146887558, 2.7785885195502673, 2.7789616712656611,
        2.779282367214889, 2.7795511048490793, 2.7797683815715239,
        2.7799346948099943, 2.7800505419910868, 2.7801164205323046,
        2.7801328278664452, 2.7801002614045447, 2.780019218576574,
        2.7798901968059737, 2.7797136935080782, 2.7794902061139206,
        2.7792202320469617, 2.778904268722481, 2.7785428135639476,
        2.7781363640084065, 2.7776854174598773, 2.7771904713525073,
        2.7766520231064638, 2.7760705701450759, 2.7754466098909272,
        2.7747806397639527, 2.7740731571944997, 2.7733246595969958,
        2.7725356444012386, 2.7717066090245726, 2.7708380508923809,
        2.7699304674308229, 2.76898435606003, 2.7680002141956903,
        2.7669785392712947, 2.7659198287133884, 2.7648245799251563,
        2.76369329035025, 2.7625264574038604, 2.7613245785020362,
        2.7600881510753887, 2.7588176725474605, 2.7575136403424851,
        2.7561765518740962, 2.7548069045708274, 2.7534051958632317,
        2.75197192315967, 2.7505075838872903, 2.7490126754792952,
        2.7474876953523539, 2.745933140923543, 2.7443495096187323,
        2.7427372988671648, 2.7410970060866822, 2.7394291286983754,
        2.7377341641290895, 2.7360126098000417, 2.7342649631317695,
        2.732491721552508, 2.7306933824813595, 2.7288704433413806,
        2.72702340155722, 2.7251527545525449, 2.7232589997450161,
        2.7213426345645524, 2.7194041564270774, 2.7174440627613903,
        2.7154628509862735, 2.713461018526405, 2.7114390628091658,
        2.7093974812473323, 2.7073367712718897, 2.705257430301637,
        2.7031599557669592, 2.7010448450771145, 2.6989125956674003,
        2.6967637049585931, 2.6945986703630807, 2.6924179893192877,
        2.6902221592379854, 2.6880116775539262, 2.68578704167692,
        2.6835487490337959, 2.6812972970569242, 2.67903318315514,
        2.6767569047620214, 2.6744689592971143, 2.6721698441796957,
        2.6698600568379316, 2.6675400946907089, 2.6652104551649018,
        2.6628716356833495, 2.6605241336608914, 2.6581684465282067,
        2.6558050717106512, 2.6534345066216956, 2.6510572486928035,
        2.6486737953434503, 2.6462846439939107, 2.643890292076696,
        2.6414912370018673, 2.6390879761958268, 2.636681007090214,
        2.6342708271030189, 2.6318579336486132, 2.629442824159574,
        2.6270259960600741, 2.6246079467704591, 2.6221891737048169,
        2.619770174295502, 2.6173514459716167, 2.6149334861416578,
        2.6125167922354908, 2.610101861676962, 2.607689191889043,
        2.6052792802912794, 2.6028726243057938, 2.6004697213629133,
        2.5980710688820805, 2.5956771642769247, 2.5932885049870005,
        2.5909055884202994, 2.5885289120110246, 2.5861589731770271,
        2.5837962693338761, 2.5814412979234107, 2.5790945563512317,
        2.5767565420391758, 2.5744277524259749, 2.5721086849273709,
        2.5697998369600317, 2.5675017059504328, 2.5652147893246533,
        2.5629395845068643, 2.560676588907655, 2.5584262999638443,
        2.5561892150944563, 2.5539658317220759, 2.5517566472638809,
        2.5495621591491515, 2.5473828648049857, 2.5452192616418841,
        2.5430718470926639, 2.5409411185754305, 2.5388275735131054,
        2.5367317093377024, 2.5346540234552157, 2.5325950133115325,
        2.5305551762819012, 2.5285350099250068, 2.5265350112624425,
        2.5245556788660553, 2.5225975069519979, 2.5206610018149882,
        2.5187466364381539, 2.5168549693920368, 2.5149859889617061,
        2.5131382672123963, 2.5113098058874863, 2.5094986923489855,
        2.5077029806351003, 2.5059207368835685, 2.5041500228329459,
        2.5023889018197853, 2.5006354366074905, 2.4988876901637913,
        2.497143725386771, 2.4954016051911325, 2.4936593924916739,
        2.4919151502009118, 2.4901669412338641, 2.488412828502975,
        2.4866508749215743, 2.4848791434076656, 2.4830956968698903,
        2.481298598220727, 2.4794859103804878, 2.4776556962582039,
        2.4758060187716535, 2.4739349408221742, 2.4720405253401738,
        2.470120835235087, 2.4681739334084711, 2.4661978827890341,
        2.4641907462834478, 2.4621505868011528, 2.4600754672665288,
        2.4579634505867443, 2.4558125996775444, 2.453620977447875,
        2.4513866468182592, 2.4491076706975523, 2.446782112000133,
        2.4444080336448279, 2.44198349853308, 2.4395065695964155,
        2.4369753097352107, 2.434387781864062, 2.4317420488998116,
        2.4290361737571771, 2.4262682193464169, 2.4234362485845025,
        2.4205383243859209, 2.4175725096559559, 2.4145368673217233,
        2.4114294602864468, 2.4082483514626225, 2.4049916037764589,
        2.4016572801292511, 2.3982434434356397, 2.3947481566190678,
        2.3911694825848682, 2.3875054842455423, 2.383754224518055,
        2.37991376632185, 2.3759821725596004, 2.3719575061481897,
        2.3678378300081255, 2.3636212070478004, 2.3593057001800677,
        2.3548893723172077, 2.35037028637893, 2.3457465052704833,
        2.3410160919206064, 2.3361771092191597, 2.331227620105607,
        2.3261656874659025, 2.320989374293374, 2.3156967432230764,
        2.3102858579885579, 2.3047547792330016, 2.2991015761474016,
        2.2933242942552647, 2.2874210446131729, 2.2813897569367549,
        2.2752288613903717, 2.2689355018254251, 2.2625153922707737,
        2.2559955283349131, 2.2494114757939561, 2.2427975141319365,
        2.2361884232469791, 2.2296188017359149, 2.2231233137031738,
        2.2167365995907748, 2.2104933083909279, 2.2044280860011769,
        2.1985755794388142, 2.1929704353220716, 2.1876473004091754,
        2.1826408214018658, 2.177985645027992, 2.1737164180190645,
        2.1698677870814329, 2.1664743989325133, 2.1635709003156265,
        2.16119193793294, 2.1593721585084578, 2.1581462087659173,
        2.1575487354242719, 2.1576143852100569, 2.158377804841515,
        2.1598736410351451, 2.1621365405127766, 2.1652011500092958,
        2.1691021162060662, 2.1738740859202896, 2.1795517056271119,
        2.1861696227192833, 2.1937624820486916, 2.2023649355157775,
        2.2120116155261971, 2.2227372083903485, 2.2345762512780962,
        2.2475636940493877, 2.2617333475358845, 2.2771219503568307,
        2.293746734253749, 2.3115764911423486, 2.3305605060580246,
        2.3506509918268659, 2.3717990222488163, 2.3939560838071552,
        2.4170735138549593, 2.4411027036307504, 2.46599502490347,
        2.491701856471813, 2.5181745746016722, 2.5453645564700178,
        2.5732231789271789, 2.6017018189386665, 2.6307518534262995,
        2.6603246593383734, 2.6903716135977542, 2.7208440931540161,
        2.7516934749298829, 2.7828711358634544, 2.8143284529066759,
        2.8460168029751367, 2.8778875630106979, 2.9098921099538342,
        2.9419818207360948, 2.9741080722939257, 3.0062222415708328,
        3.038275705484569, 3.0702198409954184, 3.102006024986681,
        3.133585634552249, 3.1649100461751045, 3.1959306380473085,
        3.226598783638527, 3.2568658694662025, 3.2866832459473243,
        3.3160023634123581, 3.344774395695906, 3.3729512798029515,
        3.4004829910149734, 3.4273325750020476, 3.4534955340358966,
        3.4789804407659966, 3.503793906122822, 3.5279433042325934,
        3.5514357326736743, 3.5742783889945233, 3.5964784345999345,
        3.6180430439534841, 3.6389793868162368, 3.6592946346391275,
        3.6789959582569196, 3.6980905287370325, 3.7165855170594271,
        3.7344880942333596, 3.751805431259061, 3.7685446991374913,
        3.784713068873423, 3.8003177114722497, 3.8153657979357236,
        3.8298644992555042, 3.843820986451262, 3.8572424305148378,
        3.8701360024519933, 3.8825088732618549, 3.8943682139530629,
        3.9057211955241438, 3.9165749889768224, 3.9269367653185712,
        3.9368136955452915, 3.9462129506660686, 3.9551417016777726,
        3.9636071195877696, 3.9716163753947433, 3.9791766401047917,
        3.9862950847169381, 3.992978880237513, 3.9992351976659108,
        4.0050712080078741, 4.0104940822622659, 4.0155109914325582,
        4.0201291065237728, 4.0243555985396959, 4.0281976384756355,
        4.0316623973402432, 4.034757046136014, 4.0374887558638584,
        4.0398646975245809, 4.0418920421225319, 4.0435779606663171,
        4.0449296241446335, 4.0459542035741523, 4.0466588699520312,
        4.047050794275755, 4.0471371475536513, 4.0469251007910954,
        4.04642182498241, 4.0456344911348339, 4.0445702702527484,
        4.0432363333351269, 4.04163985138549, 4.0397879954114133,
        4.0376879364049252, 4.0353468453760222, 4.0327718933298868,
        4.02997025125946, 4.0269490901835443, 4.0237155810826328,
        4.020276894971194, 4.0166402028613577, 4.0128126757327385,
        4.0088014846290285, 4.0046138004439733, 4.0002567944843248,
        3.9957376369164974, 3.9910635010493509, 3.9862415515104925,
        3.9812789769364718, 3.9761828996943303, 3.9709606125270107,
        3.9656182728725309, 3.9601592189797876, 3.9545856537330066,
        3.9488999504618882, 3.9431044161900384, 3.9372013819509557,
        3.9311931701062925, 3.9250821061606409, 3.9188705144712235,
        3.9125607198104815, 3.9061550468103539, 3.8996558201421681,
        3.8930653644721249, 3.8863860044573966, 3.8796200647728685,
        3.8727698700809183, 3.8658377450418695, 3.8588260143225108,
        3.8517370025842594, 3.8445730344998146, 3.8373364347306196,
        3.8300295279322283, 3.822654638782637, 3.8152140919395947,
        3.8077102120733843, 3.8001453238357876, 3.7925217519046113,
        3.7848418209433734, 3.7771078556050242, 3.7693221805684045,
        3.7614871204932, 3.7536050000442067, 3.7456781438789384,
        3.7377088766722433, 3.72969952308506, 3.721652407780677,
        3.7135698554244052, 3.7054541906840428, 3.6973077382224733,
        3.6891328226977032, 3.6809317687847884, 3.6727069011463782,
        3.6644605444404474, 3.6561950233339133, 3.6479126625004978,
        3.6396157865932111, 3.6313067202816636, 3.6229877882321739,
        3.6146613151035383, 3.6063296255710413, 3.5979950442879551,
        3.5896598959237327, 3.5813265051450522, 3.5729971966117939,
        3.5646742949937598, 3.5563601249543551, 3.5480570111539862,
        3.5397672782645726, 3.5314932509443313, 3.523237253859075,
        3.5150016116767939, 3.5067886490613618, 3.4986006906748104,
        3.490440061182015, 3.4823090852508312, 3.4742100875443347,
        3.46614539272449, 3.4581173254624651, 3.4501282104138791,
        3.4421803722522979, 3.434276135635534, 3.4264178252468511,
        3.4186077656791887, 3.4108482818107411, 3.4031416977140441,
        3.3954903396949123, 3.3878965278705118, 3.3803625994991724,
        3.3728908444977508, 3.3654836745092855, 3.3581426901682669,
        3.3508674782648149, 3.3436568146258363, 3.3365095967630278,
        3.3294246748432692, 3.3224009162126733, 3.3154371819920478,
        3.308532335553247, 3.301685239452659, 3.294894756553937,
        3.2881597495938553, 3.2814790813668693, 3.2748516146379352,
        3.2682762121834932, 3.2617517367885891, 3.2552770512200531,
        3.2488510182581476, 3.2424725006680188, 3.236140361238669,
        3.229853462744634, 3.2236106679521903, 3.2174108396427852,
        3.2112528405912766, 3.2051355335703815, 3.1990577813666161,
        3.193018446740763, 3.1870163924790136, 3.181050481350443,
        3.1751195761292381, 3.1692225396059839, 3.1633582345429265,
        3.1575255237092676, 3.1517232698961206, 3.1459503358750789,
        3.1402055844142658, 3.1344878782942964, 3.1287960802923842,
        3.1231290531866036, 3.1174856597397103, 3.111864762738962,
        3.1062652249621987, 3.1006859091761791, 3.0951256781592456,
        3.0895833946868154, 3.0840579215408, 3.0785481214890433,
        3.0730528573040883, 3.06757099177559, 3.0621013876662784,
        3.056642907756197, 3.0511944148228518, 3.0457547716396549,
        3.0403228409822787, 3.0348974856239641, 3.0294775683448085,
        3.0240619519178393, 3.0186494991224637, 3.0132390727286729,
        3.0078295355122768, 3.0024197502558052, 2.9970085797310633,
        2.9915948867031368, 2.9861775339651615, 2.9807553842864576,
        2.9753273004346328, 2.9698921451979445, 2.9644487813423006,
        2.9589960716470363, 2.9535328788878656, 2.9480580658359803,
        2.9425704952793161, 2.9370690299765516, 2.9315525327347265,
        2.9260198662309085, 2.9204698935296025, 2.9149014766132888,
        2.9093134804308383, 2.9037047617559537, 2.8980741999717914,
        2.8924206120617844, 2.886742975399653, 2.8810391987010577,
        2.8753045369226879, 2.8695331763711662, 2.863719463730495,
        2.8578576833012312, 2.8519421419777262, 2.8459671384894558,
        2.8399269745285505, 2.8338159507009757, 2.827628368014413,
        2.8213585273237838, 2.8150007295524051, 2.8085492755862442,
        2.8019984663305739, 2.7953426026849559, 2.7885759855413736,
        2.7816929158064845, 2.7746876943778829, 2.7675546221553753,
        2.7602880000318835, 2.7528821289146372, 2.7453313097037997,
        2.7376298432898882, 2.7297720305785527, 2.721752172469011,
        2.7135645698529047, 2.705203523643692, 2.6966633347257409,
        2.6879383040201503, 2.6790227323587357, 2.6699109208434582,
        2.6605971698093254, 2.6510757817161, 2.6413410531405166,
        2.6313872969451655, 2.6212087809162674, 2.6107998976065314,
        2.600154694294825, 2.5892681712147483, 2.5781328791288769,
        2.5667576887493389, 2.5551919968179515, 2.5435015200307824,
        2.531749525611485, 2.5200002337267837, 2.5083175192869849,
        2.4967653819676165, 2.4854077763462818, 2.4743086733154849,
        2.4635320378661705, 2.4531418371158806, 2.443202037418915,
        2.4337766054042649, 2.424929507605853, 2.4167247105854339,
        2.4092261808935591, 2.4024978850945269, 2.3966037897388675,
        2.3916078613862863, 2.3875740665905107, 2.3845663719103003,
        2.3826487438982782, 2.3818851491146136, 2.3823395541143793,
        2.3840759254486148, 2.3871582296906912, 2.3916504333700739,
        2.397616503073793, 2.4051204053043129, 2.4142261067791586,
        2.4249975735697316, 2.4374987735947768, 2.4517936696501237,
        2.4679462386836675, 2.4860204185022812, 2.5060802552179751,
        2.5281894952236117, 2.552412714290333, 2.5788121990828339,
        2.6074561202697883, 2.6383734453303305, 2.6714957916724837,
        2.7067155735092765, 2.7439310890902666, 2.7830383475162557,
        2.8239341872868216, 2.8665151471804253, 2.9106778742777708,
        2.9563189765298166, 3.0033350760164859, 3.0516227897157608,
        3.1010787364522785, 3.15159953438341, 3.2030818019066363,
        3.2554221573322781, 3.3085172190006373, 3.362263605238645,
        3.4165579343905761, 3.4712968247778981, 3.52637689474179,
        3.5816947626136622, 3.6371470467244604, 3.6926303654167931,
        3.7480413370095205, 3.8032765798441011, 3.8582327122616569,
        3.9128063525770114, 3.9668941191377658, 4.0203926302666666,
        4.0731985043388192, 4.1252083595415456, 4.1763188146336443,
        4.2264264867555186, 4.2754279975484906, 4.32321995620154,
        4.3696990063471821, 4.4147616962936258, 4.4583048381809895,
        4.5002245140571588, 4.5404188210244358, 4.5787806766724568,
        4.615237507830237, 4.6498024352950882, 4.6825230890956906,
        4.713441919763385, 4.7426033928687374, 4.7700512438934837,
        4.7958294721669876, 4.81998198166858, 4.8425527108349353,
        4.8635855856558576, 4.8831245366157807, 4.9012134925695774,
        4.9178963829697979, 4.9332171370495264, 4.9472196841216105,
        4.95994795346822, 4.9714458743812759, 4.9817573761583516,
        4.99092638808069, 4.9989968394513848, 5.0060126595521961,
        5.0120177776733286, 5.0170561231173458, 5.0211716251650831,
        5.02440821310941, 5.02680981624442, 5.0284203638613967,
        5.0292837852512076, 5.0294440096928525, 5.0289449665236292,
        5.0278305849098484, 5.0261447944701194, 5.023931523599896,
        5.0212347040897, 5.0180982582825049, 5.0145661347145358,
        5.0106822094166468, 5.0064905590470552, 5.0020347065510808,
        4.9973595981434826, 4.992500697233222, 4.9874699193552452,
        4.9822696972278946, 4.9769038868494029, 4.9713757905035081,
        4.9656889111012026, 4.9598466790392282, 4.9538525509280449,
        4.9477099738978429, 4.9414223985103041, 4.9349932740838947,
        4.9284260503867463, 4.92172417702763, 4.9148911036667311,
        4.9079302799485625, 4.9008451555286916, 4.8936391800512524,
        4.8863158031650888, 4.8788784745231686, 4.8713306437723594,
        4.8636757605591319, 4.8559172745332244, 4.8480586353522686,
        4.8401032926546419, 4.8320546960942012, 4.8239162953202381,
        4.81569153997997, 4.80738387972469, 4.7989967642014655,
        4.790533643059625, 4.7819979659517822, 4.773393182520369,
        4.7647227424224008, 4.7559900952998051, 4.7471986908052894,
        4.7383519785921937, 4.7294534082968687, 4.7205064295838861,
        4.7115144920912524, 4.7024810454740118, 4.6934095393781421,
        4.6843034234500927, 4.6751661473485964, 4.6660011607139964,
        4.6568119131984167, 4.647601854450504, 4.6383744341204549,
        4.6291331018563957, 4.6198813073064438, 4.610622500122826,
        4.6013601299509022, 4.5920976464402088, 4.5828384992451463,
        4.5735861380063536, 4.56434401238304, 4.5551155720156542,
        4.5459042665534071, 4.5367135456530221, 4.527546858958992,
        4.5184076561184732, 4.5092993867821365, 4.5002255006024576,
        4.4911894472203278, 4.4821946762934486, 4.4732446374696258,
        4.4643427803940972, 4.4554925547144757, 4.4466974100895893,
        4.4379607961568333, 4.4292861625759841, 4.420676958979942,
        4.4121366350586406, 4.4036686403417411, 4.3952764248230594,
        4.3869634371916151, 4.3787331297520975, 4.37058894480202,
        4.3625343523317932, 4.3545727459076318, 4.3467077155304175,
        4.3389415424282527, 4.3312732578706576, 4.3237005843630456,
        4.3162214408263209, 4.3088336697871048, 4.3015351414378156,
        4.2943237159706387, 4.28719725720104, 4.2801536276345971,
        4.2731906902411572, 4.2663063078298764, 4.2594983432681994,
        4.2527646593964379, 4.2461031190789065, 4.2395115851508924,
        4.2329879204700056, 4.2265299878886537, 4.2201356502543321,
        4.2138027704192, 4.2075292112359923, 4.2013128355520575,
        4.1951515062160949, 4.18904308608601, 4.1829854380064475,
        4.1769764248297188, 4.1710139094064553, 4.1650957545869245,
        4.1592198232257633, 4.1533839781674935, 4.14758608226706,
        4.141823998372784, 4.13609558933505, 4.13039871800939,
        4.1247312472383113, 4.1190910398808107, 4.1134759587823577,
        4.1078838667948041, 4.1023126267704484, 4.0967601015561357,
        4.0912241540067882, 4.0857026469722122, 4.08019344329913,
        4.0746944058425081, 4.0692033974530419, 4.0637182809797769,
        4.05823691927223, 4.0527571751832632, 4.0472769115656186,
        4.0417939912624519, 4.036306277132617, 4.0308116320219165,
        4.025307918781432, 4.0197930002673408, 4.01426473932264,
        4.0087209987995971, 4.0031596415533564, 3.9975785304307161,
        3.9919755282828286, 3.9863484979623829, 3.9806953023188751,
        3.975013804201029, 3.9693018664610777, 3.963557351951748,
        3.9577781235203697, 3.9519620440194108, 3.9461069762996646,
        3.9402107832088369, 3.9342713276049173, 3.928286472331505,
        3.9222540802392514, 3.9161720141819112, 3.9100381370143431,
        3.9038503115613481, 3.8976064007632107, 3.8913042672095908,
        3.88494177447016, 3.8785167834130356, 3.8720271623721616,
        3.8654707590331112, 3.8588454780559105, 3.8521490776552181,
        3.8453802917920088, 3.8385402774164983, 3.8316311672257459,
        3.8246549474662355, 3.8176136613635907, 3.810509331498229,
        3.8033439879106825, 3.7961196579398457, 3.7888383699103274,
        3.7815021517873251, 3.7741130316596041, 3.7666730375822719,
        3.7591841976141569, 3.7516485398156885, 3.7440680922472147,
        3.7364448829664871, 3.7287809400376024, 3.7210782915175691,
        3.7133389654637785, 3.7055649899402572, 3.6977583930047113,
        3.6899212027170152, 3.6820554471397831, 3.6741631543310578,
        3.6662463523471276, 3.6583070692527908, 3.6503473331081695,
        3.642369171966235, 3.6343746138972479, 3.6263656869530254,
        3.6183444191947527, 3.6103128386853816, 3.6022729734830348,
        3.5942268516452054, 3.5861765012367508, 3.5781239503124462,
        3.5700712269351969, 3.5620203591644208, 3.55397337505787,
        3.5459323026810639, 3.5378991700841138, 3.5298760053372962,
        3.521864836496218, 3.5138676916186449, 3.5058865987659513,
        3.4979235859957507, 3.4899806813759886, 3.4820599129590062,
        3.4741633088027446, 3.4662928969749842, 3.4584507055297422,
        3.4506387625277406, 3.442859096032493, 3.4351137340995352,
        3.42740470478919, 3.4197340361632911, 3.4121037562810326,
        3.4045158931996076, 3.3969724749844037, 3.3894755296912029,
        3.3820270853796743, 3.374629170113256, 3.3672838119459758,
        3.3599930389427994, 3.352758879161589, 3.3455833606624221,
        3.3384685115039847, 3.3314163597498747, 3.3244289334555508,
        3.317508260683351, 3.3106563694857263, 3.30387528795616,
        3.2971670440490524, 3.2905336661158273, 3.2839771814222978,
        3.2774996202167364, 3.2711030044966884, 3.2647893791087323,
        3.258560725817266, 3.2524191885385028, 3.2463658308771342,
        3.2403990338184387, 3.2345160980447041, 3.2287144863804955,
        3.2229915985725377, 3.2173448572151195, 3.211771676656912,
        3.2062694742156106, 3.2008356661428978, 3.1954676690741555,
        3.1901628995013822, 3.1849187739749238, 3.1797327090219509,
        3.1746021211751456, 3.1695244269683025, 3.1644970429348018,
        3.1595173856045817, 3.1545828715172259, 3.1496909171989635,
        3.1448389391839524, 3.1400243540055479, 3.1352445781995475,
        3.1304970282946547, 3.125779120825853, 3.1210882723247728,
        3.1164218993250374, 3.1117774183624847, 3.1071522459623053,
        3.1025437986644739, 3.0979494930014195, 3.0933667454995484,
        3.088792972699296, 3.0842255911301768, 3.0796620173255587,
        3.0750996678161071, 3.0705359591364423, 3.0659683078237796,
        3.0613941304025709, 3.056810843411327, 3.0522158633820489,
        3.0476066068444347, 3.0429804903377944, 3.0383349303901217,
        3.0336673435314516, 3.0289751463038184, 3.0242557552322347,
        3.0195065868509712, 3.0147250576967775, 3.0099085842972304,
        3.0050545831893585, 3.0001604709037424, 2.9952236639725065,
        2.9902415789319807, 2.9852116323127977, 2.9801312406459055,
        2.9749978204675078, 2.9698087883108362, 2.9645615607051292,
        2.9592535541839395, 2.9538821852840176, 2.9484448705346926,
        2.9429390264679425, 2.9373620696212539, 2.9317114165235991,
        2.9259844837077589, 2.9201786877081717, 2.9142914450597739,
        2.9083201722853369, 2.902262285941394, 2.8961152025073624,
        2.8898763386534463, 2.8835431105324085, 2.8771129357407585,
        2.8705832278555308, 2.8639514115881073, 2.8572148808466231,
        2.8503711147715625, 2.8434173566617895, 2.8363515007219853,
        2.8291697680748475, 2.8218795269881767, 2.8145158265517374,
        2.8071248629882217, 2.7997511594556324, 2.7924398899983833,
        2.7852359928341577, 2.7781844914137213, 2.7713303783748593,
        2.7647186574910791, 2.7583943285214212, 2.7524023926705885,
        2.7467878506187731, 2.7415957032360629, 2.7368709513279494,
        2.7326585957211762, 2.72900363722783, 2.7259510766789403,
        2.7235459148913939, 2.7218331526826529, 2.7208577908775826,
        2.72066483029316, 2.7212992717549889, 2.72280611607975,
        2.7252303640874955, 2.728617016603057, 2.7330110744454048,
        2.7384575384307275, 2.7450014093942539, 2.7526876881147562,
        2.7615613755360231, 2.7716674721318335, 2.783050979676541,
        2.7957568963669051, 2.8098302302779423, 2.82531596215871,
        2.8422591483673636, 2.8607046360245638, 2.8806978512756345,
        2.9022826221125464, 2.9255068844750123, 2.95039120449113,
        2.9768881831277, 3.00492305155046, 3.0344251488647642,
        3.0653222160137483, 3.0975425729879009, 3.1310143305256166,
        3.1656656749685146, 3.2014247653473364, 3.238219770557119,
        3.2759788559300596, 3.3146301880882647, 3.3541019331811261,
        3.3943222575350651, 3.4352193274091412, 3.4767213090895779,
        3.5187563688492802, 3.561252672967401, 3.6041383877264361,
        3.6473416793940805, 3.690790714254951, 3.7344136585904972,
        3.7781386786689217, 3.8218939407745625, 3.865607611182917,
        3.9092078561726717, 3.9526228420243128, 3.9957807350012651,
        4.0386097014195492, 4.0810379074608578, 4.1229935196613576,
        4.1644047035831795, 4.205199627491508, 4.245306452165897,
        4.2846533590808509, 4.3231684724675867, 4.3607800749827712,
        4.3974160108520923, 4.4330053344008675, 4.467473989484982,
        4.5007686439131147, 4.5328874275346518, 4.5638491941469,
        4.5936696870897631, 4.6223658597982791, 4.6499542272551579,
        4.6764514629048906, 4.7018741829342892, 4.7262390242038972,
        4.7495626161146021, 4.7718615907650985, 4.7931525792774252,
        4.8134522131187625, 4.8327771236413177, 4.8511439422358791,
        4.868569300277299, 4.8850698291491845, 4.9006621602257212,
        4.91536292489215, 4.9291887545302693, 4.9421562805131289,
        4.9542821342247585, 4.965582947050569, 4.9760753503632289,
        4.9857759755384059, 4.994701453966619, 5.0028684170282594,
        5.0102934960919612, 5.0169933225495775, 5.0229845277749474,
        5.0282837431488288, 5.03290760005402, 5.0368727298657605,
        5.04019576396776, 5.0428933337416035, 5.0449820705613293,
        5.0464786058115543, 5.0473995708763786, 5.0477615971255174,
        5.0475813159431429, 5.0468753587129163, 5.0456603568144374,
        5.0439529416229911, 5.0417697445232657, 5.0391273968918533,
        5.0360425301131473, 5.0325317755653085, 5.0286117646211572,
        5.0242991286715633, 5.01961049909502, 5.0145625072632338,
        5.0091717845641215, 5.0034549623796964, 4.9974286720816137,
        4.9911095450518106, 4.9845142126742887, 4.9776593063298806,
        4.970561457394755, 4.9632372972504584, 4.9557034572763055,
        4.9479765688557693, 4.9400732633652193, 4.9320101721818821,
        4.9238039266925338, 4.915471158275464, 4.9070284983063761,
        4.8984925781710622, 4.8898800292482116, 4.8812074829152232,
        4.8724915705478757, 4.8637489235482763, 4.854996173239428,
        4.8462499511566284, 4.8375268882462974, 4.8288436170934554,
        4.8202167657354433, 4.8116629748164819, 4.8031988500673641,
        4.7948410938372072, 4.7866061418310748, 4.7785111151122708,
        4.7705685684501233, 4.76277971748517, 4.755141211552651,
        4.7476503853515739, 4.7403043069471957, 4.7331001410088849,
        4.7260350172969421, 4.7191060781824365, 4.7123104614811355,
        4.7056453066609043, 4.699107752583064, 4.6926949383291143,
        4.6864040029083922, 4.6802320853564785, 4.674176324687104,
        4.6682338599320987, 4.6624018301203876, 4.6566773742692877,
        4.6510576314066414, 4.6455397405601193, 4.6401208407557544,
        4.6347980710173626, 4.6295685703699316, 4.6244294778400548,
        4.6193779324534177, 4.6144110732325982, 4.6095260392044883,
        4.6047199693959477, 4.5999900028334721, 4.5953332785372929,
        4.59074693553798, 4.5862281128600166, 4.5817739495271645,
        4.5773815845636463, 4.5730481569965731, 4.5687708058561327,
        4.5645466701597037, 4.5603728889362252, 4.5562466012117238,
        4.5521649460123266, 4.548125062362363, 4.5441240892840478,
        4.54015916580676, 4.5362274309571919, 4.5323260237569745,
        4.52845208323332, 4.52460274841257, 4.5207751583197, 4.5169664519760993,
        4.5131737684116358, 4.5093942466537511, 4.505625025722793,
        4.5018632446452749, 4.4981060424486126, 4.4943505581597814,
        4.490593930799041, 4.48683329939253, 4.4830658029699837,
        4.4792885805554334, 4.4754987711715382, 4.4716935138458611,
        4.4678699476045347, 4.4640252114710828, 4.4601564444704573,
        4.4562607856324092, 4.4523353739782712, 4.4483773485324241,
        4.4443838483251712, 4.4403520123791935, 4.4362789797162376,
        4.4321618893742452, 4.42799788034751, 4.4237840917447286,
        4.4195176623631758, 4.41519573185891, 4.4108154375099753,
        4.406373923177175, 4.4018683145124244, 4.3972957874122489,
        4.3926533886277319, 4.3879390253958315, 4.3831527417497771,
        4.3782954422201739, 4.3733679021816547, 4.3683709472504146,
        4.3633053848497072, 4.3581720289683243, 4.3529716912267409,
        4.3477051841031162, 4.3423733197580532, 4.336976910477766,
        4.33151676849381, 4.3259937060606246, 4.3204085354267114,
        4.3147620688397579, 4.3090551185499146, 4.3032884968024909,
        4.2974630158485985, 4.2915794879391971, 4.285638725316768,
        4.279641540235736, 4.2735887449414509, 4.2674811516821105,
        4.2613195727108986, 4.2551048202706605, 4.2488377066138519,
        4.24251904398889, 4.2361496446402311, 4.2297303208219006,
        4.2232618847812642, 4.2167451487639385, 4.2101809250208913,
        4.2035700258015414, 4.1969132633526236, 4.1902114499236918,
        4.18346539776425, 4.1766759191199965, 4.1698438262431088,
        4.1629699313798607, 4.1560550467796782, 4.14909998469238,
        4.14210555736406, 4.135072577044336, 4.1280018559831166,
        4.1208942064284306, 4.1137504406275287, 4.1065713708299727,
        4.0993578092849212, 4.0921105682412291, 4.0848304599459313,
        4.0775182966490053, 4.0701748905984134, 4.0628010540420876,
        4.0553975992313429, 4.0479653384128653, 4.0405050838351686,
        4.0330176477465249, 4.0255038423965877, 4.0179644800336609,
        4.0104003729070943, 4.0028123332646635, 3.9952011733531747,
        3.9875677054250236, 3.9799127417280431, 3.9722370945070238,
        3.9645415760151179, 3.9568269984996602, 3.9490941742063956,
        3.9413439153886891, 3.9335770342931711, 3.9257943431669773,
        3.9179966542605253, 3.9101847798207694, 3.9023595320985782,
        3.8945217233426219, 3.8866721657980539, 3.8788116717165568,
        3.8709410533470638, 3.8630611229356284, 3.8551726927320473,
        3.8472765749877822, 3.8393735819467381, 3.8314645258600866,
        3.8235502189776942, 3.8156314735450314, 3.8077091018127858,
        3.7997839160294862, 3.7918567284440172, 3.783928351302515,
        3.7759995968568671, 3.7680712773547773, 3.7601442050427862,
        3.7522191921722086, 3.7442970509910189, 3.7363785937477187,
        3.7284646326894983, 3.7205559800662074, 3.7126534481268085,
        3.7047578491202224, 3.6968699952947408, 3.6889906988973813,
        3.681120772178943, 3.6732610273866153, 3.6654122767693926,
        3.6575753325775424, 3.6497510070566617, 3.64194011245761,
        3.6341434610283705, 3.626361865018346, 3.6185961366745074,
        3.61084708824544, 3.6031155319831556, 3.5954022801314292,
        3.5877081449429018, 3.5800339386650206, 3.5723804735436238,
        3.5647485618324488, 3.5571390157762, 3.5495526476238828,
        3.5419902696265853, 3.534452694030513, 3.5269407330841633,
        3.5194551990386116, 3.5119969041408967, 3.5045666606384649,
        3.4971652807820703, 3.489793576819793, 3.4824523609985549,
        3.47514244556951, 3.4678646427804156, 3.4606197648786896,
        3.4534086241143966, 3.4462320327354981, 3.4390908029912635,
        3.4319857471290947, 3.4249176773993408, 3.4178874060488282,
        3.4108957453264384, 3.4039435074839988, 3.3970315047650281,
        3.3901605494212004, 3.3833314537014774, 3.3765450298518149,
        3.369802090124784, 3.3631034467660594, 3.3564499120247562,
        3.3498422981499494, 3.3432814173897754, 3.3367680819940251,
        3.3303031042106386, 3.3238872962883792, 3.3175214704731792,
        3.311206439022397, 3.3049430141607354, 3.2987320082052691,
        3.2925742332128696, 3.2864705019609444, 3.2804216252412322,
        3.2744284193263371, 3.2684916853719619, 3.2626122633918988,
        3.2567907344599192, 3.2510270366283209, 3.2453208490026579,
        3.2396718895553214, 3.2340798611389561, 3.2285444720808583,
        3.2230654287315268, 3.21764243815917, 3.2122752071659946,
        3.2069634426562557, 3.2017068514960321, 3.1965051405610656,
        3.1913580167283171, 3.1862651868749325, 3.1812263578724216,
        3.1762412365997359, 3.1713095299309946, 3.1664309447404531,
        3.161605187907627, 3.1568319663037743, 3.1521109868058708,
        3.1474419562899967, 3.1428245816311962, 3.138258569705521,
        3.133743627388224, 3.1292794615533133, 3.1248657790782697,
        3.1205022868387946, 3.1161886917088633, 3.1119247005640767,
        3.10771002028243, 3.1035443577365078, 3.0994274198019696,
        3.09535891335667, 3.0913385452743403, 3.0873660224303534,
        3.0834410517014716, 3.0795633399623608, 3.075732594088493,
        3.0719485209557158, 3.0682108274389028, 3.0645192204145366,
        3.0608734067571417, 3.0572730933435177, 3.0537179870481412,
        3.0502077947472372, 3.0467422233143968, 3.04332097962755,
        3.039943770562116, 3.0366103029916771, 3.0333202837928459,
        3.030073419840976, 3.02686941801271, 3.0237079851823228,
        3.0205888282252991, 3.017511654017198, 3.0144761694339519,
        3.01148208135112, 3.0085290966442297, 3.0056169221885116,
        3.0027452648590951, 2.9999138315323464, 2.9971223290831368,
        2.9943704643874831, 2.9916579443212497, 2.9889844757581132,
        2.9863497655750662, 2.983753520648242, 2.9811954478518108,
        2.9786752540616557, 2.9761926461536974, 2.9737473310027775,
        2.9713390154848378, 2.9689674064762435, 2.9666322108509471,
        2.964333135485135, 2.9620698872541467, 2.959842173033695,
        2.9576496997001418, 2.9554921741275133, 2.9533693031914834,
        2.9512807937688392, 2.9492263527344522, 2.9472056869628633,
        2.9452185033308265, 2.9432645087137033, 2.9413434099860916,
        2.9394549140244721, 2.9375987277044593, 2.9357745579010586,
        2.933982111489807, 2.9322210953457093, 2.9304912163453047,
        2.9287921813641531, 2.92712369727705, 2.9254854709591944,
        2.9238772092869292, 2.9222986191364493, 2.9207494073808524,
        2.9192292808980711, 2.917737946563038, 2.9162751112493353,
        2.9148404818351077, 2.9134337651950784, 2.9120546682042057,
        2.9107028977381222, 2.9093781606727869, 2.9080801638831413,
        2.9068086142451635, 2.9055632186347515, 2.9043436839262236,
        2.9031497169958178, 2.9019810247194924, 2.9008373139715316,
        2.8997182916284054, 2.8986236645660246, 2.8975531396583989,
        2.8965064237824905, 2.8954832238133084, 2.8944832466258061,
        2.8935061990962776, 2.8925517881003087, 2.8916197205125003,
        2.8907097032092066, 2.8898214430663649, 2.8889546469576888,
        2.888109021760024, 2.8872842743493856, 2.8864801115996146,
        2.8856962403880786, 2.8849323675894611, 2.8841882000777854,
        2.8834634447314671, 2.8827578084242869, 2.8820709980316366,
        2.8814027204299633, 2.8807526824935468, 2.880120591099093,
        2.8795061531213628, 2.8789090754365949, 2.8783290649192841,
        2.877765828445928, 2.877219072891636, 2.8766885051310545,
        2.8761738320419115, 2.8756747604975583, 2.8751909973743204,
        2.8747222495481268, 2.8742682238934667, 2.8738286272867488,
        2.8734031666035218, 2.8729915487187228, 2.8725934805074562,
        2.8722086688487076, 2.8718368206077329, 2.8714776426873647,
        2.8711308418883577, 2.8707961252932015, 2.8704731992053758,
        2.8701617720760875, 2.8698615468389161, 2.8695722631828855,
        2.8692937520695545, 2.8690258812146161, 2.8687685128180553,
        2.8685215112258056, 2.8682847400054015, 2.8680580630069268,
        2.8678413439777808, 2.8676344467024975, 2.8674372349513408,
        2.8672495725010259, 2.8670713231266678, 2.8669023506006481,
        2.8667425186994389, 2.8665916911978879, 2.8664497318692441,
        2.8663165044889807, 2.8661918728323519, 2.8660757006728033,
        2.8659678517854941, 2.8658681899453886, 2.8657765789264236,
        2.8656928825040375, 2.8656169644530736, 2.8655486885467449,
        2.8654879185607571, 2.8654345182703302, 2.8653883514490919,
        2.8653492818718718, 2.86531717331367, 2.8652918895484545,
        2.8652732943512844, 2.8652612514982136, 2.8652556247615313,
        2.8652562779168922, 2.8652630747395547, 2.8652758790029633,
        2.8652945544831989, 2.8653189649537456, 2.8653489741893621,
        2.8653844459657365, 2.8654252440561643, 2.8654712322359539,
        2.8655222742801967, 2.8655782339626086, 2.8656389750586566,
        2.8657043613423068, 2.8657742565887614, 2.8658485245729337,
        2.8659270290680783, 2.86600963385069, 2.86609620269394,
        2.866186599373036, 2.8662806876634059, 2.866378331337935,
        2.86647939417305, 2.8665837399425556, 2.86669123242087,
        2.8668017353834387, 2.8669151126040489, 2.8670312278578116,
        2.8671499449198268, 2.8672711275638689, 2.8673946395649574,
        2.867520344697978, 2.8676481067371933, 2.8677777894574779,
        2.8679092566336086, 2.8680423720400072, 2.8681769994514004,
        2.868313002642402, 2.86845024538785, 2.8685885914621823,
        2.8687279046401528, 2.8688680486964158, 2.8690088874054793,
        2.86915028454225, 2.8692921038811452, 2.8694342091970655,
        2.8695764642648149, 2.869718732858106, 2.8698608787525806,
        2.8700027657228455, 2.87014425754302, 2.870285217987838,
        2.870425510832336, 2.8705649998510205, 2.8707035488181174,
        2.8708410215093765, 2.8709772816980381, 2.8711121931596919,
        2.8712456196687706, 2.8713774249999249, 2.8715074729278416,
        2.8716356272271466, 2.8717617516721208, 2.8718857100377861,
        2.8720073660994232, 2.8721265836306054, 2.8722432264062161,
        2.8723571582013911, 2.8724682427905193, 2.8725763439483143,
        2.8726813254491432, 2.8727830510679486, 2.8728813845795043,
        2.8729761897578996, 2.8730673303785395, 2.8731546702153996,
        2.8732380730438143, 2.8733174026377806, 2.8733925227722557,
        2.8734632972220258, 2.8735295897615645, 2.8735912641652042,
        2.8736481842084065, 2.8737002136652734, 2.8737472163103286,
        2.8737890559188028, 2.8738255962645494, 2.8738567011231777,
        2.8738822342684465, 2.8739020594756828, 2.8739160405190809,
        2.8739240411735327, 2.8739259252137326, 2.8739215564139187,
        2.8739107985493959, 2.87389351539452, 2.8738695707236874,
        2.8738388283117651, 2.8738011519336237, 2.8737564053636286,
        2.8737044523765087, 2.8736451567469667, 2.8735783822495233,
        2.8735039926588084, 2.8734218517500647, 2.8733318232971312,
        2.8732337710749785, 2.8731275588585228, 2.8730130504219358,
        2.8728901095404438, 2.87275859998815, 2.8726183855397798,
        2.8724693299705248, 2.8723112970545266, 2.8721441505663936,
        2.8719677542813238, 2.8717819719731414, 2.8715866674176356,
        2.8713817043871894, 2.8711669466627026, 2.8709422580014716,
        2.8707075022258355, 2.8704625429785851, 2.8702072443997535,
        2.8699414692524394, 2.8696650841009705, 2.8693779457381452,
        2.8690799760639862, 2.8687712586540433, 2.8684519421916326,
        2.8681221655880083, 2.8677820715559621, 2.8674318014309423,
        2.867071497046251, 2.8667013000552446, 2.8663213521762381,
        2.8659317951042467, 2.865532770542492, 2.86512442019123,
        2.8647068857520273, 2.8642803089259687, 2.8638448314139087,
        2.8634005949170049, 2.8629477411366806, 2.8624864117736526,
        2.862016748529304, 2.861538893104715, 2.8610529872008148,
        2.8605591725189323, 2.8600575907600754, 2.8595483836253579,
        2.8590316928159711, 2.8585076600329211, 2.857976426977336,
        2.8574381353504394, 2.856892926853233, 2.8563409431868179,
        2.8557823260524007, 2.8552172171510146, 2.8546457581838531,
        2.8540680908519396, 2.8534843568563648, 2.8528946978984528,
        2.8522992556790614, 2.8516981718994, 2.8510915882606955,
        2.8504796464638971, 2.8498624882101335, 2.8492402552005682,
        2.8486130891364203, 2.8479811317185768, 2.8473445246482885,
        2.8467034096266937, 2.8460579283547931, 2.8454082225337864,
        2.8447544338647703, 2.8440967040488743, 2.8434351747871687,
        2.8427699877807533, 2.8421012847308207, 2.8414292073384155,
        2.8407538973046966, 2.8400754963307238, 2.8393941461176277,
        2.838709988366618, 2.8380231647785914, 2.8373338170549411,
        2.8366420868960045, 2.8359481160050635, 2.8352520460803747,
        2.8345540188243454, 2.8338541759387175, 2.8331526591232832,
        2.8324496100811283, 2.83174517051093, 2.8310394821147957,
        2.83033268659512, 2.829624925651419, 2.8289163409851761,
        2.8282070742976608, 2.8274972672895586, 2.82678706166364,
        2.8260765991182311, 2.8253660213568286, 2.8246554700796951,
        2.823945086987528, 2.8232350137824431, 2.8225253921645996,
        2.821816363835572, 2.8211080704962308, 2.8204006538476958,
        2.8196942555911653, 2.8189890174284908, 2.8182850810596269,
        2.8175825881861702, 2.8168816805089412, 2.8161824997303953,
        2.8154851875495064, 2.814789885667845, 2.8140967357893567,
        2.8134058796112709, 2.8127174588363872, 2.8120316151664251,
        2.8113484903015595, 2.8106682259433478, 2.8099909637921123,
        2.8093168455508839, 2.8086460129190551, 2.8079786075978195,
        2.807314771289112, 2.8066546456922374, 2.8059983725111515,
        2.8053460934464378, 2.8046979501956852, 2.8040540844641084,
        2.8034146379516511, 2.8027797523577216, 2.8021495693855192,
        2.8015242307361166, 2.8009038781087079, 2.8002886532064419,
        2.7996786977298211, 2.7990741533782453, 2.7984751618561665,
        2.7978818648623185, 2.7972944040975, 2.7967129212645432,
        2.7961375580643417, 2.7955684561957885, 2.795005757362409,
        2.7944496032643356, 2.7939001356033129, 2.7933574960801066,
        2.7928218263943916, 2.7922932682494461, 2.791771963346092,
        2.7912580533837272, 2.7907516800657275, 2.7902529850914521,
        2.7897621101624464, 2.7892791969801625, 2.7888043872459458,
        2.7883378226614481, 2.7878796449244772, 2.7874299957399065,
        2.7869890168085165, 2.7865568498289375, 2.7861336365038372,
        2.7857195185349259, 2.785314637622204, 2.7849191354671725,
        2.7845331537706297, 2.7841568342357181, 2.7837903185605608,
        2.7834337484464782, 2.7830872655976928, 2.7827510117125258,
        2.7824251284928381, 2.7821097576375386, 2.7818050408608248,
        2.7815111198212139, 2.7812281363362512, 2.7809562317912087,
        2.7806955487582168, 2.780446226529, 2.7802084134675171, 2.77998223290544,
        2.7797678725154351, 2.77956509128547, 2.779372583679339,
        2.7791886154775964, 2.779011516802163, 2.7788395927407974,
        2.7786711574537541, 2.7785045218232844, 2.7783379979123284,
        2.7781698973622215, 2.7779985319639691, 2.7778222134525472,
        2.7776392535875964, 2.7774479641172127, 2.7772466567930607,
        2.7770336433683043, 2.77680723559156, 2.7765657452168653,
        2.7763074839948478, 2.7760307636753567, 2.7757338960127789,
        2.7754151927573276, 2.7750729656590383, 2.7747055264709863,
        2.7743111869445314, 2.7738882588295763, 2.7734350538802754,
        2.7729498838459041, 2.772431060479239, 2.7718768955296462,
        2.7712857007508336, 2.7706557878944968, 2.7699854687096321,
        2.7692730549495792, 2.7685168583657678, 2.7677151907086066,
        2.7668663637298549, 2.7659686891817232, 2.7650204788157162,
        2.7640200443824781, 2.7629656976333048, 2.7618557503214882,
        2.7606885141952424, 2.7594623010085644, 2.7581754225137853,
        2.7568261904581721, 2.7554129165975727, 2.7539339126819478,
        2.7523874904611731, 2.7507719616893156, 2.7490856381152922,
        2.7473268314931327, 2.745493853572059, 2.7435850161050697,
        2.741598630842935, 2.7395330095364274, 2.7373864639391878,
        2.7351573057999823, 2.7328438468716416, 2.7304443989052909,
        2.7279572736540949, 2.7253807828673087, 2.722713238295623,
        2.7199529516938292, 2.7170982348104156, 2.7141473993987191,
        2.7110987572093848, 2.7079506199927517, 2.7047012995031405,
        2.7013491074879505, 2.6978923557043482, 2.6943293558962953,
        2.6906584198289574, 2.686877859212391, 2.68298598592348,
        2.6789811113604012, 2.6748615482504827, 2.6706256056508595,
        2.6662716027623117, 2.6617978307262375, 2.6572026581273813,
        2.652484254483328, 2.6476421156365513, 2.6426790309805073,
        2.6375991162304966, 2.6324062880361336, 2.6271045404934528,
        2.6216978396341588, 2.6161901616343806, 2.6105854790062333,
        2.604887765582554, 2.5991009947199006, 2.5932291399473391,
        2.5872761747332711, 2.5812460725644719, 2.5751428069222024,
        2.5689703512929944, 2.5627326791589837, 2.5564337640007344,
        2.550077579304149, 2.5436680985528559, 2.5372092952268077,
        2.5307051428139764, 2.5241596147946792, 2.5175766846522816,
        2.5109603258708657, 2.5043145119330812, 2.4976432163216,
        2.4909504125209532, 2.4842400740147461, 2.4775161742847618,
        2.4707826868148404, 2.4640435850883144, 2.4573028425889669,
        2.4505644327985752, 2.4438323292020661, 2.4371105052821376,
        2.4304029345204654, 2.42371359040261, 2.4170464464117538,
        2.4104054760290947, 2.4037946527404115, 2.3972179500263735,
        2.3906793413711895, 2.3841828002598264, 2.3777323001733217,
        2.3713318145959872, 2.3649853170113824, 2.3586967809016648,
        2.3524701797510725, 2.3463094870415149, 2.3402186762585195,
        2.3342017208836614, 2.3282625943993445, 2.3224052702914646,
        2.3166337220409883, 2.3109519231331634, 2.3053638470491684,
        2.2998734672730365, 2.2944847572885587, 2.2892016905789534,
        2.2840282406272849, 2.278968380915388, 2.2740260849294183,
        2.2692053261502667, 2.264510078061214, 2.2599443141473845,
        2.2555120078905606, 2.2512171327736885, 2.2470636622826712,
        2.2430555698953571, 2.2391968291043765, 2.2354914133730204,
        2.2319432962357904, 2.2285564510341915, 2.2253348516490141,
        2.22228247045646, 2.219403284007361, 2.2167012572973026,
        2.214180387291028, 2.2118445827226068, 2.2096979791227227,
        2.2077432009182454, 2.2059791201337338, 2.2044030976854203,
        2.2030127212957891, 2.2018054904474136, 2.2007789365922963,
        2.1999305796349775, 2.1992579436476531, 2.1987585511995178,
        2.1984299254014985, 2.1982695891678832, 2.1982750654863166,
        2.1984438773156589, 2.1987735476276713, 2.1992615993871634,
        2.1999055555621085, 2.2007029391198318, 2.20165127302724,
        2.2027480802526314, 2.2039908837625424, 2.2053772065236736,
        2.2069045715052118, 2.2085705016711223, 2.2103725199925792,
        2.2123081494348042, 2.2143749129650616, 2.2165703335510454,
        2.2188919341591542, 2.2213372377592586, 2.2239037673152531,
        2.2265890457964934, 2.2293905961695093, 2.2323059414015529,
        2.2353326044617234, 2.2384681083146321, 2.2417099759277015,
        2.2450557302709693, 2.248502894310044, 2.2520489910108852,
        2.2556915433422895, 2.2594280742719643, 2.2632561067660109,
        2.2671731637918571, 2.2711767683190338, 2.2752644433103266,
        2.27943371173691, 2.2836820965656797, 2.2880071207614421,
        2.292406307293855, 2.2968771791294795, 2.3014172592359285,
        2.30602407057936, 2.3106951361277805, 2.3154279788490806,
        2.3202201217101273, 2.3250690876773716, 2.3299723997190163,
        2.3349275808015739, 2.3399321538939528, 2.3449836419620405,
        2.3500795679731414, 2.3552174548949769, 2.3603948256939868,
        2.365609203339373, 2.3708581107967985, 2.3761390710338879,
        2.3814496070169966, 2.3867872417162244, 2.3921494980960452,
        2.3975338991243507, 2.4029379677699567, 2.4083592269974465,
        2.4137951997767875, 2.4192434090743458, 2.4247013778562887,
        2.430166629091457, 2.4356366857467489, 2.4411090707881735,
        2.446581307184831, 2.4520509179037142, 2.4575154259097922,
        2.4629723541737629, 2.4684192256614113, 2.4738535633384373,
        2.4792728901751282, 2.48467472913719, 2.4900566031900864,
        2.4954160353048556, 2.5007505484468973, 2.5060576655825972,
        2.5113349096811706, 2.5165798037082183, 2.5217898706315789,
        2.5269626334190369, 2.5320956150379197, 2.5371863384533375,
        2.5422323266364981, 2.5472311025506604, 2.5521801891661973,
        2.5570771094496014, 2.5619193863663114, 2.5667045428859803,
        2.571430101973645, 2.57609358660069, 2.5806925197290105,
        2.5852244243288696, 2.5896868233680395, 2.5940772398115959,
        2.5983931966295444, 2.6026322167875904, 2.6067918232522218,
        2.6108695389926195, 2.6148628869751573, 2.6187693901660016,
        2.6225865715353258, 2.6263119540472997, 2.6299430606708083,
        2.6334774143733903, 2.6369125381214613, 2.640245954882579,
        2.6434751876240585, 2.646597759313754, 2.6496111929182806,
        2.6525130114049422, 2.6553007377416158, 2.6579718948939388,
        2.6605240058305757, 2.6629545935202663, 2.6652611809270152,
        2.6674412910204874, 2.6694924467675092, 2.6714121711338827,
        2.6731979870893174, 2.6748474175998069, 2.6763579856320812,
        2.6777272141549484, 2.6789526261336878, 2.6800317445366244,
        2.6809620923326896, 2.6817411924873995, 2.68236656796641,
        2.6828357417390838, 2.6831462367745487, 2.6832955760379473,
        2.6832812824938528, 2.6831008791137569, 2.6827518888640123,
        2.6822318347108682, 2.6815382396213172, 2.6806686265692141,
        2.679620518500037, 2.678391438437282, 2.6769789091987621,
        2.6753804541677506, 2.6735935951499017, 2.6716158583361507,
        2.6694447577693219, 2.6670778411159022, 2.6645125630027291,
        2.6617466348298096, 2.6587771079993487, 2.6556054312607418,
        2.6522439729494556, 2.6487094987550162, 2.6450181143662963,
        2.641186182238866, 2.6372299717968941, 2.6331657860841626,
        2.6290099159973987, 2.6247786568215314, 2.6204883022560623,
        2.6161551465734978, 2.6117954838379482, 2.6074256081921781,
        2.6030618137479427, 2.59872039462696, 2.5944176449493819,
        2.5901698588383106, 2.5859933304103722, 2.5819043537874786,
        2.577919223091234, 2.5740542324382987, 2.5703256759521786,
        2.5667498477540325, 2.5633430419612013, 2.5601215526942522,
        2.5571016740764816, 2.55429970022598, 2.5517319252624628,
        2.5494146433085572, 2.5473641484813649, 2.54559673490317,
        2.5441286966976175, 2.5429763279781792, 2.5421559228696111,
        2.5416837754909194, 2.5415761799626666, 2.541849430406995,
        2.5425198209412012, 2.5436036456869142, 2.5451171987641903,
        2.5470767742936906, 2.549498666396131, 2.5523991691915278,
        2.555794576800094, 2.5597011833413523, 2.5641352829372792,
        2.5691131697068563, 2.5746511377705432, 2.5807654812506597,
        2.5874724942635083, 2.5947884709347195, 2.602729705378783,
        2.6113124917192567, 2.6205531240796547, 2.63046789657219,
        2.6410731033233765, 2.6523850384535086, 2.6644199960799622,
        2.6771942703235312, 2.6907241553082524, 2.705025945149881,
        2.7201159339693213, 2.7360104158908496, 2.7527256850307253,
        2.7702780355095378, 2.7886837614510847, 2.8079591569699307,
        2.8281205161959559, 2.8491841332249845, 2.8711663022538683,
        2.8940833171774116, 2.9179514727426716, 2.9427870613492022,
        2.9686063818611439, 2.9954257152808479, 3.0232613920337381,
        3.0521296057589, 3.082046928624175, 3.1130288880484067,
        3.1450936968935763, 3.1782416758203986, 3.2124287152273125,
        3.2475928133066803, 3.2836746536981867, 3.3206138752910248,
        3.35835049550544, 3.3968243949698906, 3.4359755037410631,
        3.4757437340144044, 3.5160690044434593, 3.5568912313410905,
        3.5981503318711563, 3.6397862228903364, 3.6817388213630742,
        3.7239480442159616, 3.7663538083902681, 3.8088960308226141,
        3.85151462845007, 3.8941495182092445, 3.9367406170347086,
        3.979227841866932, 4.0215511096440055, 4.063650337298065,
        4.1054654417696179, 4.146936339996369, 4.1880029489105164,
        4.2286051854572975, 4.2686829665610695, 4.30817620918801,
        4.3470248301882286, 4.3851687467512885, 4.42254787511007,
        4.4591021341466623, 4.4947714354197625, 4.5294957107528768,
        4.5632148358777656, 4.5958688417673512, 4.6273973298234505,
        4.6577410870327158, 4.6868378529581118, 4.7146456711271005,
        4.7411730042704283, 4.7664486190963427, 4.7904982348793395,
        4.813348756475289, 4.8350266591802882, 4.8555585735230391,
        4.8749710739417926, 4.8932907551431217, 4.9105442045087706,
        4.9267580120642807, 4.941958766886315, 4.9561730583900543,
        4.969427475867878, 4.9817486086596814, 4.993163046084768,
        5.0036973774710658, 5.013378192147413, 5.0222320794359039,
        5.0302856286637576, 5.037565429157052, 5.0440980702411471,
        5.0499101412429521, 5.0550282314899961, 5.0594789303044516,
        5.0632888270128351, 5.0664845109456227, 5.0690925714257347,
        5.07113959777154, 5.0726521793359058, 5.0736569053765832,
        5.0741803654021416, 5.07424914822955, 5.0738898446055636,
        5.0731290399167062, 5.0719933343849766, 5.0705092871910216,
        5.0687035710759156, 5.0666025453503174, 5.0642333749754194,
        5.0616178571207495, 5.0587644595516617, 5.0556762822306878,
        5.052357230779446, 5.0488108973883641, 5.0450409878077682,
        5.0410511667447455, 5.0368451137444614, 5.0324265029894475,
        5.0277990105931512, 5.0229663119797738, 5.017932082819482,
        5.0126999986907554, 5.007273735209159, 5.0016569679737506,
        4.9958533725914611, 4.989866624668168, 4.9837003998076943,
        4.9773583736101159, 4.9708442216857254, 4.9641616196377019,
        4.9573142430676231, 4.9503057675852888, 4.9431398687902215,
        4.9358202222874867, 4.928350503684725, 4.92073438858511,
        4.9129755525926111, 4.9050776713119824, 4.8970444203471741,
        4.8888794753029634, 4.8805865117875564, 4.8721692053989507,
        4.8636312317435637, 4.8549762664312954, 4.8462079850611337,
        4.8373300632392269, 4.828346176569144, 4.8192600006581889,
        4.8100752111069855, 4.8007954835222568, 4.7914244935110153,
        4.7819659166729194, 4.77242342861642, 4.7628007049437331,
        4.7531014212579681, 4.7433292531691187, 4.7334878762774766,
        4.7235809661867156, 4.7136121985060324, 4.7035852488349406,
        4.6935037927811276, 4.6833715059491059, 4.6731920639421176,
        4.6629691423624635, 4.6527064168190009, 4.6424075629169286,
        4.6320762562548952, 4.6217161724422855, 4.6113309870817467,
        4.60092437578003, 4.5905000141392032, 4.580061577763674,
        4.5696127422593582, 4.5591571832330455, 4.5486985762828773,
        4.538240597017082, 4.5277869210442905, 4.5173412239618536,
        4.5069071813772315, 4.4964884688960982, 4.486088762116732,
        4.4757117366726584, 4.4653610680817906, 4.4550404321777517,
        4.4447535039389532, 4.4345039607034256, 4.4242954732806981,
        4.4141317305425769, 4.4040163715052811, 4.3939531633421369,
        4.3839450193261822, 4.3739927322784729, 4.364096241108733,
        4.3542556129020547, 4.3444708648729584, 4.3347420323001975,
        4.3250691439417395, 4.3154522309108385, 4.3058913234619967,
        4.2963864521731345, 4.28693764749544, 4.2775449399239767,
        4.2682083599460139, 4.2589279380496707, 4.249703704721032,
        4.2405356904434761, 4.2314239257074675, 4.22236844099859,
        4.213369266804265, 4.2044264336102222, 4.1955399719011481,
        4.1867099121700315, 4.1779362848975552, 4.1692191205716318,
        4.1605584496808987, 4.1519543027135128, 4.143406710150769,
        4.1349157024822487, 4.1264813101978683, 4.1181035637782433,
        4.1097824937169314, 4.10151813049438, 4.093310504600904,
        4.085159646522567, 4.0770655867444665, 4.0690283557558038,
        4.0610479840430314, 4.0531245020921833, 4.0452579403875157,
        4.037448329421192, 4.0296956996770987, 4.022000081637775,
        4.0143615057984992, 4.00678000264049, 3.9992556026485446,
        3.991788336314924, 3.9843782341250087, 3.9770253265636808,
        3.9697296441168159, 3.9624912172746449, 3.9553100765199911,
        3.9481862523435582, 3.9411197752288856, 3.9341106756648014,
        3.9271589841356138, 3.9202647311317045, 3.913427947136495,
        3.9066486626368726, 3.8999269081253574, 3.8932627140785527,
        3.8866561109896578, 3.8801071293474503, 3.8736157996334413,
        3.8671821523370524, 3.86080621794441, 3.854488026942771,
        3.8482276098179704, 3.8420249970571212, 3.8358802191473322,
        3.8297933065739591, 3.8237642898267987, 3.8177931993954273,
        3.8118800657396918, 3.8060249194255715, 3.8002277907355078,
        3.7944887107061263, 3.7888077083063765, 3.7831848182261063,
        3.7776200593561753, 3.772113491209764, 3.766664902622221,
        3.7612734102869365, 3.7559378602200915, 3.7506571390615764,
        3.7454301176465772, 3.7402556725373506, 3.7351326782287759,
        3.7300600099578718, 3.7250365426992822, 3.7200611515194137,
        3.7151327114513597, 3.7102500975406416, 3.7054121848282215,
        3.7006178483557362, 3.69586596316718, 3.69115540430331,
        3.6864850468079915, 3.681853765723226, 3.6772604360864722,
        3.6727039329473632, 3.6681831313454127, 3.6636969063212592,
        3.659244132916561, 3.6548236861761252, 3.6504344411411411,
        3.6460752728555819, 3.6417450563593192, 3.6374426666939415,
        3.6331669789036836, 3.6289168680289809, 3.62469120911525,
        3.6204888772014407, 3.6163087473321411, 3.6121496945468246,
        3.6080105938896354, 3.6038903204056636, 3.5997877491312749,
        3.59570175511092, 3.5916312133870125, 3.5875749990048793,
        3.5835319870037163, 3.5795010524253925, 3.5754810703131321,
        3.5714709157067088, 3.5674694636534183, 3.5634755891917638,
        3.5594881673643122, 3.5555060732166, 3.5515281817838922,
        3.5475533681151172, 3.5435805072494966, 3.5396084742305933,
        3.5356361441001951, 3.5316623918987857, 3.5276860926709146,
        3.5237061214578, 3.5197213533033391, 3.515730663246559,
        3.5117329263311206, 3.5077270176000863, 3.5037118120972943,
        3.4996861848590246, 3.4956490109366047, 3.491599165362218,
        3.4875355231838365, 3.4834569594441982, 3.4793623491837469,
        3.4752505674449217, 3.4711204892688348, 3.4669709897041003,
        3.462800943782657, 3.4586092265542088, 3.4543947130580341,
        3.4501562783370545, 3.4458927974348494, 3.4416031453908311,
        3.4372861972499864, 3.4329408280543503, 3.4285659128411803,
        3.4241603266598282, 3.4197229445518476, 3.4152526415522466,
        3.4107482927093566, 3.4062087730660804, 3.4016329576619957,
        3.3970197215384061, 3.3923679397416993, 3.3876764873099408,
        3.382944239288598, 3.3781700707181974, 3.3733528566372875,
        3.3684914720963479, 3.3635847921347066, 3.3586316917887662,
        3.3536310461049768, 3.348581730130173, 3.3434826188986086,
        3.3383325874578182, 3.3331305108471843, 3.3278752641109532,
        3.3225657222896143, 3.3172007604282188, 3.311779253566312,
        3.3063000767449351, 3.300762105012248, 3.2951642133998389,
        3.28950527696407, 3.283784170739541, 3.2779997697600396,
        3.2721509490831209, 3.2662365837456693, 3.2602555487870171,
        3.2542067192488604, 3.2480889701775664, 3.2419011766145367,
        3.2356422135989287, 3.22931095617544, 3.2229062793856489,
        3.2164270582725445, 3.2098721678781925, 3.2032404832434587,
        3.1965308794127547, 3.1897422314266977, 3.1828734143263797,
        3.1759233031566785, 3.1688907729596254, 3.1617746987766484,
        3.1545739556471935, 3.1472874186213775, 3.1399139627317894,
        3.1324524630254054, 3.1249017945478403, 3.1172608323356545,
        3.10952845143212, 3.1017035268818249, 3.0937849337250958,
        3.085771547005141, 3.0776622417662765, 3.0694558930427354,
        3.06115137588702, 3.0527475653360585, 3.0442433364313333,
        3.0356375642183164, 3.0269291237344635, 3.0181168900272963,
        3.0091997381382352, 3.0001765431041356, 2.9910461799757617,
        2.9818075237842496, 2.97245944958995, 2.9630008323959491,
        2.9534305473456133, 2.9437474692056496, 2.9339504737848192,
        2.9240384339825338, 2.914010230795145, 2.9038647227900318,
        2.893600830578082, 2.8832173031196868, 2.8727133631411679,
        2.8620870155417197, 2.8513443793074176, 2.8405117224417142,
        2.8296234270481886, 2.8187126573831027, 2.8078130514848763,
        2.7969580757367964, 2.7861812585615731, 2.7755161059541327,
        2.764996132016655, 2.7546548479263464, 2.7445257659103715,
        2.7346423978217529, 2.72503825564598, 2.7157468513214678,
        2.7068016968054156, 2.6982363040447859, 2.6900841849908628,
        2.682378851596, 2.6751538158122639, 2.6684425895863946,
        2.6622786848750208, 2.6566956136275888, 2.6517268877875582,
        2.6474060193161271, 2.64376652015974, 2.6408419022677543,
        2.6386656775946724, 2.6372713580881681, 2.636692455708316,
        2.6369624823723465, 2.6381149501291921, 2.6401833706702558,
        2.643201256643009, 2.6472021180813496, 2.6522194722485186,
        2.658286816383062, 2.6654377031552117, 2.6737055318302474,
        2.6831241261698859, 2.6937261383215572, 2.7055472319800398,
        2.7186030059409294, 2.7328592335389268, 2.7482616232026236,
        2.76475889491158, 2.78229859703438, 2.8008287024284995,
        2.8202970305541624, 2.8406514563018677, 2.8618398345309655,
        2.8838100273372964, 2.9065098942060903, 2.92988729555737,
        2.9538900914867918, 2.9784661422000984, 3.0035633078558228,
        3.0291294486394196, 3.0551124247294172, 3.0814600962926351,
        3.1081203235154189, 3.1350409665692278, 3.1621698856316507,
        3.1894549408782438, 3.216843992480372, 3.2442849006250021,
        3.2717255254819246, 3.2991137272278457, 3.3263973660393864,
        3.353524302088684, 3.3804423955766403, 3.4070995066076466,
        3.433443495544799, 3.4594222220540538, 3.4849835477313897,
        3.5100753288202786, 3.5346454363782827, 3.5586417004679745,
        3.5820120646058959, 3.6047041583628792, 3.62666647778804,
        3.6478452917280149, 3.6682017081479, 3.6877336838131809,
        3.7064540146178677, 3.7243732692422054, 3.7415028828453778,
        3.7578539766499603, 3.7734377853279764, 3.7882655025514684,
        3.8023483368104549, 3.8156974912445261, 3.8283241709202773,
        3.8402395802108238, 3.8514549237444422, 3.8619814060460129,
        3.8718302316885014, 3.8810126052286322, 3.8895397312116833,
        3.897422814206299, 3.9046730587700114, 3.91130166945385,
        3.9173198508150664, 3.9227388074163358, 3.9275697438134745,
        3.9318238645605952, 3.9355123742194968, 3.938646477343053,
        3.9412373784896162, 3.9432962822187703, 3.9448343930859844,
        3.9458629156493235, 3.9463930544677037, 3.946436014091053,
        3.9460029990845169, 3.9451052140056424, 3.943753863405226,
        3.9419601518438134, 3.9397352838837745, 3.9370904640751809,
        3.9340368969750639, 3.9305857871503456, 3.9267483391439946,
        3.922535757522732, 3.9179592468472451, 3.9130300116631229,
        3.9077592565371093, 3.9021581860236818, 3.8962380046772171,
        3.890009917058491, 3.8834851277269857, 3.8766748412336209,
        3.8695902621393463, 3.862242595002646, 3.8546430443778119,
        3.8468028148240072, 3.8387331108999985, 3.8304451371594053,
        3.8219500981590535, 3.8132591984621609, 3.804383642620381,
        3.7953346351939743, 3.7861233807387822, 3.7767610838116039,
        3.7672589489692507, 3.7576281807756384, 3.747879983780376,
        3.7380255625372043, 3.7280761216181548, 3.7180428655663658,
        3.7079369989471105, 3.6977697263152756, 3.6875522522231008,
        3.6772957812386458, 3.6670115179053964, 3.65671066681278,
        3.6464044324234006, 3.6361040195585073, 3.625820632044106,
        3.61556547646135, 3.6053497517978, 3.5951846765359372,
        3.5850813393185406, 3.5750505063681532, 3.5651028140650021,
        3.55524891828347, 3.5454994673146643, 3.5358651121888927,
        3.5263565029515171, 3.5169842900057819, 3.5077591236269132,
        3.4986916541290261, 3.48979253181762, 3.4810724069995427,
        3.4725419299823042, 3.4642117510773325, 3.4560925205845487,
        3.4481948888119036, 3.4405295060706793, 3.4331070226653444,
        3.4259380889018285, 3.4190333550895722, 3.41240347153324,
        3.4060590885395605, 3.4000108564213947, 3.394269425477364,
        3.388845446016068, 3.38374956835142, 3.3789924427851292,
        3.3745847196188063, 3.3705370491721762, 3.3668600817447802,
        3.3635644676366931, 3.3606608571702972, 3.3581599006412444,
        3.356072248359911, 3.3544085506356174, 3.3531794577700929,
        3.3523956200738989, 3.3520676878532027, 3.3522063114140837,
        3.3528221410685259, 3.3539258271161048, 3.3555280198657624,
        3.3576393696293616, 3.3602705267091575, 3.3634321414131625,
        3.3671348640482388, 3.3713893449230428, 3.3762062343416233,
        3.381596182611863, 3.3875698400435752, 3.3941378569432783,
        3.4013108836110231, 3.4090995703611755, 3.4175145675019332,
        3.4265665253324094, 3.4362660941667529, 3.4466239243088226,
        3.4576506660640507, 3.4693569697442883, 3.4817534856522721,
        3.4948508640971929, 3.5086597553849344, 3.5231908098201705,
        3.5384546777144372, 3.5544620093755648, 3.5712234551058235,
        3.5887496652086912, 3.6070512900082048, 3.6261389797774095,
        3.6460233848992174, 3.6667151554852633, 3.6882249423773383,
        3.7105633943847205, 3.7337411659930888, 3.7577688959318065,
        3.7826572665305442, 3.8084168394931113, 3.835058510293611,
        3.8625922532050843, 3.8910304103846163, 3.9203695475034475,
        3.9505670537796829, 3.981564541929282, 4.0133059925640122,
        4.0457344650984082, 4.0787933527123039, 4.11242592795387,
        4.14657550698172, 4.18118539018727, 4.2161988836533917,
        4.2515592914134235, 4.2872099182438621, 4.3230940686484933,
        4.3591550472223961, 4.3953361585436648, 4.4315807071865372,
        4.4678319977225529, 4.5040333347356718, 4.5401280227973624,
        4.5760593664871072, 4.6117706703780215, 4.6472052390470839,
        4.6823063770714413, 4.7170173890306373, 4.75128157949625,
        4.7850422530431826, 4.8182427142583419, 4.8508262676955036,
        4.8827362179775244, 4.9139158695802374, 4.9443085273369389,
        4.973857495098847, 5.0025060794754683, 5.0301975794134819,
        5.0568753150580115, 5.0824825479033109, 5.1069627017605823,
        5.1302587512626623, 5.1523149107899346, 5.1730722080472074,
        5.1924929024469968, 5.210591976313399, 5.2274056436685647,
        5.2429669318807379, 5.2573101080517057, 5.2704689900927013,
        5.2824775582550565, 5.2933697341317263, 5.3031794605024647,
        5.3119406724971361, 5.3196873080059612, 5.326453303926761,
        5.332272597508255, 5.3371791258858376, 5.34120682622304,
        5.3443896356706428, 5.3467614914026518, 5.3483563305607786,
        5.349208090321576, 5.34935070783421, 5.3488181202570182,
        5.347644264751052, 5.3458630784803729, 5.3435084985987666,
        5.3406144622605067, 5.3372149066420613, 5.3333437688836831,
        5.3290349861559783, 5.3243224956088939, 5.3192402344273111,
        5.3138221396787708, 5.3081021487951725, 5.3021141981771294,
        5.2958922270909836, 5.289470166861002, 5.282881970777364,
        5.2761615313691133, 5.2693429092992243, 5.2624597011534266,
        5.2555466964187048, 5.2486307367576046, 5.2417189277315774,
        5.2348104271113405, 5.2279055855486813, 5.2210042895987927,
        5.214106593983062, 5.207212492655116, 5.2003220015116289,
        5.1934351285390266, 5.1865518845647243, 5.1796722793933068,
        5.1727963232024443, 5.1659240260379713, 5.1590553979857683,
        5.15219044912558, 5.1453291895366791, 5.1384716292866068,
        5.1316177784683674, 5.1247676471484986, 5.1179212454076275,
        5.11107858332709, 5.1042396709849909, 5.0974045184550913,
        5.090573135821141, 5.0837455331619807, 5.0769217205473707,
        5.0701017080606361, 5.06328550578071, 5.0564731237827578,
        5.0496645721470594, 5.0428598609596342, 5.0360590002812611,
        5.0292620002039552, 5.022468870804337, 5.0156796221467994,
        5.0088942643297871, 5.00211280741426, 4.9953352614896227,
        4.98856163663363, 4.9817919429222179, 4.9750261904279585,
        4.9682643892306535, 4.9615065494217472, 4.9547526810565969,
        4.948002794229323, 4.9412568990251584, 4.93451500549876,
        4.9277771237458046, 4.9210432638438517, 4.914313435856779,
        4.9075876498798356, 4.9008659159809813, 4.8941482442439614,
        4.8874346447436583, 4.8807251275590247, 4.8740197027712355,
        4.867318380449035, 4.8606211706811564, 4.8539280835424927,
        4.8472391291054988, 4.840554317456486, 4.83387365867024,
        4.8271971628288046, 4.8205248399987513, 4.8138567002649388,
        4.8071927537147907, 4.8005330104109447, 4.79387748044059,
        4.7872261738820612, 4.7805791008094349, 4.7739362713050522,
        4.7672976954420481, 4.7606633833015675, 4.7540333449591863,
        4.7474075905010054, 4.7407861299959446, 4.7341689735251355,
        4.727556131173964, 4.7209476130067269, 4.7143434291064166,
        4.7077435895569755, 4.701148104439107, 4.6945569838137029,
        4.6879702377786607, 4.68138787640269, 4.6748099097585225,
        4.6682363479364666, 4.6616672010104052, 4.6551024790501216,
        4.6485421921428252, 4.6419863503692556, 4.635434963795757,
        4.6288880425131893, 4.6223455965928082, 4.6158076361105,
        4.6092741711496092, 4.6027452117856651, 4.5962207680940317,
        4.5897008501626084, 4.5831854680649879, 4.5766746318718621,
        4.5701683516663918, 4.5636666375302477, 4.5571694995399765,
        4.5506769477697153, 4.5441889923056227, 4.5377056432100451,
        4.5312269105824132, 4.5247528044894478, 4.518283335001775,
        4.5118185122115388, 4.5053583461877444, 4.4989028470153416,
        4.4924520247692872, 4.4860058895280348, 4.4795644513631,
        4.473127720361612, 4.4666957066055177, 4.460268420157167,
        4.45384587110372, 4.4474280695323625, 4.4410150255052665,
        4.4346067491094878, 4.4282032504223094, 4.4218045395189085,
        4.41541062647797, 4.409021521376606, 4.4026372343032412,
        4.3962577753236385, 4.3898831545205113, 4.3835133819710208,
        4.377148467756137, 4.3707884219538133, 4.364433254632222,
        4.3580829758884914, 4.3517375957850355, 4.3453971244026235,
        4.3390615718229393, 4.3327309481234106, 4.3264052633857562,
        4.32008452767865, 4.3137687510854024, 4.30745794368661,
        4.3011521155602273, 4.2948512767748213, 4.2885554374248569,
        4.2822646075765771, 4.275978797308408, 4.2696980167046235,
        4.2634222758368816, 4.257151584792612, 4.2508859536352572,
        4.2446253924573742, 4.2383699113230113, 4.2321195203523638,
        4.2258742294981984, 4.2196340491628082, 4.2133989885554337,
        4.207169060165076, 4.200944267384096, 4.194724637062305,
        4.1885100396031021, 4.1822999569524626, 4.17609371458085,
        4.1698906614684423, 4.1636901374458057, 4.1574914856468315,
        4.1512940480354468, 4.1450971669707348, 4.1389001846909617,
        4.1327024434687436, 4.1265032855583463, 4.1203020532406711,
        4.1140980887603034, 4.1078907343917086, 4.1016793324016936,
        4.0954632250460339, 4.0892417545950375, 4.0830142633149453,
        4.0767800934622223, 4.0705385873054922, 4.064289087107765,
        4.0580309351333694, 4.0517634736469894, 4.0454860449192571,
        4.0391979911984057, 4.0328986547557992, 4.0265873778711709,
        4.0202635027801188, 4.0139263717668756, 4.0075753270988717,
        4.0012097110189124, 3.9948288658107387, 3.9884321337248489,
        3.9820188570419592, 3.9755883780151309, 3.9691400389004645,
        3.962673181980052, 3.956187149505376, 3.9496812837469655,
        3.9431549269646138, 3.9366074214262103, 3.9300381093923877,
        3.9234463331275395, 3.9168314349002484, 3.9101927569728816,
        3.9035296416044658, 3.8968414310638342, 3.8901274676158053,
        3.8833870935214136, 3.8766196510481254, 3.8698244824547947,
        3.8630009300092829, 3.856148335981834, 3.8492660426247305,
        3.8423533922053168, 3.8354097269969372, 3.82843438925145,
        3.8214267212406168, 3.8143860652245682, 3.8073117634692215,
        3.8002031582427898, 3.793059591798615, 3.7858804064112155,
        3.7786649443397211, 3.7714125478508564, 3.7641225592064633,
        3.7567943206714993, 3.74942717451074, 3.7420204629871794,
        3.7345735283682671, 3.7270857129045423, 3.7195563588969431,
        3.7119848085179523, 3.7043704042819288, 3.6967124877536297,
        3.68901040313547, 3.6812634873234771, 3.6734710974452822,
        3.6656325346503627, 3.6577472545507193, 3.6498143157381291,
        3.6418354221106073, 3.6338188464320238, 3.6257755067967752,
        3.6177159242455437, 3.6096507743017017, 3.6015906765051846,
        3.593546270631256, 3.5855281891492381, 3.5775470671553475,
        3.5696135388075665, 3.5617382386086827, 3.5539318009099468,
        3.5462048601427711, 3.5385680507107788, 3.531032007004451,
        3.52360736344462, 3.5163047544197985, 3.5091348143486352,
        3.5021081776298231, 3.495235478663957, 3.4885273518580138,
        3.4819944316245732, 3.475647352355244, 3.4694967484623809,
        3.4635532543481786, 3.4578275044179039, 3.4523301330780325,
        3.4470717747247557, 3.4420630637741443, 3.4373146346200567,
        3.432837121675135, 3.4286411593418542, 3.4247373820204627,
        3.4211364241187647, 3.4178489200430553, 3.4148855041947166,
        3.4122568109788518, 3.4099734747993247, 3.40804613006277,
        3.4064854111702108, 3.4053019525326871, 3.4045063885478988,
        3.40410935361877, 3.4041214821640717, 3.4045534085699951,
        3.4054157672494161, 3.4067191926095926, 3.4084743190521518,
        3.4106917809763933, 3.4133822127973432, 3.4165562489113546,
        3.4202245237239022, 3.424397671643681, 3.4290863270718033,
        3.4343011244112804, 3.4400526980715718, 3.4463516824557248,
        3.4532087119578678, 3.4606344210034008, 3.4686394439790318,
        3.4772344152937738, 3.4864299693549397, 3.4962367405693771,
        3.5066653633306983, 3.5177264720524857, 3.5294307011397006,
        3.5417886849914089, 3.5548110580197338, 3.5685084546088217,
        3.5828915092279785, 3.5979708560873918, 3.6137571301373681,
        3.6302609642840555, 3.647492997047991, 3.6654638514723721,
        3.684184193404938, 3.7036645701984958, 3.723915857166066,
        3.7449480244205824, 3.7667733688051546, 3.789388685078932,
        3.8127522730281744, 3.8368069303634171, 3.8614977815141383,
        3.8867690457163282, 3.9125652701737308, 3.938830883573889,
        3.9655103574187169, 3.9925481477585434, 4.0198887162093735,
        4.047476522383902, 4.0752560266160049, 4.1031716889834007,
        4.1311679696482377, 4.1591893287520305, 4.1871802264397644,
        4.2150851228469488, 4.2428484781266214, 4.270414752416702,
        4.2977284058638023, 4.3247338986077972, 4.3513756907936765,
        4.3775982425663793, 4.4033460140720475, 4.4285634654420232,
        4.45319505683279, 4.4771852483817431, 4.50047850023864,
        4.5230192725288667, 4.5447520254369156, 4.5656212190007146,
        4.5855713136649836, 4.6045467687185093, 4.6224920466858377,
        4.6393516011210938, 4.65506991038925, 4.6695913842165186,
        4.6828606222781666, 4.6948216986112827, 4.7054201380419665,
        4.7145977362141238, 4.7223211351581753, 4.7286186758637738,
        4.7335435457131458, 4.7371452028835614, 4.739474556373688,
        4.7405819895227115, 4.740518075617131, 4.7393333193303828,
        4.7370782501066548, 4.7338033884529578, 4.7295592581047625,
        4.7243963816191075, 4.7183652820023125, 4.7115164820691113,
        4.7039005047283622, 4.6955678728386552, 4.6865691092893522,
        4.6769547369464464, 4.6667752786897925, 4.6560812573962789,
        4.6449231959469346, 4.6333516172099607, 4.6214170440665079,
        4.6091699993954682, 4.5966610060648128, 4.583940586960769,
        4.5710592649473831, 4.5580675629261744, 4.5450160037281568,
        4.5319551103121158, 4.5189354053647568, 4.5060074122646032,
        4.4932216525210249, 4.4806286527698154, 4.4682789254882183,
        4.456223022357908, 4.4445113864987142, 4.4331947615020662,
        4.4223230616444758, 4.4119483328758609, 4.4021084184616486,
        4.3928058932564449, 4.3840291294028129, 4.3757686307366175,
        4.3680140717754519, 4.3607554275185629, 4.35398256437805,
        4.3476853879931578, 4.3418537898501146, 4.3364776665195475,
        4.3315469127533888, 4.3270514239640052, 4.32298109531876,
        4.319325822069878, 4.3160754994425288, 4.313220022678351,
        4.310749286994997, 4.3086531876449046, 4.3069216198414795,
        4.3055444788288382, 4.3045116598359838, 4.3038130580986209,
        4.3034385688441574, 4.303378087305858, 4.3036215087207941,
        4.30415872831483, 4.3049796413296955, 4.3060741429911218,
        4.3074321285333435, 4.3090434931884722, 4.3108981321850912,
        4.3129859407693836, 4.3152968141593933, 4.3178206475879959,
        4.3205473363081905, 4.3234667755190133, 4.326568860482161,
        4.3298434864189286, 4.3332805485554511, 4.3368699421414183,
        4.3406015623876808, 4.3444653045440145, 4.3484510638399581,
        4.3525487354994379, 4.3567482147642336, 4.3610393968627816,
        4.3654121770286025, 4.3698564504928834, 4.3743621124896617,
        4.3789190582514053, 4.3835171830111941, 4.38814638200417,
        4.3927965504500373, 4.3974575836006666, 4.40211937667476,
        4.4067718249107042, 4.41140482353724, 4.4160082677892758,
        4.4205720529019237, 4.4250860741006939, 4.42954022662953,
        4.4339244057042251, 4.43822850657682, 4.4424424244626, 4.446556054605125,
        4.4505592922350834, 4.4544420325836516, 4.4581941708799109,
        4.46180560236281, 4.4652662222509907, 4.468565925822527,
        4.4716946081862163, 4.4746421649128587, 4.4773984902781434,
        4.4799534821940474, 4.4822970284662205, 4.484419044874338,
        4.4863093698102059, 4.4879580552344525, 4.4893546041301038,
        4.4904921773025643, 4.491373018728928, 4.492003030236253,
        4.4923875646463323, 4.4925321883521256, 4.4924423903785176,
        4.4921236877048791, 4.4915815872090787, 4.4908215994194185,
        4.4898492335415039, 4.4886699992638111, 4.4872894060971618,
        4.485712963626816, 4.48394618139623, 4.4819945689663072,
        4.4798636359012818, 4.4775588917570772, 4.4750858460910168,
        4.4724500084710126, 4.4696568884461252, 4.466711995579141,
        4.4636208394306518, 4.4603889295585706, 4.4570217755242725,
        4.4535248868828106, 4.4499037731999138, 4.4461639440262815,
        4.4423109089224866, 4.4383501774581395, 4.4342872591797731,
        4.4301276636515254, 4.4258769004416241, 4.42154047908712,
        4.4171239091686427, 4.4126327002357089, 4.4080723618475508,
        4.4034484035642416, 4.3987663349522608, 4.3940316655568363,
        4.389249904945701, 4.3844265626804635, 4.3795671483110254,
        4.3746771714095587, 4.3697621415200718, 4.3648275682129816,
        4.3598789610474817, 4.3549218295781751, 4.3499616833566677,
        4.3450040319632759, 4.3400543849417854, 4.335118251847196,
        4.3302011422541362, 4.3253085657111408, 4.3204460317738977,
        4.3156190500183165, 4.3108331299868361, 4.30609378124474,
        4.3014065133542907, 4.2967768358650877, 4.2922102583516564,
        4.2877122903562963, 4.2832884414515329, 4.2789442211886186,
        4.2746851391310274, 4.2705167048376262, 4.26644442786109,
        4.2624738177693837, 4.2586103841212557, 4.2548596364641895,
        4.2512270843789253, 4.2477182373897122, 4.2443386051326755,
        4.2410936969490445, 4.2379890230339941, 4.2350300911409056,
        4.2322224158676187, 4.22957149281954, 4.22708287014596,
        4.224761951003595, 4.2226145112519911, 4.2206438434727556,
        4.2188470737811636, 4.2172188449855685, 4.2157541726193681,
        4.2144479272293456, 4.21329503188269, 4.2122903906657339,
        4.2114289145296695, 4.2107055119398682, 4.210115092265168,
        4.209652564538259, 4.2093128379234557, 4.2090908215343665,
        4.2089814245024018, 4.208979555948912, 4.2090801249995105,
        4.2092780407905908, 4.2095682124426839, 4.2099455490770312,
        4.21040495982921, 4.210941353822947, 4.21154964017405,
        4.2122247280312735, 4.2129615264983942, 4.2137549447177651,
        4.2145998918085485, 4.21549127689822, 4.2164240091153955,
        4.2173929975788038, 4.2183931514300088, 4.2194193797817805,
        4.2204665917669635, 4.2215296965102223, 4.2226036031328738,
        4.2236832207759081, 4.2247634585544258, 4.2258392255938571,
        4.2269054310326029, 4.2279569839764779, 4.22898879357216,
        4.229995768937207, 4.2309728192022424, 4.2319148534901068,
        4.2328167809255488, 4.2336735106354446, 4.2344799517533813,
        4.235231013397458, 4.2359216046979595, 4.2365466347859835,
        4.237101012780836, 4.2375796478077552, 4.2379774489975564,
        4.23828932548468, 4.2385101863766614, 4.238634940811548,
        4.2386584979199906, 4.238575766821346, 4.2383816566429369,
        4.2380710765165759, 4.2376389355571851, 4.2370801429026548,
        4.2363896076751679, 4.23556223900385, 4.2345929460084006,
        4.2334766378238449, 4.2322082235726652, 4.2307826123774568,
        4.2291947133752483, 4.2274394356786607, 4.225511688423623,
        4.2234063807462778, 4.2211184217115507, 4.2186427206515269,
        4.215974186117684, 4.2131077287918739, 4.2100382535196248,
        4.2067606812899507, 4.2032698883881361, 4.1995608744996,
        4.1956283221200232, 4.1914690270509505, 4.1870850328524583,
        4.1824804963979529, 4.1776592573639775, 4.17262527883463,
        4.1673824791793113, 4.1619347929209445, 4.15628614874619,
        4.1504404774608847, 4.1444017090982976, 4.138173773961431,
        4.131760602272867, 4.125166124273111, 4.1183942702034351,
        4.1114489702997075, 4.1043341547992194, 4.097053753941343,
        4.0896116979589454, 4.0820119171015214, 4.0742583415980285,
        4.0663549016873457, 4.0583055276164481, 4.0501141496047115,
        4.0417846979133092, 4.0333211027680687, 4.0247272944039345,
        4.0160072030630332, 4.0071647589883561, 3.99820389241039,
        3.9891285335746534, 3.9799426127181858, 3.970650060065898,
        3.9612548058744812, 3.9517607803681227, 3.9421719138007822,
        3.9324921363900729, 3.9227253783940226, 3.9128755700395339,
        3.9029466415630316, 3.8929425232110488, 3.8828671452121353,
        3.8727244378163723, 3.862518331252669, 3.8522527557587933,
        3.8419316415830149, 3.831558918947775, 3.8211385181044868,
        3.8106743692912426, 3.8001704027326526, 3.789630548685436,
        3.7790587373681661, 3.7684588990391292, 3.7578349639218467,
        3.7471908622559442, 3.7365305242916582, 3.7258578802506985,
        3.7151768603834139, 3.7044913949229725, 3.6938054141058316,
        3.6831228481734755, 3.6724476273627498, 3.6617836819147289,
        3.6511349420653065, 3.6405053380443748, 3.6298988001065555,
        3.6193192584793641, 3.6087706434057942, 3.5982568851152892,
        3.5877819138564693, 3.5773496598645318, 3.5669640533727329,
        3.5566290246199617, 3.5463485038697824, 3.5361264212661205,
        3.525966707323485, 3.5158732915125457, 3.5058501061731064,
        3.495901075767688, 3.4860301464981647, 3.4762412045542579,
        3.46653829045543, 3.4569244162813861, 3.4474000402610434,
        3.4379645921864466, 3.4286176561902511, 3.4193587563765946,
        3.4101874385797419, 3.4011032408014672, 3.3921057038629043,
        3.3831943675579947, 3.3743687720735376, 3.3656284574362827,
        3.3569729637433996, 3.3484018310436392, 3.339914599435664,
        3.3315108089790435, 3.3231899997540348, 3.3149517118348006,
        3.3067954852797072, 3.29872086018253, 3.2907273766097864,
        3.2828145746316029, 3.27498199432374, 3.2672291757610497,
        3.2595556590183858, 3.2519609841639046, 3.2444446912800973,
        3.237006320429745, 3.2296454116900226, 3.2223615051453565,
        3.2151541408556734, 3.2080228589061717, 3.2009671993517346,
        3.1939867022837056, 3.187080907777581, 3.1802493558928933,
        3.1734915867113531, 3.1668071403043134, 3.1601955567489548,
        3.15365637611663, 3.1471891384810404, 3.1407933839140165,
        3.1344686524942929, 3.1282144842846029, 3.1220304193770976,
        3.1159159978290916, 3.1098707597204931, 3.1038942451280849,
        3.0979859941097474, 3.0921455467685788, 3.0863724431473805,
        3.0806662233335453, 3.0750264274118395, 3.0694525954323884,
        3.0639442674909176, 3.0585009836432206, 3.0531222839700645,
        3.0478077085594397, 3.0425567974542349, 3.0373690907548454,
        3.0322441285238972, 3.0271814508372481, 3.0221805977720555,
        3.0172411093896416, 3.0123625257749369, 3.0075443870034606,
        3.0027862331386252, 2.9980876042615381, 2.9934480404433752,
        2.9888670817582428, 2.9843442682787216, 2.9798791400829847,
        2.9754712372446166, 2.9711200998185272, 2.9668252679095444,
        2.9625862815779249, 2.9584026808592059, 2.9542740059699208,
        2.950199796582841, 2.9461795938670474, 2.9422129351718671,
        2.9382993831648889, 2.9344385633249344, 2.9306301264238392,
        2.9268737194568457, 2.9231689908766496, 2.9195155886052033,
        2.9159131607664235, 2.9123613554007988, 2.9088598205919145,
        2.9054082043934946, 2.9020061548843206, 2.8986533201160256,
        2.8953493481626724, 2.8920938870924084, 2.8888865849708494,
        2.8857270898693561, 2.8826150498430576, 2.8795501129610237,
        2.8765319272982057, 2.8735601409157527, 2.8706344018749728,
        2.8677543582537934, 2.8649196581102081, 2.862129949512846,
        2.8593848805262554, 2.8566840992263431, 2.854027253666688,
        2.8514139919167136, 2.8488439620509909, 2.8463168121273359,
        2.8438321902189947, 2.8413897443875431, 2.8389891226964017,
        2.8366299732263669, 2.8343119440278572, 2.832034683168648,
        2.8297978387334881, 2.8276010587617244, 2.8254439913414968,
        2.8233262845316576, 2.8212475863995503, 2.8192075450055505,
        2.8172058084227203, 2.8152420247211292, 2.8133158419588966,
        2.8114269082038135, 2.8095748715257143, 2.80775937999417,
        2.8059800816639853, 2.8042366246154589, 2.8025286569057952,
        2.8008558266054626, 2.79921778177632, 2.7976141704932611,
        2.7960446408174353, 2.7945088408129384, 2.7930064185508128,
        2.7915370220941309, 2.7901002995136173, 2.78869589887467,
        2.787323468233661, 2.7859826556752041, 2.7846731092550114,
        2.7833944770332923, 2.7821464070948907, 2.7809285474885552,
        2.7797405462913631, 2.7785820515647321, 2.7774527113755911,
        2.7763521737910608, 2.7752800868852034, 2.7742360987110093,
        2.7732198573363376, 2.7722310108435777, 2.7712692072843281,
        2.7703340947248711, 2.7694253212424309, 2.7685425348915667,
        2.7676853837480992, 2.7668535158735676, 2.76604657932834,
        2.7652642222024584, 2.7645060925313034, 2.7637718384037004,
        2.7630611078733431, 2.7623735490207775, 2.7617088098939653,
        2.761066538570752, 2.7604463831192096, 2.7598479915999765,
        2.7592710120866362, 2.7587150926370683, 2.7581798813260296,
        2.7576650262064253, 2.7571701753624969, 2.7566949768542566,
        2.7562390787414768, 2.755802129088762, 2.7553837759891264,
        2.7549836674726751, 2.7546014516294393, 2.7542367765165823,
        2.7538892901989844, 2.7535586407582859, 2.7532444762417678,
        2.752946444723797, 2.7526641942756456, 2.7523973729576094,
        2.7521456288332189, 2.7519086099785426, 2.7516859644517457,
        2.7514773403311463, 2.7512823856623423, 2.7511007485312398,
        2.7509320770010564, 2.750776019126973, 2.75063222298768,
        2.7505003366437606, 2.7503800081627152, 2.7502708856112572,
        2.75017261705942, 2.7500848505661466, 2.7500072342002966,
        2.7499394160363408, 2.749881044126155, 2.7498317665496121,
        2.7497912313739445, 2.7497590866517494, 2.7497349804634617,
        2.7497185608649453, 2.7497094759262306, 2.7497073737175879,
        2.7497119023045755, 2.749722709753236, 2.7497394441269809,
        2.7497617534937784, 2.7497892859214232, 2.7498216894701382,
        2.7498586122211468, 2.7498997022289515, 2.7499446075564786,
        2.7499929762891178, 2.7500444564677764, 2.7500986961718219,
        2.7501553434795731, 2.7502140464417688, 2.75027445312701,
        2.75033621160167, 2.7503989699337126, 2.7504623761935107,
        2.7505260784461569, 2.7505897247572761, 2.7506529631776786,
        2.7507154418049633, 2.7507768086730313, 2.7508367119119153,
        2.75089479940436, 2.7509507197397047, 2.7510041194990325,
        2.7510546508951665, 2.7511019515969997, 2.7511457563536621,
        2.7511860411241931, 2.7512228789934809, 2.7512563284236577,
        2.7512864536138935, 2.7513133166540058, 2.7513369803887646,
        2.7513575074250483, 2.7513749604321589, 2.7513894020385767,
        2.7514008949402573, 2.7514095017552473, 2.7514152851603373,
        2.7514183078014613, 2.7514186323376175, 2.7514163214279712,
        2.7514114377109213, 2.7514040438641056, 2.7513942025248785,
        2.7513819763572087, 2.7513674280184977, 2.7513506201486275,
        2.7513316154200642, 2.7513104764853011, 2.7512872659878127,
        2.7512620465979047, 2.7512348809527003, 2.7512058317293819,
        2.7511749615625511, 2.7511423331191618, 2.7511080090565692,
        2.7510720520156937, 2.7510345246642274, 2.7509954896569653,
        2.7509550096408968, 2.7509131472765129, 2.7508699652201112,
        2.7508255261249697, 2.7507798926502294, 2.7507331274345646,
        2.75068529315546, 2.7506364524578544, 2.7505866679960143,
        2.750536002423583, 2.7504845184038227, 2.7504322785793787,
        2.7503793456141565, 2.7503257821650142, 2.7502716508817358,
        2.7502170144158247, 2.7501619354290838, 2.7501064765841363,
        2.7500507005184143, 2.7499946698966724, 2.7499384473759179,
        2.7498820956033296, 2.7498256772415739, 2.7497692549537294,
        2.7497128913639162, 2.7496566491625383, 2.7496005909855623,
        2.7495447794925951, 2.7494892773397481, 2.7494341471812831,
        2.7493794516653978, 2.7493252534615351, 2.7492716152161876,
        2.7492185995807863, 2.7491662692186267, 2.7491146867753149,
        2.7490639149226483, 2.749014016307946, 2.7489650535340684,
        2.7489170894463624, 2.7488701861881495, 2.7488244078021653,
        2.748779813109167, 2.7487364753859258, 2.7486944278785472,
        2.7486538143535246, 2.7486144944694288, 2.7485782207793861,
        2.7485514462294489, 2.7485425166469351, 2.7485594937835387,
        2.7486105498819566, 2.7487038171599778, 2.7488474422951339,
        2.7490495667494095, 2.7493183338833158, 2.7496618863163289,
        2.7500883669909513, 2.7506059187183993, 2.7512226843461605,
        2.7519468067198152, 2.7527864286696655, 2.7537496930526615,
        2.7548447427106533, 2.756079720485026, 2.7574627692076494,
        2.7590020317388215, 2.7607056509070396, 2.7625817695605135,
        2.7646385305496612, 2.7668840767034024, 2.769326550877151,
        2.771974095903269, 2.7748348546288581, 2.777916969904497,
        2.7812285845527085, 2.7847778414416693, 2.7885728833947874,
        2.7926218532648575, 2.7969328938953271, 2.8015141481165697,
        2.80637375878878, 2.8115198687401248, 2.8169606208280333,
        2.8227041578839431, 2.8287586227544028, 2.835132158274932,
        2.8418329073007889, 2.8488690126721594, 2.8562486172194772,
        2.8639798638073648, 2.8720708952580618, 2.8805298544203448,
        2.889364884152692, 2.8985841272768069, 2.9081957266368494,
        2.9182078250926913, 2.9286285654678967, 2.9394660906227088,
        2.9507285433857011, 2.9624240666069479, 2.9745608031327118,
        2.9871468957945098, 3.0001904874483216, 3.0136997209188308,
        3.0276827390651508, 3.042147684731193, 3.0571027007506362,
        3.07255592997112, 3.0885155152316175, 3.10498959937238,
        3.1219863252523159, 3.1395138356952725, 3.1575802735531515,
        3.1761937816728536, 3.1953625028750139, 3.2150945800672788,
        3.2353981559250458, 3.2562813737204337, 3.2777523751296029,
        3.2998193062577839, 3.3224903008656992, 3.345773526929511,
        3.3696770577922095, 3.3942092285498915, 3.4193776519444254,
        3.4451917974593953, 3.4716487636171856, 3.4987149290870949,
        3.5263443015694294, 3.5544927455185986, 3.5831154030249075,
        3.6121676779297749, 3.6416048794676157, 3.6713823510551391,
        3.7014554237617125, 3.7317794331101686, 3.7623097130483538,
        3.7930015980498437, 3.8238104224073806, 3.8546915204989372,
        3.8856002266636924, 3.9164918752507432, 3.9473218005973889,
        3.9780453370608604, 4.0086178189859636, 4.0389945807142329,
        4.0691309565990785, 4.09898228098241, 4.1285038882155538,
        4.15765111264618, 4.1863792886173341, 4.214643750477399,
        4.2423998325761687, 4.2696028692561976, 4.2962081948826345,
        4.3221711437250168, 4.3474470503629394, 4.3719912484970056,
        4.3957590742766408, 4.4187058570519513, 4.440786944999715,
        4.461957634165441, 4.482173364924483, 4.5013891782500153,
        4.5195612174294277, 4.5366427923370356, 4.5526060910352673,
        4.5674701801910551, 4.5812730046800123, 4.5940496759530216,
        4.6058364077558593, 4.6166690144490659, 4.6265834547532441,
        4.6356156352015336, 4.6438014811921766, 4.6511769113003707,
        4.657777846569803, 4.6636402071632137, 4.6687999135383329,
        4.6732928860721215, 4.6771550451620474, 4.6804223111732846,
        4.68313060452159, 4.6853158455880539, 4.6870139547605838,
        4.6882608524260876, 4.6890924589754786, 4.6895446948056243,
        4.68965348029561, 4.6894547358456693, 4.6889843818253247,
        4.6882783386544524, 4.6873725267034239, 4.68630286635742,
        4.6851052780078284, 4.6838156820873245, 4.6824699988450122,
        4.6811041490372851, 4.6797540520817007, 4.6784556310226746,
        4.67724479890504, 4.676157496475625, 4.6752295877430123,
        4.6744971491475491, 4.6739956707263133, 4.6737621499363629,
        4.6738235404089847, 4.6741818547335425, 4.6748290616203452,
        4.6757586372558713, 4.6769634713982615, 4.6784366662608958,
        4.6801712472587722, 4.6821602675884044, 4.68439677038921,
        4.6868738024498073, 4.6895844092240875, 4.6925216366559157,
        4.6956785305085544, 4.6990481366237979, 4.7026235007905486,
        4.7063976688430573, 4.7103636865838014, 4.7145145998298164,
        4.7188434543948983, 4.7233432960762851, 4.7280071707194633,
        4.7328281241181918, 4.7377992020787669, 4.7429134504284534,
        4.74816391497798, 4.7535436415319845, 4.7590456759175135,
        4.7646630639372827, 4.7703888514104609, 4.7762160841464913,
        4.7821378079603125, 4.788147068669903, 4.7942369120769852,
        4.800400384002887, 4.8066305302701089, 4.812920396673686,
        4.8192630290375762, 4.8256514731804048, 4.8320787748952547,
        4.838537980017553, 4.8450221343475475, 4.8515242837095442,
        4.858037473902943, 4.8645547507613021, 4.87106916007152,
        4.8775737476529732, 4.8840615593577335, 4.8905256409446389,
        4.8969590382510422, 4.9033547971021054, 4.9097059632831321,
        4.9160055826406381, 4.92224670096888, 4.9284223640672264,
        4.9345256177805883, 4.9405495079040476, 4.94648708024404,
        4.9523313806379621, 4.9580754548838062, 4.963712348786931,
        4.9692351081791823, 4.9746367788618766, 4.9799104066461544,
        4.985049037354834, 4.9900457167995063, 4.9948934907904894,
        4.9995854051390509, 5.00411450566137, 5.0084738381766423,
        5.01265644848238, 5.0166553824219369, 5.0204636857398155,
        5.0240744044408752, 5.0274805837879244, 5.030675271145661,
        5.0336515080155682, 5.0364023521600751, 5.0389208162981109,
        5.0412000375197179, 5.0432328332207668, 5.0450141506723183,
        5.04654422596943, 5.0478254250835191, 5.0488597943004638,
        5.0496495042764815, 5.0501966805999405, 5.050503465150193,
        5.0505719939264111, 5.05040440503681, 5.0500028358509326,
        5.0493694239805382, 5.0485063069638887, 5.047415622361985,
        5.0460995077278907, 5.0445601006187983, 5.0427995385891968,
        5.0408199591971616, 5.03862349998265, 5.0362122985416855,
        5.0335884923884615, 5.0307542190819685, 5.0277116162079549,
        5.0244628212829294, 5.0210099718959063, 5.0173552055855373,
        5.0135006599017391, 5.0094484724170494, 5.0052007806782592,
        5.000759722231364, 4.9961274346479811, 4.9913060554800373,
        4.9862977222678841, 4.9811045725898833, 4.9757287439824429,
        4.9701723740126873, 4.9644376002315633, 4.9585265601985,
        4.9524413914603951, 4.9461842315816336, 4.9397572181104241,
        4.9331624886095291, 4.9264021806318645, 4.9194784317265272,
        4.91239337945799, 4.9051491613734246, 4.8977479150345653,
        4.8901917779987105, 4.8824828878138788, 4.8746233820500535,
        4.8666153982396869, 4.8584610739418768, 4.8501625467417391,
        4.8417219541608558, 4.8331414337795975, 4.8244231231306687,
        4.8155691597789625, 4.8065816812883249, 4.7974628252052538,
        4.78821472907916, 4.7788395304871987, 4.7693393669637851,
        4.7597163760687753, 4.74997269536557, 4.7401104624008843,
        4.730131814743622, 4.7200388899228827, 4.7098338255277854,
        4.6995187590809984, 4.689095828164457, 4.6785671703260086,
        4.667934923104279, 4.6572012240843268, 4.6463682107942867,
        4.6354380208015806, 4.624412791671034, 4.61329466093813,
        4.6020857661657768, 4.59078824493042, 4.5794042347542554,
        4.5679358732040205, 4.5563852978521222, 4.5447546462326143,
        4.5330460559146344, 4.52126166444576, 4.509403609375549,
        4.4974740282809185, 4.4854750586961663, 4.4734088381921628,
        4.4612775043025055, 4.4490831946080833, 4.4368280466561352,
        4.4245141979931937, 4.41214378618522, 4.3997189487802055,
        4.3872418233385311, 4.3747145474057758, 4.3621392585638539,
        4.3495180943342193, 4.3368531922953917, 4.3241466899904006,
        4.3114007249761928, 4.2986174348310477, 4.2857989570715516,
        4.2729474292756766, 4.2600649890111271, 4.247153773798253,
        4.2342159212225914, 4.2212535688276187, 4.2082688541739,
        4.195263914814392, 4.1822408882894164, 4.1692019121920412,
        4.1561491240388522, 4.1430846614090457, 4.1300106618482459,
        4.1169292629019552, 4.1038426021506984, 4.0907528171373118,
        4.077662045414475, 4.0645724245345427, 4.0514860920644757,
        4.0384051855511931, 4.0253318425487956, 4.01226820061877,
        3.9992163973217525, 3.9861785701944004, 3.9731568568107027,
        3.9601533947173384, 3.947170321471829, 3.9342097746265492,
        3.9212738917439536, 3.9083648103743323, 3.8954846680655293,
        3.8826356023910775, 3.8698197508929932, 3.8570392511304581,
        3.844296240661512, 3.8315928570353512, 3.818931237813747,
        3.8063135205522203, 3.7937418427943923, 3.7812183421169516,
        3.7687451560526473, 3.756324422174564, 3.7439582780312382,
        3.7316488611719105, 3.7193983091670884, 3.7072087595599235,
        3.695082349903732, 3.6830212177636277, 3.6710275006953763,
        3.6591033362447285, 3.6472508619827639, 3.6354722154099695,
        3.6237695342744622, 3.612144955568334, 3.6006006183916814,
        3.5891386560806429, 3.5777612177921765, 3.5664704091626054,
        3.5552684475904313, 3.5441568059441293, 3.5331351082003928,
        3.5222022338532937, 3.5113571740850218, 3.5005988766289633,
        3.4899263049848126, 3.4793384169286972, 3.4688341723205181,
        3.4584125302655857, 3.448072450125065, 3.4378128911933792,
        3.4276328127723623, 3.4175311741755996, 3.4075069346781324,
        3.3975590536046045, 3.3876864902709309, 3.3778882039476135,
        3.3681631539579868, 3.3585102996005527, 3.3489286001876684,
        3.3394170150051523, 3.3299745033572972, 3.320600024569178,
        3.3112925379247735, 3.3020510027262633, 3.2928743782863124,
        3.2837616238986453, 3.2747116988820717, 3.2657235625157246,
        3.2567961741231315, 3.2479284930077834, 3.239119478443377,
        3.2303680897674232, 3.2216732862767885, 3.2130340272532063,
        3.2044492720288225, 3.1959179798788324, 3.1874391101261104,
        3.1790116220711542, 3.1706344750041593, 3.1623066282473338,
        3.1540270410794333, 3.1457946728380617, 3.1376084827864132,
        3.1294674302548486, 3.1213704745425375, 3.1133165749399785,
        3.1053046907651081, 3.0973337813133011, 3.0894028058868606,
        3.0815107237928863, 3.073656494333147, 3.0658390768083845,
        3.0580574305342316, 3.0503105147814709, 3.0425972888909896,
        3.034916712143982, 3.0272677438455782, 3.0196493433114955,
        3.0120604698283353, 3.0045000827057255, 2.9969671412532075,
        2.9894606047670811, 2.981979432545208, 2.974522583898596,
        2.9670890181333842, 2.9596776945397254, 2.9522875724372968,
        2.9449176111252831, 2.9375667698850036, 2.9302340080221665,
        2.9229182849551547, 2.9156185596038093, 2.9083337923131194,
        2.9010629395376122, 2.8938049684436118, 2.88655881657347,
        2.8793235034700126, 2.8720978217206148, 2.8648811903441134,
        2.8576714182106349, 2.8504770417464282, 2.8433332362832271,
        2.8362859046882192, 2.8293793396739728, 2.822658460438157,
        2.816167959157764, 2.8099526100473442, 2.8040571576945719,
        2.7985263573748576, 2.793404960517146, 2.7887377199339078,
        2.7845693879316062, 2.7809447170054629, 2.777908459588597,
        2.7755053681252106, 2.7737801950667929, 2.7727776928560544,
        2.7725426139253795, 2.7731197107291723, 2.7745537357289427,
        2.776889441333513, 2.7801715800186058, 2.7844449042154817,
        2.7897541663746654, 2.7961441189395755, 2.8036595143467649,
        2.8123451050553605, 2.822245643511093, 2.833405882109564,
        2.8458705734483849, 2.8596844695891681, 2.8748923239769155,
        2.8915388862993261, 2.909668916696754, 2.9293271462307042,
        2.9505583865867142, 2.973407226229722, 2.9979188713601008,
        3.0241368232226349, 3.0521089655038445, 3.0818539833286329,
        3.1133180555799096, 3.1464181625885614, 3.1810756671042055,
        3.2172102269552787, 3.2547421176785893, 3.2935913915889,
        3.3336781816654866, 3.3749225917254524, 3.417244736147528,
        3.4605647254643639, 3.5048026716250504, 3.5498786860579945,
        3.5957128803580485, 3.6422253661032338, 3.689336254845724,
        3.7369656581419255, 3.7850336875790109, 3.8334604547133515,
        3.8821660710966626, 3.93107064831182, 3.9800942979120313,
        4.029157131468196, 4.0781792605435969, 4.1270807967035079,
        4.1757818515134471, 4.2242025365419034, 4.2722629633376306,
        4.31988324350424, 4.3669834885057659, 4.4134838101981471,
        4.4593043193298643, 4.5043651297113465, 4.548586346724524,
        4.5918880990503776, 4.6341904508284566, 4.6754136449161789,
        4.7154774295996429, 4.75430291810771, 4.7918077152880905,
        4.827932801037317, 4.862677200706214, 4.8960633146784573,
        4.9281100349711293, 4.95883761852404, 4.9882658277480365,
        5.0164146037346322, 5.0433038230288956, 5.0689533854909312,
        5.0933831825791449, 5.1166131087555371, 5.1386630574178316,
        5.1595529223524741, 5.1793025971621081, 5.1979319755843552,
        5.2154609512641317, 5.2319094178853227, 5.2472972691082047,
        5.2616443986253154, 5.2749707001236406, 5.28729606724239,
        5.2986403936927751, 5.3090235731393225, 5.3184654992428024,
        5.3269860656985744, 5.3346051661771385, 5.341342694349847,
        5.3472185438985491, 5.3522526084948572, 5.3564647818193256,
        5.3598749575382021, 5.3625030293356621, 5.3643688908883727,
        5.3654924358714089, 5.3658935579510674, 5.3655921508075917,
        5.3646081081327273, 5.3629613235842157, 5.3606716908475454,
        5.3577591035863961, 5.3542434554856859, 5.3501446402243342,
        5.34548255147972, 5.3402770829049429, 5.3345481282062064,
        5.3283155810490506, 5.3215993350937971, 5.314419284041783,
        5.3067953215523813, 5.2987473412977506, 5.2902952369660907,
        5.28145890223035, 5.2722582307709391, 5.2627131162442851,
        5.2528434523366272, 5.2426691327441475, 5.2322100511134328,
        5.2214861011371765, 5.2105171764780218, 5.199323170828734,
        5.1879239778446395, 5.1763394912292959, 5.1645896046351671,
        5.1526942117351107, 5.1406732062367837, 5.1285464817739248,
        5.1163339320529095, 5.1040554507415754, 5.091730931508712,
        5.0793802680356537, 5.0670233540153529, 5.0546800830554544,
        5.0423703490276059, 5.0301140450456989, 5.0179310663475363,
        5.0058413023402863, 4.9938646584385618, 4.9820209958689681,
        4.9703302980866919, 4.9588122112074711, 4.9474872484717709,
        4.9363701455981577, 4.9254612914805973, 4.9147552974396858,
        4.9042476419562968, 4.893933466193972, 4.8838080334949581,
        4.8738665630607363, 4.8641042900354092, 4.8545164438195947,
        4.8450982558780291, 4.8358449569235828, 4.826751777947341,
        4.8178139498473165, 4.8090267035467233, 4.8003852699528222,
        4.791884879989313, 4.7835207645813913, 4.7752881546361312,
        4.7671822810721167, 4.7591983748053091, 4.7513316667513754,
        4.7435773878538052, 4.7359307689864272, 4.7283870410918007,
        4.7209414350927537, 4.7135891818953874, 4.7063255124176822,
        4.6991456575769961, 4.6920448482912214, 4.6850183154818659,
        4.6780612900664469, 4.6711690029571242, 4.6643366850720893,
        4.6575595673307806, 4.650832880655515, 4.6441518559472224,
        4.6375117241447317, 4.6309077161379539, 4.6243350628824835,
        4.6177889952651014, 4.6112647442086834, 4.6047575406421322,
        4.5982626154662167, 4.5917751996156877, 4.58529052399006,
        4.5788038195191785, 4.5723103171183981, 4.5658052477104256,
        4.5592838421923183, 4.5527413315038334, 4.54617294655667,
        4.5395739182534509, 4.5329394775317127, 4.5262648553021068,
        4.5195452824757343, 4.512775989974366, 4.5059522087121122,
        4.4990691696243061, 4.4921221035994616, 4.485106241564238,
        4.4780168144585684, 4.4708490531684264, 4.4635981886301632,
        4.456259451760336, 4.4488280734529839, 4.4412992846724846,
        4.4336683162934767, 4.4259303992353578, 4.4180807644493783,
        4.4101146428177795, 4.4020272652989219, 4.3938138626982211,
        4.385469666224787, 4.3769899060673714, 4.36836981508744,
        4.359604618825478, 4.3506895631177747, 4.3416198375770056,
        4.3323907869994445, 4.32299735733778, 4.313437151822118,
        4.3037143723153122, 4.2938358779487453, 4.2838081290525745,
        4.2736377410703739, 4.2633312732814339, 4.2528953052337837,
        4.2423364091617719, 4.2316611599625551, 4.220876131536313,
        4.2099878981548047, 4.199003033969352, 4.1879281131724868,
        4.1767697099255479, 4.1655343984092008, 4.1542287528007069,
        4.142859347283, 4.13143275601605, 4.1199555531958483, 4.1084343129935927,
        4.0968756095812937, 4.0852860171257, 4.0736721098267923,
        4.0620404618523942, 4.050397647365517, 4.038750240559013,
        4.0271048156011151, 4.01546794667096, 4.0038462079551138,
        3.9922461736055008, 3.9806744178253473, 3.9691375147820329,
        3.9576420386270281, 3.9461945635889855, 3.9348016637975607,
        3.9234699134535114, 3.9122058867294833, 3.9010161577898921,
        3.889907300835461, 3.878885890018581, 3.8679584995217975,
        3.8571317035375383, 3.8464120762217795, 3.8358061917595574,
        3.8253206243421944, 3.8149619481167503, 3.804736737278068,
        3.7946515660035467, 3.7847130084644354, 3.7749276388391291,
        3.7653020313105872, 3.755842760037889, 3.7465563992148527,
        3.7374495230159952, 3.7285287056043606, 3.7198005211674694,
        3.7112715438925621, 3.7029483479364673, 3.6948375074746633,
        3.6869455967100571, 3.6792791897965627, 3.6718448609136076,
        3.6646491842363038, 3.6576987339518734, 3.651000084232694,
        3.644559809250504, 3.638384483180698, 3.6324806802104637,
        3.6268549745058754, 3.6215139402549132, 3.6164641515923455,
        3.6117121828363277, 3.6072646077890185, 3.6031280015783107,
        3.5993089358448587, 3.5958139917400356, 3.5926497241142337,
        3.5898227606905704, 3.5873395279048457, 3.5852069697458431,
        3.5834285815657942, 3.581999294915009, 3.5809105926730469,
        3.5801544753807493, 3.5797227421389879, 3.5796072650397095,
        3.5797998898279362, 3.5802924717396225, 3.5810768625971687,
        3.5821449154444713, 3.5834884828893623, 3.5850994177052029,
        3.5869695725962143, 3.5890908002851845, 3.5914549535064975,
        3.5940538849818569, 3.5968794474209478, 3.5999234935731921,
        3.603177876144442, 3.6066344478578847, 3.6102850614503548,
        3.6141215696470312, 3.6181358251460973, 3.6223196806964992,
        3.6266649890248743, 3.6311636028312155, 3.6358073748680617,
        3.6405881578380273, 3.6454978044717712, 3.6505281675050885,
        3.6556710996492416, 3.6609184536223771, 3.6662620821653347,
        3.6716938379921826, 3.6772055738268916, 3.6827891423978407,
        3.6884363964327669, 3.69413918864257, 3.6998893717619876,
        3.7056787985181545, 3.7114993216191268, 3.7173427938058272,
        3.7232010677911922, 3.7290659963053439, 3.7349294320778079,
        3.740783227821058, 3.7466192362677142, 3.7524293101335013,
        3.7582053021486326, 3.7639390650402058, 3.7696224515203274,
        3.775247314332725, 3.7808055061797696, 3.7862888797974414,
        3.7916892879153274, 3.7969985832424795, 3.802208618521417,
        3.807311246449145, 3.8122983197800839, 3.8171616912282897,
        3.8218932135012604, 3.8264847393420065, 3.8309281214752793,
        3.8352152126143624, 3.8393378654828378, 3.8432879328170397,
        3.8470572673324512, 3.8506377217531149, 3.8540211488068774,
        3.8571994012119917, 3.8601643317223449, 3.8629077929626456,
        3.8654216379162039, 3.8676977186024772, 3.8697278897443943,
        3.8715039984491622, 3.87301791305462, 3.8742614430532,
        3.8752265603654905, 3.8759048194246826, 3.8762905562144629,
        3.8763850138388349, 3.8761922169914409, 3.87571577283176,
        3.8749594509752754, 3.873926962184207, 3.8726220384543186,
        3.8710484041386524, 3.8692097863601305, 3.8671099112121614,
        3.8647525051718268, 3.8621412945761757, 3.85928000582141,
        3.8561723652746305, 3.8528220993167448, 3.8492329343216327,
        3.8454085966646243, 3.8413528127257783, 3.8370693088689172,
        3.8325618114886644, 3.82783404695449, 3.8228897416336296,
        3.817732621905356, 3.8123664141560183, 3.80679484474823,
        3.8010216400766765, 3.7950505264972318, 3.7888852303872893,
        3.7825294781384873, 3.7759869961247734, 3.7692615106837488,
        3.7623567482954821, 3.7552764350941454, 3.7480242980856455,
        3.7406040619230643, 3.7330194577306473, 3.7252741987425617,
        3.7173720477617986, 3.7093166302682588, 3.7011119509324026,
        3.6927610396716886, 3.6842734207296237, 3.6756747451576244,
        3.6669971583155054, 3.6582718308055053, 3.6495303124837211,
        3.6408040157801418, 3.6321244027740396, 3.6235229176460062,
        3.6150310109980941, 3.6066801311277445, 3.598501727183598,
        3.590527247964812, 3.5827881424544925, 3.5753158595293186,
        3.5681418481028317, 3.5612975570922174, 3.5548144354031903,
        3.548723931955482, 3.5430574956524437, 3.5378465754139174,
        3.5331226201454649, 3.5289170787483086, 3.5252614001611571,
        3.5221870332782017, 3.519725427005036, 3.5179080302683632,
        3.51676629195287, 3.5163316610183957, 3.5166355863187886,
        3.5177095168321912, 3.5195849013163554, 3.5222931891442726,
        3.5258658279168968, 3.5303342701533196, 3.5357299548229428,
        3.5420843583532382, 3.5494288534388718, 3.5577950999800723,
        3.567213965057999, 3.5777183536467563, 3.5893275930820265,
        3.6020272944370286, 3.6157894911840676, 3.6305882546438095,
        3.6463968633240964, 3.6631888829828423, 3.6809377755769956,
        3.6996170405806312, 3.7192001638886634, 3.7396606363207847,
        3.7609719469173104, 3.7831075853463516, 3.8060410410535224,
        3.8297458035840086, 3.8541953624305645, 3.879363207104455,
        3.9052228271035614, 3.9317477119512581, 3.9589113511488732,
        3.9866872342144619, 4.0150488506309552, 4.0439696899306723,
        4.0734232416193272, 4.1033829952005121, 4.1338224401772417,
        4.16471506605897, 4.1960343623750775, 4.2277538186087593,
        4.2598469242696106, 4.292287168878306, 4.32504804194703,
        4.3581030329643085, 4.3914256314554656, 4.4249893269187224,
        4.4587676088727193, 4.4927339668153037, 4.5268618902618227,
        4.5611248687113566, 4.5954963916907348, 4.6299499486893012,
        4.6644590292173564, 4.6989971228033243, 4.7335377189259225,
        4.768054307116345, 4.80252037687779, 4.8369094177108183,
        4.8711949191237984, 4.9053503706455439, 4.9393492617582035,
        4.9731650819860631, 5.0067713208301523, 5.0401414677990548,
        5.073249012398481, 5.10606744415228, 5.138570252561764,
        5.1707309271132926, 5.2025229573453409, 5.2339198327579473,
        5.2648950428445733, 5.2954220771313656, 5.3254744251302517,
        5.3550255763187131, 5.3840490202403473, 5.4125182463869121,
        5.4404067442741244, 5.4676880033947981, 5.4943355132717207,
        5.5203227634095891, 5.5456232433256449, 5.5702104424878094,
        5.5940578505247611, 5.6171389565960066, 5.6394272511999342,
        5.6608962211132781, 5.6815193633479346, 5.7012701467031128,
        5.720122118000404, 5.7380486080709279, 5.7550235441055735,
        5.7710193201744922, 5.7860185453075212, 5.8000291944634439,
        5.8130694575855655, 5.8251559914146718, 5.836306049172725,
        5.8465366679841093, 5.8558649630687167, 5.8643080214005892,
        5.8718829401703729, 5.878606812886332, 5.8844967343848227,
        5.8895697990138824, 5.8938431013085344, 5.897333735727698,
        5.9000587967669755, 5.9020353788993054, 5.9032805765994363,
        5.9038114843618477, 5.9036451966752592, 5.9027988079954747,
        5.901289412829323, 5.89913410564131, 5.8963499809224169,
        5.8929541331563025, 5.8889636568103549, 5.8843956463932994,
        5.8792671963493452, 5.8735954011912437, 5.8673973553903132,
        5.86069015343354, 5.8534908897812317, 5.8458166589458695,
        5.8376845553903092, 5.8291116735907993, 5.8201151080513576,
        5.8107119532382727, 5.8009193036337079, 5.7907542537153942,
        5.7802338979809864, 5.76937533089411, 5.7581956469524727,
        5.7467119406191074, 5.7349413064019519, 5.7229008387518849,
        5.7106076321743249, 5.6980787811424243, 5.6853313801266507,
        5.6723825236437646, 5.6592493061319642, 5.6459488220938434,
        5.6324981660212, 5.6189144323754192, 5.6052147156475121,
        5.5914161103247988, 5.5775357108847929, 5.5635906118003042,
        5.5495979075541868, 5.5355746926527054, 5.5215380615422287,
        5.5075051087330387, 5.4934929286905083, 5.4795186158969056,
        5.4655992648452925, 5.4517519700056285, 5.4379938258663323,
        5.4243419269014019, 5.4108133676061128, 5.397425242447107,
        5.3841946459221859, 5.3711386724931955, 5.3582744166666973,
        5.3456189728557355, 5.3331894357632805, 5.3210028992684792,
        5.3090764594573843, 5.2974272064365433, 5.2860722467315737,
        5.275028641545541, 5.264313577447318, 5.2539438949717061,
        5.2439373241179288, 5.2343056686477274, 5.2250460160044687,
        5.2161495274231955, 5.2076082535758541, 5.1994138990989223,
        5.1915582939950324, 5.1840332229829826, 5.1768304871485906,
        5.1699418816414964, 5.1633592037522993, 5.1570742500247224,
        5.1510788172698083, 5.1453647021778162, 5.1399237014935375,
        5.1347476119593765, 5.1298282302914862, 5.1251573532322041,
        5.1207267775098666, 5.11652829985995, 5.112553717017998,
        5.108794825698511, 5.1052434226683587, 5.1018913046233871,
        5.0987302683172278, 5.0957521104746526, 5.0929486278377629,
        5.0903116171305074, 5.0878328750740618, 5.0855041984197458,
        5.0833173838991321, 5.0812642282408449, 5.0793365281655811,
        5.0775260804254474, 5.0758246817364183, 5.0742241288427028,
        5.0727162184639942, 5.0712927473525831, 5.06994551222088,
        5.0686663098087328, 5.0674469368612591, 5.0662791900876449,
        5.065154866236913, 5.0640657620327927, 5.0630036742057563,
        5.0619603995062024, 5.060927734653017, 5.0598974763723747,
        5.05886142141048, 5.0578113664848026, 5.0567391083444182,
        5.0556364437107968, 5.0544951693202043, 5.0533070819074695,
        5.052063978206399, 5.0507576549282556, 5.0493799088313223,
        5.0479225366478557, 5.0463773350933065, 5.0447361009058884,
        5.0429906308256305, 5.0411327215753259, 5.03915416989184,
        5.0370467725126762, 5.03480232615914, 5.0324126275870213,
        5.0298694734907246, 5.0271646606214979, 5.0242899857285632,
        5.0212372455370566, 5.0179982367444884, 5.0145647561512527,
        5.0109286003540641, 5.0070815664620847, 5.00301545012217,
        4.9987220511122752, 4.9941931577300016, 4.9894205900112354,
        4.9843960802596134, 4.9791116028583078, 4.9735585099016326,
        4.9677322999274907, 4.9616387680072265, 4.9552878557055031,
        4.9486888822667092, 4.9418514089917194, 4.9347849095001761,
        4.9274988891081826, 4.920002841661808, 4.9123062651708285,
        4.90441865611984, 4.8963495115624012, 4.8881083283325495,
        4.879704603348884, 4.871147833500082, 4.8624475156903051,
        4.8536131468059507, 4.8446542237558088, 4.835580243420023,
        4.8264007027189466, 4.8171250985260814, 4.8077629277352569,
        4.7983236872672439, 4.7888168740178712, 4.7792519848461863,
        4.7696385166828259, 4.7599859664239341, 4.7503038309412053,
        4.74060160716791, 4.7308887919622782, 4.7211748822475874,
        4.7114693748998873, 4.7017817668364037, 4.6921215549454764,
        4.6824982361145651, 4.6729213072557938, 4.663400265254988,
        4.6539446070045347, 4.6445638294205143, 4.6352674293723846,
        4.6260649037756849, 4.6169657495280116, 4.6079794635130735,
        4.5991155426365316, 4.5903834837875257, 4.581792783878214,
        4.5733529397818646, 4.5650734484175342, 4.5569638066666265,
        4.5490335114316984, 4.5412920596053086, 4.5337489480902438,
        4.5264136737764895, 4.5192957335627852, 4.5124046243503892,
        4.5057498430346925, 4.4993408865007893, 4.4931872516528832,
        4.4872984353966032, 4.4816839346157584, 4.476353246209861,
        4.471315867070353, 4.4665812941162164, 4.462159024210786,
        4.4580585542862572, 4.4542893812034032, 4.4508610018732888,
        4.4477829132104008, 4.4450646120793325, 4.4427155954095836,
        4.4407453600696583, 4.4391634029713991, 4.4379792209755591,
        4.4372023111283845, 4.4368421698851925, 4.43690829534347,
        4.43741018110041, 4.4383573331472261, 4.4397592232544412,
        4.4416254178427508, 4.443965222068762, 4.4467886126968326,
        4.4501010917737149, 4.45389704971451, 4.458166402215018,
        4.4628997366118082, 4.4680873788988551, 4.4737197497940491,
        4.47978723579962, 4.4862802357253591, 4.493189143997542,
        4.5005043565683884, 4.5082162688801519, 4.5163152765553942,
        4.5247917751406153, 4.5336361601842912, 4.5428388272762019,
        4.5523901719896029, 4.5622805898810022, 4.572500476518055,
        4.583040227457337, 4.5938902382892746, 4.6050409045734986,
        4.6164826218761368, 4.6282057857520229, 4.6402007917888808,
        4.65245803554025, 4.6649679125813535, 4.6777208184782353,
        4.6907071488020184, 4.7039172991066769, 4.7173416649715536,
        4.7309706419609858, 4.7447946256496669, 4.7588040115871033,
        4.7729891953562262, 4.7873405725290983, 4.8018485386469116,
        4.8165034893167551, 4.8312958200700518, 4.8462159264864555,
        4.8612542041384819, 4.8764010485979083, 4.8916468554199426,
        4.9069820201822871, 4.9223969384399409, 4.9378820057701818,
        4.9534276177431709, 4.969024169910873, 4.984662057863253,
        5.00033167715408, 5.01602342334266, 5.0317276920254992,
        5.047434878738744, 5.0631353790586475, 5.078819588577689,
        5.0944779028136669, 5.1101007173905248, 5.1256784278317777,
        5.1412014297188948, 5.1566601186370056, 5.17204489013488,
        5.1873461397787786, 5.2025542631322832, 5.2176596557866732,
        5.2326527132955416, 5.2475238312133632, 5.2622634051339077,
        5.2768618306058253, 5.29130950319219, 5.3055968184909883,
        5.3197141720183749, 5.3336519594057723, 5.3474005761306467,
        5.3609504180004137, 5.3742918798667549, 5.3874153592744847,
        5.4003112463835263, 5.4129699516954117, 5.4253818294301119,
        5.4375373891915642, 5.4494267411609592, 5.4610426566549721,
        5.4723845150784722, 5.4834543569179948, 5.49425382333579,
        5.504784710788476, 5.51504875949038, 5.5250477299794731,
        5.5347833754440527, 5.5442574517312151, 5.5534717137249912,
        5.5624279166749311, 5.5711278156775954, 5.5795731658735876,
        5.5877657224300465, 5.5957072404801504, 5.6033994751414786,
        5.6108441815820953, 5.6180431149468637, 5.6249980303667471,
        5.63171068296834, 5.6381828279198194, 5.64441622034714,
        5.6504126154164629, 5.6561737682261413, 5.6617014339584086,
        5.66699736774446, 5.6720633247028145, 5.6769010600070429,
        5.6815123287857388, 5.6858988861832334, 5.6900624873281238,
        5.6940048873779183, 5.6977278414753529, 5.7012331047657723,
        5.7045224323529204, 5.7075975794350544, 5.7104603011198245,
        5.7131123525529128, 5.7155554888798834, 5.7177914652451332,
        5.71982203679085, 5.7216489586429775, 5.7232739859709314,
        5.7246988738914748, 5.7259253775668437, 5.7269552521203746,
        5.7277902527114639, 5.7284321344625928, 5.7288826525348124,
        5.7291435620533537, 5.7292166181770483, 5.7291035760330127,
        5.7288061907730006, 5.7283262175346916, 5.72766541145566,
        5.7268255276898534, 5.7258083213549389, 5.724615547631271,
        5.7232489616277658, 5.7217103184951119, 5.7200013733839272,
        5.7181238814321178, 5.716079597774459, 5.7138702775572536,
        5.7114976759111729, 5.7089635480100007, 5.7062696489768054,
        5.7034177339300536, 5.7004095580498655, 5.6972468764629651,
        5.6939314442851376, 5.6904650167382833, 5.6868493487899521,
        5.6830861960228454, 5.6791773123803, 5.6751244562640251,
        5.6709293738632764, 5.6665938450786175, 5.6621195567515636,
        5.657508434994992, 5.6527608113165728, 5.6478730571915481,
        5.6428399494739772, 5.6376565043261673, 5.6323176447779852,
        5.6268183276504757, 5.6211534975300745, 5.6153181034301607,
        5.6093070927595079, 5.6031154135008858, 5.5967380134492091,
        5.5901698404516491, 5.5834058423222617, 5.576440966936735,
        5.5692701620930736, 5.5618883756541342, 5.554290555446169,
        5.5464716492987218, 5.5384266050598354, 5.5301503705467736,
        5.5216378936182711, 5.512884122099452, 5.50388400382711,
        5.4946324866388627, 5.4851245183612738, 5.4753550468494714,
        5.4653190199156656, 5.455011385423834, 5.444427091188734,
        5.433561085038451, 5.4224083148556872, 5.4109637284120433,
        5.3992222735980224, 5.3871788982210322, 5.374828550117229,
        5.3621661771400957, 5.3491867270998332, 5.3358851478534506,
        5.3222563872360062, 5.3082953930812407, 5.2939971132183832,
        5.2793564954811751, 5.26436848772034, 5.249028037761863,
        5.2333300934371367, 5.2172696026044862, 5.20084151306294,
        5.1840407726957167, 5.166862329293699, 5.1493011307151884,
        5.1313521248018361, 5.1130102593723654, 5.0942704822782812,
        5.0751277413450122, 5.0555769844222844, 5.0356131593246927,
        5.0152312139058006, 4.9944260959901507, 4.9731927534387559,
        4.9515261340426715, 4.9294211856783372, 4.9068728561738633,
        4.8838760933490812, 4.8604258450598632, 4.8365170591203013,
        4.8121446833901755, 4.7873036656896284, 4.761988953849027,
        4.7361954957379941, 4.7099182391293413, 4.6831521319972866,
        4.6558921218029585, 4.6281331574618321, 4.5998701837445282,
        4.57109815702776, 4.5418120014854741, 4.5120067304579941, 4.481677110529,
        4.4508185892963326, 4.4194248639220417, 4.3875012940186062,
        4.3550821996150884, 4.32221356318571, 4.2889396167932095,
        4.2553052734848427, 4.2213551995460019, 4.1871341504648534,
        4.1526868495036835, 4.1180580315540487, 4.0832924272846389,
        4.0484347689392806, 4.0135297881451493, 3.9786222167929788,
        3.9437567866474876, 3.9089782295217028, 3.8743312772291874,
        3.8398606615622608, 3.8056111143575362, 3.7716273673878238,
        3.7379541524670392, 3.7046362014186909, 3.671718246028699,
        3.6392450181055418, 3.6072612494755023, 3.5758116719077737,
        3.5449410172469191, 3.514694017279151, 3.4851154038106515,
        3.456249908653271, 3.4281422635916137, 3.4008372005164293,
        3.3743794509818912, 3.3488137474719926, 3.3241848199391315,
        3.3005374053156897, 3.2779162212650537, 3.2563660386364717,
        3.235931481251312, 3.2166575797398647, 3.1985882418311387,
        3.1817702617296133, 3.16623120196501, 3.1519508688363445,
        3.1388898369784441, 3.1270115675057757, 3.1162783985912048,
        3.1066530752810793, 3.0980981955598694, 3.0905764105727918,
        3.0840503522279734, 3.0784826594047465, 3.0738359684941914,
        3.0700729167225984, 3.0671561410463664, 3.0650482785352748,
        3.0637119662036674, 3.0631098410527544, 3.0632045401646586,
        3.0639587005184103, 3.0653349591616283, 3.067295953133812,
        3.0698043194263716, 3.0728226951076074, 3.0763137171876989,
        3.0802400226896487, 3.0845642486538321, 3.0892490320945121,
        3.0942570100503257, 3.0995508195415122, 3.105093097601713,
        3.1108464812391743, 3.1167736075809378, 3.1228371133444424,
        3.1289996364108048, 3.1352238114321564, 3.1414722820714767,
        3.1477076669404993, 3.1538926540261176, 3.1599897394192689,
        3.1659619487174981, 3.1717709464907911, 3.1773874657083949,
        3.1828047582584804, 3.1880251445360597, 3.1930495837733637,
        3.1978795647752958, 3.2025163844769309, 3.2069614091332577,
        3.2112159799645688, 3.215281447233989, 3.21915915793746,
        3.2228504602405357, 3.2263567019088466, 3.2296792308384128,
        3.23281939486695, 3.2357785418700105, 3.2385580197131669,
        3.2411591762582215, 3.2435833593640733, 3.2458319168815395,
        3.2479061966882381, 3.2498075466519678, 3.2515373146267548,
        3.2530968484542693, 3.25448749602743, 3.2557106052040021,
        3.2567675238319613, 3.2576595997827771, 3.2583881809070268,
        3.2589546150994129, 3.2593602501875272, 3.2596064340366206,
        3.2596945145342531, 3.2596258395220592, 3.2594017568564646,
        3.2590236144194282, 3.2584927600657094, 3.2578105416662773,
        3.2569783070545237, 3.255997404112712, 3.2548691807128169,
        3.2535949847001864, 3.2521761639431261, 3.2506140663095167,
        3.2489100396512138, 3.2470654318297747, 3.2450815907282822,
        3.2429598641813291, 3.2407016000652566, 3.2383081462475825,
        3.2357808505778167, 3.2331210609239429, 3.2303301251556831,
        3.2274093911160535, 3.2243602066940427, 3.2211839197296803,
        3.2178818780954241, 3.2144554296536807, 3.2109059222526479,
        3.20723470377946, 3.2034431220698112, 3.1995325250175144,
        3.1955042604478874, 3.1913596762592644, 3.1871001202832185,
        3.1827269403958534, 3.1782414844848046, 3.173645100348101,
        3.16893913590854, 3.1641249390011361, 3.1592038574962462,
        3.154177239269381, 3.1490464321206364, 3.1438127840216281,
        3.1384776426505216, 3.1330423563533221, 3.127508271618602,
        3.1218767400737235, 3.1161490992255452, 3.1103267254806295,
        3.1044108952349481, 3.0984035516996222, 3.0923082943013216,
        3.0861293893357877, 3.0798710029916139, 3.0735373404797497,
        3.067132592794227, 3.0606609560942708, 3.0541266247105945,
        3.0475337935678537, 3.0408866574188971, 3.0341894110655989,
        3.0274462492919993, 3.0206613669084734, 3.0138389586774936,
        3.0069832194205333, 3.0000983439000679, 2.9931885269265637,
        2.9862579632853694, 2.9793108477625472, 2.9723513751564479,
        2.9653837402554024, 2.9584121378354928, 2.9514407627196246,
        2.9444738096772025, 2.9375154734884994, 2.93056994898431,
        2.9236414309002852, 2.9167341140846021, 2.9098521932891668,
        2.9029998633049297, 2.8961813189532171, 2.8894007549997323,
        2.8826623662453659, 2.8759703474568923, 2.8693288934679586,
        2.8627421990428767, 2.8562144589738581, 2.8497498680547491,
        2.8433526210830191, 2.837026912834181, 2.8307769381107493,
        2.82460689171332, 2.81852096840857, 2.8125233630061905,
        2.8066182702822462, 2.8008098850468519, 2.7951024020756607,
        2.7895000161605936, 2.7840069220985861, 2.7786273146782192,
        2.7733653887021048, 2.7682253389283469, 2.7632113601809896,
        2.7583276472395264, 2.7535783948965222, 2.748967797937794,
        2.7445000511566913, 2.7401793493475317, 2.7360098872861252,
        2.7319958597941962, 2.72814146162901, 2.724450887620363,
        2.7209283325070168, 2.7175779911265026, 2.7144040582552384,
        2.7114107286626536, 2.7086021971758476, 2.7059826585605271,
        2.7035563076123852, 2.7013273391305979, 2.6992999478905344,
        2.6974783286925796, 2.6958666763825665, 2.6944691855234324,
        2.6932900515296927, 2.6923334675744637, 2.6916036327818498,
        2.6911047300285995, 2.690840987066291, 2.6908165078083464,
        2.691035714435241, 2.6915009089067721, 2.6922091282242748,
        2.6931552891941495, 2.6943346268732542, 2.6957422524877588,
        2.6973733221196134, 2.6992229756806427, 2.7012863588555116,
        2.7035586152877076, 2.7060348893971331, 2.7087103252455584,
        2.7115800670561776, 2.7146392590166477, 2.7178830452815941,
        2.7213065700721852, 2.7249049775422081, 2.7286734118716671,
        2.7326070172689887, 2.7367009378897356, 2.7409503179313925,
        2.7453503015751717, 2.7498960330089415, 2.7545826563831555,
        2.7594053159177117, 2.7643591557869396, 2.7694393201561383,
        2.77464095323146, 2.7799591991580446, 2.7853892021679081,
        2.7909261064117343, 2.7965650560595803, 2.802301195343321,
        2.8081296683963504, 2.8140456194192782, 2.8200441926142088,
        2.8261205321143916, 2.8322697821467036, 2.8384870868758649,
        2.8447675904869185, 2.8511064371626329, 2.8574987710834425,
        2.8639397364373975, 2.8704244773992329, 2.8769481381540847,
        2.8835058628865298, 2.8900927957711757, 2.8967040810027327,
        2.903334862750123, 2.9099802852099774, 2.9166354925637941,
        2.923295628975306, 2.9299558386436058, 2.9366112657421368,
        2.9432570544704011, 2.9498883489779675, 2.9565002934874078,
        2.96308803215446, 2.9696467091537486, 2.976171468701621,
        2.9826574549494742, 2.9890998120974923, 2.9954936843212039,
        3.0018342157865452, 3.0081165507126011, 3.0143358332631753,
        3.0204872076019322, 3.0265658179439496, 3.0325668084377138,
        3.0384853233007405, 3.0443165066918172, 3.0500555028107077,
        3.0556974557902854, 3.0612375099542146, 3.06667080910994,
        3.0719924984834863, 3.0771977193296931, 3.082281623922289,
        3.0872393341123305, 3.0920660556643385, 3.0967567785807746,
        3.1013079310014637, 3.1057195125791348, 3.1099929612235786,
        3.1141294989206, 3.1181304316825091, 3.1219970350860735,
        3.1257305956865697, 3.1293323960764146, 3.1328037202894521,
        3.1361458518363965, 3.1393600744014458, 3.1424476716328447,
        3.1454099271793381, 3.14824812467773, 3.1509635477741336,
        3.1535574801256785, 3.1560312053815043, 3.1583860071770031,
        3.1606231691738205, 3.1627439750039361, 3.1647497083421716,
        3.1666416527988255, 3.1684210920458988, 3.1700893097145797,
        3.1716475894610983, 3.1730972149636854, 3.1744394698130334,
        3.175675637674725, 3.1768070022339532, 3.1778348470870941,
        3.1787604559130176, 3.1795851123440246, 3.1803101000317757,
        3.180936702636572, 3.1814662037781, 3.1818998871370776,
        3.182239036328848, 3.1824849350268312, 3.1826388668551817,
        3.1827021154839166, 3.1826759645604272, 3.1825616977019808,
        3.1823605985729522, 3.1820739508393485, 3.1817030381223463,
        3.181249144084116, 3.1807135523738697, 3.180097546622378,
        3.1794024104929979, 3.1786294276248142, 3.1777798816766576,
        3.1768550562762345, 3.1758562350848183, 3.1747847017525666,
        3.1736417399180206, 3.1724286332305622, 3.1711466653481,
        3.1697971199001911, 3.16838128053629, 3.1669004309257458,
        3.1653558547079914, 3.1637488354988461, 3.1620806569729218,
        3.1603526027959141, 3.1585659565794133, 3.1567220019936193,
        3.1548220226701407, 3.1528673022740974, 3.1508591244346835,
        3.1487987728056424, 3.1466875310480815, 3.1445266827885741,
        3.1423175116765036, 3.1400613013827097, 3.1377593355449216,
        3.1354128977753613, 3.1330232717719122, 3.1305917411561888,
        3.1281195895791658, 3.12560810069223, 3.1230585581302641,
        3.1204722455616789, 3.1178504466101837, 3.1151944449468152,
        3.112505524196258, 3.1097849680221037, 3.1070340600679809,
        3.1042540839774428, 3.1014463234120258, 3.098612061985583,
        3.0957525833745669, 3.092869171239633, 3.0899631091764652,
        3.0870356808866872, 3.0840881699956841, 3.0811218601313066,
        3.0781380349731133, 3.0751379781554147, 3.0721229733311941,
        3.0690943041287579, 3.0660532542076488, 3.0630011072215844,
        3.0599391468315993, 3.0568686566323633, 3.0537909203236384,
        3.0507072215502582, 3.0476188439123364, 3.0445270711096617,
        3.04143318676227, 3.0383384745147413, 3.0352442180538279,
        3.0321517009642949, 3.0290622069294977, 3.02597701961392,
        3.0228974226394754, 3.0198246996529874, 3.0167601343082433,
        3.0137050102488, 3.0106606111275918, 3.0076282206061511,
        3.00460912228931, 3.0016045998600323, 2.9986159369624077,
        2.9956444172411105, 2.9926913243260205, 2.9897579418969271,
        2.98684555357113, 2.9839554430011637, 2.9810888938595346,
        2.9782471897713503, 2.9754316143765451, 2.9726434513455313,
        2.96988398431675, 2.9671544969172015, 2.9644562728271047,
        2.9617905956765873, 2.9591587491346161, 2.9565620167901492,
        2.9540016823569104, 2.9514790294603817, 2.9489953417183945,
        2.9465519028290283, 2.944149996399855, 2.9417909060930723,
        2.9394759155613559, 2.9372063084414428, 2.934983368375466,
        2.9328083790356159, 2.9306826240628214, 2.9286073870692304,
        2.9265839517485643, 2.9246136017326996, 2.9226976206493971,
        2.9208372921525387, 2.9190338999695764, 2.9172887274636827,
        2.9156030591384745, 2.9139781761339618, 2.912415369157936,
        2.9109159024003137, 2.9094811080995218, 2.9081118657964082,
        2.9068079310752868, 2.9055686069055984, 2.9043932642049213,
        2.90328124741381, 2.9022319105995789, 2.9012446043253695,
        2.9003186804682906, 2.8994534903741447, 2.8986483855918967,
        2.8979027176257657, 2.8972158379842714, 2.8965870981589115,
        2.8960158496666333, 2.8955014439982079, 2.8950432326596118,
        2.8946405671739455, 2.8942927990162786, 2.8939992797161773,
        2.8937593607522643, 2.8935723936557842, 2.89343772990442,
        2.8933547210205033, 2.8933227185008286, 2.8933410738525103,
        2.8934091385843135, 2.893526264183941, 2.8936918021573854,
        2.893905104020384, 2.8941655212887922, 2.894472405439239,
        2.8948251079770642, 2.89522298042622, 2.8956653742794689,
        2.8961516410444204, 2.8966811322119637, 2.8972531992992847,
        2.897867193814081, 2.8985224672342449, 2.8992183711105191,
        2.8999542569037087, 2.9007294761321361, 2.9015433803054202,
        2.9023953209175573, 2.9032846494837607, 2.904210717489613,
        2.9051728764640448, 2.9061704778939887, 2.9072028732903252,
        2.9082694141525653, 2.9093694519787729, 2.9105023382819413,
        2.9116674245731806, 2.9128640623402604, 2.9140916031084845,
        2.9153493983514496, 2.9166367995870481, 2.9179531583351332,
        2.919297826076483, 2.9206701543332754, 2.9220694945870447,
        2.9234951983726161, 2.9249466171731586, 2.9264231024803666,
        2.9279240058215339, 2.9294486787064482, 2.9309964726050395,
        2.9325667390627466, 2.9341588295439558, 2.9357720955736895,
        2.9374058886755692, 2.9390595603038849, 2.9407324620025825,
        2.9424239452642347, 2.94413336157265, 2.9458600624763669,
        2.9476033994483934, 2.9493627239839069, 2.9511373876261873,
        2.952926741825221, 2.9547301381213482, 2.956546928015261,
        2.9583764630111267, 2.9602180945954237, 2.9620711743052532,
        2.96393505360129, 2.9658090840129865, 2.9676926170546647,
        2.9695850042091632, 2.9714855969835248, 2.9733937468984366,
        2.9753088054411858, 2.9772301241070931, 2.9791570544285677,
        2.9810889478899143, 2.9830251560083112, 2.9849650302556658,
        2.9869079221779362, 2.9888531832491632, 2.9908001649989133,
        2.9927482189063443, 2.99469669647091, 2.9966449492455927,
        2.9985923286714495, 3.0005381862886371, 3.0024818736018521,
        3.0044227420898584, 3.0063601432815839, 3.0082934286704988,
        3.010221949773666, 3.0121450580819014, 3.0140621050699496,
        3.0159724423058063, 3.0178754212623704, 3.0197703934245874,
        3.0216567103188687, 3.0235337234322248, 3.0254007843042112,
        3.0272572444021804, 3.0291024552411323, 3.0309357683281792,
        3.0327565351538817, 3.0345641072611631, 3.0363578360984484,
        3.0381370732061037, 3.0399011700863849, 3.0416494782325776,
        3.0433813491639445, 3.0450961343303566, 3.0467931853367007,
        3.0484718535969835, 3.050131490640982, 3.0517714480045619,
        3.0533910771305592, 3.0549897295814397, 3.0565667568369603,
        3.0581215103869712, 3.0596533417572012, 3.0611616024585295,
        3.0626456439528096, 3.0641048177995645, 3.0655384754631596,
        3.0669459684537692, 3.0683266482873419, 3.0696798664596945,
        3.0710049744849552, 3.0723013238437926, 3.0735682660545258,
        3.0748051526348132, 3.0760113350646257, 3.0771861648799943,
        3.0783289934734954, 3.0794391726915977, 3.0805160531452782,
        3.0815589887299994, 3.0825673243128255, 3.0835404298054381,
        3.084477605752, 3.0853783442646119, 3.0862416448683367, 3.08706978939093,
        3.0878732103129258, 3.0886656224726368, 3.0894602480094213,
        3.0902705007718709, 3.0911097251383151, 3.0919912905627083,
        3.0929285574930163, 3.0939348896066203, 3.0950236493824312,
        3.096208199796783, 3.0975019035774358, 3.0989181235968153,
        3.1004702226522456, 3.10217156355311, 3.1040355091334311,
        3.1060754221945417, 3.1083046655848166, 3.1107366020880169,
        3.113384594537187, 3.1162620057515866, 3.1193821985385388,
        3.1227585357322951, 3.1264043801528278, 3.1303330945897656,
        3.1345580418967325, 3.1390925848603839, 3.1439500863164982,
        3.1491439090813622, 3.1546874159758, 3.1605939698217709,
        3.1668769333620288, 3.1735496696754621, 3.180625540881262,
        3.1881179116438965, 3.1960401397418936, 3.2044056019452376,
        3.2132276224428273, 3.2225196709275092, 3.2322948154949374,
        3.2425671565900549, 3.2533439165459463, 3.2646152383079952,
        3.2763643866949144, 3.2885756589182371, 3.3012329505562996,
        3.314320302666022, 3.3278217037905238, 3.3417211613892075,
        3.3560026761376127, 3.3706502511507135, 3.3856478886172927,
        3.4009795911073382, 3.4166293610758718, 3.4325812009338632,
        3.4488191131798605, 3.4653271002577513, 3.4820891646355028,
        3.4990893087765427, 3.5163115351309329, 3.5337398461538347,
        3.5513582443090881, 3.5691507320832625, 3.58710131189702,
        3.6051939862215847, 3.6234127575399033, 3.6417416282694086,
        3.6601646009143654, 3.6786656779024405, 3.6972288617011597,
        3.7158381547777362, 3.7344775596104691, 3.7531310785519665,
        3.771782714386982, 3.7904164686675257, 3.8090163463745039,
        3.82756634298723, 3.8460504803174489, 3.8644527071549164,
        3.8827571740617444, 3.9009475128710935, 3.9190108111721242,
        3.9369427377658286, 3.9547424171461047, 3.9724084551331758,
        3.9899396593594152, 4.007334764300456, 4.0245925309091808,
        4.0417117105361875, 4.058691058025687, 4.0755293269448423,
        4.092225271310161, 4.1087776450302931, 4.1251852019842739,
        4.1414466960972129, 4.1575608813013316, 4.1735265114693725,
        4.189342340529433, 4.2050071223862329, 4.2205196109547947,
        4.2358785601602635, 4.2510827238683744, 4.2661308560159181,
        4.2810217105315074, 4.2957540413146891, 4.3103266022383417,
        4.3247381472479036, 4.3389874302542681, 4.3530732051678083,
        4.3669942258769012, 4.3807492462988415, 4.3943370203578178,
        4.407756301962177, 4.4210058450089589, 4.4340844034040252,
        4.4469907310838153, 4.4597235819310814, 4.4722817098684775,
        4.4846638688042635, 4.4968688126524761, 4.50889529531452,
        4.5207420706973256, 4.5324078927278766, 4.5438915152999915,
        4.5551916923308973, 4.5663071777374613, 4.5772367254059292,
        4.5879790892681651, 4.5985330232286747, 4.60889728119438,
        4.6190706170845273, 4.6290517847811241, 4.6388395382410694,
        4.6484326313364246, 4.6578298179748314, 4.667029852091825,
        4.676031487593785, 4.6848334783644505, 4.6934345783399154,
        4.7018335414175825, 4.7100291215143226, 4.7180200725371986,
        4.7258051483956205, 4.7333831029939537, 4.7407526902582333,
        4.7479126640829037, 4.7548617783675793, 4.76159878707012,
        4.76812244403273, 4.7744315032230862, 4.780524718523596,
        4.7864008438442234, 4.792058633101834, 4.7974968402287637,
        4.8027142189967149, 4.8077095237807184, 4.8124815071694984,
        4.8170289266507531, 4.8213505263844683, 4.8254450870933754,
        4.829311295544958, 4.8329484650493244, 4.8363574651638679,
        4.8395397920585275, 4.8424968478170118, 4.8452300711967142,
        4.8477408876411854, 4.8500307273715668, 4.8521010189236247,
        4.8539531914338507, 4.8555886738193887, 4.8570088950633927,
        4.8582152841402779, 4.8592092700277965, 4.85999228168864,
        4.8605657481226867, 4.8609310982452971, 4.8610897611087234,
        4.8610431656210054, 4.8607927407927969, 4.8603399155895692,
        4.859686118953582, 4.858832779934863, 4.857781327419481,
        4.8565331904315112, 4.8550897979456584, 4.8534525789054248,
        4.8516229623078209, 4.8496023771224994, 4.8473922523140711,
        4.8449940168684167, 4.84240909976229, 4.839638929944269,
        4.8366849364096529, 4.8335485481391771, 4.8302311940706666,
        4.82673430321859, 4.8230593045387495, 4.8192076269951558,
        4.8151806995747988, 4.8109799512462, 4.8066068109927214,
        4.8020627077689566, 4.7973490705640272, 4.7924673283373549,
        4.7874189101060844, 4.7822052447580949, 4.7768277613422692,
        4.7712878888212309, 4.7655870561419658, 4.7597266922888952,
        4.7537082262394952, 4.7475330869765608, 4.7412027034471151,
        4.7347185046466782, 4.7280819195441941, 4.7212943771222466,
        4.714357306319938, 4.707272136161774, 4.70004029557762,
        4.6926632135538551, 4.6851423190812582, 4.6774790411063618,
        4.6696748086339426, 4.6617310506164538, 4.6536491960112976,
        4.6454306738458744, 4.6370769130316134, 4.6285893425725311,
        4.6199693914554061, 4.6112184886349663, 4.60233806307,
        4.5933295437684531, 4.5841943596864114, 4.57493393981591,
        4.5655497130266616, 4.55604310866637, 4.5464155545732883,
        4.53666848290333, 4.5268033138928638, 4.5168215005462322,
        4.5067244116900875, 4.4965139771133229, 4.486193519764015,
        4.4757669235505366, 4.4652379881834445, 4.4546105461249805,
        4.4438884179758666, 4.4330754286597021, 4.4221754014926189,
        4.4111921604202555, 4.4001295290867324, 4.3889913312740969,
        4.377781390784639, 4.3665035313265621, 4.3551615766639848,
        4.343759350564544, 4.3323006767672076, 4.32078937904183,
        4.3092292811330024, 4.2976242068116077, 4.2859779798049713,
        4.2742944238820559, 4.2625773628126069, 4.2508306203404373,
        4.2390580202207584, 4.2272633861814946, 4.2154505420358221,
        4.203623311496389, 4.1917855183260437, 4.1799409862804158,
        4.1680935391227321, 4.1562470006061378, 4.1444051944754046,
        4.1325719444907563, 4.1207510744163969, 4.10894640798711,
        4.0971617689725743, 4.0854009811489274, 4.0736678682300642,
        4.0619662539928036, 4.0502999621888316, 4.038672816583083,
        4.0270886409025231, 4.0155512589304569, 4.0040644944109047,
        3.9926321711067341, 3.9812581127737583, 3.9699461431477214,
        3.958700085993184, 3.947523765076693, 3.9364210041453984,
        3.9253956269628638, 3.9144514572633056, 3.9035923188047486,
        3.8928220353812737, 3.8821444307018966, 3.8715633285371007,
        3.8610825526428703, 3.8507059267819779, 3.8404372746974209,
        3.8302804201516287, 3.8202391868917158, 3.81031739870011,
        3.8005188792823956, 3.7908474524312674, 3.7813069419112844,
        3.771901171432078, 3.762633964791676, 3.7535091457317549,
        3.7445305379991489, 3.73570196534411, 3.7270272515342975,
        3.7185102203802174, 3.7101546954091016, 3.701964500996517,
        3.6939434592122886, 3.6860953984858429, 3.6784241295830311,
        3.6709335122359432, 3.6636272710357023, 3.6565094778465781,
        3.64958189100164, 3.6428405238774468, 3.6362790763085622,
        3.6298915954037945, 3.623671993101957, 3.6176142303990888,
        3.6117122505551511, 3.6059600031486787, 3.6003514355780504,
        3.5948804959853842, 3.5895411322307686, 3.5843272923193483,
        3.5792329241568313, 3.5742519756989468, 3.5693783949095956,
        3.5646061297391793, 3.5599291281108125, 3.5553413379788772,
        3.5508367073223752, 3.5464091840203804, 3.5420527161015669,
        3.5377612514523276, 3.5335287380420306, 3.5293491238436574,
        3.525216356742491, 3.5211243847554967, 3.5170671557920352,
        3.5130386177772133, 3.5090327187151478, 3.5050434065145315,
        3.5010646291323009, 3.4970903345214883, 3.4931144706319932,
        3.4891309853867623, 3.4851338267374934, 3.4811169426877591,
        3.477074281114013, 3.4729997899988403, 3.4688874172758473,
        3.464731110905658, 3.4605248188240814, 3.4562624889954043,
        3.4519380693353527, 3.4475455078242421, 3.4430787524051532,
        3.4385317509947377, 3.4338984515815518, 3.4291728020843766,
        3.4243487504582015, 3.419420244670015, 3.4143812326508183,
        3.4092256623129598, 3.4039474816775339, 3.3985406386304571,
        3.3929990811437509, 3.3873167571764431, 3.3814876146458563,
        3.3755056015378457, 3.3693646657481011, 3.3630587552722835,
        3.3565818180367786, 3.3499278019850234, 3.3430906550948465,
        3.3360643252622606, 3.3288427604672037, 3.3214199086784446,
        3.3137897177991094, 3.3059461357850655, 3.2978831106041349,
        3.2895945902071326, 3.2810745224602029, 3.272316855558274,
        3.2633155368503428, 3.2540645157832131, 3.2445577362381495,
        3.2347891574278682, 3.2247526961284887, 3.2144423865436935,
        3.203851937846534, 3.1929759565328477, 3.1818067423993335,
        3.1703519643786486, 3.1586574564611123, 3.1467844218153918,
        3.1347917568355763, 3.1227392553874367, 3.1106863860867042,
        3.0986927351451841, 3.0868178462980089, 3.0751212785816384,
        3.06366258553108, 3.0525013226491953, 3.041697044747198,
        3.0313093068938493, 3.0213976640161451, 3.0120216711633954,
        3.0032408832844819, 2.9951148553903453, 2.9877031424762643,
        2.9810652994933933, 2.9752608814693371, 2.9703494433750692,
        2.9663905402353214, 2.9634437269858704, 2.9615685586380187,
        2.9608245902074066, 2.9612713766570931, 2.9629684729611485,
        2.9659754341750504, 2.9703518151926165, 2.9761571711577828,
        2.9834510566594017, 2.9922930278376252, 3.0027426364916487,
        3.01485944640038, 3.0287029882675705, 3.0443328842784778,
        3.061808503452315, 3.0811899154405231, 3.1025352559693777,
        3.1259076319546253, 3.1513370287653344, 3.1787711838667252,
        3.2081247132077921, 3.2393172038779192, 3.2722663090748023,
        3.3068903826208209, 3.3431075251825373, 3.380835928888736,
        3.4199937528090438, 3.4604991679703381, 3.5022703410792961,
        3.5452254403899164, 3.5892826336360084, 3.63436008868667,
        3.68037597339423, 3.7272484556243759, 3.77489570318943,
        3.8232358840069893, 3.8721871658391565, 3.9216677166275988,
        3.9715957041669578, 4.021889296336961, 4.0724666609725633,
        4.123245965912921, 4.1741453790524714, 4.2250830682100391,
        4.27597720124932, 4.3267459460064437, 4.3773074703741708,
        4.4275799421026489, 4.4774815293197587, 4.52693039907861,
        4.5758447214899478, 4.6241426582041587, 4.6717423940620924,
        4.7185620500354153, 4.7645199236742037, 4.809533823982779,
        4.8535229082854343, 4.89640286811784, 4.9381124871153723,
        4.978647891974715, 5.0180283016004523, 5.0562694689671908,
        5.0933884954666615, 5.1294019939043585, 5.1643267536537349,
        5.1981795003195339, 5.2309769825024048, 5.2627359405255687,
        5.2934731176946368, 5.3232052562407013, 5.3519490987706044,
        5.3797213877489227, 5.4065388657397992, 5.4324182752282262,
        5.4573763587238036, 5.4814298587455141, 5.5045955178175214,
        5.5268900784343646, 5.5483302831185064, 5.5689328743856565,
        5.5887145947451531, 5.6076921867110014, 5.6258823927940389,
        5.6433019555273081, 5.65996761739307, 5.6758961209094725,
        5.6911042086204278, 5.7056086229922842, 5.7194261065993555,
        5.7325734018950367, 5.7450672514044916, 5.7569243976769684,
        5.7681615832010351, 5.7787955504808943, 5.788843042048387,
        5.798320800408086, 5.8072455680571213, 5.8156340875492543,
        5.8235031013531193, 5.8308693520295529, 5.8377495820419476,
        5.8441605339182461, 5.8501189502052622, 5.8556415733797493,
        5.8607451459626745, 5.8654464104880129, 5.8697621094245074,
        5.8737089853325637, 5.87730378069876, 5.8805632380384543,
        5.8835040998899739, 5.8861431087209439, 5.8884970070876381,
        5.89058253747517, 5.8924164424264616, 5.8940154644068325,
        5.89539634596944, 5.89657582963168, 5.897570657882178,
        5.8983975732264833, 5.8990733182232251, 5.8996146353489625,
        5.9000382671086786, 5.9003609560611441, 5.9005994446554926,
        5.9007704754733927, 5.9008907909778507, 5.9009771337125025,
        5.9010462461377973, 5.9011148709292458, 5.901199750134082,
        5.9013176275429471, 5.9014852421873849, 5.9017193461755628,
        5.90203665543404, 5.902453986127421, 5.9029878775704763,
        5.9036555809130489, 5.9044696042298543, 5.9054306777337988,
        5.9065347885259447, 5.9077786355918773, 5.909158641051917,
        5.910671327244116, 5.9123131803463531, 5.9140806995689834,
        5.9159703794476375, 5.9179787162220343, 5.9201022054322827,
        5.9223373429324324, 5.924680624467892, 5.92712854584631,
        5.9296776027645661, 5.9323242910165961, 5.9350651063838278,
        5.9378965446128973, 5.9408151014575354, 5.943817272713467,
        5.9468995541259684, 5.9500584414611968, 5.953290430490342,
        5.9565920169483455, 5.95995969666154, 5.9633899653409213,
        5.9668793187645806, 5.9704242527262439, 5.9740212629217924,
        5.9776668452035713, 5.9813574952702284, 5.9850897089084594,
        5.9888599819100845, 5.9926648099861906, 5.9965006889426178,
        6.0003641145190842, 6.0042515825162832, 6.008159588651556,
        6.0120846287060781, 6.0160231984838086, 6.0199717937023536,
        6.0239269101379609, 6.027885043571632, 6.0318426897285935,
        6.0357963444355409, 6.0397425034120173, 6.0436776624177808,
        6.0475983172480658, 6.0515009636467054, 6.0553820973838892,
        6.0592382142359913, 6.0630658099583732, 6.0668613802953129,
        6.0706214210485756, 6.0743424279681246, 6.0780208967983542,
        6.0816533233414454, 6.0852362033388276, 6.0887660325437221,
        6.0922393067548963, 6.0956525217006341, 6.0990021731877331,
        6.1022847569415806, 6.105496768733639, 6.1086347043612621,
        6.1116950595515727, 6.1146743300822433, 6.1175690117237123,
        6.1203756002456728, 6.1230905913893636, 6.1257104809423089,
        6.1282317646416757, 6.1306509383496106, 6.1329644975342674,
        6.1351689387487287, 6.1372607556346939, 6.1392364498239917,
        6.1410925008286306, 6.142825449216307, 6.1444316786312552,
        6.1459086184841611, 6.1472562951613865, 6.1484757809011477,
        6.1495679909403771, 6.1505339016087337, 6.1513744670866481,
        6.1520906495697094, 6.1526834083708275, 6.1531537037794184,
        6.1535024958099118, 6.1537307445554692, 6.153839410047576,
        6.1538294523730137, 6.1537018315981289, 6.1534575077566318,
        6.1530974409616919, 6.152622591265402, 6.1520339186787361,
        6.1513323833349833, 6.15051894526672, 6.1495945645519567,
        6.1485602012393006, 6.1474168154286648, 6.1461653671179945,
        6.1448068164077743, 6.1433421234227072, 6.141772248136868,
        6.1400981506487646, 6.1383207910337685, 6.1364411293376637,
        6.1344601256565694, 6.1323787400338361, 6.1301979325099243,
        6.1279186631914158, 6.1255418921412925, 6.1230685793838475,
        6.1204996850202713, 6.1178361691293679, 6.1150789917054658,
        6.1122291128908213, 6.10928749270932, 6.1062550912364983,
        6.1031328685420139, 6.0999217846804994, 6.0966227997193272,
        6.09323687373565, 6.0897649667786453, 6.0862080389004172,
        6.0825670502151654, 6.0788429607390437, 6.07503673055361,
        6.0711493197309974, 6.0671816883382128, 6.0631347964129274,
        6.0590096040371479, 6.0548070713004059, 6.0505281582121624,
        6.0461738248940886, 6.0417450313810033, 6.0372427377295246,
        6.0326679040419684, 6.0280214903492562, 6.0233044567197478,
        6.018517763235125, 6.0136623699322769, 6.0087392369029269,
        6.0037493242139419, 5.9986935919110929, 5.9935730000513,
        5.98838850874641, 5.9831410780067369, 5.9778316679174823, 5.972461238566,
        5.9670307499763888, 5.9615411622442513, 5.9559934354178488,
        5.95038852959744, 5.9447274047896324, 5.9390110210863689,
        5.9332403385729826, 5.9274163172967853, 5.9215399172970775,
        5.9156120987069531, 5.9096338215098649, 5.9036060458267574,
        5.8975297317037505, 5.8914058391816377, 5.8852353283950967,
        5.8790191593263188, 5.8727582920924783, 5.8664536867454693,
        5.8601063033424277, 5.8537171019641185, 5.84728704266458,
        5.8408170855090225, 5.8343081905555518, 5.8277613178595296,
        5.8211774275529153, 5.81455747961198, 5.8079024341401926,
        5.8012132512169243, 5.794490890877146, 5.787736313230802,
        5.78095047827155, 5.7741343461197161, 5.7672888768487516,
        5.7604150304740722, 5.7535137670966128, 5.7465860467772671,
        5.7396328295716321, 5.7326550755288679, 5.7256537447532354,
        5.7186297972991547, 5.7115841931939029, 5.7045178925607516,
        5.6974318554172161, 5.6903270418373761, 5.683204411907246,
        5.6760649256674327, 5.6689095432081524, 5.6617392245623357,
        5.6545549298181523, 5.6473576190378045, 5.6401482522721453,
        5.6329277896221823, 5.6256971911049005, 5.6184574167938237,
        5.6112094268083172, 5.6039541811460065, 5.596692639899091,
        5.5894257631321018, 5.5821545108969062, 5.5748798432915994,
        5.56760272037007, 5.560324102161645, 5.5530449487512685,
        5.5457662202322222, 5.53848887662422, 5.5312138780310525,
        5.5239421844964527, 5.516674756091196, 5.5094125528636759,
        5.5021565349078658, 5.4949076622725821, 5.4876668950354279,
        5.4804351932294786, 5.4732135169182365, 5.4660028262438463,
        5.45880408119149, 5.45161824184421, 5.4444462682713208,
        5.4372891205444986, 5.4301477587360409, 5.4230231428653255,
        5.4159162330771951, 5.4088279892791453, 5.4017593719930925,
        5.39471133988679, 5.38768485706159, 5.3806808723787247,
        5.3737003738831568, 5.366744088941525, 5.359812097754161,
        5.352904219860644, 5.3460203138731721, 5.3391602232946624,
        5.3323237970425348, 5.3255108820987349, 5.3187213261522883,
        5.3119549765973968, 5.3052116809908174, 5.2984912867638485,
        5.2917936414735172, 5.28511859255875, 5.2784659875742515,
        5.2718356739321726, 5.2652274991551629, 5.2586413107935321,
        5.2520769562659275, 5.2455342830848775, 5.2390131387677181,
        5.23251337077043, 5.2260348266036685, 5.21957735377229,
        5.2131407997381407, 5.2067250120451156, 5.2003298381066276,
        5.1939551254807483, 5.1876007216443067, 5.18126647409684,
        5.1749522302981976, 5.1686578377466326, 5.1623831439881283,
        5.1561279964427467, 5.1498922426546363, 5.143675730084718,
        5.1374783062503635, 5.131299818631601, 5.1251401147161646,
        5.1189990419984444, 5.1128764480011588, 5.1067721801511707,
        5.1006860859842886, 5.0946180130186658, 5.0885678086954647,
        5.0825353205181063, 5.0765203960066492, 5.0705228826305433,
        5.0645426278850953, 5.05857947926785, 5.0526332842644495,
        5.0467038903700949, 5.0407911450664553, 5.0348948958966853,
        5.0290149902720351, 5.0231512757506884, 5.0173035998005853,
        5.0114718098965616, 5.0056557535506174, 4.9998552782983863,
        4.9940702315338807, 4.9883004608118862, 4.9825458136543643,
        4.9768061374758341, 4.9710812798281232, 4.96537108817364,
        4.9596754100211564, 4.9539940928651305, 4.9483269841769539,
        4.94267393145442, 4.93703478224193, 4.9314093839275728, 4.92579758409481,
        4.92019923022065, 4.91461416975336, 4.9090422502211988, 4.90348331910162,
        4.8979372239280314, 4.89240381213946, 4.88688293123267,
        4.8813744287344409, 4.8758781521188768, 4.8703939488734669,
        4.8649216664912034, 4.859461152469315, 4.854012254320013,
        4.8485748194770926, 4.8431486954989831, 4.83773372984952,
        4.8323297700206993, 4.8269366634886408, 4.8215542577705062,
        4.8161824003700682, 4.8108209387510836, 4.8054697204029981,
        4.8001285928378827, 4.79479740355535, 4.7894760000052985,
        4.78416422974823, 4.7788619401981594, 4.7735689789191005,
        4.7682851933369044, 4.7630104309868964, 4.7577445393816209,
        4.7524873659777178, 4.7472387582450022, 4.7419985637035333,
        4.736766629889626, 4.7315428042272476, 4.726326934222965,
        4.7211188673951767, 4.7159184512277177, 4.7107255331931563,
        4.7055399608202269, 4.7003615815528672, 4.6951902429286561,
        4.6900257924099069, 4.6848680775051443, 4.679716945716434,
        4.6745722445060656, 4.6694338213734952, 4.6643015238314511,
        4.65917519936537, 4.6540546954534365, 4.6489398596104579,
        4.6438305393085759, 4.6387265820386627, 4.6336278353328986,
        4.6285341466214662, 4.6234453634165167, 4.6183613332788607,
        4.6132819035970769, 4.6082069219151256, 4.6031362357346906,
        4.5980696925384121, 4.59300713979606, 4.58794842503108,
        4.5828933957227225, 4.5778418993707755, 4.5727937834584953,
        4.5677488954630743, 4.5627070829135343, 4.557668193286708,
        4.5526320740582538, 4.5475985727340644, 4.5425675368287726,
        4.5375388137847539, 4.5325122511410552, 4.527487696372682,
        4.522464996966213, 4.5174440004015253, 4.5124245542126307,
        4.50740650586202, 4.5023897028466306, 4.4973739926429044,
        4.4923592228086031, 4.4873452406595185, 4.4823318941995627,
        4.4773190294422394, 4.4723064980138787, 4.4672941360282419,
        4.4622818194328477, 4.4572691585615543, 4.4522551040317122,
        4.4472383407778464, 4.4422175935396186, 4.4371915716623231,
        4.4321589900519234, 4.4271185615788919, 4.42206899985731,
        4.4170090182255395, 4.4119373301544229, 4.4068526490258764,
        4.4017536882745736, 4.3966391613254912, 4.3915077816030337,
        4.386358262470182, 4.3811893174231864, 4.375999659836399,
        4.3707880031547566, 4.3655530607660946, 4.36029354608028,
        4.355008172595265, 4.34969565362696, 4.344354702673062,
        4.3389840331056808, 4.3335823583572912, 4.3281483918798713,
        4.322680847014654, 4.3171784372504245, 4.3116398760051231,
        4.3060638766359167, 4.3004491526241191, 4.2947944173628114,
        4.2890983842799768, 4.2833597667815146, 4.2775772782791286,
        4.2717496322283468, 4.2658755420331191, 4.2599537210843428,
        4.2539828828158779, 4.2479617406862049, 4.2418890080427092,
        4.2357633983632592, 4.2295836250313394, 4.22334840150674,
        4.2170564411503468, 4.2107064574121829, 4.204297163746193,
        4.1978272735145961, 4.1912955001544336, 4.1847005570771767,
        4.1780411577338077, 4.1713160155323967, 4.1645238438560943,
        4.1576633561360659, 4.1507332658478377, 4.1437322863385173,
        4.1366591310636585, 4.1295125134291295, 4.1222911468637387,
        4.1149937447757887, 4.1076190205924012, 4.100165687732237,
        4.0926324596257642, 4.0850180496382817, 4.0773211712439865,
        4.06954053787389, 4.0616748628805723, 4.0537228597480075,
        4.0456832418503748, 4.0375547226430211, 4.029336015502472,
        4.0210258338962159, 4.0126228911926809, 4.0041259008984031,
        3.9955335761843211, 3.9868446311423091, 3.978057777401061,
        3.9691717333285759, 3.9601851985573071, 3.9510969245168468,
        3.9419055296009793, 3.9326105183522642, 3.9232135959173395,
        3.9137173536414012, 3.9041242498919386, 3.8944367946806007,
        3.8846574793851167, 3.8747888020725547, 3.8648332584320522,
        3.8547933449718816, 3.8446715579139106, 3.834470393629112,
        3.8241923483825011, 3.8138399184462215, 3.8034156001725057,
        3.7929218898808119, 3.7823612838288843, 3.771736278322916,
        3.7610493696988136, 3.7503030542423881, 3.7394998282445662,
        3.7286421880513343, 3.7177326299178062, 3.7067736501947781,
        3.6957677451117363, 3.6847174110504541, 3.6736251443100856,
        3.6624934411501582, 3.6513247978691421, 3.6401217108168247,
        3.6288866763008669, 3.6176221905656534, 3.6063307499742931,
        3.595014850793556, 3.5836769893413685, 3.5723196619306883,
        3.5609453648458218, 3.549556594401571, 3.5381558469215908,
        3.5267456186586337, 3.5153284059537242, 3.5039067051131276,
        3.49248301245662, 3.4810598242274957, 3.4696396367683362,
        3.4582249463921269, 3.4468182493932065, 3.4354220420653334,
        3.42403882071116, 3.4126710816835022, 3.4013213212057858,
        3.3899920356455242, 3.3786857212553478, 3.3674048744099978,
        3.3561519913522435, 3.344929568389392, 3.3337401018527277,
        3.3225860880642935, 3.3114700232552687, 3.3003944037749511,
        3.2893617259577854, 3.2783744860510002, 3.267435180386085,
        3.2565463052532562, 3.2457103570004984, 3.2349298318614692,
        3.2242072261943044, 3.2135450362728775, 3.2029457584362633,
        3.1924118889434592, 3.1819459241180859, 3.1715503602880926,
        3.1612276937254937, 3.1509804207277257, 3.1408110376392924,
        3.1307220407024818, 3.1207159264068149, 3.1107951904631967,
        3.1009623308115781, 3.0912198393179859, 3.0815702233951141,
        3.0720158865071605, 3.0625589737528349, 3.053201526227908,
        3.0439456006129286, 3.0347932475301307, 3.0257465198303262,
        3.0168074695168277, 3.0079781488872022, 2.999260610196361,
        2.990656905664034, 2.9821690875279536, 2.9737992080407838,
        2.9655493194251936, 2.9574214739448692, 2.9494177238072736,
        2.941540121272888, 2.9337907185639538, 2.9261715679582005,
        2.9186847216314256, 2.9113322318659547, 2.9041161509031159,
        2.8970385309823148, 2.8901014242986833, 2.8833068831216715,
        2.8766569597270042, 2.8701537062974118, 2.8637991750886203,
        2.8575954183453911, 2.8515444883304921, 2.8456484372114468,
        2.8399093172972023, 2.8343291808164603, 2.8289100799718527,
        2.8236540670274364, 2.81856319422744, 2.8136395138077486,
        2.8088850780005665, 2.8043019390174506, 2.7998921491507827,
        2.7956577606189521, 2.7916008256557681, 2.7877233964996311,
        2.7840275253717719, 2.7805152645528755, 2.7771886662483247,
        2.7740497827240209, 2.7711006661861313, 2.76834336889819,
        2.7657799430809811, 2.7634124409925191, 2.7612429148746251,
        2.7592734169399247, 2.757505999439986, 2.7559427146102573,
        2.754585614711345, 2.7534367519701566, 2.7524981786016549,
        2.7517719468666404, 2.75126010901603, 2.7509647172764455,
        2.7508878238676147, 2.7510314810519674, 2.751397741075853,
        2.7519886561362439, 2.7528062785325864, 2.7538526604507703,
        2.7551298541643985, 2.7566399118865959, 2.7583848858673714,
        2.7603668283533431, 2.7625877915886705, 2.7650498277727018,
        2.767754989181924, 2.7707053280825442, 2.7739028966149837,
        2.7773497471061317, 2.7810479317639794, 2.7849995028611056,
        2.7892065125780494, 2.7936710131707918, 2.7983950569119123,
        2.803380696037, 2.8086299827259089, 2.8141449692599156,
        2.8199277079163423, 2.8259802508656557, 2.8323046503686764,
        2.8389029586882408, 2.8457772280386942, 2.8529295106646821,
        2.8603618588069195, 2.8680763246916472, 2.876074960604456,
        2.8843598187166162, 2.8929329513072366, 2.9017964106185787,
        2.9109522488818054, 2.9204025183134785, 2.9301492711822377,
        2.9401945597267654, 2.9505404361704479, 2.9611889527436412,
        2.9721421617050456, 2.9834021153118666, 2.9949708657611964,
        3.0068504653018708, 3.0190429661818547, 3.0315504206747432,
        3.0443748809408917, 3.0575183992729555, 3.0709830279116836,
        3.0847708190755423, 3.0988838250116073, 3.1133240979547896,
        3.1280936901579266, 3.1431946538483682, 3.158629041258691,
        3.1743989046256234, 3.1905062962337407, 3.2069532682541841,
        3.2237418729542111, 3.2408741626001025, 3.2583521894071295,
        3.2761780056022345, 3.2943536634289, 3.3128812151415188,
        3.3317627129974245, 3.3510002091629163, 3.3705957559417228,
        3.3905514055717041, 3.4108692102558229, 3.4315512222563052,
        3.4525994937992746, 3.4740160771471511, 3.4958030245088616,
        3.5179623881467781, 3.5404962202850312, 3.5634065731767763,
        3.5866954990351654, 3.6103650501368016, 3.6344172786909934,
        3.6588542369567478, 3.6836779771486796, 3.7088905515031358,
        3.734494012302263, 3.7604904117630018, 3.7868818021012727,
        3.8136702355731655, 3.8408577644295359, 3.868446440896482,
        3.8964383171613, 3.92483544569423, 3.9536398780732811,
        3.9828536684150579, 4.0124788635622473, 4.0425175308928853,
        4.0729716807655709, 4.1038434811253257, 4.1351346652834442,
        4.1668480835955961, 4.1989791440789128, 4.2315047737695384,
        4.2643944573905843, 4.2976187967223964, 4.3311479589055466,
        4.3649522685679631, 4.3990019934484348, 4.43326742184306,
        4.4677188345897054, 4.5023265152326335, 4.5370607463679509,
        4.5718918109235887, 4.6067899916339057, 4.6417255714030548,
        4.6766688330681205, 4.711590059404962, 4.7464595332603174,
        4.78124753751215, 4.8159243549407984, 4.8504602683740643,
        4.884825560687891, 4.9189905146675121, 4.9529254131816893,
        4.9866005390089683, 5.0199861750351831, 5.0530526040560932,
        5.0857701089449261, 5.1181089724715338, 5.1500394774680558,
        5.1815319068561934, 5.2125565433659871, 5.2430836699274392,
        5.273083569059132, 5.3025265245543416, 5.3313828166182757,
        5.3596227352954466, 5.3872165434239037, 5.4141345791625533,
        5.4403469723133764, 5.4658244275500989, 5.4905361721267516,
        5.51446127630481, 5.5376032525262193, 5.55997545625625,
        5.581589765639424, 5.6024586334715076, 5.6225943043956477,
        5.6420090983022275, 5.6607153078307473, 5.6787252354859,
        5.6960511802508496, 5.7127054423089518, 5.728700321506393,
        5.7440481176885223, 5.7587611307996056, 5.7728516607324636,
        5.7863320073605573, 5.799214470559356, 5.811511350243677,
        5.8232349463487578, 5.8343975587000028, 5.84501148722209,
        5.8550890318048685, 5.8646424923574445, 5.87368416873976,
        5.8822263608730019, 5.8902813686377211, 5.8978614919482233,
        5.9049790306524095, 5.9116462847001641, 5.9178755539700418,
        5.9236791382939167, 5.92906933766825, 5.9340584518497312,
        5.9386587810630536, 5.9428826242656445, 5.9467422839100532,
        5.9502500528313194, 5.9534182505522741, 5.956259122716129,
        5.9587851050075189, 5.96100736696334, 5.9629339341631074,
        5.96457156597892, 5.9659272119037938, 5.9670077474340291,
        5.9678200749365855, 5.9683710869611977, 5.968667679672051,
        5.9687167480009693, 5.9685251871667173, 5.96809989238322,
        5.9674477588514758, 5.9665756817820013, 5.9654905563676248,
        5.9641992778114021, 5.96270874128953, 5.9610258420566007,
        5.9591574752715033, 5.957110536128063, 5.9548919198798647,
        5.9525085216945, 5.94996723674461, 5.9472749602723782,
        5.9444385874817653, 5.94146501354983, 5.9383611336844471,
        5.9351338430794547, 5.9317900369617664, 5.9283366105160544,
        5.9247804589529611, 5.9211284774222674, 5.9173875612349729,
        5.91356460545619, 5.9096665054342, 5.90570015606794, 5.9016724535520391,
        5.8975902901134489, 5.8934605691845636, 5.8892901634087504,
        5.8850860244034013, 5.8808545773882557, 5.8766009402420645,
        5.8723297043554128, 5.8680455401305478, 5.86375308727303,
        5.8594569966014172, 5.85516191491663, 5.8508724904043543,
        5.8465933708691908, 5.8423292041862993, 5.8380846381931644,
        5.8338643207619363, 5.8296728997271208, 5.8255150230279629,
        5.8213953384440984, 5.8173184938681031, 5.813289137170492,
        5.8093119161981566, 5.8053914788197627, 5.8015324728904449,
        5.7977395462577634, 5.7940173468484266, 5.7903705224273452,
        5.7868037209101679, 5.7833215901829735, 5.779928778054563,
        5.7766299324114527, 5.77342970108977, 5.7703327320161142,
        5.7673436729961223, 5.764467171862691, 5.7617078766044925,
        5.7590704347897477, 5.756559495047755, 5.7541797031413173,
        5.7519357127197859, 5.7498321556135936, 5.747873724077861,
        5.7460649435178057, 5.7444107682791925, 5.7429132946961579,
        5.7415675218542779, 5.74036559084061, 5.7393000717486284,
        5.7383633677011838, 5.7375479423126725, 5.7368462373750333,
        5.7362507025834892, 5.7357537847003632, 5.73534793161191,
        5.735025590795046, 5.7347792098506192, 5.7346012363243224,
        5.7344841178291688, 5.7344203019301707, 5.7344022361743088,
        5.734422368166717, 5.7344731454516005, 5.7345470156606089,
        5.7346364263161282, 5.7347338249838185, 5.7348316592967006,
        5.734922376782988, 5.7349984250310779, 5.7350522516211981,
        5.7350763040988237, 5.735063030105314, 5.735004877136971,
        5.7348942928141344, 5.7347237247094105, 5.7344856203885444,
        5.734172427418442, 5.7337765933734337, 5.7332905658679731,
        5.73270679242701, 5.7320177206533582, 5.7312157981067875,
        5.7302934723778218, 5.7292431910278614, 5.7280574016298109,
        5.7267285517673781, 5.7252490890392842, 5.7236114609583968,
        5.7218081151429825, 5.7198314991704144, 5.717674060611265,
        5.715328247030099, 5.7127865059802891, 5.7100412850975282,
        5.707085031913449, 5.7039101939876629, 5.7005092189481532,
        5.6968745543219876, 5.6929986477107946, 5.6888739466789335,
        5.6844928987864, 5.6798479516673561, 5.6749315528154343,
        5.6697361498352405, 5.6642541903333381, 5.6584781218824434,
        5.6524003919975394, 5.646013448271499, 5.6393097383662978,
        5.6322817097679634, 5.6249218100363745, 5.6172224868189273,
        5.6091761876539739, 5.600775360100311, 5.5920124517610814,
        5.582879910171493, 5.5733701830220408, 5.56347571765399,
        5.5531889620460726, 5.542502362843817, 5.5314083702310644,
        5.5198994244375656, 5.5079679933211487, 5.4956064684974342,
        5.4828074518608707, 5.4695630046038932, 5.4558687912335868,
        5.4417294241390293, 5.427153119006964, 5.4121475507614942,
        5.3967206046945915, 5.3808800898818649, 5.3646338429022489,
        5.34798969045671, 5.3309554627889009, 5.3135389888632094,
        5.29574809808936, 5.27759061971809, 5.2590743831169142,
        5.2402072175535146, 5.2209969523032447, 5.20145141671053,
        5.1815784400505995, 5.1613858516534563, 5.1408814808199574,
        5.1200731568108315, 5.0989687089856721, 5.0775759665972977,
        5.0559027589969912, 5.033956915447976, 5.0117462652692915,
        4.989278637753209, 4.9665618622052019, 4.9436037679637774,
        4.9204121842845518, 4.8969949404817337, 4.8733598658754573,
        4.8495147897486879, 4.82546754142142, 4.8012259501737553,
        4.7767978453320081, 4.7521910562112177, 4.7274134120460936,
        4.70247274221199, 4.67737687597293, 4.6521336426754321,
        4.6267508715638144, 4.6012363919559212, 4.5755980332044617,
        4.5498436245547573, 4.5239809953201524, 4.49801797481532,
        4.4719623923703979, 4.4458220772219015, 4.4196048587168546,
        4.393318566163388, 4.3669710288526, 4.340570076060267, 4.314123537129551,
        4.2876392413605613, 4.26112501804227, 4.2345886964493413,
        4.2080381059374954, 4.1814810758052889, 4.1549254353050831,
        4.128379013787451, 4.1018496405208111, 4.0753451448447491,
        4.0488733560437451, 4.02244210341167, 3.99605921624238,
        3.9697325239076857, 3.9434698556210583, 3.9172790407124305,
        3.8911679085311461, 3.8651442883335676, 3.8392160094119254,
        3.81339090112321, 3.7876767926322561, 3.7620815137030807,
        3.7366128922931727, 3.7112787613800493, 3.6860869401765761,
        3.66104528590065, 3.6361615508044172, 3.6114437570766706,
        3.5868981284191821, 3.5625264224503468, 3.5383285984525146,
        3.5143048856023009, 3.4904554080077004, 3.4667803278239031,
        3.4432797935280384, 3.41995395855346, 3.39680297443919,
        3.3738269935005576, 3.351026167743056, 3.3284006492818805,
        3.3059505901920843, 3.2836761425631904, 3.2615774584777104,
        3.2396546900459113, 3.21790798931759, 3.1963375083992798,
        3.1749433994123648, 3.1537258143891682, 3.1326849054348158,
        3.1118208246502075, 3.0911337241121335, 3.0706237559345069,
        3.0502910721475276, 3.0301358248592791, 3.0101581662199988,
        2.9903582482399962, 2.9707362230281333, 2.9512922426787784,
        2.9320264593002912, 2.9129390249253198, 2.8940300916929598,
        2.8752998116734534, 2.8567483369427347, 2.8383758196184905,
        2.8201824117264587, 2.8021682654250029, 2.7843335327815573,
        2.7666783658499225, 2.7492029167552285, 2.7319073375562692,
        2.7147917803588446, 2.6978563972706895, 2.6811013403308843,
        2.6645267616705346, 2.6481328133293638, 2.6319196474390316,
        2.6158874160732979, 2.6000362713252563, 2.5843663652738518,
        2.5688778499660092, 2.5535708775856087, 2.5384456001563045,
        2.5235021697387308, 2.5087407384824489, 2.4941614584587679,
        2.4797644817172744, 2.4655499603860358, 2.45151804652899,
        2.4376688922853544, 2.4240026496449332, 2.4105194707650717,
        2.3972195077506995, 2.3841029126370974, 2.3711698375352941,
        2.3584204345021891, 2.3458548557019143, 2.3334732531608862,
        2.3212757789516876, 2.309262585208522, 2.2974338239973768,
        2.2857896474020771, 2.2743302075215706, 2.2630556564283824,
        2.2519661462436908, 2.2410618289944875, 2.2303428568174128,
        2.2198093818044575, 2.2094615560046242, 2.19929953153352,
        2.1893234604624316, 2.1795334948948883, 2.1699297869138556,
        2.160512488594839, 2.1512817520422631, 2.1422377293194219,
        2.1333805725475785, 2.1247104337941396, 2.1162274651211792,
        2.107931818691736, 2.0998236465035407, 2.0919031006924373,
        2.0841703333647756, 2.0766254965755206, 2.0692687423872957,
        2.0621002229319654, 2.0551200902955875, 2.0483284965589168,
        2.0417255938020138, 2.0353115340961843, 2.0290864695570274,
        2.0230505522592543, 2.0172039343207326, 2.0115467677636638,
        2.0060792047470462, 2.0008013972901497, 1.9957134975244715,
        1.9908156575435991, 1.98610802943473, 1.9815907652280322,
        1.9772640170624525, 1.9731279370407542, 1.9691826771947236,
        1.9654283896735021, 1.9618652265083076, 1.9584933398213553,
        1.9553128817023397, 1.9523240042160628, 1.9495268594356623,
        1.946921599529027, 1.9445083764859377, 1.942287342434553,
        1.940258649491523, 1.9384224497086906, 1.936778895168668,
        1.9353281379683946, 1.9340703302262725, 1.9330056239830324,
        1.9321341713321194, 1.9314561244032218, 1.9309716352465174,
        1.9306808559455066, 1.9305839386084176, 1.9306810353048958,
        1.9309722981536699, 1.9314578791974488, 1.9321379305444304,
        1.9330126042983933, 1.9340820525560285, 1.9353464273377026,
        1.9368058807649322, 1.9384605649787947, 1.9403106320114132,
        1.9423562339491385, 1.9445975228807941, 1.9470346509228573,
        1.9496677701491576, 1.9524970326286404, 1.9555225904680322,
        1.9587445957528864, 1.9621632005724372, 1.9657785569166657,
        1.9695908172931533, 1.9736001325986974, 1.9778066581210825,
        1.9822105370391487, 1.9868119462590972, 1.9916109692330193,
        1.9966079470976166, 2.0018025589333335, 2.0071988949012116,
        2.0128119991816056, 2.01866132712465, 2.0247656719614748,
        2.0311440846084436, 2.0378155224790682, 2.0447989768729169,
        2.0521134268328738, 2.0597778558456761, 2.0678112457164617,
        2.076232578918793, 2.08506083771537, 2.0943150043611496,
        2.1040140611840545, 2.1141769904374743, 2.124822774441582,
        2.1359703954555176, 2.1476388357969802, 2.1598470777329966,
        2.1726141035845759, 2.1859588955813836, 2.1999004360544285,
        2.2144577073068854, 2.229649691602261, 2.2454953712217116,
        2.2620137284559925, 2.2792237456290794, 2.2971444050065557,
        2.3157946888460357, 2.335193579513926, 2.3553600591397164,
        2.3763131104613491, 2.3980717144332435, 2.4206548571156956,
        2.4440815103533584, 2.4683706853618284, 2.4935412842902354,
        2.519612511255223, 2.5466027368823423, 2.57453247415158,
        2.603407962266286, 2.633199995345048, 2.6638650936631967,
        2.6953619198550376, 2.7276483031303682, 2.7606823746796372,
        2.7944221564808185, 2.8288257100109186, 2.8638510825010957,
        2.8994563263760553, 2.9355994920204997, 2.9722386307055113,
        3.0093317933655928, 3.0468370309895056, 3.0847123945770694,
        3.1229159351242219, 3.1614057037049879, 3.2001397512130052,
        3.2390761286807472, 3.2781728871802658, 3.3173880776288884,
        3.3566797510614332, 3.396005958483733, 3.435324750857204,
        3.4745941792619304, 3.5137722946226466, 3.5528171479745945,
        3.5916867903078429, 3.6303392726515789, 3.6687326459480896,
        3.7068249613098065, 3.7445742694093056, 3.7819386223003257,
        3.8188760679340721, 3.8553446657217409, 3.891302443517787,
        3.9267075163502949, 3.9615177579610124, 3.9956917082598449,
        4.0291861946054288, 4.0619694550335872, 4.0940380627800534,
        4.1254000019142349, 4.1560615437647348, 4.186029625952334,
        4.2153109447402137, 4.2439122836163969, 4.2718403945228687,
        4.2991020407801663, 4.3257039816480161, 4.3516529778473991,
        4.3769557895207569, 4.401619177063453, 4.425649900786258,
        4.44905472101547, 4.4718403980368713, 4.4940136922054394,
        4.5155813638544195, 4.5365501732674511, 4.556926880764161,
        4.576718246698845, 4.5959310313915474, 4.6145719951463731,
        4.6326478982667236, 4.6501655011153993, 4.6671315639659285,
        4.6835528471860943, 4.6994361110848226, 4.7147881159492533,
        4.7296156221603729, 4.7439253899618041, 4.7577241797505438,
        4.7710187518041858, 4.783815866469193, 4.796122284033749,
        4.8079447648177487, 4.8192900692038121, 4.830164957467515,
        4.84057618990469, 4.8505305268824976, 4.860034728704929,
        4.8690955556875659, 4.8777197681597029, 4.8859141264308814,
        4.8936853908710525, 4.9010403217006733, 4.90798567933387,
        4.9145282240684267, 4.9206747161977145, 4.92643191605836,
        4.9318065839906868, 4.9368054802975783, 4.9414353652983971,
        4.9457029993038031, 4.9496151426538368, 4.9531785556905552,
        4.9563999986734286, 4.959286231990772, 4.9618440159047568,
        4.9640801107804124, 4.9660012769136621, 4.9676142746363476,
        4.968925864278904, 4.9699428061489455, 4.9706718605523275,
        4.9711197878267637, 4.9712933483160509, 4.97119930231605,
        4.9708444101279508, 4.97023543211263, 4.969379128569, 4.9682822598410317,
        4.9669515861784648, 4.9653938681014864, 4.9636158653729971,
        4.9616243399139464, 4.9594260475133787, 4.9570277611147961,
        4.9544362061941438, 4.9516582303966281, 4.948699866335601,
        4.9455651227279374, 4.9422571929244414, 4.9387793927505852,
        4.935134990512819, 4.93132727160478, 4.9273595152722987,
        4.923235003044363, 4.9189570155338167, 4.9145288337640576,
        4.9099537385586869, 4.905235010836912, 4.9003759314747892,
        4.8953797813963247, 4.8902498414174929, 4.88498939247449,
        4.87960171550461, 4.8740900913043772, 4.8684578008265946,
        4.86270812494133, 4.8568443445524307, 4.8508697405109968,
        4.8447875937348677, 4.8386011851247375, 4.8323137955751747,
        4.8259287059083134, 4.8194491970863247, 4.8128785499984525,
        4.8062200455106048, 4.7994769644835467, 4.7926525878352679,
        4.7857501965017786, 4.7787730713039078, 4.771724493148553,
        4.76460774293378, 4.7574261015750885, 4.750182849910364,
        4.7428812688634814, 4.7355246393119055, 4.7281162421799028,
        4.72065935828153, 4.7131572685742409, 4.705613253923798,
        4.6980305952332353, 4.6904125733776221, 4.6827624692111964,
        4.67508356372638, 4.6673791377143061, 4.6596524720884744,
        4.651906847753704, 4.64414554562926, 4.636371846558359, 4.62858903141066,
        4.6208003811316969, 4.61300917662243, 4.6052186986939558,
        4.59743222826825, 4.5896530462870553, 4.5818844336133786,
        4.5741296710905477, 4.5663920396340307, 4.5586748201778411,
        4.5509812935677525, 4.5433147406798184, 4.5356784424278844,
        4.5280756797301489, 4.5205097334312443, 4.51298388442997,
        4.5055014135955416, 4.4980656019014, 4.4906797301513413,
        4.4833470792374763, 4.4760709301344725, 4.4688545635881427,
        4.4617012607951905, 4.4546143017819579, 4.447596969798334,
        4.4406525393904843, 4.4337843088243654, 4.42699551106419,
        4.4202895470185037, 4.4136686984691371, 4.4071324681585011,
        4.4006792397461361, 4.3943075649003251, 4.3880159298231138,
        4.3818028445111636, 4.375666810349335, 4.369606331814424,
        4.3636199122642312, 4.3577060554941847, 4.3518632651450719,
        4.3460900448509339, 4.3403848983044524, 4.3347463291243926,
        4.3291728410524692, 4.3236629377040368, 4.3182151227332648,
        4.3128278998228691, 4.3074997726532782, 4.3022292448654715,
        4.2970148201420049, 4.2918550021269395, 4.286748294510967,
        4.281693200928018, 4.2766882250356595, 4.2717318705806928,
        4.2668226411458541, 4.2619590403856042, 4.2571395720205931,
        4.2523627397077748, 4.24762704708168, 4.24293099782683,
        4.2382730955867771, 4.2336518440740756, 4.2290657468904227,
        4.2245133077447212, 4.2199930302854067, 4.2155034182027071,
        4.2110429751176319, 4.2066102047099019, 4.2022036106800682,
        4.197821696634378, 4.1934629662758978, 4.1891259232725542,
        4.1848090712732722, 4.1805109139287842, 4.1762299549425785,
        4.1719646979474128, 4.1677136466397107, 4.1634753046213211,
        4.1592481756206761, 4.1550307633010357, 4.1508215712860759,
        4.1466191032586375, 4.1424218628866312, 4.1382283538551183,
        4.1340370797856947, 4.1298465443591166, 4.1256552512821054,
        4.1214617041570465, 4.1172644066661457, 4.1130618625008193,
        4.1088525753026275, 4.1046350487686381, 4.1004077864736406,
        4.0961692921979935, 4.0919180695434711, 4.0876526221683775,
        4.0833714537884616, 4.0790730679749947, 4.0747559685145678,
        4.0704186589961671, 4.0660596430644516, 4.0616774244289129,
        4.0572705068064048, 4.0528373935911537, 4.0483765892374075,
        4.0438865951954206, 4.0393659212214095, 4.0348130558706545,
        4.03022662759064, 4.02560561170057, 4.0209491231047947,
        4.0162562558160317, 4.0115261120397063, 4.0067577910007071,
        4.0019503929279834, 3.9971030177383091, 3.9922147654682916,
        3.9872847361073238, 3.982312029636093, 3.9772957460660869,
        3.9722349853774968, 3.9671288476108391, 3.9619764326696107,
        3.9567768406456891, 3.9515291715297267, 3.9462325252415345,
        3.9408860018251444, 3.93548870126588, 3.9300397235999971,
        3.9245381687983052, 3.9189831368136012, 3.9133737276817127,
        3.9077090414260645, 3.9019881779921062, 3.896210237403543,
        3.8903743196373695, 3.8844795247220234, 3.8785249526356953,
        3.8725097033464686, 3.8664328768931417, 3.8602935732805319,
        3.8540908924213859, 3.8478239344091567, 3.8414917992243152,
        3.8350935868017992, 3.8286283971856983, 3.8220953303638141,
        3.8154934863410266, 3.808821965084455, 3.8020798666302378,
        3.795266290933347, 3.7883803380382024, 3.7814211078896247,
        3.77438770052232, 3.7672792159016542, 3.7600947540637426,
        3.7528334149622804, 3.745494298601459, 3.7380765050448113,
        3.7305791341614056, 3.7230012860359061, 3.7153420606895797,
        3.7076005580358573, 3.6997758781194547, 3.6918671209290412,
        3.6838733864571069, 3.6757937747329397, 3.6676273856722488,
        3.6593733193379232, 3.6510306757254147, 3.6425985548195796,
        3.63407605658316, 3.6254622810612553, 3.6167563282204771,
        3.6079572980802475, 3.5990642906050021, 3.5900764058273831,
        3.5809927437259716, 3.5718124042607919, 3.5625344875361127,
        3.5531580933252287, 3.5436823221598472, 3.534106272649518,
        3.5244290484404357, 3.5146497396307841, 3.5047674735294576,
        3.4947812745937594, 3.4846904509916374, 3.47449358182542,
        3.4641941037702293, 3.4538075160784643, 3.44335417562624,
        3.4328537102380823, 3.4223260312927422, 3.41179094754819,
        3.4012683047881445, 3.3907779353891394, 3.38033967659283,
        3.3699733639046281, 3.3596988334406137, 3.3495359211430533,
        3.3395044628773451, 3.3296242946532359, 3.3199152525010036,
        3.3103971722635146, 3.3010898899262449, 3.2920132414770134,
        3.28318706283913, 3.27463118999352, 3.2663654588809683,
        3.258409705429389, 3.2507837656611924, 3.2435074754680162,
        3.2366006708282651, 3.2300831876869989, 3.2239748620635722,
        3.218295529796475, 3.2130650269015315, 3.208303189402566,
        3.2040298531317708, 3.2002648541186254, 3.1970280283028223,
        3.1943392116344809, 3.1922182400818544, 3.1906849495561254,
        3.1897591760642676, 3.1894607555807504, 3.1898095239608164,
        3.1908253172552183, 3.1925279713880137, 3.1949373223264486,
        3.1980732059820531, 3.2019554583172676, 3.206603915374274,
        3.2120384129994424, 3.2182787872058047, 3.2253448739335187,
        3.2332565091273331, 3.2420335287734416, 3.2516957688035326,
        3.2622630651456719, 3.2737552538355388, 3.2861921707465775,
        3.2995936518587619, 3.3139795331642214, 3.3293696505900612,
        3.3457838400538793, 3.3632419375557, 3.3817637790692503,
        3.401369200504587, 3.4220780378399112, 3.4439101270058825,
        3.4668853040125893, 3.491023404751743, 3.5163442652047592,
        3.5428677213238102, 3.5706136091073137, 3.5996017644146012,
        3.6298520233482656, 3.6613842215403372, 3.694218195849305,
        3.7283737797142269, 3.7638708159775449, 3.8007291216220791,
        3.8389685851296758, 3.8786088969851038, 3.9196702957477889,
        3.9621715072226138, 4.0061351456185337, 4.0515579176806993,
        4.0983721959278148, 4.1464844453955276, 4.1958050195140322,
        4.246242758993402, 4.2977070526243324, 4.3501070911439781,
        4.4033521368129245, 4.457351426062341, 4.5120142047605363,
        4.5672497151961426, 4.6229672010436884, 4.6790759054768589,
        4.735485071854713, 4.7921039434222763, 4.8488417634877923,
        4.90560777538922, 4.9623112224069938, 5.0188613477927229,
        5.0751673949004248, 5.1311386070158642, 5.1866842274177616,
        5.2417134994141481, 5.2961356662986949, 5.3498599714157322,
        5.4027956579921641, 5.4548519693361586, 5.5059381487990464,
        5.5559634396187239, 5.6048370851803568, 5.652468328570686,
        5.6987664135913283, 5.7436405822833487, 5.7870000813338045,
        5.8287541445955178, 5.8688120417048522, 5.9070829428997982,
        5.9434762937437506, 5.9779007777228541, 6.0102671819554887,
        6.0404808861836718, 6.0684832970397187, 6.0943052835318507,
        6.1180137414026392, 6.1396701591593761, 6.1593381289713989,
        6.177080480741485, 6.1929603199657848, 6.2070406524225934,
        6.21938451999904, 6.2300549515301746, 6.2391149805375195,
        6.2466276388757809, 6.252655958958683, 6.2572629730784994,
        6.2605117135207422, 6.2624652125422395, 6.2631865024528119,
        6.2627386155568132, 6.2611845841219367, 6.2585874404396673,
        6.2550102168140871, 6.25051594555623, 6.2451676588847462,
        6.2390283891285572, 6.2321611686138843, 6.2246290296013926,
        6.2164950043315015, 6.2078221251480423, 6.19867342434218,
        6.1891119342555374, 6.1792006869314458, 6.16900271523029,
        6.1585810500508718, 6.1479987275307479, 6.137318769145911,
        6.1266042372236695, 6.115918080918167, 6.1053235625726412,
        6.0948830800381755, 6.0846612534880338, 6.07470789557959,
        6.0650360483454193, 6.0556439462104708, 6.0465320461349688,
        6.0376999403265978, 6.029147534345241, 6.0208746205463974,
        6.0128810321860557, 6.0051665876576887, 5.9977311108095517,
        5.9905744235170619, 5.9836963483709118, 5.9770967076142441,
        5.9707753237101322, 5.9647320190531135, 5.9589666159803354,
        5.9534789368800558, 5.9482688041619127, 5.94333604022166,
        5.93868046740026, 5.9343019081117445, 5.9302001847230033,
        5.92637511964822, 5.9228265352318319, 5.9195542538765071,
        5.9165580979696655, 5.9138378899095994, 5.9113934520461,
        5.9092246067723231, 5.9073311765066476, 5.9057129836030242,
        5.9043698504156135, 5.9033015993769817, 5.9025080529063505,
        5.9019890332890208, 5.901744362961221, 5.9017738643175992,
        5.9020773597120337, 5.9026546715710078, 5.90350562224577,
        5.9046300341079467, 5.9060277296044461, 5.9076985310389318,
        5.9096422608357635, 5.9118587413954744, 5.9143477950909311,
        5.9171092442662676, 5.9201429113519835, 5.9234486187356881,
        5.9270261887617934, 5.9308754438600007, 5.9349962063617836,
        5.9393882987094164, 5.9440515432494374, 5.9489857623677107,
        5.9541907784575256, 5.9596664139257216, 5.965412491098891,
        5.9714288323993578, 5.9777152602248949, 5.98427159692818,
        5.9910976649213286, 5.9981932865438221, 6.0055582842357946,
        6.0131924803557109, 6.0210956972799226, 6.0292677573858651,
        6.0377084830866252, 6.0464176967560777, 6.0553952207620565,
        6.06464087748661, 6.0741544894059381, 6.0839358786347129,
        6.09398486821539, 6.1043012787740327, 6.1148849376066075,
        6.1257356534718319, 6.13685328632882, 6.1482375551213417,
        6.1598885411101527, 6.1718039112088254, 6.1839753369506605,
        6.19639207548801, 6.2090437462749914, 6.2219198278861425,
        6.2350098499231814, 6.2483033235282424, 6.2617897664784365,
        6.2754586942597923, 6.2892996231108329, 6.3033020690248636,
        6.3174555480010568, 6.3317495761532445, 6.3461736695714555,
        6.3607173441876332, 6.3753701161134426, 6.3901215013932982,
        6.4049610161148411, 6.4198781762640786, 6.4348624979087674,
        6.449903497100145, 6.4649906899176894, 6.4801135923638924,
        6.4952617204959635, 6.5104245903946971, 6.5255917180951641,
        6.5407526196087167, 6.5558968110264839, 6.5710138084108651,
        6.5860931277583852, 6.6011242851538547, 6.6160967966043112,
        6.6310001782590033, 6.6458239460568231, 6.66055761606553,
        6.6751907044150212, 6.6897127270662535, 6.7041132000774528,
        6.7183816395502838, 6.7325075614900234, 6.7464804819727275,
        6.7602899169976975, 6.7739253826561709, 6.7873763950250448,
        6.8006324700775416, 6.8136831239194962, 6.8265178725630022,
        6.8391262321106838, 6.851497718535839, 6.8636218479510482,
        6.8754881363792775, 6.8870860998831569, 6.8984052544691616,
        6.909435116229603, 6.9201652012232007, 6.9305850254748407,
        6.9406841050036832, 6.9504519558807605, 6.95987809422765,
        6.9689520359927224, 6.97766329724101, 6.9860013940698131,
        6.9939558425107613, 7.0015161585665124, 7.008671858349703,
        7.0154124578544668, 7.0217274731940744, 7.0276064203530453,
        7.033038815402449, 7.0380141744066158, 7.0425220134543824,
        7.0465518483815721, 7.0500931957216038, 7.0531355703801539,
        7.0556684915734014, 7.057681466471, 7.0591640356772452,
        7.0601056473647841, 7.0604960048146133, 7.0603241553165681,
        7.0595835181439988, 7.0582783696841576, 7.0564173584118706,
        7.0540084765705364, 7.0510599718075557, 7.04757999917254,
        7.04357674711443, 7.0390583921269929, 7.0340331150041342,
        7.0285090948950346, 7.0224945116330142, 7.0159975447883234,
        7.0090263740121044, 7.0015891788829983, 6.9936941390993166,
        6.9853494342641387, 6.9765632440186662, 6.9673437479557352,
        6.9576991257519136, 6.9476375570339055, 6.9371672214135147,
        6.9262962985201071, 6.9150329679773082, 6.9033854094738292,
        6.8913618025937211, 6.87897032695172, 6.8662191622239668,
        6.853116487999058, 6.8396704839360876, 6.8258893296845278,
        6.8117812047926325, 6.7973542890399212, 6.7826167618431459,
        6.7675768031705639, 6.7522425918704956, 6.7366223098261937,
        6.7207241302539931, 6.7045562505490208, 6.6881268015303315,
        6.6714440848759766, 6.6545152639892713, 6.6373446755587917,
        6.61993551792855, 6.602291160251526, 6.5844149053300205,
        6.5663100799224114, 6.5479800022162689, 6.5294279933891515,
        6.5106573736018278, 6.4916714633762478, 6.4724735831137918,
        6.4530670532392573, 6.4334551941252878, 6.413641326305493,
        6.3936287701084131, 6.3734208459782611, 6.3530208743471785,
        6.3324321756608208, 6.31165807030756, 6.2907018787010527,
        6.2695669212908527, 6.2482565185266772, 6.226773990766894,
        6.2051226584726242, 6.1833058420767424, 6.1613268619741994,
        6.1391890385890111, 6.11689569237356, 6.0944501437374541,
        6.0718557130788051, 6.04911572085383, 6.0262334874476355,
        6.0032123333826073, 5.9800555788584084, 5.9567665448103684,
        5.9333485502457846, 5.9098049195902869, 5.8861389622660658,
        5.8623540287905476, 5.8384533643726044, 5.8144409157451165,
        5.7903223714589931, 5.7661041217185067, 5.7417924514521994,
        5.7173936863578785, 5.6929141375014645, 5.6683601212200792,
        5.6437379518885189, 5.6190539446529666, 5.5943144143508317,
        5.5695256759001346, 5.5446940443067341, 5.5198258343701028,
        5.4949273610786378, 5.470004939391325, 5.445064884182278,
        5.4201135103456171, 5.3951571328586869, 5.3702020666417463,
        5.3452546265745209, 5.3203211276228162, 5.2954078846714729,
        5.2705212126812748, 5.2456674265573886, 5.2208528412010988,
        5.1960837715602, 5.1713665325714189, 5.1467074390962546,
        5.1221128061164221, 5.097588948552608, 5.0731421813051449,
        5.0487788192634708, 5.0245051774091891, 5.0003275706612138,
        4.9762523138984118, 4.9522857220713847, 4.9284341101059708,
        4.9047037929212491, 4.8811010854234445, 4.8576323025457517,
        4.8343037592159259, 4.81112177037645, 4.7880926508675037,
        4.7652227157058329, 4.7425182797844831, 4.7199856580141644,
        4.6976311652977412, 4.6754611165916655, 4.6534818268322589,
        4.6316996108952591, 4.6101207837132279, 4.5887516602305389,
        4.5675985553794973, 4.5466677840313574, 4.5259656611571968,
        4.5054985016313287, 4.485272620449682, 4.4652943324791305,
        4.4455699526200476, 4.4261057958540944, 4.4069081771045493,
        4.3879834112151483, 4.3693378131680705, 4.350977697910305,
        4.33290938031793, 4.3151391753199695, 4.297673397820212,
        4.2805183628192571, 4.26368038515417, 4.2471657797660587,
        4.2309808616200293, 4.2151319455876806, 4.1996253466732654,
        4.1844673795108127, 4.1696643598768226, 4.1552226004133468,
        4.1411484222653341, 4.1274481230188345, 4.11412806577816,
        4.1011944325567669, 4.0886538707963087, 4.07650992742719,
        4.0647584500652671, 4.0533921858580388, 4.0424043473608275,
        4.0317879659675, 4.0215361387919755, 4.0116419390982312,
        4.002098448904845, 3.9928987469654471, 3.9840359132317897,
        3.9755030272363103, 3.9672931686517789, 3.9593994171541955,
        3.9518148523172569, 3.944532553826372, 3.9375456012947421,
        3.9308470744065076, 3.9244300527733742, 3.9182876160181559,
        3.9124128438108241, 3.906798815811825, 3.9014386116033584,
        3.8963253108837894, 3.891451993251672, 3.8868117383981846,
        3.8823976259001887, 3.8782027354133048, 3.8742201466569579,
        3.8704429391846515, 3.8668641926421592, 3.8634769867152912,
        3.8602744010449945, 3.8572495152147548, 3.8543954089212975,
        3.85170516177525, 3.8491718534605175, 3.8467885635546084,
        3.8445483717320155, 3.842444357651869, 3.8404696009582726,
        3.8386171812270975, 3.8368801781407367, 3.8352516713902007,
        3.8337247405299792, 3.8322924652549548, 3.8309479251985623,
        3.8296842000012563, 3.8284943692882272, 3.8273715127135013,
        3.8263087098998527, 3.8252990405352412, 3.8243355842043476,
        3.8234114205968748, 3.8225196293234092, 3.82165329004391,
        3.8208054823859108, 3.8199692859540635, 3.81913778049331,
        3.8183040455552275, 3.8174611607846067, 3.8166022058805322,
        3.8157202604374909, 3.81480840409758, 3.8138597165247066,
        3.8128672773341812, 3.8118241662113355, 3.8107234627262181,
        3.8095582465911906, 3.8083215974155582, 3.8070065948421625,
        3.8056063184996729, 3.8041138479774652, 3.8025222632270341,
        3.8008246430848693, 3.7990140691329821, 3.7970836158834964,
        3.7950263771054455, 3.7928353933411811, 3.7905038520089835,
        3.7880245630012883, 3.7853928520884872, 3.7826102924871705,
        3.7796809735174923, 3.7766086068237388, 3.7733970508787409,
        3.7700501109706148, 3.766571611703144, 3.7629653706500208,
        3.7592352078825555, 3.7553849426144827, 3.7514183943425663,
        3.7473393825167252, 3.7431517265316234, 3.7388592458272973,
        3.7344657598664814, 3.7299750880802294, 3.7253910498400598,
        3.720717464654554, 3.715958151926408, 3.7111169311038976,
        3.7061976215835517, 3.7012040427944712, 3.696140014277483,
        3.6910093553475316, 3.6858158854473495, 3.6805634240778291,
        3.675255790617193, 3.6698968045332925, 3.6644902852334775,
        3.6590400521421143, 3.6535499247701906, 3.6480237224341057,
        3.64246526461862, 3.6368783708076848, 3.6312668604076843,
        3.6256345527731688, 3.6199852674112467, 3.6143228238045455,
        3.6086510412786716, 3.6029737392980259, 3.5972947373483706,
        3.5916178548279603, 3.5859469111503666, 3.5802857257752296,
        3.5746381181315372, 3.5690079076763488, 3.5633989137713384,
        3.5578149559199366, 3.5522598535272651, 3.5467374260651807,
        3.5412514928974876, 3.5358058734624254, 3.5304043872974078,
        3.5250508537346454, 3.5197490922214882, 3.5145029222091635,
        3.5093161631416883, 3.5041926344439274, 3.4991361555370148,
        3.4941505458500774, 3.4892396248494837, 3.4844072119266221,
        3.4796571265456167, 3.4749931881406164, 3.4704192161457126,
        3.4659390299536419, 3.4615564490237305, 3.4572752928173669,
        3.4530993807526613, 3.4490325322250288, 3.445078566719372,
        3.441241303655056, 3.4375245624555917, 3.433932162459195,
        3.4304679235310518, 3.4271356638403661, 3.4239392063848224,
        3.4208823605702365, 3.4179689737371546, 3.4152027883472389,
        3.4125878164132022, 3.4101262740153788, 3.4078159175559866,
        3.4056527072946881, 3.403632873216702, 3.4017525404153828,
        3.4000078718690716, 3.3983950169703352, 3.39691012997674,
        3.3955493633833616, 3.3943088703749624, 3.3931848037432775,
        3.3921733166143904, 3.391270561895364, 3.3904726925757909,
        3.3897758615670672, 3.3891762218720851, 3.3886699265307882,
        3.3882531284016237, 3.3879219805257117, 3.3876726358645981,
        3.3875012473584367, 3.3874039680082086, 3.3873769507670004,
        3.3874163486047388, 3.3875183145305447, 3.3876790014432987,
        3.3878945623593086, 3.3881611502554096, 3.3884749180953331,
        3.3888320188294125, 3.3892286054208189, 3.38966083088996,
        3.3901248481813231, 3.3906168102288148, 3.3911328700496819,
        3.3916691806166961, 3.3922218948693734, 3.3927871657831128,
        3.3933611463334548, 3.3939399895250442, 3.3945198482564742,
        3.3950968755433446, 3.3956672243801136, 3.3962270477072058,
        3.3967724984602454, 3.3972997296449639, 3.3978048942716104,
        3.3982841452584522, 3.39873363555723, 3.3991495181913738,
        3.3995279461143273, 3.3998650722801016, 3.4001570496604017,
        3.400400031237234, 3.4005901700227734, 3.4007236188717029,
        3.4007965308522023, 3.4008050589271139, 3.4007453560406229,
        3.4006135751426361, 3.4004058692387473, 3.4001183913218624,
        3.399747294308193, 3.3992887311872066, 3.3987388549255626,
        3.3980938185244525, 3.3973497749169983, 3.3965028770789525,
        3.3955492779997853, 3.3944851306501733, 3.3933065879810829,
        3.3920098028536363, 3.3905909287994658, 3.3890461171371773,
        3.3873715253242911, 3.3855632939352431, 3.3836176102968585,
        3.3815305323672171, 3.3792984750247035, 3.3769169357564772,
        3.3743875245969468, 3.3717270302437692, 3.3689583540363577,
        3.3661034798734382, 3.3631847485522037, 3.3602243714277784,
        3.3572446068359612, 3.3542676960737814, 3.351315886530577,
        3.3484114233938325, 3.3455765527231085, 3.3428335202723964,
        3.3402045717767019, 3.3377119531048054, 3.3353779100268128,
        3.3332246883678809, 3.3312745339013135, 3.3295496924320855,
        3.3280724097639096, 3.3268649317064272, 3.3259495040647158,
        3.3253483726209185, 3.3250837831786617, 3.3251779815727818,
        3.3256532135415049, 3.3265317249143793, 3.3278357615396943,
        3.3295875691548562, 3.3318093935761568, 3.334523480589529,
        3.337752076100152, 3.3415174256374507, 3.3458417756650292,
        3.3507473701751893, 3.3562564600207385, 3.3623912770618989,
        3.3691741054502224, 3.3766270850863456, 3.3847727548090649,
        3.3936325522382131, 3.40323074575179, 3.4135727432275913,
        3.4246171175174323, 3.43630358095551, 3.4485746766764933,
        3.4613718464834071, 3.4746369312145959, 3.4883116275533754,
        3.50233768421047, 3.5166568311591453, 3.5312108050440894,
        3.5459413402115874, 3.5607901718051922, 3.5756990346091513,
        3.5906096636104889, 3.6054637937161131, 3.6202031598708833,
        3.6347694969662427, 3.6491045399133508, 3.663150023662701,
        3.6768476831204597, 3.690139253205079, 3.702966468843059,
        3.7152710649225842, 3.7269947764648212, 3.738079338245063,
        3.7484664852632537, 3.7580979524717582, 3.7669154747577576,
        3.7748607869973267, 3.7818756240926641, 3.7879017212673509,
        3.7928808126309117, 3.7967546351156796, 3.7994649181647611,
        3.8009534119193273, 3.8011618091925241, 3.8000319613221505,
        3.7975052811652472, 3.7935243916559016, 3.7880288055106148,
        3.7809787565331225, 3.7723859330653537, 3.76228274431054,
        3.7506984894883959, 3.7376636777658985, 3.7232083799723528,
        3.70736282521755, 3.690157185508474, 3.6716216534550496,
        3.6517864142665077, 3.630681655746391, 3.608337564805935,
        3.5847843287133658, 3.5600521345451739, 3.5341711694523417,
        3.5071716206268375, 3.4790836751680665, 3.449937520298473,
        3.4197633430335097, 3.388591330587392, 3.356451670198938,
        3.3233745488837427, 3.2893901538276178, 3.2545286721851423,
        3.2188202911248913, 3.1822951977599589, 3.1449835792600824,
        3.1069156227385011, 3.068121515398003, 3.0286314443232047,
        2.9884755966878775, 2.947684159672769, 2.9062873203953838,
        2.8643152659739086, 2.8217981835568957, 2.7787662603766621,
        2.7352496835105122, 2.691278640069255, 2.6468833172746509,
        2.6020939022428085, 2.5569405821282771, 2.5114535440741879,
        2.4656629751834322, 2.419599062697364, 2.3732919936699108,
        2.3267719553077963, 2.2800691347327739, 2.2332137191057106,
        2.1862358955592893, 2.1391658512046452, 2.09203377329596,
        2.044869848883764, 1.9977042651164396, 1.9505672092134854,
        1.9034888682686208, 1.8564994293969657, 1.8096290798212955,
        1.762908006650477, 1.7163663970415257, 1.6700344381069789,
        1.6239423170227922, 1.5781202209596179, 1.5325983370287861,
        1.4874068523514696, 1.4425759541365821, 1.3981358295117567,
        1.3541166656223016, 1.3105486495714211, 1.2674619685501787,
        1.224886809754981, 1.1828533601680586, 1.1413918072631715,
        1.1005323373283498, 1.0603051399763366, 1.02074039552897,
        0.98186830982940465, 0.94371901850376128, 0.906322851198588,
        0.86970960210255621, 0.83391044166801109, 0.79894737179588793,
        0.764819627036948, 0.7315172735666271, 0.69903175361984116,
        0.66735397412666064, 0.63647503584338283, 0.60638596959979607,
        0.57707783158569992, 0.54854166832716, 0.52076853120354849 } ;

      carPlant_smoothUp_followerSt_DW.FromWorkspace_PWORK.TimePtr = (void *)
        pTimeValues0;
      carPlant_smoothUp_followerSt_DW.FromWorkspace_PWORK.DataPtr = (void *)
        pDataValues0;
      carPlant_smoothUp_followerSt_DW.FromWorkspace_IWORK.PrevIndex = 0;
    }

    /* Start for Probe: '<S23>/Probe' */
    carPlant_smoothUp_followerSto_B.Probe_g[0] = 0.01;
    carPlant_smoothUp_followerSto_B.Probe_g[1] = 0.0;

    /* Start for Atomic SubSystem: '<Root>/Subscribe' */
    /* Start for MATLABSystem: '<S15>/SourceBlock' */
    carPlant_smoothUp_followerSt_DW.obj_j.matlabCodegenIsDeleted = true;
    carPlant_smoothUp_followerSt_DW.obj_j.isInitialized = 0;
    carPlant_smoothUp_followerSt_DW.obj_j.matlabCodegenIsDeleted = false;
    carPlant_smoothUp_followerSt_DW.objisempty = true;
    carPlant_smoothUp_followerSt_DW.obj_j.isSetupComplete = false;
    carPlant_smoothUp_followerSt_DW.obj_j.isInitialized = 1;
    for (i = 0; i < 7; i++) {
      tmp_3[i] = tmp[i];
    }

    tmp_3[7] = '\x00';
    Sub_carPlant_smoothUp_followerStopper_335.createSubscriber(tmp_3, 1);
    carPlant_smoothUp_followerSt_DW.obj_j.isSetupComplete = true;

    /* End of Start for MATLABSystem: '<S15>/SourceBlock' */
    /* End of Start for SubSystem: '<Root>/Subscribe' */

    /* Start for Scope: '<Root>/Acceleration' */
    {
      RTWLogSignalInfo rt_ScopeSignalInfo;
      static int_T rt_ScopeSignalWidths[] = { 1, 1 };

      static int_T rt_ScopeSignalNumDimensions[] = { 1, 1 };

      static int_T rt_ScopeSignalDimensions[] = { 1, 1 };

      static void *rt_ScopeCurrSigDims[] = { (NULL), (NULL) };

      static int_T rt_ScopeCurrSigDimsSize[] = { 4, 4 };

      static const char_T *rt_ScopeSignalLabels[] = { "filtered",
        "follower vehicle acceleration" };

      static char_T rt_ScopeSignalTitles[] =
        "filteredfollower vehicle acceleration";
      static int_T rt_ScopeSignalTitleLengths[] = { 8, 29 };

      static boolean_T rt_ScopeSignalIsVarDims[] = { 0, 0 };

      static int_T rt_ScopeSignalPlotStyles[] = { 0, 0 };

      BuiltInDTypeId dTypes[2] = { SS_DOUBLE, SS_DOUBLE };

      static char_T rt_ScopeBlockName[] =
        "carPlant_smoothUp_followerStopper/Acceleration";
      static int_T rt_ScopeFrameData[] = { 0, 0 };

      static RTWPreprocessingFcnPtr rt_ScopeSignalLoggingPreprocessingFcnPtrs[] =
      {
        (NULL), (NULL)
      };

      rt_ScopeSignalInfo.numSignals = 2;
      rt_ScopeSignalInfo.numCols = rt_ScopeSignalWidths;
      rt_ScopeSignalInfo.numDims = rt_ScopeSignalNumDimensions;
      rt_ScopeSignalInfo.dims = rt_ScopeSignalDimensions;
      rt_ScopeSignalInfo.isVarDims = rt_ScopeSignalIsVarDims;
      rt_ScopeSignalInfo.currSigDims = rt_ScopeCurrSigDims;
      rt_ScopeSignalInfo.currSigDimsSize = rt_ScopeCurrSigDimsSize;
      rt_ScopeSignalInfo.dataTypes = dTypes;
      rt_ScopeSignalInfo.complexSignals = (NULL);
      rt_ScopeSignalInfo.frameData = rt_ScopeFrameData;
      rt_ScopeSignalInfo.preprocessingPtrs =
        rt_ScopeSignalLoggingPreprocessingFcnPtrs;
      rt_ScopeSignalInfo.labels.cptr = rt_ScopeSignalLabels;
      rt_ScopeSignalInfo.titles = rt_ScopeSignalTitles;
      rt_ScopeSignalInfo.titleLengths = rt_ScopeSignalTitleLengths;
      rt_ScopeSignalInfo.plotStyles = rt_ScopeSignalPlotStyles;
      rt_ScopeSignalInfo.blockNames.cptr = (NULL);
      rt_ScopeSignalInfo.stateNames.cptr = (NULL);
      rt_ScopeSignalInfo.crossMdlRef = (NULL);
      rt_ScopeSignalInfo.dataTypeConvert = (NULL);
      carPlant_smoothUp_followerSt_DW.Acceleration_PWORK.LoggedData[0] =
        rt_CreateStructLogVar(
        carPlant_smoothUp_followerSt_M->rtwLogInfo,
        0.0,
        rtmGetTFinal(carPlant_smoothUp_followerSt_M),
        carPlant_smoothUp_followerSt_M->Timing.stepSize0,
        (&rtmGetErrorStatus(carPlant_smoothUp_followerSt_M)),
        "leaderConstant_followerStopper_acc",
        1,
        0,
        1,
        0.01,
        &rt_ScopeSignalInfo,
        rt_ScopeBlockName);
      if (carPlant_smoothUp_followerSt_DW.Acceleration_PWORK.LoggedData[0] ==
          (NULL))
        return;
    }

    /* Start for MATLABSystem: '<S13>/Get Parameter1' */
    carPlant_smoothUp_followerSt_DW.obj.matlabCodegenIsDeleted = true;
    carPlant_smoothUp_followerSt_DW.obj.isInitialized = 0;
    carPlant_smoothUp_followerSt_DW.obj.ticksUntilNextHit = 0.0;
    carPlant_smoothUp_followerSt_DW.obj.matlabCodegenIsDeleted = false;
    carPlant_smoothUp_followerSt_DW.objisempty_d = true;
    carPlant_smoothUp_followerSt_DW.obj.isSetupComplete = false;
    carPlant_smoothUp_followerSt_DW.obj.isInitialized = 1;
    for (i = 0; i < 7; i++) {
      tmp_3[i] = tmp_1[i];
    }

    tmp_3[7] = '\x00';
    ParamGet_carPlant_smoothUp_followerStopper_99.initialize(tmp_3);
    ParamGet_carPlant_smoothUp_followerStopper_99.initialize_error_codes(0, 1, 2,
      3);
    ParamGet_carPlant_smoothUp_followerStopper_99.set_initial_value(4.5);
    carPlant_smoothUp_followerSt_DW.obj.isSetupComplete = true;

    /* End of Start for MATLABSystem: '<S13>/Get Parameter1' */

    /* Start for MATLABSystem: '<S13>/Get Parameter2' */
    carPlant_smoothUp_followerSt_DW.obj_m.matlabCodegenIsDeleted = true;
    carPlant_smoothUp_followerSt_DW.obj_m.isInitialized = 0;
    carPlant_smoothUp_followerSt_DW.obj_m.ticksUntilNextHit = 0.0;
    carPlant_smoothUp_followerSt_DW.obj_m.matlabCodegenIsDeleted = false;
    carPlant_smoothUp_followerSt_DW.objisempty_f = true;
    carPlant_smoothUp_followerSt_DW.obj_m.isSetupComplete = false;
    carPlant_smoothUp_followerSt_DW.obj_m.isInitialized = 1;
    for (i = 0; i < 12; i++) {
      tmp_2[i] = tmp_0[i];
    }

    tmp_2[12] = '\x00';
    ParamGet_carPlant_smoothUp_followerStopper_100.initialize(tmp_2);
    ParamGet_carPlant_smoothUp_followerStopper_100.initialize_error_codes(0, 1,
      2, 3);
    ParamGet_carPlant_smoothUp_followerStopper_100.set_initial_value(6.0);
    carPlant_smoothUp_followerSt_DW.obj_m.isSetupComplete = true;

    /* End of Start for MATLABSystem: '<S13>/Get Parameter2' */

    /* Start for Atomic SubSystem: '<Root>/Publish' */
    /* Start for MATLABSystem: '<S14>/SinkBlock' */
    carPlant_smoothUp_followerSt_DW.obj_p.matlabCodegenIsDeleted = true;
    carPlant_smoothUp_followerSt_DW.obj_p.isInitialized = 0;
    carPlant_smoothUp_followerSt_DW.obj_p.matlabCodegenIsDeleted = false;
    carPlant_smoothUp_followerSt_DW.objisempty_c = true;
    carPlant_smoothUp_followerSt_DW.obj_p.isSetupComplete = false;
    carPlant_smoothUp_followerSt_DW.obj_p.isInitialized = 1;
    for (i = 0; i < 7; i++) {
      tmp_3[i] = tmp[i];
    }

    tmp_3[7] = '\x00';
    Pub_carPlant_smoothUp_followerStopper_334.createPublisher(tmp_3, 1);
    carPlant_smoothUp_followerSt_DW.obj_p.isSetupComplete = true;

    /* End of Start for MATLABSystem: '<S14>/SinkBlock' */
    /* End of Start for SubSystem: '<Root>/Publish' */

    /* Start for Scope: '<Root>/Velocity' */
    {
      RTWLogSignalInfo rt_ScopeSignalInfo;
      static int_T rt_ScopeSignalWidths[] = { 1, 1 };

      static int_T rt_ScopeSignalNumDimensions[] = { 1, 1 };

      static int_T rt_ScopeSignalDimensions[] = { 1, 1 };

      static void *rt_ScopeCurrSigDims[] = { (NULL), (NULL) };

      static int_T rt_ScopeCurrSigDimsSize[] = { 4, 4 };

      static const char_T *rt_ScopeSignalLabels[] = { "lead vehicle velocity",
        "follower vehicle velocity" };

      static char_T rt_ScopeSignalTitles[] =
        "lead vehicle velocityfollower vehicle velocity";
      static int_T rt_ScopeSignalTitleLengths[] = { 21, 25 };

      static boolean_T rt_ScopeSignalIsVarDims[] = { 0, 0 };

      static int_T rt_ScopeSignalPlotStyles[] = { 0, 0 };

      BuiltInDTypeId dTypes[2] = { SS_DOUBLE, SS_DOUBLE };

      static char_T rt_ScopeBlockName[] =
        "carPlant_smoothUp_followerStopper/Velocity";
      static int_T rt_ScopeFrameData[] = { 0, 0 };

      static RTWPreprocessingFcnPtr rt_ScopeSignalLoggingPreprocessingFcnPtrs[] =
      {
        (NULL), (NULL)
      };

      rt_ScopeSignalInfo.numSignals = 2;
      rt_ScopeSignalInfo.numCols = rt_ScopeSignalWidths;
      rt_ScopeSignalInfo.numDims = rt_ScopeSignalNumDimensions;
      rt_ScopeSignalInfo.dims = rt_ScopeSignalDimensions;
      rt_ScopeSignalInfo.isVarDims = rt_ScopeSignalIsVarDims;
      rt_ScopeSignalInfo.currSigDims = rt_ScopeCurrSigDims;
      rt_ScopeSignalInfo.currSigDimsSize = rt_ScopeCurrSigDimsSize;
      rt_ScopeSignalInfo.dataTypes = dTypes;
      rt_ScopeSignalInfo.complexSignals = (NULL);
      rt_ScopeSignalInfo.frameData = rt_ScopeFrameData;
      rt_ScopeSignalInfo.preprocessingPtrs =
        rt_ScopeSignalLoggingPreprocessingFcnPtrs;
      rt_ScopeSignalInfo.labels.cptr = rt_ScopeSignalLabels;
      rt_ScopeSignalInfo.titles = rt_ScopeSignalTitles;
      rt_ScopeSignalInfo.titleLengths = rt_ScopeSignalTitleLengths;
      rt_ScopeSignalInfo.plotStyles = rt_ScopeSignalPlotStyles;
      rt_ScopeSignalInfo.blockNames.cptr = (NULL);
      rt_ScopeSignalInfo.stateNames.cptr = (NULL);
      rt_ScopeSignalInfo.crossMdlRef = (NULL);
      rt_ScopeSignalInfo.dataTypeConvert = (NULL);
      carPlant_smoothUp_followerSt_DW.Velocity_PWORK.LoggedData[0] =
        rt_CreateStructLogVar(
        carPlant_smoothUp_followerSt_M->rtwLogInfo,
        0.0,
        rtmGetTFinal(carPlant_smoothUp_followerSt_M),
        carPlant_smoothUp_followerSt_M->Timing.stepSize0,
        (&rtmGetErrorStatus(carPlant_smoothUp_followerSt_M)),
        "leaderConstant_followerStopper_vel",
        1,
        0,
        1,
        0.01,
        &rt_ScopeSignalInfo,
        rt_ScopeBlockName);
      if (carPlant_smoothUp_followerSt_DW.Velocity_PWORK.LoggedData[0] == (NULL))
        return;
    }

    /* Start for Scope: '<Root>/X-Position' */
    {
      RTWLogSignalInfo rt_ScopeSignalInfo;
      static int_T rt_ScopeSignalWidths[] = { 1, 1 };

      static int_T rt_ScopeSignalNumDimensions[] = { 1, 1 };

      static int_T rt_ScopeSignalDimensions[] = { 1, 1 };

      static void *rt_ScopeCurrSigDims[] = { (NULL), (NULL) };

      static int_T rt_ScopeCurrSigDimsSize[] = { 4, 4 };

      static const char_T *rt_ScopeSignalLabels[] = { "<X>",
        "<X>" };

      static char_T rt_ScopeSignalTitles[] = "<X><X>";
      static int_T rt_ScopeSignalTitleLengths[] = { 3, 3 };

      static boolean_T rt_ScopeSignalIsVarDims[] = { 0, 0 };

      static int_T rt_ScopeSignalPlotStyles[] = { 0, 0 };

      BuiltInDTypeId dTypes[2] = { SS_DOUBLE, SS_DOUBLE };

      static char_T rt_ScopeBlockName[] =
        "carPlant_smoothUp_followerStopper/X-Position";
      static int_T rt_ScopeFrameData[] = { 0, 0 };

      static RTWPreprocessingFcnPtr rt_ScopeSignalLoggingPreprocessingFcnPtrs[] =
      {
        (NULL), (NULL)
      };

      rt_ScopeSignalInfo.numSignals = 2;
      rt_ScopeSignalInfo.numCols = rt_ScopeSignalWidths;
      rt_ScopeSignalInfo.numDims = rt_ScopeSignalNumDimensions;
      rt_ScopeSignalInfo.dims = rt_ScopeSignalDimensions;
      rt_ScopeSignalInfo.isVarDims = rt_ScopeSignalIsVarDims;
      rt_ScopeSignalInfo.currSigDims = rt_ScopeCurrSigDims;
      rt_ScopeSignalInfo.currSigDimsSize = rt_ScopeCurrSigDimsSize;
      rt_ScopeSignalInfo.dataTypes = dTypes;
      rt_ScopeSignalInfo.complexSignals = (NULL);
      rt_ScopeSignalInfo.frameData = rt_ScopeFrameData;
      rt_ScopeSignalInfo.preprocessingPtrs =
        rt_ScopeSignalLoggingPreprocessingFcnPtrs;
      rt_ScopeSignalInfo.labels.cptr = rt_ScopeSignalLabels;
      rt_ScopeSignalInfo.titles = rt_ScopeSignalTitles;
      rt_ScopeSignalInfo.titleLengths = rt_ScopeSignalTitleLengths;
      rt_ScopeSignalInfo.plotStyles = rt_ScopeSignalPlotStyles;
      rt_ScopeSignalInfo.blockNames.cptr = (NULL);
      rt_ScopeSignalInfo.stateNames.cptr = (NULL);
      rt_ScopeSignalInfo.crossMdlRef = (NULL);
      rt_ScopeSignalInfo.dataTypeConvert = (NULL);
      carPlant_smoothUp_followerSt_DW.XPosition_PWORK.LoggedData[0] =
        rt_CreateStructLogVar(
        carPlant_smoothUp_followerSt_M->rtwLogInfo,
        0.0,
        rtmGetTFinal(carPlant_smoothUp_followerSt_M),
        carPlant_smoothUp_followerSt_M->Timing.stepSize0,
        (&rtmGetErrorStatus(carPlant_smoothUp_followerSt_M)),
        "leaderConstant_followerStopper_xpos",
        1,
        0,
        1,
        0.01,
        &rt_ScopeSignalInfo,
        rt_ScopeBlockName);
      if (carPlant_smoothUp_followerSt_DW.XPosition_PWORK.LoggedData[0] == (NULL))
        return;
    }

    /* Start for Scope: '<Root>/Y-Position' */
    {
      RTWLogSignalInfo rt_ScopeSignalInfo;
      static int_T rt_ScopeSignalWidths[] = { 1, 1 };

      static int_T rt_ScopeSignalNumDimensions[] = { 1, 1 };

      static int_T rt_ScopeSignalDimensions[] = { 1, 1 };

      static void *rt_ScopeCurrSigDims[] = { (NULL), (NULL) };

      static int_T rt_ScopeCurrSigDimsSize[] = { 4, 4 };

      static const char_T *rt_ScopeSignalLabels[] = { "<Y>",
        "<Y>" };

      static char_T rt_ScopeSignalTitles[] = "<Y><Y>";
      static int_T rt_ScopeSignalTitleLengths[] = { 3, 3 };

      static boolean_T rt_ScopeSignalIsVarDims[] = { 0, 0 };

      static int_T rt_ScopeSignalPlotStyles[] = { 0, 0 };

      BuiltInDTypeId dTypes[2] = { SS_DOUBLE, SS_DOUBLE };

      static char_T rt_ScopeBlockName[] =
        "carPlant_smoothUp_followerStopper/Y-Position";
      static int_T rt_ScopeFrameData[] = { 0, 0 };

      static RTWPreprocessingFcnPtr rt_ScopeSignalLoggingPreprocessingFcnPtrs[] =
      {
        (NULL), (NULL)
      };

      rt_ScopeSignalInfo.numSignals = 2;
      rt_ScopeSignalInfo.numCols = rt_ScopeSignalWidths;
      rt_ScopeSignalInfo.numDims = rt_ScopeSignalNumDimensions;
      rt_ScopeSignalInfo.dims = rt_ScopeSignalDimensions;
      rt_ScopeSignalInfo.isVarDims = rt_ScopeSignalIsVarDims;
      rt_ScopeSignalInfo.currSigDims = rt_ScopeCurrSigDims;
      rt_ScopeSignalInfo.currSigDimsSize = rt_ScopeCurrSigDimsSize;
      rt_ScopeSignalInfo.dataTypes = dTypes;
      rt_ScopeSignalInfo.complexSignals = (NULL);
      rt_ScopeSignalInfo.frameData = rt_ScopeFrameData;
      rt_ScopeSignalInfo.preprocessingPtrs =
        rt_ScopeSignalLoggingPreprocessingFcnPtrs;
      rt_ScopeSignalInfo.labels.cptr = rt_ScopeSignalLabels;
      rt_ScopeSignalInfo.titles = rt_ScopeSignalTitles;
      rt_ScopeSignalInfo.titleLengths = rt_ScopeSignalTitleLengths;
      rt_ScopeSignalInfo.plotStyles = rt_ScopeSignalPlotStyles;
      rt_ScopeSignalInfo.blockNames.cptr = (NULL);
      rt_ScopeSignalInfo.stateNames.cptr = (NULL);
      rt_ScopeSignalInfo.crossMdlRef = (NULL);
      rt_ScopeSignalInfo.dataTypeConvert = (NULL);
      carPlant_smoothUp_followerSt_DW.YPosition_PWORK.LoggedData[0] =
        rt_CreateStructLogVar(
        carPlant_smoothUp_followerSt_M->rtwLogInfo,
        0.0,
        rtmGetTFinal(carPlant_smoothUp_followerSt_M),
        carPlant_smoothUp_followerSt_M->Timing.stepSize0,
        (&rtmGetErrorStatus(carPlant_smoothUp_followerSt_M)),
        "leaderConstant_followerStopper_ypos",
        1,
        0,
        1,
        0.01,
        &rt_ScopeSignalInfo,
        rt_ScopeBlockName);
      if (carPlant_smoothUp_followerSt_DW.YPosition_PWORK.LoggedData[0] == (NULL))
        return;
    }
  }

  /* InitializeConditions for Integrator: '<S378>/Integrator' incorporates:
   *  Integrator: '<S182>/Integrator'
   */
  if (rtmIsFirstInitCond(carPlant_smoothUp_followerSt_M)) {
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE[0] = 0.0;
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE[1] = 0.0;
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE[2] = 0.0;
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE[3] = 0.0;
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE_j[0] = 0.0;
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE_j[1] = 0.0;
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE_j[2] = 0.0;
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE_j[3] = 0.0;
  }

  carPlant_smoothUp_followerSt_DW.Integrator_IWORK = 1;

  /* End of InitializeConditions for Integrator: '<S378>/Integrator' */

  /* InitializeConditions for DiscreteIntegrator: '<S22>/Integrator' */
  carPlant_smoothUp_followerSt_DW.Integrator_IC_LOADING = 1U;
  carPlant_smoothUp_followerSt_DW.Integrator_PrevResetState = 0;

  /* InitializeConditions for Integrator: '<S182>/Integrator' */
  carPlant_smoothUp_followerSt_DW.Integrator_IWORK_h = 1;

  /* InitializeConditions for DiscreteIntegrator: '<S26>/Integrator' */
  carPlant_smoothUp_followerSt_DW.Integrator_IC_LOADING_g = 1U;
  carPlant_smoothUp_followerSt_DW.Integrator_PrevResetState_c = 0;

  /* InitializeConditions for Integrator: '<S290>/Integrator' incorporates:
   *  Integrator: '<S94>/Integrator'
   */
  if (rtmIsFirstInitCond(carPlant_smoothUp_followerSt_M)) {
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE_o[0] = 0.0;
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE_o[1] = 0.0;
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE_a[0] = 0.0;
    carPlant_smoothUp_followerSto_X.Integrator_CSTATE_a[1] = 0.0;
  }

  carPlant_smoothUp_followerSt_DW.Integrator_IWORK_p = 1;

  /* End of InitializeConditions for Integrator: '<S290>/Integrator' */

  /* InitializeConditions for Integrator: '<S94>/Integrator' */
  carPlant_smoothUp_followerSt_DW.Integrator_IWORK_l = 1;

  /* InitializeConditions for Backlash: '<S205>/Backlash' */
  carPlant_smoothUp_followerSt_DW.PrevY =
    carPlant_smoothUp_followerSto_P.Backlash_InitialOutput;

  /* InitializeConditions for Backlash: '<S186>/Backlash' */
  carPlant_smoothUp_followerSt_DW.PrevY_e =
    carPlant_smoothUp_followerSto_P.Backlash_InitialOutput_h;

  /* InitializeConditions for Integrator: '<S181>/lateral' */
  carPlant_smoothUp_followerSto_X.lateral_CSTATE[0] =
    carPlant_smoothUp_followerSto_P.lateral_IC;

  /* InitializeConditions for Integrator: '<S377>/lateral' */
  carPlant_smoothUp_followerSto_X.lateral_CSTATE_f[0] =
    carPlant_smoothUp_followerSto_P.lateral_IC_g;

  /* InitializeConditions for Integrator: '<S181>/lateral' */
  carPlant_smoothUp_followerSto_X.lateral_CSTATE[1] =
    carPlant_smoothUp_followerSto_P.lateral_IC;

  /* InitializeConditions for Integrator: '<S377>/lateral' */
  carPlant_smoothUp_followerSto_X.lateral_CSTATE_f[1] =
    carPlant_smoothUp_followerSto_P.lateral_IC_g;

  /* InitializeConditions for Integrator: '<S181>/lateral' */
  carPlant_smoothUp_followerSto_X.lateral_CSTATE[2] =
    carPlant_smoothUp_followerSto_P.lateral_IC;

  /* InitializeConditions for Integrator: '<S377>/lateral' */
  carPlant_smoothUp_followerSto_X.lateral_CSTATE_f[2] =
    carPlant_smoothUp_followerSto_P.lateral_IC_g;

  /* InitializeConditions for Integrator: '<S181>/lateral' */
  carPlant_smoothUp_followerSto_X.lateral_CSTATE[3] =
    carPlant_smoothUp_followerSto_P.lateral_IC;

  /* InitializeConditions for Integrator: '<S377>/lateral' */
  carPlant_smoothUp_followerSto_X.lateral_CSTATE_f[3] =
    carPlant_smoothUp_followerSto_P.lateral_IC_g;

  /* SystemInitialize for Atomic SubSystem: '<Root>/Subscribe' */
  /* SystemInitialize for Enabled SubSystem: '<S15>/Enabled Subsystem' */
  /* SystemInitialize for Outport: '<S382>/Out1' */
  carPlant_smoothUp_followerSto_B.In1 = carPlant_smoothUp_followerSto_P.Out1_Y0;

  /* End of SystemInitialize for SubSystem: '<S15>/Enabled Subsystem' */
  /* End of SystemInitialize for SubSystem: '<Root>/Subscribe' */

  /* SystemInitialize for MATLAB Function: '<Root>/MATLAB Function' */
  carPlant_smoothUp_followerSt_DW.y_not_empty = false;

  /* SystemInitialize for Atomic SubSystem: '<Root>/Dead Man's Switch' */
  /* SystemInitialize for MATLAB Function: '<S4>/timeout set to 0 output' */
  carPlant_smoothUp_followerSt_DW.sinceLastMsg_not_empty = false;

  /* End of SystemInitialize for SubSystem: '<Root>/Dead Man's Switch' */

  /* set "at time zero" to false */
  if (rtmIsFirstInitCond(carPlant_smoothUp_followerSt_M)) {
    rtmSetFirstInitCond(carPlant_smoothUp_followerSt_M, 0);
  }
}

/* Model terminate function */
void carPlant_smoothUp_followerStopper_terminate(void)
{
  /* Terminate for Atomic SubSystem: '<Root>/Subscribe' */
  /* Terminate for MATLABSystem: '<S15>/SourceBlock' */
  matlabCodegenHandle_matlabC_em4(&carPlant_smoothUp_followerSt_DW.obj_j);

  /* End of Terminate for SubSystem: '<Root>/Subscribe' */

  /* Terminate for MATLABSystem: '<S13>/Get Parameter1' */
  matlabCodegenHandle_matlabCodeg(&carPlant_smoothUp_followerSt_DW.obj);

  /* Terminate for MATLABSystem: '<S13>/Get Parameter2' */
  matlabCodegenHandle_matlabCodeg(&carPlant_smoothUp_followerSt_DW.obj_m);

  /* Terminate for Atomic SubSystem: '<Root>/Publish' */
  /* Terminate for MATLABSystem: '<S14>/SinkBlock' */
  matlabCodegenHandle_matlabCo_em(&carPlant_smoothUp_followerSt_DW.obj_p);

  /* End of Terminate for SubSystem: '<Root>/Publish' */
}
