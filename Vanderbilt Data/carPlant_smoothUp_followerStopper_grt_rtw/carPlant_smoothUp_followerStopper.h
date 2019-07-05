/*
 * carPlant_smoothUp_followerStopper.h
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

#ifndef RTW_HEADER_carPlant_smoothUp_followerStopper_h_
#define RTW_HEADER_carPlant_smoothUp_followerStopper_h_
#include <math.h>
#include <stddef.h>
#include <float.h>
#include <string.h>
#ifndef carPlant_smoothUp_followerStopper_COMMON_INCLUDES_
# define carPlant_smoothUp_followerStopper_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "rt_logging.h"
#include "slros_initialize.h"
#endif                                 /* carPlant_smoothUp_followerStopper_COMMON_INCLUDES_ */

#include "carPlant_smoothUp_followerStopper_types.h"

/* Shared type includes */
#include "multiword_types.h"
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtGetInf.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetContStateDisabled
# define rtmGetContStateDisabled(rtm)  ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
# define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
# define rtmGetContStates(rtm)         ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
# define rtmSetContStates(rtm, val)    ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
# define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
# define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
# define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
# define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetFinalTime
# define rtmGetFinalTime(rtm)          ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetIntgData
# define rtmGetIntgData(rtm)           ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
# define rtmSetIntgData(rtm, val)      ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeF
# define rtmGetOdeF(rtm)               ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
# define rtmSetOdeF(rtm, val)          ((rtm)->odeF = (val))
#endif

#ifndef rtmGetOdeY
# define rtmGetOdeY(rtm)               ((rtm)->odeY)
#endif

#ifndef rtmSetOdeY
# define rtmSetOdeY(rtm, val)          ((rtm)->odeY = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
# define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
# define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
# define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
# define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetRTWLogInfo
# define rtmGetRTWLogInfo(rtm)         ((rtm)->rtwLogInfo)
#endif

#ifndef rtmGetZCCacheNeedsReset
# define rtmGetZCCacheNeedsReset(rtm)  ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
# define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
# define rtmGetdX(rtm)                 ((rtm)->derivs)
#endif

#ifndef rtmSetdX
# define rtmSetdX(rtm, val)            ((rtm)->derivs = (val))
#endif

#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
# define rtmGetStopRequested(rtm)      ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
# define rtmSetStopRequested(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
# define rtmGetStopRequestedPtr(rtm)   (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
# define rtmGetT(rtm)                  (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTFinal
# define rtmGetTFinal(rtm)             ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetTPtr
# define rtmGetTPtr(rtm)               ((rtm)->Timing.t)
#endif

#define carPlant_smoothUp_followerStopper_M (carPlant_smoothUp_followerSt_M)

/* Block signals for system '<S94>/COMB2I' */
typedef struct {
  real_T y[2];                         /* '<S94>/COMB2I' */
} B_COMB2I_carPlant_smoothUp_fo_T;

/* Block signals (default storage) */
typedef struct {
  SL_Bus_carPlant_smoothUp_followe_Twist_domkw2 BusAssignment7;/* '<Root>/Bus Assignment7' */
  SL_Bus_carPlant_smoothUp_followe_Twist_domkw2 In1;/* '<S382>/In1' */
  real_T Probe[2];                     /* '<S19>/Probe' */
  real_T AvoidDividebyZero;            /* '<S19>/Avoid Divide by Zero' */
  real_T VectorConcatenate[4];         /* '<S378>/Vector Concatenate' */
  real_T VectorConcatenate1[4];        /* '<S378>/Vector Concatenate1' */
  real_T Saturation;                   /* '<S22>/Saturation' */
  real_T uT;                           /* '<S5>/1//T' */
  real_T AB;                           /* '<S5>/[A,B]' */
  real_T Probe_g[2];                   /* '<S23>/Probe' */
  real_T AvoidDividebyZero_l;          /* '<S23>/Avoid Divide by Zero' */
  real_T UnitConversion6;              /* '<S182>/Unit Conversion6' */
  real_T VectorConcatenate_g[4];       /* '<S182>/Vector Concatenate' */
  real_T VectorConcatenate1_i[4];      /* '<S182>/Vector Concatenate1' */
  real_T Saturation_l;                 /* '<S26>/Saturation' */
  real_T uT_k;                         /* '<S6>/1//T' */
  real_T AB_f;                         /* '<S6>/[A,B]' */
  real_T VectorConcatenate3[2];        /* '<S290>/Vector Concatenate3' */
  real_T VectorConcatenate_l[3];       /* '<S293>/Vector Concatenate' */
  real_T Add[3];                       /* '<S312>/Add' */
  real_T VectorConcatenate3_a[2];      /* '<S94>/Vector Concatenate3' */
  real_T VectorConcatenate_n[3];       /* '<S97>/Vector Concatenate' */
  real_T Add_b[3];                     /* '<S116>/Add' */
  real_T relativeposition;             /* '<Root>/Subtract2' */
  real_T Subtract;                     /* '<Root>/Subtract' */
  real_T VectorConcatenate1_b[3];      /* '<S185>/Vector Concatenate1' */
  real_T UnaryMinus[3];                /* '<S45>/Unary Minus' */
  real_T VectorConcatenate_a[6];       /* '<S45>/Vector Concatenate' */
  real_T VectorConcatenate2[2];        /* '<S41>/Vector Concatenate2' */
  real_T VectorConcatenate1_bi[2];     /* '<S41>/Vector Concatenate1' */
  real_T Backlash;                     /* '<S205>/Backlash' */
  real_T VectorConcatenate1_a[2];      /* '<Root>/Vector Concatenate1' */
  real_T VectorConcatenate4[2];        /* '<S174>/Vector Concatenate4' */
  real_T VectorConcatenate2_j[2];      /* '<S160>/Vector Concatenate2' */
  real_T VectorConcatenate3_b[2];      /* '<S160>/Vector Concatenate3' */
  real_T VectorConcatenate2_o[2];      /* '<S170>/Vector Concatenate2' */
  real_T VectorConcatenate3_b5[2];     /* '<S170>/Vector Concatenate3' */
  real_T VectorConcatenate_o[4];       /* '<S178>/Vector Concatenate' */
  real_T Product1[4];                  /* '<S181>/Product1' */
  real_T Backlash_a;                   /* '<S186>/Backlash' */
  real_T VectorConcatenate1_iz[3];     /* '<S381>/Vector Concatenate1' */
  real_T UnaryMinus_h[3];              /* '<S241>/Unary Minus' */
  real_T VectorConcatenate_c[6];       /* '<S241>/Vector Concatenate' */
  real_T VectorConcatenate2_i[2];      /* '<S237>/Vector Concatenate2' */
  real_T VectorConcatenate1_m[2];      /* '<S237>/Vector Concatenate1' */
  real_T VectorConcatenate_i[2];       /* '<Root>/Vector Concatenate' */
  real_T VectorConcatenate4_c[2];      /* '<S370>/Vector Concatenate4' */
  real_T VectorConcatenate4_b[4];      /* '<S350>/Vector Concatenate4' */
  real_T VectorConcatenate2_h[2];      /* '<S356>/Vector Concatenate2' */
  real_T VectorConcatenate3_m[2];      /* '<S356>/Vector Concatenate3' */
  real_T VectorConcatenate2_c[2];      /* '<S366>/Vector Concatenate2' */
  real_T VectorConcatenate3_f[2];      /* '<S366>/Vector Concatenate3' */
  real_T VectorConcatenate_e[4];       /* '<S374>/Vector Concatenate' */
  real_T Product1_l[4];                /* '<S377>/Product1' */
  real_T stateDer[4];                  /* '<S11>/vehicle model' */
  real_T stateDer_p[4];                /* '<S8>/vehicle model' */
  real_T safeValue;                    /* '<S4>/timeout set to 0 output' */
  B_COMB2I_carPlant_smoothUp_fo_T sf_COMB2I_f;/* '<S290>/COMB2I' */
  B_COMB2I_carPlant_smoothUp_fo_T sf_COMB2I;/* '<S94>/COMB2I' */
} B_carPlant_smoothUp_followerS_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  robotics_slros_internal_block_T obj; /* '<S13>/Get Parameter1' */
  robotics_slros_internal_block_T obj_m;/* '<S13>/Get Parameter2' */
  robotics_slros_internal_blo_e_T obj_p;/* '<S14>/SinkBlock' */
  robotics_slros_internal_bl_em_T obj_j;/* '<S15>/SourceBlock' */
  real_T Integrator_DSTATE;            /* '<S22>/Integrator' */
  real_T Integrator_DSTATE_m;          /* '<S26>/Integrator' */
  real_T PrevY;                        /* '<S205>/Backlash' */
  real_T PrevY_e;                      /* '<S186>/Backlash' */
  real_T y;                            /* '<Root>/MATLAB Function' */
  real_T sinceLastMsg;                 /* '<S4>/timeout set to 0 output' */
  struct {
    void *TimePtr;
    void *DataPtr;
    void *RSimInfoPtr;
  } FromWorkspace_PWORK;               /* '<Root>/From Workspace' */

  struct {
    void *LoggedData[2];
  } Acceleration_PWORK;                /* '<Root>/Acceleration' */

  struct {
    void *LoggedData[2];
  } Velocity_PWORK;                    /* '<Root>/Velocity' */

  struct {
    void *LoggedData[2];
  } XPosition_PWORK;                   /* '<Root>/X-Position' */

  struct {
    void *LoggedData[2];
  } YPosition_PWORK;                   /* '<Root>/Y-Position' */

  struct {
    void *LoggedData[2];
  } Scope1_PWORK;                      /* '<S7>/Scope1' */

  struct {
    int_T PrevIndex;
  } FromWorkspace_IWORK;               /* '<Root>/From Workspace' */

  int_T Integrator_IWORK;              /* '<S378>/Integrator' */
  int_T Integrator_IWORK_h;            /* '<S182>/Integrator' */
  int_T Integrator_IWORK_p;            /* '<S290>/Integrator' */
  int_T Integrator_IWORK_l;            /* '<S94>/Integrator' */
  int8_T Integrator_PrevResetState;    /* '<S22>/Integrator' */
  int8_T Integrator_PrevResetState_c;  /* '<S26>/Integrator' */
  uint8_T Integrator_IC_LOADING;       /* '<S22>/Integrator' */
  uint8_T Integrator_IC_LOADING_g;     /* '<S26>/Integrator' */
  boolean_T objisempty;                /* '<S15>/SourceBlock' */
  boolean_T objisempty_c;              /* '<S14>/SinkBlock' */
  boolean_T objisempty_d;              /* '<S13>/Get Parameter1' */
  boolean_T objisempty_f;              /* '<S13>/Get Parameter2' */
  boolean_T y_not_empty;               /* '<Root>/MATLAB Function' */
  boolean_T sinceLastMsg_not_empty;    /* '<S4>/timeout set to 0 output' */
} DW_carPlant_smoothUp_follower_T;

/* Continuous states (default storage) */
typedef struct {
  real_T Integrator_CSTATE[4];         /* '<S378>/Integrator' */
  real_T Integrator_CSTATE_j[4];       /* '<S182>/Integrator' */
  real_T Integrator_CSTATE_o[2];       /* '<S290>/Integrator' */
  real_T Integrator_CSTATE_a[2];       /* '<S94>/Integrator' */
  real_T lateral_CSTATE[4];            /* '<S181>/lateral' */
  real_T lateral_CSTATE_f[4];          /* '<S377>/lateral' */
} X_carPlant_smoothUp_followerS_T;

/* State derivatives (default storage) */
typedef struct {
  real_T Integrator_CSTATE[4];         /* '<S378>/Integrator' */
  real_T Integrator_CSTATE_j[4];       /* '<S182>/Integrator' */
  real_T Integrator_CSTATE_o[2];       /* '<S290>/Integrator' */
  real_T Integrator_CSTATE_a[2];       /* '<S94>/Integrator' */
  real_T lateral_CSTATE[4];            /* '<S181>/lateral' */
  real_T lateral_CSTATE_f[4];          /* '<S377>/lateral' */
} XDot_carPlant_smoothUp_follow_T;

/* State disabled  */
typedef struct {
  boolean_T Integrator_CSTATE[4];      /* '<S378>/Integrator' */
  boolean_T Integrator_CSTATE_j[4];    /* '<S182>/Integrator' */
  boolean_T Integrator_CSTATE_o[2];    /* '<S290>/Integrator' */
  boolean_T Integrator_CSTATE_a[2];    /* '<S94>/Integrator' */
  boolean_T lateral_CSTATE[4];         /* '<S181>/lateral' */
  boolean_T lateral_CSTATE_f[4];       /* '<S377>/lateral' */
} XDis_carPlant_smoothUp_follow_T;

#ifndef ODE3_INTG
#define ODE3_INTG

/* ODE3 Integration Data */
typedef struct {
  real_T *y;                           /* output */
  real_T *f[3];                        /* derivatives */
} ODE3_IntgData;

#endif

/* Parameters (default storage) */
struct P_carPlant_smoothUp_followerS_T_ {
  real_T FilteredDerivativeDiscreteorCon;/* Mask Parameter: FilteredDerivativeDiscreteorCon
                                          * Referenced by: '<S5>/[A,B]'
                                          */
  real_T FilteredDerivativeDiscreteorC_f;/* Mask Parameter: FilteredDerivativeDiscreteorC_f
                                          * Referenced by: '<S6>/[A,B]'
                                          */
  real_T FollowerVehicle_Af;           /* Mask Parameter: FollowerVehicle_Af
                                        * Referenced by: '<S45>/.5.*A.*Pabs.//R.//T'
                                        */
  real_T LeadVehicle_Af;               /* Mask Parameter: LeadVehicle_Af
                                        * Referenced by: '<S241>/.5.*A.*Pabs.//R.//T'
                                        */
  real_T FilteredDerivativeDiscreteorC_p;/* Mask Parameter: FilteredDerivativeDiscreteorC_p
                                          * Referenced by: '<S5>/[A,B]'
                                          */
  real_T FilteredDerivativeDiscreteorC_m;/* Mask Parameter: FilteredDerivativeDiscreteorC_m
                                          * Referenced by: '<S6>/[A,B]'
                                          */
  real_T FollowerVehicle_Cd;           /* Mask Parameter: FollowerVehicle_Cd
                                        * Referenced by: '<S45>/Constant'
                                        */
  real_T LeadVehicle_Cd;               /* Mask Parameter: LeadVehicle_Cd
                                        * Referenced by: '<S241>/Constant'
                                        */
  real_T FollowerVehicle_Cl;           /* Mask Parameter: FollowerVehicle_Cl
                                        * Referenced by: '<S45>/Constant1'
                                        */
  real_T LeadVehicle_Cl;               /* Mask Parameter: LeadVehicle_Cl
                                        * Referenced by: '<S241>/Constant1'
                                        */
  real_T FollowerVehicle_Cpm;          /* Mask Parameter: FollowerVehicle_Cpm
                                        * Referenced by: '<S45>/Constant2'
                                        */
  real_T LeadVehicle_Cpm;              /* Mask Parameter: LeadVehicle_Cpm
                                        * Referenced by: '<S241>/Constant2'
                                        */
  real_T FollowerVehicle_Cs[31];       /* Mask Parameter: FollowerVehicle_Cs
                                        * Referenced by: '<S45>/Cs'
                                        */
  real_T LeadVehicle_Cs[31];           /* Mask Parameter: LeadVehicle_Cs
                                        * Referenced by: '<S241>/Cs'
                                        */
  real_T FollowerVehicle_Cy_f;         /* Mask Parameter: FollowerVehicle_Cy_f
                                        * Referenced by: '<S41>/Cyf'
                                        */
  real_T LeadVehicle_Cy_f;             /* Mask Parameter: LeadVehicle_Cy_f
                                        * Referenced by: '<S237>/Cyf'
                                        */
  real_T FollowerVehicle_Cy_r;         /* Mask Parameter: FollowerVehicle_Cy_r
                                        * Referenced by: '<S41>/Cyr'
                                        */
  real_T LeadVehicle_Cy_r;             /* Mask Parameter: LeadVehicle_Cy_r
                                        * Referenced by: '<S237>/Cyr'
                                        */
  real_T FollowerVehicle_Cym[31];      /* Mask Parameter: FollowerVehicle_Cym
                                        * Referenced by: '<S45>/Cym'
                                        */
  real_T LeadVehicle_Cym[31];          /* Mask Parameter: LeadVehicle_Cym
                                        * Referenced by: '<S241>/Cym'
                                        */
  real_T KinematicSteering1_Db;        /* Mask Parameter: KinematicSteering1_Db
                                        * Referenced by: '<S205>/Backlash'
                                        */
  real_T KinematicSteering_Db;         /* Mask Parameter: KinematicSteering_Db
                                        * Referenced by: '<S186>/Backlash'
                                        */
  real_T FollowerVehicle_Fznom;        /* Mask Parameter: FollowerVehicle_Fznom
                                        * Referenced by: '<S8>/vehicle model'
                                        */
  real_T LeadVehicle_Fznom;            /* Mask Parameter: LeadVehicle_Fznom
                                        * Referenced by: '<S11>/vehicle model'
                                        */
  real_T FollowerVehicle_Izz;          /* Mask Parameter: FollowerVehicle_Izz
                                        * Referenced by: '<S8>/vehicle model'
                                        */
  real_T LeadVehicle_Izz;              /* Mask Parameter: LeadVehicle_Izz
                                        * Referenced by: '<S11>/vehicle model'
                                        */
  real_T FilteredDerivativeDiscreteor_me;/* Mask Parameter: FilteredDerivativeDiscreteor_me
                                          * Referenced by: '<S5>/Gain'
                                          */
  real_T FilteredDerivativeDiscreteorC_h;/* Mask Parameter: FilteredDerivativeDiscreteorC_h
                                          * Referenced by: '<S6>/Gain'
                                          */
  real_T FollowerVehicle_NF;           /* Mask Parameter: FollowerVehicle_NF
                                        * Referenced by: '<S8>/vehicle model'
                                        */
  real_T LeadVehicle_NF;               /* Mask Parameter: LeadVehicle_NF
                                        * Referenced by: '<S11>/vehicle model'
                                        */
  real_T FollowerVehicle_NR;           /* Mask Parameter: FollowerVehicle_NR
                                        * Referenced by: '<S8>/vehicle model'
                                        */
  real_T LeadVehicle_NR;               /* Mask Parameter: LeadVehicle_NR
                                        * Referenced by: '<S11>/vehicle model'
                                        */
  real_T FollowerVehicle_Pabs;         /* Mask Parameter: FollowerVehicle_Pabs
                                        * Referenced by: '<S45>/.5.*A.*Pabs.//R.//T'
                                        */
  real_T LeadVehicle_Pabs;             /* Mask Parameter: LeadVehicle_Pabs
                                        * Referenced by: '<S241>/.5.*A.*Pabs.//R.//T'
                                        */
  real_T DragForce_R;                  /* Mask Parameter: DragForce_R
                                        * Referenced by: '<S45>/.5.*A.*Pabs.//R.//T'
                                        */
  real_T DragForce_R_e;                /* Mask Parameter: DragForce_R_e
                                        * Referenced by: '<S241>/.5.*A.*Pabs.//R.//T'
                                        */
  real_T KinematicSteering1_StrgRatio; /* Mask Parameter: KinematicSteering1_StrgRatio
                                        * Referenced by: '<S206>/Gain'
                                        */
  real_T KinematicSteering_StrgRatio;  /* Mask Parameter: KinematicSteering_StrgRatio
                                        * Referenced by: '<S187>/Gain'
                                        */
  real_T KinematicSteering1_StrgRng;   /* Mask Parameter: KinematicSteering1_StrgRng
                                        * Referenced by: '<S205>/Saturation'
                                        */
  real_T KinematicSteering_StrgRng;    /* Mask Parameter: KinematicSteering_StrgRng
                                        * Referenced by: '<S186>/Saturation'
                                        */
  real_T FilteredDerivativeDiscreteor_hf;/* Mask Parameter: FilteredDerivativeDiscreteor_hf
                                          * Referenced by: '<S19>/Time constant'
                                          */
  real_T FilteredDerivativeDiscreteor_md;/* Mask Parameter: FilteredDerivativeDiscreteor_md
                                          * Referenced by: '<S23>/Time constant'
                                          */
  real_T FollowerVehicle_Tair;         /* Mask Parameter: FollowerVehicle_Tair
                                        * Referenced by: '<S45>/.5.*A.*Pabs.//R.//T'
                                        */
  real_T LeadVehicle_Tair;             /* Mask Parameter: LeadVehicle_Tair
                                        * Referenced by: '<S241>/.5.*A.*Pabs.//R.//T'
                                        */
  real_T KinematicSteering1_TrckWdth;  /* Mask Parameter: KinematicSteering1_TrckWdth
                                        * Referenced by: '<S207>/Constant1'
                                        */
  real_T KinematicSteering_TrckWdth;   /* Mask Parameter: KinematicSteering_TrckWdth
                                        * Referenced by: '<S188>/Constant1'
                                        */
  real_T KinematicSteering1_WhlBase;   /* Mask Parameter: KinematicSteering1_WhlBase
                                        * Referenced by: '<S207>/Constant'
                                        */
  real_T KinematicSteering_WhlBase;    /* Mask Parameter: KinematicSteering_WhlBase
                                        * Referenced by: '<S188>/Constant'
                                        */
  real_T LeadVehicle_X_o;              /* Mask Parameter: LeadVehicle_X_o
                                        * Referenced by: '<S290>/X_o'
                                        */
  real_T FollowerVehicle_X_o;          /* Mask Parameter: FollowerVehicle_X_o
                                        * Referenced by: '<S94>/X_o'
                                        */
  real_T LeadVehicle_Y_o;              /* Mask Parameter: LeadVehicle_Y_o
                                        * Referenced by: '<S290>/Y_o'
                                        */
  real_T FollowerVehicle_Y_o;          /* Mask Parameter: FollowerVehicle_Y_o
                                        * Referenced by: '<S94>/Y_o'
                                        */
  real_T FollowerVehicle_a;            /* Mask Parameter: FollowerVehicle_a
                                        * Referenced by:
                                        *   '<S8>/vehicle model'
                                        *   '<S45>/Constant3'
                                        *   '<S95>/a'
                                        *   '<S96>/a'
                                        */
  real_T LeadVehicle_a;                /* Mask Parameter: LeadVehicle_a
                                        * Referenced by:
                                        *   '<S11>/vehicle model'
                                        *   '<S241>/Constant3'
                                        *   '<S291>/a'
                                        *   '<S292>/a'
                                        */
  real_T FollowerVehicle_b;            /* Mask Parameter: FollowerVehicle_b
                                        * Referenced by:
                                        *   '<S8>/vehicle model'
                                        *   '<S45>/Constant3'
                                        *   '<S98>/b'
                                        *   '<S99>/b'
                                        */
  real_T LeadVehicle_b;                /* Mask Parameter: LeadVehicle_b
                                        * Referenced by:
                                        *   '<S11>/vehicle model'
                                        *   '<S241>/Constant3'
                                        *   '<S294>/b'
                                        *   '<S295>/b'
                                        */
  real_T FollowerVehicle_beta_w[31];   /* Mask Parameter: FollowerVehicle_beta_w
                                        * Referenced by:
                                        *   '<S45>/Cs'
                                        *   '<S45>/Cym'
                                        */
  real_T LeadVehicle_beta_w[31];       /* Mask Parameter: LeadVehicle_beta_w
                                        * Referenced by:
                                        *   '<S241>/Cs'
                                        *   '<S241>/Cym'
                                        */
  real_T FollowerVehicle_d;            /* Mask Parameter: FollowerVehicle_d
                                        * Referenced by:
                                        *   '<S8>/vehicle model'
                                        *   '<S95>/d'
                                        *   '<S96>/d'
                                        *   '<S98>/d'
                                        *   '<S99>/d'
                                        */
  real_T LeadVehicle_d;                /* Mask Parameter: LeadVehicle_d
                                        * Referenced by:
                                        *   '<S11>/vehicle model'
                                        *   '<S291>/d'
                                        *   '<S292>/d'
                                        *   '<S294>/d'
                                        *   '<S295>/d'
                                        */
  real_T FollowerVehicle_g;            /* Mask Parameter: FollowerVehicle_g
                                        * Referenced by: '<S8>/vehicle model'
                                        */
  real_T LeadVehicle_g;                /* Mask Parameter: LeadVehicle_g
                                        * Referenced by: '<S11>/vehicle model'
                                        */
  real_T FollowerVehicle_h;            /* Mask Parameter: FollowerVehicle_h
                                        * Referenced by:
                                        *   '<S8>/vehicle model'
                                        *   '<S95>/h'
                                        *   '<S96>/h'
                                        *   '<S98>/h'
                                        *   '<S99>/h'
                                        */
  real_T LeadVehicle_h;                /* Mask Parameter: LeadVehicle_h
                                        * Referenced by:
                                        *   '<S11>/vehicle model'
                                        *   '<S291>/h'
                                        *   '<S292>/h'
                                        *   '<S294>/h'
                                        *   '<S295>/h'
                                        */
  real_T LeadVehicle_latOff;           /* Mask Parameter: LeadVehicle_latOff
                                        * Referenced by: '<S293>/latOff'
                                        */
  real_T FollowerVehicle_latOff;       /* Mask Parameter: FollowerVehicle_latOff
                                        * Referenced by: '<S97>/latOff'
                                        */
  real_T LeadVehicle_longOff;          /* Mask Parameter: LeadVehicle_longOff
                                        * Referenced by: '<S293>/longOff'
                                        */
  real_T FollowerVehicle_longOff;      /* Mask Parameter: FollowerVehicle_longOff
                                        * Referenced by: '<S97>/longOff'
                                        */
  real_T FollowerVehicle_m;            /* Mask Parameter: FollowerVehicle_m
                                        * Referenced by: '<S8>/vehicle model'
                                        */
  real_T LeadVehicle_m;                /* Mask Parameter: LeadVehicle_m
                                        * Referenced by: '<S11>/vehicle model'
                                        */
  real_T FilteredDerivativeDiscreteorC_g;/* Mask Parameter: FilteredDerivativeDiscreteorC_g
                                          * Referenced by: '<S19>/Minimum sampling to time constant ratio'
                                          */
  real_T FilteredDerivativeDiscreteorC_d;/* Mask Parameter: FilteredDerivativeDiscreteorC_d
                                          * Referenced by: '<S23>/Minimum sampling to time constant ratio'
                                          */
  real_T LeadVehicle_psi_o;            /* Mask Parameter: LeadVehicle_psi_o
                                        * Referenced by: '<S378>/psi_o'
                                        */
  real_T FollowerVehicle_psi_o;        /* Mask Parameter: FollowerVehicle_psi_o
                                        * Referenced by: '<S182>/psi_o'
                                        */
  real_T LeadVehicle_r_o;              /* Mask Parameter: LeadVehicle_r_o
                                        * Referenced by: '<S378>/r_o'
                                        */
  real_T FollowerVehicle_r_o;          /* Mask Parameter: FollowerVehicle_r_o
                                        * Referenced by: '<S182>/r_o'
                                        */
  real_T FollowerVehicle_sigma_f;      /* Mask Parameter: FollowerVehicle_sigma_f
                                        * Referenced by: '<S178>/Constant1'
                                        */
  real_T LeadVehicle_sigma_f;          /* Mask Parameter: LeadVehicle_sigma_f
                                        * Referenced by: '<S374>/Constant1'
                                        */
  real_T FollowerVehicle_sigma_r;      /* Mask Parameter: FollowerVehicle_sigma_r
                                        * Referenced by: '<S178>/Constant2'
                                        */
  real_T LeadVehicle_sigma_r;          /* Mask Parameter: LeadVehicle_sigma_r
                                        * Referenced by: '<S374>/Constant2'
                                        */
  real_T DeadMansSwitch_stepSize;      /* Mask Parameter: DeadMansSwitch_stepSize
                                        * Referenced by: '<S4>/Simulate step size'
                                        */
  real_T div0protectpoly_thresh;       /* Mask Parameter: div0protectpoly_thresh
                                        * Referenced by:
                                        *   '<S383>/Constant'
                                        *   '<S384>/Constant'
                                        */
  real_T div0protectabspoly3_thresh;   /* Mask Parameter: div0protectabspoly3_thresh
                                        * Referenced by:
                                        *   '<S219>/Constant'
                                        *   '<S220>/Constant'
                                        */
  real_T div0protectabspoly_thresh;    /* Mask Parameter: div0protectabspoly_thresh
                                        * Referenced by:
                                        *   '<S213>/Constant'
                                        *   '<S214>/Constant'
                                        */
  real_T div0protectabspoly3_thresh_g; /* Mask Parameter: div0protectabspoly3_thresh_g
                                        * Referenced by:
                                        *   '<S200>/Constant'
                                        *   '<S201>/Constant'
                                        */
  real_T div0protectabspoly_thresh_b;  /* Mask Parameter: div0protectabspoly_thresh_b
                                        * Referenced by:
                                        *   '<S194>/Constant'
                                        *   '<S195>/Constant'
                                        */
  real_T DeadMansSwitch_timeout;       /* Mask Parameter: DeadMansSwitch_timeout
                                        * Referenced by: '<S4>/Timeout in seconds'
                                        */
  real_T LeadVehicle_vertOff;          /* Mask Parameter: LeadVehicle_vertOff
                                        * Referenced by: '<S293>/vertOff'
                                        */
  real_T FollowerVehicle_vertOff;      /* Mask Parameter: FollowerVehicle_vertOff
                                        * Referenced by: '<S97>/vertOff'
                                        */
  real_T FollowerVehicle_w[2];         /* Mask Parameter: FollowerVehicle_w
                                        * Referenced by:
                                        *   '<S8>/vehicle model'
                                        *   '<S95>/w'
                                        *   '<S96>/w'
                                        *   '<S98>/w'
                                        *   '<S99>/w'
                                        */
  real_T LeadVehicle_w[2];             /* Mask Parameter: LeadVehicle_w
                                        * Referenced by:
                                        *   '<S11>/vehicle model'
                                        *   '<S291>/w'
                                        *   '<S292>/w'
                                        *   '<S294>/w'
                                        *   '<S295>/w'
                                        */
  real_T FollowerVehicle_xdot_tol;     /* Mask Parameter: FollowerVehicle_xdot_tol
                                        * Referenced by: '<S8>/vehicle model'
                                        */
  real_T LeadVehicle_xdot_tol;         /* Mask Parameter: LeadVehicle_xdot_tol
                                        * Referenced by: '<S11>/vehicle model'
                                        */
  real_T LeadVehicle_ydot_o;           /* Mask Parameter: LeadVehicle_ydot_o
                                        * Referenced by: '<S378>/ydot_o'
                                        */
  real_T FollowerVehicle_ydot_o;       /* Mask Parameter: FollowerVehicle_ydot_o
                                        * Referenced by: '<S182>/ydot_o'
                                        */
  SL_Bus_carPlant_smoothUp_followe_Twist_domkw2 Out1_Y0;/* Computed Parameter: Out1_Y0
                                                         * Referenced by: '<S382>/Out1'
                                                         */
  SL_Bus_carPlant_smoothUp_followe_Twist_domkw2 Constant_Value;/* Computed Parameter: Constant_Value
                                                                * Referenced by: '<S15>/Constant'
                                                                */
  SL_Bus_carPlant_smoothUp_followe_Twist_domkw2 Constant_Value_b;/* Computed Parameter: Constant_Value_b
                                                                  * Referenced by: '<S3>/Constant'
                                                                  */
  SL_Bus_carPlant_smoothUp_followe_Twist_domkw2 Constant_Value_k;/* Computed Parameter: Constant_Value_k
                                                                  * Referenced by: '<S2>/Constant'
                                                                  */
  real_T Switch1_Threshold;            /* Expression: 0
                                        * Referenced by: '<S16>/Switch1'
                                        */
  real_T Constant1_Value;              /* Expression: 0
                                        * Referenced by: '<Root>/Constant1'
                                        */
  real_T vehiclemodel_Fxtire_sat;      /* Expression: Fxtire_sat
                                        * Referenced by: '<S8>/vehicle model'
                                        */
  real_T vehiclemodel_Fytire_sat;      /* Expression: Fytire_sat
                                        * Referenced by: '<S8>/vehicle model'
                                        */
  real_T vehiclemodel_Fxtire_sat_d;    /* Expression: Fxtire_sat
                                        * Referenced by: '<S11>/vehicle model'
                                        */
  real_T vehiclemodel_Fytire_sat_e;    /* Expression: Fytire_sat
                                        * Referenced by: '<S11>/vehicle model'
                                        */
  real_T Constant_Value_f;             /* Expression: 0
                                        * Referenced by: '<S5>/Constant'
                                        */
  real_T Integrator_gainval;           /* Computed Parameter: Integrator_gainval
                                        * Referenced by: '<S22>/Integrator'
                                        */
  real_T Integrator_UpperSat;          /* Expression: antiwindupUpperLimit
                                        * Referenced by: '<S22>/Integrator'
                                        */
  real_T Integrator_LowerSat;          /* Expression: antiwindupLowerLimit
                                        * Referenced by: '<S22>/Integrator'
                                        */
  real_T Saturation_UpperSat;          /* Expression: windupUpperLimit
                                        * Referenced by: '<S22>/Saturation'
                                        */
  real_T Saturation_LowerSat;          /* Expression: windupLowerLimit
                                        * Referenced by: '<S22>/Saturation'
                                        */
  real_T Constant_Value_i;             /* Expression: 0
                                        * Referenced by: '<S6>/Constant'
                                        */
  real_T Constant_Value_g;             /* Expression: 1
                                        * Referenced by: '<S16>/Constant'
                                        */
  real_T Switch_Threshold;             /* Expression: 0
                                        * Referenced by: '<Root>/Switch'
                                        */
  real_T Integrator_gainval_p;         /* Computed Parameter: Integrator_gainval_p
                                        * Referenced by: '<S26>/Integrator'
                                        */
  real_T Integrator_UpperSat_h;        /* Expression: antiwindupUpperLimit
                                        * Referenced by: '<S26>/Integrator'
                                        */
  real_T Integrator_LowerSat_j;        /* Expression: antiwindupLowerLimit
                                        * Referenced by: '<S26>/Integrator'
                                        */
  real_T Saturation_UpperSat_h;        /* Expression: windupUpperLimit
                                        * Referenced by: '<S26>/Saturation'
                                        */
  real_T Saturation_LowerSat_i;        /* Expression: windupLowerLimit
                                        * Referenced by: '<S26>/Saturation'
                                        */
  real_T Constant12_Value;             /* Expression: 1
                                        * Referenced by: '<Root>/Constant12'
                                        */
  real_T reference_vel1_Value;         /* Expression: 3.748509500812156
                                        * Referenced by: '<S17>/reference_vel1'
                                        */
  real_T max_accel1_Value;             /* Expression: 1.5
                                        * Referenced by: '<S17>/max_accel1'
                                        */
  real_T max_decel1_Value;             /* Expression: -1.5
                                        * Referenced by: '<S17>/max_decel1'
                                        */
  real_T Steeringangle_Value;          /* Expression: 0.0
                                        * Referenced by: '<Root>/Steering angle'
                                        */
  real_T Constant_Value_o;             /* Expression: 0
                                        * Referenced by: '<S290>/Constant'
                                        */
  real_T Constant7_Value;              /* Expression: 0
                                        * Referenced by: '<S290>/Constant7'
                                        */
  real_T Constant2_Value;              /* Expression: 0
                                        * Referenced by: '<S290>/Constant2'
                                        */
  real_T Constant_Value_e;             /* Expression: 0
                                        * Referenced by: '<S94>/Constant'
                                        */
  real_T Constant7_Value_h;            /* Expression: 0
                                        * Referenced by: '<S94>/Constant7'
                                        */
  real_T Constant2_Value_l;            /* Expression: 0
                                        * Referenced by: '<S94>/Constant2'
                                        */
  real_T decel_Value[3];               /* Expression: [1.5 1.0 0.5]
                                        * Referenced by: '<S13>/decel'
                                        */
  real_T CoefficientofFrictionDry_Value[4];/* Expression: [0.7 0.7; 0.7 0.7]
                                            * Referenced by: '<Root>/Coefficient of Friction: Dry'
                                            */
  real_T Constant3_Value;              /* Expression: 0
                                        * Referenced by: '<Root>/Constant3'
                                        */
  real_T Constant4_Value;              /* Expression: 0
                                        * Referenced by: '<Root>/Constant4'
                                        */
  real_T Crm_tableData[2];             /* Expression: [0 0]
                                        * Referenced by: '<S45>/Crm'
                                        */
  real_T Crm_bp01Data[2];              /* Expression: [-1 1]
                                        * Referenced by: '<S45>/Crm'
                                        */
  real_T Constant4_Value_d[3];         /* Expression: ones(1,3)
                                        * Referenced by: '<S45>/Constant4'
                                        */
  real_T Switch_Threshold_f;           /* Expression: 0
                                        * Referenced by: '<S45>/Switch'
                                        */
  real_T Gain_Gain;                    /* Expression: -1/2
                                        * Referenced by: '<S95>/Gain'
                                        */
  real_T Gain_Gain_k;                  /* Expression: 1/2
                                        * Referenced by: '<S96>/Gain'
                                        */
  real_T Gain_Gain_i;                  /* Expression: -1/2
                                        * Referenced by: '<S98>/Gain'
                                        */
  real_T Gain_Gain_e;                  /* Expression: 1/2
                                        * Referenced by: '<S99>/Gain'
                                        */
  real_T Backlash_InitialOutput;       /* Expression: 0
                                        * Referenced by: '<S205>/Backlash'
                                        */
  real_T index_Value;                  /* Expression: 1
                                        * Referenced by: '<S10>/index'
                                        */
  real_T Switch_Threshold_j;           /* Expression: 0
                                        * Referenced by: '<S10>/Switch'
                                        */
  real_T Switch1_Threshold_h;          /* Expression: 0
                                        * Referenced by: '<S10>/Switch1'
                                        */
  real_T Constant_Value_a;             /* Expression: 0
                                        * Referenced by: '<S174>/Constant'
                                        */
  real_T Constant4_Value_k;            /* Expression: 0
                                        * Referenced by: '<S94>/Constant4'
                                        */
  real_T Constant5_Value;              /* Expression: 0
                                        * Referenced by: '<S94>/Constant5'
                                        */
  real_T Constant6_Value;              /* Expression: 0
                                        * Referenced by: '<S94>/Constant6'
                                        */
  real_T lateral_IC;                   /* Expression: 0
                                        * Referenced by: '<S181>/lateral'
                                        */
  real_T Backlash_InitialOutput_h;     /* Expression: 0
                                        * Referenced by: '<S186>/Backlash'
                                        */
  real_T index_Value_d;                /* Expression: 1
                                        * Referenced by: '<S9>/index'
                                        */
  real_T Switch_Threshold_d;           /* Expression: 0
                                        * Referenced by: '<S9>/Switch'
                                        */
  real_T Switch1_Threshold_p;          /* Expression: 0
                                        * Referenced by: '<S9>/Switch1'
                                        */
  real_T Crm_tableData_n[2];           /* Expression: [0 0]
                                        * Referenced by: '<S241>/Crm'
                                        */
  real_T Crm_bp01Data_p[2];            /* Expression: [-1 1]
                                        * Referenced by: '<S241>/Crm'
                                        */
  real_T Constant4_Value_l[3];         /* Expression: ones(1,3)
                                        * Referenced by: '<S241>/Constant4'
                                        */
  real_T Switch_Threshold_p;           /* Expression: 0
                                        * Referenced by: '<S241>/Switch'
                                        */
  real_T Gain_Gain_b;                  /* Expression: -1/2
                                        * Referenced by: '<S291>/Gain'
                                        */
  real_T Gain_Gain_a;                  /* Expression: 1/2
                                        * Referenced by: '<S292>/Gain'
                                        */
  real_T Gain_Gain_be;                 /* Expression: -1/2
                                        * Referenced by: '<S294>/Gain'
                                        */
  real_T Gain_Gain_p;                  /* Expression: 1/2
                                        * Referenced by: '<S295>/Gain'
                                        */
  real_T Constant_Value_bd;            /* Expression: 0
                                        * Referenced by: '<S370>/Constant'
                                        */
  real_T Constant_Value_j;             /* Expression: mu
                                        * Referenced by: '<S350>/Constant'
                                        */
  real_T Constant4_Value_o;            /* Expression: 0
                                        * Referenced by: '<S290>/Constant4'
                                        */
  real_T Constant5_Value_j;            /* Expression: 0
                                        * Referenced by: '<S290>/Constant5'
                                        */
  real_T Constant6_Value_d;            /* Expression: 0
                                        * Referenced by: '<S290>/Constant6'
                                        */
  real_T lateral_IC_g;                 /* Expression: 0
                                        * Referenced by: '<S377>/lateral'
                                        */
};

/* Real-time Model Data Structure */
struct tag_RTM_carPlant_smoothUp_fol_T {
  const char_T *errorStatus;
  RTWLogInfo *rtwLogInfo;
  RTWSolverInfo solverInfo;
  X_carPlant_smoothUp_followerS_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  boolean_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeY[20];
  real_T odeF[3][20];
  ODE3_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    boolean_T firstInitCondFlag;
    time_T tFinal;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

/* Block parameters (default storage) */
extern P_carPlant_smoothUp_followerS_T carPlant_smoothUp_followerSto_P;

/* Block signals (default storage) */
extern B_carPlant_smoothUp_followerS_T carPlant_smoothUp_followerSto_B;

/* Continuous states (default storage) */
extern X_carPlant_smoothUp_followerS_T carPlant_smoothUp_followerSto_X;

/* Block states (default storage) */
extern DW_carPlant_smoothUp_follower_T carPlant_smoothUp_followerSt_DW;

/* Model entry point functions */
extern void carPlant_smoothUp_followerStopper_initialize(void);
extern void carPlant_smoothUp_followerStopper_step(void);
extern void carPlant_smoothUp_followerStopper_terminate(void);

/* Real-time Model object */
extern RT_MODEL_carPlant_smoothUp_fo_T *const carPlant_smoothUp_followerSt_M;

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'carPlant_smoothUp_followerStopper'
 * '<S1>'   : 'carPlant_smoothUp_followerStopper/Blank Message'
 * '<S2>'   : 'carPlant_smoothUp_followerStopper/Blank Message Leader'
 * '<S3>'   : 'carPlant_smoothUp_followerStopper/Blank Message7'
 * '<S4>'   : 'carPlant_smoothUp_followerStopper/Dead Man's Switch'
 * '<S5>'   : 'carPlant_smoothUp_followerStopper/Filtered Derivative (Discrete or Continuous)'
 * '<S6>'   : 'carPlant_smoothUp_followerStopper/Filtered Derivative (Discrete or Continuous)1'
 * '<S7>'   : 'carPlant_smoothUp_followerStopper/Follower Stopper Controller'
 * '<S8>'   : 'carPlant_smoothUp_followerStopper/Follower Vehicle'
 * '<S9>'   : 'carPlant_smoothUp_followerStopper/Kinematic Steering'
 * '<S10>'  : 'carPlant_smoothUp_followerStopper/Kinematic Steering1'
 * '<S11>'  : 'carPlant_smoothUp_followerStopper/Lead Vehicle'
 * '<S12>'  : 'carPlant_smoothUp_followerStopper/MATLAB Function'
 * '<S13>'  : 'carPlant_smoothUp_followerStopper/Parameters'
 * '<S14>'  : 'carPlant_smoothUp_followerStopper/Publish'
 * '<S15>'  : 'carPlant_smoothUp_followerStopper/Subscribe'
 * '<S16>'  : 'carPlant_smoothUp_followerStopper/div0protect - poly'
 * '<S17>'  : 'carPlant_smoothUp_followerStopper/smoothUp Parameters'
 * '<S18>'  : 'carPlant_smoothUp_followerStopper/Dead Man's Switch/timeout set to 0 output'
 * '<S19>'  : 'carPlant_smoothUp_followerStopper/Filtered Derivative (Discrete or Continuous)/Enable//disable time constant'
 * '<S20>'  : 'carPlant_smoothUp_followerStopper/Filtered Derivative (Discrete or Continuous)/Integrator (Discrete or Continuous)'
 * '<S21>'  : 'carPlant_smoothUp_followerStopper/Filtered Derivative (Discrete or Continuous)/Integrator (Discrete or Continuous)/Continuous'
 * '<S22>'  : 'carPlant_smoothUp_followerStopper/Filtered Derivative (Discrete or Continuous)/Integrator (Discrete or Continuous)/Discrete'
 * '<S23>'  : 'carPlant_smoothUp_followerStopper/Filtered Derivative (Discrete or Continuous)1/Enable//disable time constant'
 * '<S24>'  : 'carPlant_smoothUp_followerStopper/Filtered Derivative (Discrete or Continuous)1/Integrator (Discrete or Continuous)'
 * '<S25>'  : 'carPlant_smoothUp_followerStopper/Filtered Derivative (Discrete or Continuous)1/Integrator (Discrete or Continuous)/Continuous'
 * '<S26>'  : 'carPlant_smoothUp_followerStopper/Filtered Derivative (Discrete or Continuous)1/Integrator (Discrete or Continuous)/Discrete'
 * '<S27>'  : 'carPlant_smoothUp_followerStopper/Follower Stopper Controller/MATLAB Function'
 * '<S28>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Cy'
 * '<S29>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Drag'
 * '<S30>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing'
 * '<S31>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/friction'
 * '<S32>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/front forces'
 * '<S33>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/front steer'
 * '<S34>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/rear forces'
 * '<S35>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/rear steer'
 * '<S36>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/sigma'
 * '<S37>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/state'
 * '<S38>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/vehicle model'
 * '<S39>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/wind'
 * '<S40>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Cy/Cy const'
 * '<S41>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Cy/Cy const dual'
 * '<S42>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Cy/Cy table'
 * '<S43>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Cy/Cy table dual'
 * '<S44>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Cy/Cy table dual/For Each Subsystem'
 * '<S45>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Drag/Drag Force'
 * '<S46>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Drag/inertial2body'
 * '<S47>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing'
 * '<S48>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual'
 * '<S49>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Forces 3DOF'
 * '<S50>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF'
 * '<S51>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Moments'
 * '<S52>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Power'
 * '<S53>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/state2bus'
 * '<S54>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Front'
 * '<S55>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric'
 * '<S56>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Rear'
 * '<S57>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Front/Rotation Angles to Direction Cosine Matrix'
 * '<S58>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Front/transform to Inertial axes'
 * '<S59>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Front/transform to Inertial axes1'
 * '<S60>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Front/wxR'
 * '<S61>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Front/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S62>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Front/wxR/Subsystem'
 * '<S63>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Front/wxR/Subsystem1'
 * '<S64>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta'
 * '<S65>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Body Slip'
 * '<S66>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Rotation Angles to Direction Cosine Matrix'
 * '<S67>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/transform to Inertial axes'
 * '<S68>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/transform to Inertial axes1'
 * '<S69>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/wxR'
 * '<S70>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Body Slip/div0protect - abs poly'
 * '<S71>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Body Slip/div0protect - abs poly/Compare To Constant'
 * '<S72>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Body Slip/div0protect - abs poly/Compare To Constant1'
 * '<S73>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S74>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/wxR/Subsystem'
 * '<S75>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/wxR/Subsystem1'
 * '<S76>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Rear/Rotation Angles to Direction Cosine Matrix'
 * '<S77>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Rear/transform to Inertial axes'
 * '<S78>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Rear/transform to Inertial axes1'
 * '<S79>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Rear/wxR'
 * '<S80>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Rear/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S81>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Rear/wxR/Subsystem'
 * '<S82>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Rear/wxR/Subsystem1'
 * '<S83>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/state2bus/Body Slip'
 * '<S84>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/state2bus/COMB2I'
 * '<S85>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/state2bus/xddot2ax'
 * '<S86>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/state2bus/Body Slip/div0protect - abs poly'
 * '<S87>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/state2bus/Body Slip/div0protect - abs poly/Compare To Constant'
 * '<S88>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/state2bus/Body Slip/div0protect - abs poly/Compare To Constant1'
 * '<S89>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing/state2bus/xddot2ax/m^22gn'
 * '<S90>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Forces 3DOF'
 * '<S91>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF'
 * '<S92>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Moments'
 * '<S93>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Power'
 * '<S94>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/state2bus'
 * '<S95>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Left'
 * '<S96>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Right'
 * '<S97>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric'
 * '<S98>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Left'
 * '<S99>'  : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Right'
 * '<S100>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Left/Hard Point Coordinate Transform External Displacement'
 * '<S101>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Left/Hard Point Coordinate Transform External Displacement/Rotation Angles to Direction Cosine Matrix'
 * '<S102>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Left/Hard Point Coordinate Transform External Displacement/transform to Inertial axes'
 * '<S103>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Left/Hard Point Coordinate Transform External Displacement/transform to Inertial axes1'
 * '<S104>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Left/Hard Point Coordinate Transform External Displacement/wxR'
 * '<S105>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Left/Hard Point Coordinate Transform External Displacement/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S106>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Left/Hard Point Coordinate Transform External Displacement/wxR/Subsystem'
 * '<S107>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Left/Hard Point Coordinate Transform External Displacement/wxR/Subsystem1'
 * '<S108>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Right/Hard Point Coordinate Transform External Displacement'
 * '<S109>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Right/Hard Point Coordinate Transform External Displacement/Rotation Angles to Direction Cosine Matrix'
 * '<S110>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Right/Hard Point Coordinate Transform External Displacement/transform to Inertial axes'
 * '<S111>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Right/Hard Point Coordinate Transform External Displacement/transform to Inertial axes1'
 * '<S112>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Right/Hard Point Coordinate Transform External Displacement/wxR'
 * '<S113>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Right/Hard Point Coordinate Transform External Displacement/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S114>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Right/Hard Point Coordinate Transform External Displacement/wxR/Subsystem'
 * '<S115>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Right/Hard Point Coordinate Transform External Displacement/wxR/Subsystem1'
 * '<S116>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta'
 * '<S117>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Body Slip'
 * '<S118>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Rotation Angles to Direction Cosine Matrix'
 * '<S119>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/transform to Inertial axes'
 * '<S120>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/transform to Inertial axes1'
 * '<S121>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/wxR'
 * '<S122>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Body Slip/div0protect - abs poly'
 * '<S123>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Body Slip/div0protect - abs poly/Compare To Constant'
 * '<S124>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Body Slip/div0protect - abs poly/Compare To Constant1'
 * '<S125>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S126>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/wxR/Subsystem'
 * '<S127>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/wxR/Subsystem1'
 * '<S128>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Left/Hard Point Coordinate Transform External Displacement'
 * '<S129>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Left/Hard Point Coordinate Transform External Displacement/Rotation Angles to Direction Cosine Matrix'
 * '<S130>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Left/Hard Point Coordinate Transform External Displacement/transform to Inertial axes'
 * '<S131>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Left/Hard Point Coordinate Transform External Displacement/transform to Inertial axes1'
 * '<S132>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Left/Hard Point Coordinate Transform External Displacement/wxR'
 * '<S133>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Left/Hard Point Coordinate Transform External Displacement/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S134>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Left/Hard Point Coordinate Transform External Displacement/wxR/Subsystem'
 * '<S135>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Left/Hard Point Coordinate Transform External Displacement/wxR/Subsystem1'
 * '<S136>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Right/Hard Point Coordinate Transform External Displacement'
 * '<S137>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Right/Hard Point Coordinate Transform External Displacement/Rotation Angles to Direction Cosine Matrix'
 * '<S138>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Right/Hard Point Coordinate Transform External Displacement/transform to Inertial axes'
 * '<S139>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Right/Hard Point Coordinate Transform External Displacement/transform to Inertial axes1'
 * '<S140>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Right/Hard Point Coordinate Transform External Displacement/wxR'
 * '<S141>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Right/Hard Point Coordinate Transform External Displacement/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S142>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Right/Hard Point Coordinate Transform External Displacement/wxR/Subsystem'
 * '<S143>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Right/Hard Point Coordinate Transform External Displacement/wxR/Subsystem1'
 * '<S144>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/state2bus/Body Slip'
 * '<S145>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/state2bus/COMB2I'
 * '<S146>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/state2bus/xddot2ax'
 * '<S147>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/state2bus/Body Slip/div0protect - abs poly'
 * '<S148>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/state2bus/Body Slip/div0protect - abs poly/Compare To Constant'
 * '<S149>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/state2bus/Body Slip/div0protect - abs poly/Compare To Constant1'
 * '<S150>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/Signal Routing/Signal Routing Dual/state2bus/xddot2ax/m^22gn'
 * '<S151>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/friction/mu ext'
 * '<S152>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/friction/mu ext dual'
 * '<S153>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/friction/mu int'
 * '<S154>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/friction/mu int dual'
 * '<S155>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/front forces/ext'
 * '<S156>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/front forces/ext dual'
 * '<S157>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/front forces/ext long'
 * '<S158>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/front forces/ext long dual'
 * '<S159>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/front forces/int'
 * '<S160>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/front forces/int dual'
 * '<S161>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/front steer/delta ext'
 * '<S162>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/front steer/delta ext dual'
 * '<S163>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/front steer/delta int'
 * '<S164>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/front steer/delta int dual'
 * '<S165>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/rear forces/ext'
 * '<S166>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/rear forces/ext dual'
 * '<S167>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/rear forces/ext long'
 * '<S168>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/rear forces/ext long dual'
 * '<S169>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/rear forces/int'
 * '<S170>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/rear forces/int dual'
 * '<S171>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/rear steer/delta ext'
 * '<S172>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/rear steer/delta ext dual'
 * '<S173>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/rear steer/delta int'
 * '<S174>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/rear steer/delta int dual'
 * '<S175>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/sigma/no sigma'
 * '<S176>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/sigma/no sigma dual'
 * '<S177>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/sigma/sigma'
 * '<S178>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/sigma/sigma dual'
 * '<S179>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/sigma/sigma/relaxation front'
 * '<S180>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/sigma/sigma/relaxation rear'
 * '<S181>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/sigma/sigma dual/relaxation'
 * '<S182>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/state/xdot ext'
 * '<S183>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/state/xdot int'
 * '<S184>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/wind/wind ext'
 * '<S185>' : 'carPlant_smoothUp_followerStopper/Follower Vehicle/wind/wind int'
 * '<S186>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering/AngInput'
 * '<S187>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering/AngInput/AckermanConstantRatio'
 * '<S188>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering/AngInput/AckermanConstantRatio/Ackerman'
 * '<S189>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering/AngInput/AckermanConstantRatio/AckermanJacobian'
 * '<S190>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering/AngInput/AckermanConstantRatio/Ackerman/div0protect - abs poly'
 * '<S191>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering/AngInput/AckermanConstantRatio/Ackerman/div0protect - abs poly1'
 * '<S192>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering/AngInput/AckermanConstantRatio/Ackerman/div0protect - abs poly2'
 * '<S193>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering/AngInput/AckermanConstantRatio/Ackerman/div0protect - abs poly3'
 * '<S194>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering/AngInput/AckermanConstantRatio/Ackerman/div0protect - abs poly/Compare To Constant'
 * '<S195>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering/AngInput/AckermanConstantRatio/Ackerman/div0protect - abs poly/Compare To Constant1'
 * '<S196>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering/AngInput/AckermanConstantRatio/Ackerman/div0protect - abs poly1/Compare To Constant'
 * '<S197>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering/AngInput/AckermanConstantRatio/Ackerman/div0protect - abs poly1/Compare To Constant1'
 * '<S198>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering/AngInput/AckermanConstantRatio/Ackerman/div0protect - abs poly2/Compare To Constant'
 * '<S199>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering/AngInput/AckermanConstantRatio/Ackerman/div0protect - abs poly2/Compare To Constant1'
 * '<S200>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering/AngInput/AckermanConstantRatio/Ackerman/div0protect - abs poly3/Compare To Constant'
 * '<S201>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering/AngInput/AckermanConstantRatio/Ackerman/div0protect - abs poly3/Compare To Constant1'
 * '<S202>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering/AngInput/AckermanConstantRatio/AckermanJacobian/div0protect - abs poly1'
 * '<S203>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering/AngInput/AckermanConstantRatio/AckermanJacobian/div0protect - abs poly1/Compare To Constant'
 * '<S204>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering/AngInput/AckermanConstantRatio/AckermanJacobian/div0protect - abs poly1/Compare To Constant1'
 * '<S205>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering1/AngInput'
 * '<S206>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering1/AngInput/AckermanConstantRatio'
 * '<S207>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering1/AngInput/AckermanConstantRatio/Ackerman'
 * '<S208>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering1/AngInput/AckermanConstantRatio/AckermanJacobian'
 * '<S209>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering1/AngInput/AckermanConstantRatio/Ackerman/div0protect - abs poly'
 * '<S210>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering1/AngInput/AckermanConstantRatio/Ackerman/div0protect - abs poly1'
 * '<S211>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering1/AngInput/AckermanConstantRatio/Ackerman/div0protect - abs poly2'
 * '<S212>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering1/AngInput/AckermanConstantRatio/Ackerman/div0protect - abs poly3'
 * '<S213>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering1/AngInput/AckermanConstantRatio/Ackerman/div0protect - abs poly/Compare To Constant'
 * '<S214>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering1/AngInput/AckermanConstantRatio/Ackerman/div0protect - abs poly/Compare To Constant1'
 * '<S215>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering1/AngInput/AckermanConstantRatio/Ackerman/div0protect - abs poly1/Compare To Constant'
 * '<S216>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering1/AngInput/AckermanConstantRatio/Ackerman/div0protect - abs poly1/Compare To Constant1'
 * '<S217>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering1/AngInput/AckermanConstantRatio/Ackerman/div0protect - abs poly2/Compare To Constant'
 * '<S218>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering1/AngInput/AckermanConstantRatio/Ackerman/div0protect - abs poly2/Compare To Constant1'
 * '<S219>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering1/AngInput/AckermanConstantRatio/Ackerman/div0protect - abs poly3/Compare To Constant'
 * '<S220>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering1/AngInput/AckermanConstantRatio/Ackerman/div0protect - abs poly3/Compare To Constant1'
 * '<S221>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering1/AngInput/AckermanConstantRatio/AckermanJacobian/div0protect - abs poly1'
 * '<S222>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering1/AngInput/AckermanConstantRatio/AckermanJacobian/div0protect - abs poly1/Compare To Constant'
 * '<S223>' : 'carPlant_smoothUp_followerStopper/Kinematic Steering1/AngInput/AckermanConstantRatio/AckermanJacobian/div0protect - abs poly1/Compare To Constant1'
 * '<S224>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Cy'
 * '<S225>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Drag'
 * '<S226>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing'
 * '<S227>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/friction'
 * '<S228>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/front forces'
 * '<S229>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/front steer'
 * '<S230>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/rear forces'
 * '<S231>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/rear steer'
 * '<S232>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/sigma'
 * '<S233>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/state'
 * '<S234>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/vehicle model'
 * '<S235>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/wind'
 * '<S236>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Cy/Cy const'
 * '<S237>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Cy/Cy const dual'
 * '<S238>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Cy/Cy table'
 * '<S239>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Cy/Cy table dual'
 * '<S240>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Cy/Cy table dual/For Each Subsystem'
 * '<S241>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Drag/Drag Force'
 * '<S242>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Drag/inertial2body'
 * '<S243>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing'
 * '<S244>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual'
 * '<S245>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Forces 3DOF'
 * '<S246>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF'
 * '<S247>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Moments'
 * '<S248>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Power'
 * '<S249>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/state2bus'
 * '<S250>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Front'
 * '<S251>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric'
 * '<S252>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Rear'
 * '<S253>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Front/Rotation Angles to Direction Cosine Matrix'
 * '<S254>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Front/transform to Inertial axes'
 * '<S255>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Front/transform to Inertial axes1'
 * '<S256>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Front/wxR'
 * '<S257>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Front/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S258>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Front/wxR/Subsystem'
 * '<S259>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Front/wxR/Subsystem1'
 * '<S260>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta'
 * '<S261>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Body Slip'
 * '<S262>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Rotation Angles to Direction Cosine Matrix'
 * '<S263>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/transform to Inertial axes'
 * '<S264>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/transform to Inertial axes1'
 * '<S265>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/wxR'
 * '<S266>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Body Slip/div0protect - abs poly'
 * '<S267>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Body Slip/div0protect - abs poly/Compare To Constant'
 * '<S268>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Body Slip/div0protect - abs poly/Compare To Constant1'
 * '<S269>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S270>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/wxR/Subsystem'
 * '<S271>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/wxR/Subsystem1'
 * '<S272>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Rear/Rotation Angles to Direction Cosine Matrix'
 * '<S273>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Rear/transform to Inertial axes'
 * '<S274>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Rear/transform to Inertial axes1'
 * '<S275>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Rear/wxR'
 * '<S276>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Rear/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S277>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Rear/wxR/Subsystem'
 * '<S278>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/Lateral 3DOF/Hard Point Coordinate Transform Rear/wxR/Subsystem1'
 * '<S279>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/state2bus/Body Slip'
 * '<S280>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/state2bus/COMB2I'
 * '<S281>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/state2bus/xddot2ax'
 * '<S282>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/state2bus/Body Slip/div0protect - abs poly'
 * '<S283>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/state2bus/Body Slip/div0protect - abs poly/Compare To Constant'
 * '<S284>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/state2bus/Body Slip/div0protect - abs poly/Compare To Constant1'
 * '<S285>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing/state2bus/xddot2ax/m^22gn'
 * '<S286>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Forces 3DOF'
 * '<S287>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF'
 * '<S288>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Moments'
 * '<S289>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Power'
 * '<S290>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/state2bus'
 * '<S291>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Left'
 * '<S292>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Right'
 * '<S293>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric'
 * '<S294>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Left'
 * '<S295>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Right'
 * '<S296>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Left/Hard Point Coordinate Transform External Displacement'
 * '<S297>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Left/Hard Point Coordinate Transform External Displacement/Rotation Angles to Direction Cosine Matrix'
 * '<S298>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Left/Hard Point Coordinate Transform External Displacement/transform to Inertial axes'
 * '<S299>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Left/Hard Point Coordinate Transform External Displacement/transform to Inertial axes1'
 * '<S300>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Left/Hard Point Coordinate Transform External Displacement/wxR'
 * '<S301>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Left/Hard Point Coordinate Transform External Displacement/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S302>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Left/Hard Point Coordinate Transform External Displacement/wxR/Subsystem'
 * '<S303>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Left/Hard Point Coordinate Transform External Displacement/wxR/Subsystem1'
 * '<S304>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Right/Hard Point Coordinate Transform External Displacement'
 * '<S305>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Right/Hard Point Coordinate Transform External Displacement/Rotation Angles to Direction Cosine Matrix'
 * '<S306>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Right/Hard Point Coordinate Transform External Displacement/transform to Inertial axes'
 * '<S307>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Right/Hard Point Coordinate Transform External Displacement/transform to Inertial axes1'
 * '<S308>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Right/Hard Point Coordinate Transform External Displacement/wxR'
 * '<S309>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Right/Hard Point Coordinate Transform External Displacement/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S310>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Right/Hard Point Coordinate Transform External Displacement/wxR/Subsystem'
 * '<S311>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Front Right/Hard Point Coordinate Transform External Displacement/wxR/Subsystem1'
 * '<S312>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta'
 * '<S313>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Body Slip'
 * '<S314>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Rotation Angles to Direction Cosine Matrix'
 * '<S315>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/transform to Inertial axes'
 * '<S316>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/transform to Inertial axes1'
 * '<S317>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/wxR'
 * '<S318>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Body Slip/div0protect - abs poly'
 * '<S319>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Body Slip/div0protect - abs poly/Compare To Constant'
 * '<S320>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Body Slip/div0protect - abs poly/Compare To Constant1'
 * '<S321>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S322>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/wxR/Subsystem'
 * '<S323>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Geometric/Hard Point Coordinate Transform External Displacement Beta/wxR/Subsystem1'
 * '<S324>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Left/Hard Point Coordinate Transform External Displacement'
 * '<S325>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Left/Hard Point Coordinate Transform External Displacement/Rotation Angles to Direction Cosine Matrix'
 * '<S326>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Left/Hard Point Coordinate Transform External Displacement/transform to Inertial axes'
 * '<S327>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Left/Hard Point Coordinate Transform External Displacement/transform to Inertial axes1'
 * '<S328>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Left/Hard Point Coordinate Transform External Displacement/wxR'
 * '<S329>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Left/Hard Point Coordinate Transform External Displacement/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S330>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Left/Hard Point Coordinate Transform External Displacement/wxR/Subsystem'
 * '<S331>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Left/Hard Point Coordinate Transform External Displacement/wxR/Subsystem1'
 * '<S332>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Right/Hard Point Coordinate Transform External Displacement'
 * '<S333>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Right/Hard Point Coordinate Transform External Displacement/Rotation Angles to Direction Cosine Matrix'
 * '<S334>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Right/Hard Point Coordinate Transform External Displacement/transform to Inertial axes'
 * '<S335>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Right/Hard Point Coordinate Transform External Displacement/transform to Inertial axes1'
 * '<S336>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Right/Hard Point Coordinate Transform External Displacement/wxR'
 * '<S337>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Right/Hard Point Coordinate Transform External Displacement/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S338>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Right/Hard Point Coordinate Transform External Displacement/wxR/Subsystem'
 * '<S339>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/Lateral 3DOF/Hard Point Coordinate Transform Rear Right/Hard Point Coordinate Transform External Displacement/wxR/Subsystem1'
 * '<S340>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/state2bus/Body Slip'
 * '<S341>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/state2bus/COMB2I'
 * '<S342>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/state2bus/xddot2ax'
 * '<S343>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/state2bus/Body Slip/div0protect - abs poly'
 * '<S344>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/state2bus/Body Slip/div0protect - abs poly/Compare To Constant'
 * '<S345>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/state2bus/Body Slip/div0protect - abs poly/Compare To Constant1'
 * '<S346>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/Signal Routing/Signal Routing Dual/state2bus/xddot2ax/m^22gn'
 * '<S347>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/friction/mu ext'
 * '<S348>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/friction/mu ext dual'
 * '<S349>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/friction/mu int'
 * '<S350>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/friction/mu int dual'
 * '<S351>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/front forces/ext'
 * '<S352>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/front forces/ext dual'
 * '<S353>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/front forces/ext long'
 * '<S354>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/front forces/ext long dual'
 * '<S355>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/front forces/int'
 * '<S356>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/front forces/int dual'
 * '<S357>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/front steer/delta ext'
 * '<S358>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/front steer/delta ext dual'
 * '<S359>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/front steer/delta int'
 * '<S360>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/front steer/delta int dual'
 * '<S361>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/rear forces/ext'
 * '<S362>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/rear forces/ext dual'
 * '<S363>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/rear forces/ext long'
 * '<S364>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/rear forces/ext long dual'
 * '<S365>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/rear forces/int'
 * '<S366>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/rear forces/int dual'
 * '<S367>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/rear steer/delta ext'
 * '<S368>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/rear steer/delta ext dual'
 * '<S369>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/rear steer/delta int'
 * '<S370>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/rear steer/delta int dual'
 * '<S371>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/sigma/no sigma'
 * '<S372>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/sigma/no sigma dual'
 * '<S373>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/sigma/sigma'
 * '<S374>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/sigma/sigma dual'
 * '<S375>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/sigma/sigma/relaxation front'
 * '<S376>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/sigma/sigma/relaxation rear'
 * '<S377>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/sigma/sigma dual/relaxation'
 * '<S378>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/state/xdot ext'
 * '<S379>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/state/xdot int'
 * '<S380>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/wind/wind ext'
 * '<S381>' : 'carPlant_smoothUp_followerStopper/Lead Vehicle/wind/wind int'
 * '<S382>' : 'carPlant_smoothUp_followerStopper/Subscribe/Enabled Subsystem'
 * '<S383>' : 'carPlant_smoothUp_followerStopper/div0protect - poly/Compare To Constant'
 * '<S384>' : 'carPlant_smoothUp_followerStopper/div0protect - poly/Compare To Constant1'
 */
#endif                                 /* RTW_HEADER_carPlant_smoothUp_followerStopper_h_ */
