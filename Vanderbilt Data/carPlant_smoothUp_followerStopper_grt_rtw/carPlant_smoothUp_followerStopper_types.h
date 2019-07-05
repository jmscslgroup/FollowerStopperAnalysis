/*
 * carPlant_smoothUp_followerStopper_types.h
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

#ifndef RTW_HEADER_carPlant_smoothUp_followerStopper_types_h_
#define RTW_HEADER_carPlant_smoothUp_followerStopper_types_h_
#include "rtwtypes.h"
#include "builtin_typeid_types.h"
#include "multiword_types.h"
#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_carPlant_smoothUp_followe_Point_dmaai3_
#define DEFINED_TYPEDEF_FOR_SL_Bus_carPlant_smoothUp_followe_Point_dmaai3_

typedef struct {
  real_T X;
  real_T Y;
  real_T Z;
} SL_Bus_carPlant_smoothUp_followe_Point_dmaai3;

#endif

#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_carPlant_smoothUp_followe_Vector3_ufmq3v_
#define DEFINED_TYPEDEF_FOR_SL_Bus_carPlant_smoothUp_followe_Vector3_ufmq3v_

typedef struct {
  real_T X;
  real_T Y;
  real_T Z;
} SL_Bus_carPlant_smoothUp_followe_Vector3_ufmq3v;

#endif

#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_carPlant_smoothUp_followe_Twist_domkw2_
#define DEFINED_TYPEDEF_FOR_SL_Bus_carPlant_smoothUp_followe_Twist_domkw2_

typedef struct {
  SL_Bus_carPlant_smoothUp_followe_Vector3_ufmq3v Linear;
  SL_Bus_carPlant_smoothUp_followe_Vector3_ufmq3v Angular;
} SL_Bus_carPlant_smoothUp_followe_Twist_domkw2;

#endif

#ifndef typedef_robotics_slcore_internal_bloc_T
#define typedef_robotics_slcore_internal_bloc_T

typedef struct {
  int32_T __dummy;
} robotics_slcore_internal_bloc_T;

#endif                                 /*typedef_robotics_slcore_internal_bloc_T*/

#ifndef typedef_robotics_slros_internal_block_T
#define typedef_robotics_slros_internal_block_T

typedef struct {
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T isSetupComplete;
  real_T ticksUntilNextHit;
  robotics_slcore_internal_bloc_T SampleTimeHandler;
} robotics_slros_internal_block_T;

#endif                                 /*typedef_robotics_slros_internal_block_T*/

#ifndef typedef_robotics_slros_internal_blo_e_T
#define typedef_robotics_slros_internal_blo_e_T

typedef struct {
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T isSetupComplete;
} robotics_slros_internal_blo_e_T;

#endif                                 /*typedef_robotics_slros_internal_blo_e_T*/

#ifndef typedef_robotics_slros_internal_bl_em_T
#define typedef_robotics_slros_internal_bl_em_T

typedef struct {
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T isSetupComplete;
} robotics_slros_internal_bl_em_T;

#endif                                 /*typedef_robotics_slros_internal_bl_em_T*/

#ifndef typedef_struct_T_carPlant_smoothUp_fo_T
#define typedef_struct_T_carPlant_smoothUp_fo_T

typedef struct {
  char_T Value[4];
} struct_T_carPlant_smoothUp_fo_T;

#endif                                 /*typedef_struct_T_carPlant_smoothUp_fo_T*/

#ifndef typedef_struct_T_carPlant_smoothUp__e_T
#define typedef_struct_T_carPlant_smoothUp__e_T

typedef struct {
  char_T Value[9];
} struct_T_carPlant_smoothUp__e_T;

#endif                                 /*typedef_struct_T_carPlant_smoothUp__e_T*/

#ifndef struct_tag_smnSVBYMKVOFO5RozNZjEpF
#define struct_tag_smnSVBYMKVOFO5RozNZjEpF

struct tag_smnSVBYMKVOFO5RozNZjEpF
{
  char_T Disallow[9];
  char_T Type[9];
};

#endif                                 /*struct_tag_smnSVBYMKVOFO5RozNZjEpF*/

#ifndef typedef_smnSVBYMKVOFO5RozNZjEpF_carPl_T
#define typedef_smnSVBYMKVOFO5RozNZjEpF_carPl_T

typedef struct tag_smnSVBYMKVOFO5RozNZjEpF smnSVBYMKVOFO5RozNZjEpF_carPl_T;

#endif                                 /*typedef_smnSVBYMKVOFO5RozNZjEpF_carPl_T*/

/* Parameters (default storage) */
typedef struct P_carPlant_smoothUp_followerS_T_ P_carPlant_smoothUp_followerS_T;

/* Forward declaration for rtModel */
typedef struct tag_RTM_carPlant_smoothUp_fol_T RT_MODEL_carPlant_smoothUp_fo_T;

#endif                                 /* RTW_HEADER_carPlant_smoothUp_followerStopper_types_h_ */
