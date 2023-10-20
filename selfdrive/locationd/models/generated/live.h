#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_8401443258435973941);
void live_err_fun(double *nom_x, double *delta_x, double *out_3011540060554697934);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_2329848449004316646);
void live_H_mod_fun(double *state, double *out_5953446052952728839);
void live_f_fun(double *state, double dt, double *out_576689838510283432);
void live_F_fun(double *state, double dt, double *out_4959727181402762163);
void live_h_4(double *state, double *unused, double *out_6587257671490569481);
void live_H_4(double *state, double *unused, double *out_1072245686001335260);
void live_h_9(double *state, double *unused, double *out_608170080952672459);
void live_H_9(double *state, double *unused, double *out_6214973249263112210);
void live_h_10(double *state, double *unused, double *out_6232880459400070621);
void live_H_10(double *state, double *unused, double *out_922693051656624723);
void live_h_12(double *state, double *unused, double *out_4570450692652821415);
void live_H_12(double *state, double *unused, double *out_3947210722030626535);
void live_h_31(double *state, double *unused, double *out_7993436279119541516);
void live_H_31(double *state, double *unused, double *out_6692773754355640244);
void live_h_32(double *state, double *unused, double *out_739000771032314965);
void live_H_32(double *state, double *unused, double *out_2425581156312012585);
void live_h_13(double *state, double *unused, double *out_8307942746337328291);
void live_H_13(double *state, double *unused, double *out_4318811407375332569);
void live_h_14(double *state, double *unused, double *out_608170080952672459);
void live_H_14(double *state, double *unused, double *out_6214973249263112210);
void live_h_33(double *state, double *unused, double *out_7702337563609058904);
void live_H_33(double *state, double *unused, double *out_8603413314715053768);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}