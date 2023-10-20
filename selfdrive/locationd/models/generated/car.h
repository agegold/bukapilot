#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_628289647837459206);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_1927844397472153779);
void car_H_mod_fun(double *state, double *out_8843017787680630764);
void car_f_fun(double *state, double dt, double *out_7261182143492256267);
void car_F_fun(double *state, double dt, double *out_8085554514994619377);
void car_h_25(double *state, double *unused, double *out_4515812170498523752);
void car_H_25(double *state, double *unused, double *out_1773635699508951289);
void car_h_24(double *state, double *unused, double *out_790092840729285282);
void car_H_24(double *state, double *unused, double *out_3046685805147036974);
void car_h_30(double *state, double *unused, double *out_2255023055851827184);
void car_H_30(double *state, double *unused, double *out_1644296752365711219);
void car_h_26(double *state, double *unused, double *out_7256696323512672829);
void car_H_26(double *state, double *unused, double *out_2430489763619263193);
void car_h_27(double *state, double *unused, double *out_7138688929579296085);
void car_H_27(double *state, double *unused, double *out_530466559434713692);
void car_h_29(double *state, double *unused, double *out_6462813864016865142);
void car_H_29(double *state, double *unused, double *out_2154528096680103403);
void car_h_28(double *state, double *unused, double *out_1657620534686701055);
void car_H_28(double *state, double *unused, double *out_4118158368245429654);
void car_h_31(double *state, double *unused, double *out_21260442564898839);
void car_H_31(double *state, double *unused, double *out_1804281661385911717);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}