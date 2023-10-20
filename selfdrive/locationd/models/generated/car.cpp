#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with sympy 1.9                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_628289647837459206) {
   out_628289647837459206[0] = delta_x[0] + nom_x[0];
   out_628289647837459206[1] = delta_x[1] + nom_x[1];
   out_628289647837459206[2] = delta_x[2] + nom_x[2];
   out_628289647837459206[3] = delta_x[3] + nom_x[3];
   out_628289647837459206[4] = delta_x[4] + nom_x[4];
   out_628289647837459206[5] = delta_x[5] + nom_x[5];
   out_628289647837459206[6] = delta_x[6] + nom_x[6];
   out_628289647837459206[7] = delta_x[7] + nom_x[7];
   out_628289647837459206[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_1927844397472153779) {
   out_1927844397472153779[0] = -nom_x[0] + true_x[0];
   out_1927844397472153779[1] = -nom_x[1] + true_x[1];
   out_1927844397472153779[2] = -nom_x[2] + true_x[2];
   out_1927844397472153779[3] = -nom_x[3] + true_x[3];
   out_1927844397472153779[4] = -nom_x[4] + true_x[4];
   out_1927844397472153779[5] = -nom_x[5] + true_x[5];
   out_1927844397472153779[6] = -nom_x[6] + true_x[6];
   out_1927844397472153779[7] = -nom_x[7] + true_x[7];
   out_1927844397472153779[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_8843017787680630764) {
   out_8843017787680630764[0] = 1.0;
   out_8843017787680630764[1] = 0;
   out_8843017787680630764[2] = 0;
   out_8843017787680630764[3] = 0;
   out_8843017787680630764[4] = 0;
   out_8843017787680630764[5] = 0;
   out_8843017787680630764[6] = 0;
   out_8843017787680630764[7] = 0;
   out_8843017787680630764[8] = 0;
   out_8843017787680630764[9] = 0;
   out_8843017787680630764[10] = 1.0;
   out_8843017787680630764[11] = 0;
   out_8843017787680630764[12] = 0;
   out_8843017787680630764[13] = 0;
   out_8843017787680630764[14] = 0;
   out_8843017787680630764[15] = 0;
   out_8843017787680630764[16] = 0;
   out_8843017787680630764[17] = 0;
   out_8843017787680630764[18] = 0;
   out_8843017787680630764[19] = 0;
   out_8843017787680630764[20] = 1.0;
   out_8843017787680630764[21] = 0;
   out_8843017787680630764[22] = 0;
   out_8843017787680630764[23] = 0;
   out_8843017787680630764[24] = 0;
   out_8843017787680630764[25] = 0;
   out_8843017787680630764[26] = 0;
   out_8843017787680630764[27] = 0;
   out_8843017787680630764[28] = 0;
   out_8843017787680630764[29] = 0;
   out_8843017787680630764[30] = 1.0;
   out_8843017787680630764[31] = 0;
   out_8843017787680630764[32] = 0;
   out_8843017787680630764[33] = 0;
   out_8843017787680630764[34] = 0;
   out_8843017787680630764[35] = 0;
   out_8843017787680630764[36] = 0;
   out_8843017787680630764[37] = 0;
   out_8843017787680630764[38] = 0;
   out_8843017787680630764[39] = 0;
   out_8843017787680630764[40] = 1.0;
   out_8843017787680630764[41] = 0;
   out_8843017787680630764[42] = 0;
   out_8843017787680630764[43] = 0;
   out_8843017787680630764[44] = 0;
   out_8843017787680630764[45] = 0;
   out_8843017787680630764[46] = 0;
   out_8843017787680630764[47] = 0;
   out_8843017787680630764[48] = 0;
   out_8843017787680630764[49] = 0;
   out_8843017787680630764[50] = 1.0;
   out_8843017787680630764[51] = 0;
   out_8843017787680630764[52] = 0;
   out_8843017787680630764[53] = 0;
   out_8843017787680630764[54] = 0;
   out_8843017787680630764[55] = 0;
   out_8843017787680630764[56] = 0;
   out_8843017787680630764[57] = 0;
   out_8843017787680630764[58] = 0;
   out_8843017787680630764[59] = 0;
   out_8843017787680630764[60] = 1.0;
   out_8843017787680630764[61] = 0;
   out_8843017787680630764[62] = 0;
   out_8843017787680630764[63] = 0;
   out_8843017787680630764[64] = 0;
   out_8843017787680630764[65] = 0;
   out_8843017787680630764[66] = 0;
   out_8843017787680630764[67] = 0;
   out_8843017787680630764[68] = 0;
   out_8843017787680630764[69] = 0;
   out_8843017787680630764[70] = 1.0;
   out_8843017787680630764[71] = 0;
   out_8843017787680630764[72] = 0;
   out_8843017787680630764[73] = 0;
   out_8843017787680630764[74] = 0;
   out_8843017787680630764[75] = 0;
   out_8843017787680630764[76] = 0;
   out_8843017787680630764[77] = 0;
   out_8843017787680630764[78] = 0;
   out_8843017787680630764[79] = 0;
   out_8843017787680630764[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_7261182143492256267) {
   out_7261182143492256267[0] = state[0];
   out_7261182143492256267[1] = state[1];
   out_7261182143492256267[2] = state[2];
   out_7261182143492256267[3] = state[3];
   out_7261182143492256267[4] = state[4];
   out_7261182143492256267[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_7261182143492256267[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_7261182143492256267[7] = state[7];
   out_7261182143492256267[8] = state[8];
}
void F_fun(double *state, double dt, double *out_8085554514994619377) {
   out_8085554514994619377[0] = 1;
   out_8085554514994619377[1] = 0;
   out_8085554514994619377[2] = 0;
   out_8085554514994619377[3] = 0;
   out_8085554514994619377[4] = 0;
   out_8085554514994619377[5] = 0;
   out_8085554514994619377[6] = 0;
   out_8085554514994619377[7] = 0;
   out_8085554514994619377[8] = 0;
   out_8085554514994619377[9] = 0;
   out_8085554514994619377[10] = 1;
   out_8085554514994619377[11] = 0;
   out_8085554514994619377[12] = 0;
   out_8085554514994619377[13] = 0;
   out_8085554514994619377[14] = 0;
   out_8085554514994619377[15] = 0;
   out_8085554514994619377[16] = 0;
   out_8085554514994619377[17] = 0;
   out_8085554514994619377[18] = 0;
   out_8085554514994619377[19] = 0;
   out_8085554514994619377[20] = 1;
   out_8085554514994619377[21] = 0;
   out_8085554514994619377[22] = 0;
   out_8085554514994619377[23] = 0;
   out_8085554514994619377[24] = 0;
   out_8085554514994619377[25] = 0;
   out_8085554514994619377[26] = 0;
   out_8085554514994619377[27] = 0;
   out_8085554514994619377[28] = 0;
   out_8085554514994619377[29] = 0;
   out_8085554514994619377[30] = 1;
   out_8085554514994619377[31] = 0;
   out_8085554514994619377[32] = 0;
   out_8085554514994619377[33] = 0;
   out_8085554514994619377[34] = 0;
   out_8085554514994619377[35] = 0;
   out_8085554514994619377[36] = 0;
   out_8085554514994619377[37] = 0;
   out_8085554514994619377[38] = 0;
   out_8085554514994619377[39] = 0;
   out_8085554514994619377[40] = 1;
   out_8085554514994619377[41] = 0;
   out_8085554514994619377[42] = 0;
   out_8085554514994619377[43] = 0;
   out_8085554514994619377[44] = 0;
   out_8085554514994619377[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_8085554514994619377[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_8085554514994619377[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8085554514994619377[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8085554514994619377[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_8085554514994619377[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_8085554514994619377[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_8085554514994619377[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_8085554514994619377[53] = -9.8000000000000007*dt;
   out_8085554514994619377[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_8085554514994619377[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_8085554514994619377[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8085554514994619377[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8085554514994619377[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_8085554514994619377[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_8085554514994619377[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_8085554514994619377[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8085554514994619377[62] = 0;
   out_8085554514994619377[63] = 0;
   out_8085554514994619377[64] = 0;
   out_8085554514994619377[65] = 0;
   out_8085554514994619377[66] = 0;
   out_8085554514994619377[67] = 0;
   out_8085554514994619377[68] = 0;
   out_8085554514994619377[69] = 0;
   out_8085554514994619377[70] = 1;
   out_8085554514994619377[71] = 0;
   out_8085554514994619377[72] = 0;
   out_8085554514994619377[73] = 0;
   out_8085554514994619377[74] = 0;
   out_8085554514994619377[75] = 0;
   out_8085554514994619377[76] = 0;
   out_8085554514994619377[77] = 0;
   out_8085554514994619377[78] = 0;
   out_8085554514994619377[79] = 0;
   out_8085554514994619377[80] = 1;
}
void h_25(double *state, double *unused, double *out_4515812170498523752) {
   out_4515812170498523752[0] = state[6];
}
void H_25(double *state, double *unused, double *out_1773635699508951289) {
   out_1773635699508951289[0] = 0;
   out_1773635699508951289[1] = 0;
   out_1773635699508951289[2] = 0;
   out_1773635699508951289[3] = 0;
   out_1773635699508951289[4] = 0;
   out_1773635699508951289[5] = 0;
   out_1773635699508951289[6] = 1;
   out_1773635699508951289[7] = 0;
   out_1773635699508951289[8] = 0;
}
void h_24(double *state, double *unused, double *out_790092840729285282) {
   out_790092840729285282[0] = state[4];
   out_790092840729285282[1] = state[5];
}
void H_24(double *state, double *unused, double *out_3046685805147036974) {
   out_3046685805147036974[0] = 0;
   out_3046685805147036974[1] = 0;
   out_3046685805147036974[2] = 0;
   out_3046685805147036974[3] = 0;
   out_3046685805147036974[4] = 1;
   out_3046685805147036974[5] = 0;
   out_3046685805147036974[6] = 0;
   out_3046685805147036974[7] = 0;
   out_3046685805147036974[8] = 0;
   out_3046685805147036974[9] = 0;
   out_3046685805147036974[10] = 0;
   out_3046685805147036974[11] = 0;
   out_3046685805147036974[12] = 0;
   out_3046685805147036974[13] = 0;
   out_3046685805147036974[14] = 1;
   out_3046685805147036974[15] = 0;
   out_3046685805147036974[16] = 0;
   out_3046685805147036974[17] = 0;
}
void h_30(double *state, double *unused, double *out_2255023055851827184) {
   out_2255023055851827184[0] = state[4];
}
void H_30(double *state, double *unused, double *out_1644296752365711219) {
   out_1644296752365711219[0] = 0;
   out_1644296752365711219[1] = 0;
   out_1644296752365711219[2] = 0;
   out_1644296752365711219[3] = 0;
   out_1644296752365711219[4] = 1;
   out_1644296752365711219[5] = 0;
   out_1644296752365711219[6] = 0;
   out_1644296752365711219[7] = 0;
   out_1644296752365711219[8] = 0;
}
void h_26(double *state, double *unused, double *out_7256696323512672829) {
   out_7256696323512672829[0] = state[7];
}
void H_26(double *state, double *unused, double *out_2430489763619263193) {
   out_2430489763619263193[0] = 0;
   out_2430489763619263193[1] = 0;
   out_2430489763619263193[2] = 0;
   out_2430489763619263193[3] = 0;
   out_2430489763619263193[4] = 0;
   out_2430489763619263193[5] = 0;
   out_2430489763619263193[6] = 0;
   out_2430489763619263193[7] = 1;
   out_2430489763619263193[8] = 0;
}
void h_27(double *state, double *unused, double *out_7138688929579296085) {
   out_7138688929579296085[0] = state[3];
}
void H_27(double *state, double *unused, double *out_530466559434713692) {
   out_530466559434713692[0] = 0;
   out_530466559434713692[1] = 0;
   out_530466559434713692[2] = 0;
   out_530466559434713692[3] = 1;
   out_530466559434713692[4] = 0;
   out_530466559434713692[5] = 0;
   out_530466559434713692[6] = 0;
   out_530466559434713692[7] = 0;
   out_530466559434713692[8] = 0;
}
void h_29(double *state, double *unused, double *out_6462813864016865142) {
   out_6462813864016865142[0] = state[1];
}
void H_29(double *state, double *unused, double *out_2154528096680103403) {
   out_2154528096680103403[0] = 0;
   out_2154528096680103403[1] = 1;
   out_2154528096680103403[2] = 0;
   out_2154528096680103403[3] = 0;
   out_2154528096680103403[4] = 0;
   out_2154528096680103403[5] = 0;
   out_2154528096680103403[6] = 0;
   out_2154528096680103403[7] = 0;
   out_2154528096680103403[8] = 0;
}
void h_28(double *state, double *unused, double *out_1657620534686701055) {
   out_1657620534686701055[0] = state[0];
}
void H_28(double *state, double *unused, double *out_4118158368245429654) {
   out_4118158368245429654[0] = 1;
   out_4118158368245429654[1] = 0;
   out_4118158368245429654[2] = 0;
   out_4118158368245429654[3] = 0;
   out_4118158368245429654[4] = 0;
   out_4118158368245429654[5] = 0;
   out_4118158368245429654[6] = 0;
   out_4118158368245429654[7] = 0;
   out_4118158368245429654[8] = 0;
}
void h_31(double *state, double *unused, double *out_21260442564898839) {
   out_21260442564898839[0] = state[8];
}
void H_31(double *state, double *unused, double *out_1804281661385911717) {
   out_1804281661385911717[0] = 0;
   out_1804281661385911717[1] = 0;
   out_1804281661385911717[2] = 0;
   out_1804281661385911717[3] = 0;
   out_1804281661385911717[4] = 0;
   out_1804281661385911717[5] = 0;
   out_1804281661385911717[6] = 0;
   out_1804281661385911717[7] = 0;
   out_1804281661385911717[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_628289647837459206) {
  err_fun(nom_x, delta_x, out_628289647837459206);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_1927844397472153779) {
  inv_err_fun(nom_x, true_x, out_1927844397472153779);
}
void car_H_mod_fun(double *state, double *out_8843017787680630764) {
  H_mod_fun(state, out_8843017787680630764);
}
void car_f_fun(double *state, double dt, double *out_7261182143492256267) {
  f_fun(state,  dt, out_7261182143492256267);
}
void car_F_fun(double *state, double dt, double *out_8085554514994619377) {
  F_fun(state,  dt, out_8085554514994619377);
}
void car_h_25(double *state, double *unused, double *out_4515812170498523752) {
  h_25(state, unused, out_4515812170498523752);
}
void car_H_25(double *state, double *unused, double *out_1773635699508951289) {
  H_25(state, unused, out_1773635699508951289);
}
void car_h_24(double *state, double *unused, double *out_790092840729285282) {
  h_24(state, unused, out_790092840729285282);
}
void car_H_24(double *state, double *unused, double *out_3046685805147036974) {
  H_24(state, unused, out_3046685805147036974);
}
void car_h_30(double *state, double *unused, double *out_2255023055851827184) {
  h_30(state, unused, out_2255023055851827184);
}
void car_H_30(double *state, double *unused, double *out_1644296752365711219) {
  H_30(state, unused, out_1644296752365711219);
}
void car_h_26(double *state, double *unused, double *out_7256696323512672829) {
  h_26(state, unused, out_7256696323512672829);
}
void car_H_26(double *state, double *unused, double *out_2430489763619263193) {
  H_26(state, unused, out_2430489763619263193);
}
void car_h_27(double *state, double *unused, double *out_7138688929579296085) {
  h_27(state, unused, out_7138688929579296085);
}
void car_H_27(double *state, double *unused, double *out_530466559434713692) {
  H_27(state, unused, out_530466559434713692);
}
void car_h_29(double *state, double *unused, double *out_6462813864016865142) {
  h_29(state, unused, out_6462813864016865142);
}
void car_H_29(double *state, double *unused, double *out_2154528096680103403) {
  H_29(state, unused, out_2154528096680103403);
}
void car_h_28(double *state, double *unused, double *out_1657620534686701055) {
  h_28(state, unused, out_1657620534686701055);
}
void car_H_28(double *state, double *unused, double *out_4118158368245429654) {
  H_28(state, unused, out_4118158368245429654);
}
void car_h_31(double *state, double *unused, double *out_21260442564898839) {
  h_31(state, unused, out_21260442564898839);
}
void car_H_31(double *state, double *unused, double *out_1804281661385911717) {
  H_31(state, unused, out_1804281661385911717);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);
