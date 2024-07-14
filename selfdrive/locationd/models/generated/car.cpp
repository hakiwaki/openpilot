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
 *                       Code generated with SymPy 1.12                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_1244754007084266224) {
   out_1244754007084266224[0] = delta_x[0] + nom_x[0];
   out_1244754007084266224[1] = delta_x[1] + nom_x[1];
   out_1244754007084266224[2] = delta_x[2] + nom_x[2];
   out_1244754007084266224[3] = delta_x[3] + nom_x[3];
   out_1244754007084266224[4] = delta_x[4] + nom_x[4];
   out_1244754007084266224[5] = delta_x[5] + nom_x[5];
   out_1244754007084266224[6] = delta_x[6] + nom_x[6];
   out_1244754007084266224[7] = delta_x[7] + nom_x[7];
   out_1244754007084266224[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_5303699367332946512) {
   out_5303699367332946512[0] = -nom_x[0] + true_x[0];
   out_5303699367332946512[1] = -nom_x[1] + true_x[1];
   out_5303699367332946512[2] = -nom_x[2] + true_x[2];
   out_5303699367332946512[3] = -nom_x[3] + true_x[3];
   out_5303699367332946512[4] = -nom_x[4] + true_x[4];
   out_5303699367332946512[5] = -nom_x[5] + true_x[5];
   out_5303699367332946512[6] = -nom_x[6] + true_x[6];
   out_5303699367332946512[7] = -nom_x[7] + true_x[7];
   out_5303699367332946512[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_182754799911660332) {
   out_182754799911660332[0] = 1.0;
   out_182754799911660332[1] = 0;
   out_182754799911660332[2] = 0;
   out_182754799911660332[3] = 0;
   out_182754799911660332[4] = 0;
   out_182754799911660332[5] = 0;
   out_182754799911660332[6] = 0;
   out_182754799911660332[7] = 0;
   out_182754799911660332[8] = 0;
   out_182754799911660332[9] = 0;
   out_182754799911660332[10] = 1.0;
   out_182754799911660332[11] = 0;
   out_182754799911660332[12] = 0;
   out_182754799911660332[13] = 0;
   out_182754799911660332[14] = 0;
   out_182754799911660332[15] = 0;
   out_182754799911660332[16] = 0;
   out_182754799911660332[17] = 0;
   out_182754799911660332[18] = 0;
   out_182754799911660332[19] = 0;
   out_182754799911660332[20] = 1.0;
   out_182754799911660332[21] = 0;
   out_182754799911660332[22] = 0;
   out_182754799911660332[23] = 0;
   out_182754799911660332[24] = 0;
   out_182754799911660332[25] = 0;
   out_182754799911660332[26] = 0;
   out_182754799911660332[27] = 0;
   out_182754799911660332[28] = 0;
   out_182754799911660332[29] = 0;
   out_182754799911660332[30] = 1.0;
   out_182754799911660332[31] = 0;
   out_182754799911660332[32] = 0;
   out_182754799911660332[33] = 0;
   out_182754799911660332[34] = 0;
   out_182754799911660332[35] = 0;
   out_182754799911660332[36] = 0;
   out_182754799911660332[37] = 0;
   out_182754799911660332[38] = 0;
   out_182754799911660332[39] = 0;
   out_182754799911660332[40] = 1.0;
   out_182754799911660332[41] = 0;
   out_182754799911660332[42] = 0;
   out_182754799911660332[43] = 0;
   out_182754799911660332[44] = 0;
   out_182754799911660332[45] = 0;
   out_182754799911660332[46] = 0;
   out_182754799911660332[47] = 0;
   out_182754799911660332[48] = 0;
   out_182754799911660332[49] = 0;
   out_182754799911660332[50] = 1.0;
   out_182754799911660332[51] = 0;
   out_182754799911660332[52] = 0;
   out_182754799911660332[53] = 0;
   out_182754799911660332[54] = 0;
   out_182754799911660332[55] = 0;
   out_182754799911660332[56] = 0;
   out_182754799911660332[57] = 0;
   out_182754799911660332[58] = 0;
   out_182754799911660332[59] = 0;
   out_182754799911660332[60] = 1.0;
   out_182754799911660332[61] = 0;
   out_182754799911660332[62] = 0;
   out_182754799911660332[63] = 0;
   out_182754799911660332[64] = 0;
   out_182754799911660332[65] = 0;
   out_182754799911660332[66] = 0;
   out_182754799911660332[67] = 0;
   out_182754799911660332[68] = 0;
   out_182754799911660332[69] = 0;
   out_182754799911660332[70] = 1.0;
   out_182754799911660332[71] = 0;
   out_182754799911660332[72] = 0;
   out_182754799911660332[73] = 0;
   out_182754799911660332[74] = 0;
   out_182754799911660332[75] = 0;
   out_182754799911660332[76] = 0;
   out_182754799911660332[77] = 0;
   out_182754799911660332[78] = 0;
   out_182754799911660332[79] = 0;
   out_182754799911660332[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_1416607489860240731) {
   out_1416607489860240731[0] = state[0];
   out_1416607489860240731[1] = state[1];
   out_1416607489860240731[2] = state[2];
   out_1416607489860240731[3] = state[3];
   out_1416607489860240731[4] = state[4];
   out_1416607489860240731[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_1416607489860240731[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_1416607489860240731[7] = state[7];
   out_1416607489860240731[8] = state[8];
}
void F_fun(double *state, double dt, double *out_782414791261192023) {
   out_782414791261192023[0] = 1;
   out_782414791261192023[1] = 0;
   out_782414791261192023[2] = 0;
   out_782414791261192023[3] = 0;
   out_782414791261192023[4] = 0;
   out_782414791261192023[5] = 0;
   out_782414791261192023[6] = 0;
   out_782414791261192023[7] = 0;
   out_782414791261192023[8] = 0;
   out_782414791261192023[9] = 0;
   out_782414791261192023[10] = 1;
   out_782414791261192023[11] = 0;
   out_782414791261192023[12] = 0;
   out_782414791261192023[13] = 0;
   out_782414791261192023[14] = 0;
   out_782414791261192023[15] = 0;
   out_782414791261192023[16] = 0;
   out_782414791261192023[17] = 0;
   out_782414791261192023[18] = 0;
   out_782414791261192023[19] = 0;
   out_782414791261192023[20] = 1;
   out_782414791261192023[21] = 0;
   out_782414791261192023[22] = 0;
   out_782414791261192023[23] = 0;
   out_782414791261192023[24] = 0;
   out_782414791261192023[25] = 0;
   out_782414791261192023[26] = 0;
   out_782414791261192023[27] = 0;
   out_782414791261192023[28] = 0;
   out_782414791261192023[29] = 0;
   out_782414791261192023[30] = 1;
   out_782414791261192023[31] = 0;
   out_782414791261192023[32] = 0;
   out_782414791261192023[33] = 0;
   out_782414791261192023[34] = 0;
   out_782414791261192023[35] = 0;
   out_782414791261192023[36] = 0;
   out_782414791261192023[37] = 0;
   out_782414791261192023[38] = 0;
   out_782414791261192023[39] = 0;
   out_782414791261192023[40] = 1;
   out_782414791261192023[41] = 0;
   out_782414791261192023[42] = 0;
   out_782414791261192023[43] = 0;
   out_782414791261192023[44] = 0;
   out_782414791261192023[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_782414791261192023[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_782414791261192023[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_782414791261192023[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_782414791261192023[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_782414791261192023[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_782414791261192023[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_782414791261192023[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_782414791261192023[53] = -9.8000000000000007*dt;
   out_782414791261192023[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_782414791261192023[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_782414791261192023[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_782414791261192023[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_782414791261192023[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_782414791261192023[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_782414791261192023[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_782414791261192023[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_782414791261192023[62] = 0;
   out_782414791261192023[63] = 0;
   out_782414791261192023[64] = 0;
   out_782414791261192023[65] = 0;
   out_782414791261192023[66] = 0;
   out_782414791261192023[67] = 0;
   out_782414791261192023[68] = 0;
   out_782414791261192023[69] = 0;
   out_782414791261192023[70] = 1;
   out_782414791261192023[71] = 0;
   out_782414791261192023[72] = 0;
   out_782414791261192023[73] = 0;
   out_782414791261192023[74] = 0;
   out_782414791261192023[75] = 0;
   out_782414791261192023[76] = 0;
   out_782414791261192023[77] = 0;
   out_782414791261192023[78] = 0;
   out_782414791261192023[79] = 0;
   out_782414791261192023[80] = 1;
}
void h_25(double *state, double *unused, double *out_7237448606815643908) {
   out_7237448606815643908[0] = state[6];
}
void H_25(double *state, double *unused, double *out_1406740492386300976) {
   out_1406740492386300976[0] = 0;
   out_1406740492386300976[1] = 0;
   out_1406740492386300976[2] = 0;
   out_1406740492386300976[3] = 0;
   out_1406740492386300976[4] = 0;
   out_1406740492386300976[5] = 0;
   out_1406740492386300976[6] = 1;
   out_1406740492386300976[7] = 0;
   out_1406740492386300976[8] = 0;
}
void h_24(double *state, double *unused, double *out_8901795115737982182) {
   out_8901795115737982182[0] = state[4];
   out_8901795115737982182[1] = state[5];
}
void H_24(double *state, double *unused, double *out_3579390091391800542) {
   out_3579390091391800542[0] = 0;
   out_3579390091391800542[1] = 0;
   out_3579390091391800542[2] = 0;
   out_3579390091391800542[3] = 0;
   out_3579390091391800542[4] = 1;
   out_3579390091391800542[5] = 0;
   out_3579390091391800542[6] = 0;
   out_3579390091391800542[7] = 0;
   out_3579390091391800542[8] = 0;
   out_3579390091391800542[9] = 0;
   out_3579390091391800542[10] = 0;
   out_3579390091391800542[11] = 0;
   out_3579390091391800542[12] = 0;
   out_3579390091391800542[13] = 0;
   out_3579390091391800542[14] = 1;
   out_3579390091391800542[15] = 0;
   out_3579390091391800542[16] = 0;
   out_3579390091391800542[17] = 0;
}
void h_30(double *state, double *unused, double *out_4438460240543556772) {
   out_4438460240543556772[0] = state[4];
}
void H_30(double *state, double *unused, double *out_1536079439529541046) {
   out_1536079439529541046[0] = 0;
   out_1536079439529541046[1] = 0;
   out_1536079439529541046[2] = 0;
   out_1536079439529541046[3] = 0;
   out_1536079439529541046[4] = 1;
   out_1536079439529541046[5] = 0;
   out_1536079439529541046[6] = 0;
   out_1536079439529541046[7] = 0;
   out_1536079439529541046[8] = 0;
}
void h_26(double *state, double *unused, double *out_6626747883244547030) {
   out_6626747883244547030[0] = state[7];
}
void H_26(double *state, double *unused, double *out_5148243811260357200) {
   out_5148243811260357200[0] = 0;
   out_5148243811260357200[1] = 0;
   out_5148243811260357200[2] = 0;
   out_5148243811260357200[3] = 0;
   out_5148243811260357200[4] = 0;
   out_5148243811260357200[5] = 0;
   out_5148243811260357200[6] = 0;
   out_5148243811260357200[7] = 1;
   out_5148243811260357200[8] = 0;
}
void h_27(double *state, double *unused, double *out_2235775339154798263) {
   out_2235775339154798263[0] = state[3];
}
void H_27(double *state, double *unused, double *out_3710842751329965957) {
   out_3710842751329965957[0] = 0;
   out_3710842751329965957[1] = 0;
   out_3710842751329965957[2] = 0;
   out_3710842751329965957[3] = 1;
   out_3710842751329965957[4] = 0;
   out_3710842751329965957[5] = 0;
   out_3710842751329965957[6] = 0;
   out_3710842751329965957[7] = 0;
   out_3710842751329965957[8] = 0;
}
void h_29(double *state, double *unused, double *out_7761832775705821376) {
   out_7761832775705821376[0] = state[1];
}
void H_29(double *state, double *unused, double *out_5424205478199516990) {
   out_5424205478199516990[0] = 0;
   out_5424205478199516990[1] = 1;
   out_5424205478199516990[2] = 0;
   out_5424205478199516990[3] = 0;
   out_5424205478199516990[4] = 0;
   out_5424205478199516990[5] = 0;
   out_5424205478199516990[6] = 0;
   out_5424205478199516990[7] = 0;
   out_5424205478199516990[8] = 0;
}
void h_28(double *state, double *unused, double *out_3406954531926158572) {
   out_3406954531926158572[0] = state[0];
}
void H_28(double *state, double *unused, double *out_3460575206634190739) {
   out_3460575206634190739[0] = 1;
   out_3460575206634190739[1] = 0;
   out_3460575206634190739[2] = 0;
   out_3460575206634190739[3] = 0;
   out_3460575206634190739[4] = 0;
   out_3460575206634190739[5] = 0;
   out_3460575206634190739[6] = 0;
   out_3460575206634190739[7] = 0;
   out_3460575206634190739[8] = 0;
}
void h_31(double *state, double *unused, double *out_5145760610766649551) {
   out_5145760610766649551[0] = state[8];
}
void H_31(double *state, double *unused, double *out_1376094530509340548) {
   out_1376094530509340548[0] = 0;
   out_1376094530509340548[1] = 0;
   out_1376094530509340548[2] = 0;
   out_1376094530509340548[3] = 0;
   out_1376094530509340548[4] = 0;
   out_1376094530509340548[5] = 0;
   out_1376094530509340548[6] = 0;
   out_1376094530509340548[7] = 0;
   out_1376094530509340548[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_1244754007084266224) {
  err_fun(nom_x, delta_x, out_1244754007084266224);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_5303699367332946512) {
  inv_err_fun(nom_x, true_x, out_5303699367332946512);
}
void car_H_mod_fun(double *state, double *out_182754799911660332) {
  H_mod_fun(state, out_182754799911660332);
}
void car_f_fun(double *state, double dt, double *out_1416607489860240731) {
  f_fun(state,  dt, out_1416607489860240731);
}
void car_F_fun(double *state, double dt, double *out_782414791261192023) {
  F_fun(state,  dt, out_782414791261192023);
}
void car_h_25(double *state, double *unused, double *out_7237448606815643908) {
  h_25(state, unused, out_7237448606815643908);
}
void car_H_25(double *state, double *unused, double *out_1406740492386300976) {
  H_25(state, unused, out_1406740492386300976);
}
void car_h_24(double *state, double *unused, double *out_8901795115737982182) {
  h_24(state, unused, out_8901795115737982182);
}
void car_H_24(double *state, double *unused, double *out_3579390091391800542) {
  H_24(state, unused, out_3579390091391800542);
}
void car_h_30(double *state, double *unused, double *out_4438460240543556772) {
  h_30(state, unused, out_4438460240543556772);
}
void car_H_30(double *state, double *unused, double *out_1536079439529541046) {
  H_30(state, unused, out_1536079439529541046);
}
void car_h_26(double *state, double *unused, double *out_6626747883244547030) {
  h_26(state, unused, out_6626747883244547030);
}
void car_H_26(double *state, double *unused, double *out_5148243811260357200) {
  H_26(state, unused, out_5148243811260357200);
}
void car_h_27(double *state, double *unused, double *out_2235775339154798263) {
  h_27(state, unused, out_2235775339154798263);
}
void car_H_27(double *state, double *unused, double *out_3710842751329965957) {
  H_27(state, unused, out_3710842751329965957);
}
void car_h_29(double *state, double *unused, double *out_7761832775705821376) {
  h_29(state, unused, out_7761832775705821376);
}
void car_H_29(double *state, double *unused, double *out_5424205478199516990) {
  H_29(state, unused, out_5424205478199516990);
}
void car_h_28(double *state, double *unused, double *out_3406954531926158572) {
  h_28(state, unused, out_3406954531926158572);
}
void car_H_28(double *state, double *unused, double *out_3460575206634190739) {
  H_28(state, unused, out_3460575206634190739);
}
void car_h_31(double *state, double *unused, double *out_5145760610766649551) {
  h_31(state, unused, out_5145760610766649551);
}
void car_H_31(double *state, double *unused, double *out_1376094530509340548) {
  H_31(state, unused, out_1376094530509340548);
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

ekf_lib_init(car)
