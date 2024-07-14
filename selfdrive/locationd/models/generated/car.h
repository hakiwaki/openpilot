#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_1244754007084266224);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_5303699367332946512);
void car_H_mod_fun(double *state, double *out_182754799911660332);
void car_f_fun(double *state, double dt, double *out_1416607489860240731);
void car_F_fun(double *state, double dt, double *out_782414791261192023);
void car_h_25(double *state, double *unused, double *out_7237448606815643908);
void car_H_25(double *state, double *unused, double *out_1406740492386300976);
void car_h_24(double *state, double *unused, double *out_8901795115737982182);
void car_H_24(double *state, double *unused, double *out_3579390091391800542);
void car_h_30(double *state, double *unused, double *out_4438460240543556772);
void car_H_30(double *state, double *unused, double *out_1536079439529541046);
void car_h_26(double *state, double *unused, double *out_6626747883244547030);
void car_H_26(double *state, double *unused, double *out_5148243811260357200);
void car_h_27(double *state, double *unused, double *out_2235775339154798263);
void car_H_27(double *state, double *unused, double *out_3710842751329965957);
void car_h_29(double *state, double *unused, double *out_7761832775705821376);
void car_H_29(double *state, double *unused, double *out_5424205478199516990);
void car_h_28(double *state, double *unused, double *out_3406954531926158572);
void car_H_28(double *state, double *unused, double *out_3460575206634190739);
void car_h_31(double *state, double *unused, double *out_5145760610766649551);
void car_H_31(double *state, double *unused, double *out_1376094530509340548);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}