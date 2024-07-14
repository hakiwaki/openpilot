#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_35(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_688119966555173306);
void live_err_fun(double *nom_x, double *delta_x, double *out_5924447345377752939);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_4459438890970929893);
void live_H_mod_fun(double *state, double *out_1714341149114002571);
void live_f_fun(double *state, double dt, double *out_3846009461295015815);
void live_F_fun(double *state, double dt, double *out_4941625080249828538);
void live_h_4(double *state, double *unused, double *out_3729668420870475872);
void live_H_4(double *state, double *unused, double *out_4876946195465417178);
void live_h_9(double *state, double *unused, double *out_1151493546623340019);
void live_H_9(double *state, double *unused, double *out_2410272739799030292);
void live_h_10(double *state, double *unused, double *out_6167467416513163071);
void live_H_10(double *state, double *unused, double *out_6412941887746337931);
void live_h_12(double *state, double *unused, double *out_1810615206417119516);
void live_H_12(double *state, double *unused, double *out_142510212566544617);
void live_h_35(double *state, double *unused, double *out_2945186851113316000);
void live_H_35(double *state, double *unused, double *out_2888073244891558326);
void live_h_32(double *state, double *unused, double *out_7977252132068576231);
void live_H_32(double *state, double *unused, double *out_137477440014514439);
void live_h_13(double *state, double *unused, double *out_7260199594589518418);
void live_H_13(double *state, double *unused, double *out_4332738930436728016);
void live_h_14(double *state, double *unused, double *out_1151493546623340019);
void live_H_14(double *state, double *unused, double *out_2410272739799030292);
void live_h_33(double *state, double *unused, double *out_655295692279212237);
void live_H_33(double *state, double *unused, double *out_6038630249530415930);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}