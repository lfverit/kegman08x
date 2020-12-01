/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_3473795810698563034);
void inv_err_fun(double *nom_x, double *true_x, double *out_7752315882522080299);
void H_mod_fun(double *state, double *out_5107623311201448163);
void f_fun(double *state, double dt, double *out_6775092890868093634);
void F_fun(double *state, double dt, double *out_8814329179518472778);
void h_25(double *state, double *unused, double *out_3877277127908247968);
void H_25(double *state, double *unused, double *out_6762120823040284824);
void h_24(double *state, double *unused, double *out_9022773902133125623);
void H_24(double *state, double *unused, double *out_4520271851034976936);
void h_30(double *state, double *unused, double *out_5327908520537131135);
void H_30(double *state, double *unused, double *out_4975304948914773313);
void h_26(double *state, double *unused, double *out_8117162285511574836);
void H_26(double *state, double *unused, double *out_5944829733989251112);
void h_27(double *state, double *unused, double *out_1189140390810094755);
void H_27(double *state, double *unused, double *out_6262886936751398625);
void h_29(double *state, double *unused, double *out_6271418373438793570);
void H_29(double *state, double *unused, double *out_1753074429408748441);
void h_28(double *state, double *unused, double *out_3201592284613370746);
void H_28(double *state, double *unused, double *out_6908639744438079814);
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
