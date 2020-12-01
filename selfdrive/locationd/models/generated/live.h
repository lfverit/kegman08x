/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_4874561107718521871);
void inv_err_fun(double *nom_x, double *true_x, double *out_2791063465481048667);
void H_mod_fun(double *state, double *out_4584880498007367256);
void f_fun(double *state, double dt, double *out_9112291710325674070);
void F_fun(double *state, double dt, double *out_959455092875613869);
void h_3(double *state, double *unused, double *out_6819425761516477051);
void H_3(double *state, double *unused, double *out_2295073677443666431);
void h_4(double *state, double *unused, double *out_4943826168677602330);
void H_4(double *state, double *unused, double *out_1930894265245222956);
void h_9(double *state, double *unused, double *out_7863544349885867869);
void H_9(double *state, double *unused, double *out_4789968661610001142);
void h_10(double *state, double *unused, double *out_5808934889666736191);
void H_10(double *state, double *unused, double *out_8502709149966356383);
void h_12(double *state, double *unused, double *out_1367679388951826364);
void H_12(double *state, double *unused, double *out_4574018629478354741);
void h_31(double *state, double *unused, double *out_6709474561669204312);
void H_31(double *state, double *unused, double *out_1898928854334934920);
void h_32(double *state, double *unused, double *out_3931352398272551390);
void H_32(double *state, double *unused, double *out_4401436889925466266);
void h_13(double *state, double *unused, double *out_3540073942281531726);
void H_13(double *state, double *unused, double *out_4693257929497467162);
void h_14(double *state, double *unused, double *out_7863544349885867869);
void H_14(double *state, double *unused, double *out_4789968661610001142);
void h_19(double *state, double *unused, double *out_5671454769117231308);
void H_19(double *state, double *unused, double *out_9222648512712432323);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);