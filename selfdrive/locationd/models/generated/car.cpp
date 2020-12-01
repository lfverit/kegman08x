
extern "C"{

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

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_3473795810698563034) {
   out_3473795810698563034[0] = delta_x[0] + nom_x[0];
   out_3473795810698563034[1] = delta_x[1] + nom_x[1];
   out_3473795810698563034[2] = delta_x[2] + nom_x[2];
   out_3473795810698563034[3] = delta_x[3] + nom_x[3];
   out_3473795810698563034[4] = delta_x[4] + nom_x[4];
   out_3473795810698563034[5] = delta_x[5] + nom_x[5];
   out_3473795810698563034[6] = delta_x[6] + nom_x[6];
   out_3473795810698563034[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_7752315882522080299) {
   out_7752315882522080299[0] = -nom_x[0] + true_x[0];
   out_7752315882522080299[1] = -nom_x[1] + true_x[1];
   out_7752315882522080299[2] = -nom_x[2] + true_x[2];
   out_7752315882522080299[3] = -nom_x[3] + true_x[3];
   out_7752315882522080299[4] = -nom_x[4] + true_x[4];
   out_7752315882522080299[5] = -nom_x[5] + true_x[5];
   out_7752315882522080299[6] = -nom_x[6] + true_x[6];
   out_7752315882522080299[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_5107623311201448163) {
   out_5107623311201448163[0] = 1.0;
   out_5107623311201448163[1] = 0.0;
   out_5107623311201448163[2] = 0.0;
   out_5107623311201448163[3] = 0.0;
   out_5107623311201448163[4] = 0.0;
   out_5107623311201448163[5] = 0.0;
   out_5107623311201448163[6] = 0.0;
   out_5107623311201448163[7] = 0.0;
   out_5107623311201448163[8] = 0.0;
   out_5107623311201448163[9] = 1.0;
   out_5107623311201448163[10] = 0.0;
   out_5107623311201448163[11] = 0.0;
   out_5107623311201448163[12] = 0.0;
   out_5107623311201448163[13] = 0.0;
   out_5107623311201448163[14] = 0.0;
   out_5107623311201448163[15] = 0.0;
   out_5107623311201448163[16] = 0.0;
   out_5107623311201448163[17] = 0.0;
   out_5107623311201448163[18] = 1.0;
   out_5107623311201448163[19] = 0.0;
   out_5107623311201448163[20] = 0.0;
   out_5107623311201448163[21] = 0.0;
   out_5107623311201448163[22] = 0.0;
   out_5107623311201448163[23] = 0.0;
   out_5107623311201448163[24] = 0.0;
   out_5107623311201448163[25] = 0.0;
   out_5107623311201448163[26] = 0.0;
   out_5107623311201448163[27] = 1.0;
   out_5107623311201448163[28] = 0.0;
   out_5107623311201448163[29] = 0.0;
   out_5107623311201448163[30] = 0.0;
   out_5107623311201448163[31] = 0.0;
   out_5107623311201448163[32] = 0.0;
   out_5107623311201448163[33] = 0.0;
   out_5107623311201448163[34] = 0.0;
   out_5107623311201448163[35] = 0.0;
   out_5107623311201448163[36] = 1.0;
   out_5107623311201448163[37] = 0.0;
   out_5107623311201448163[38] = 0.0;
   out_5107623311201448163[39] = 0.0;
   out_5107623311201448163[40] = 0.0;
   out_5107623311201448163[41] = 0.0;
   out_5107623311201448163[42] = 0.0;
   out_5107623311201448163[43] = 0.0;
   out_5107623311201448163[44] = 0.0;
   out_5107623311201448163[45] = 1.0;
   out_5107623311201448163[46] = 0.0;
   out_5107623311201448163[47] = 0.0;
   out_5107623311201448163[48] = 0.0;
   out_5107623311201448163[49] = 0.0;
   out_5107623311201448163[50] = 0.0;
   out_5107623311201448163[51] = 0.0;
   out_5107623311201448163[52] = 0.0;
   out_5107623311201448163[53] = 0.0;
   out_5107623311201448163[54] = 1.0;
   out_5107623311201448163[55] = 0.0;
   out_5107623311201448163[56] = 0.0;
   out_5107623311201448163[57] = 0.0;
   out_5107623311201448163[58] = 0.0;
   out_5107623311201448163[59] = 0.0;
   out_5107623311201448163[60] = 0.0;
   out_5107623311201448163[61] = 0.0;
   out_5107623311201448163[62] = 0.0;
   out_5107623311201448163[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_6775092890868093634) {
   out_6775092890868093634[0] = state[0];
   out_6775092890868093634[1] = state[1];
   out_6775092890868093634[2] = state[2];
   out_6775092890868093634[3] = state[3];
   out_6775092890868093634[4] = state[4];
   out_6775092890868093634[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_6775092890868093634[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_6775092890868093634[7] = state[7];
}
void F_fun(double *state, double dt, double *out_8814329179518472778) {
   out_8814329179518472778[0] = 1;
   out_8814329179518472778[1] = 0;
   out_8814329179518472778[2] = 0;
   out_8814329179518472778[3] = 0;
   out_8814329179518472778[4] = 0;
   out_8814329179518472778[5] = 0;
   out_8814329179518472778[6] = 0;
   out_8814329179518472778[7] = 0;
   out_8814329179518472778[8] = 0;
   out_8814329179518472778[9] = 1;
   out_8814329179518472778[10] = 0;
   out_8814329179518472778[11] = 0;
   out_8814329179518472778[12] = 0;
   out_8814329179518472778[13] = 0;
   out_8814329179518472778[14] = 0;
   out_8814329179518472778[15] = 0;
   out_8814329179518472778[16] = 0;
   out_8814329179518472778[17] = 0;
   out_8814329179518472778[18] = 1;
   out_8814329179518472778[19] = 0;
   out_8814329179518472778[20] = 0;
   out_8814329179518472778[21] = 0;
   out_8814329179518472778[22] = 0;
   out_8814329179518472778[23] = 0;
   out_8814329179518472778[24] = 0;
   out_8814329179518472778[25] = 0;
   out_8814329179518472778[26] = 0;
   out_8814329179518472778[27] = 1;
   out_8814329179518472778[28] = 0;
   out_8814329179518472778[29] = 0;
   out_8814329179518472778[30] = 0;
   out_8814329179518472778[31] = 0;
   out_8814329179518472778[32] = 0;
   out_8814329179518472778[33] = 0;
   out_8814329179518472778[34] = 0;
   out_8814329179518472778[35] = 0;
   out_8814329179518472778[36] = 1;
   out_8814329179518472778[37] = 0;
   out_8814329179518472778[38] = 0;
   out_8814329179518472778[39] = 0;
   out_8814329179518472778[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_8814329179518472778[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_8814329179518472778[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8814329179518472778[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8814329179518472778[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_8814329179518472778[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_8814329179518472778[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_8814329179518472778[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_8814329179518472778[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_8814329179518472778[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_8814329179518472778[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8814329179518472778[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8814329179518472778[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_8814329179518472778[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_8814329179518472778[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_8814329179518472778[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8814329179518472778[56] = 0;
   out_8814329179518472778[57] = 0;
   out_8814329179518472778[58] = 0;
   out_8814329179518472778[59] = 0;
   out_8814329179518472778[60] = 0;
   out_8814329179518472778[61] = 0;
   out_8814329179518472778[62] = 0;
   out_8814329179518472778[63] = 1;
}
void h_25(double *state, double *unused, double *out_3877277127908247968) {
   out_3877277127908247968[0] = state[6];
}
void H_25(double *state, double *unused, double *out_6762120823040284824) {
   out_6762120823040284824[0] = 0;
   out_6762120823040284824[1] = 0;
   out_6762120823040284824[2] = 0;
   out_6762120823040284824[3] = 0;
   out_6762120823040284824[4] = 0;
   out_6762120823040284824[5] = 0;
   out_6762120823040284824[6] = 1;
   out_6762120823040284824[7] = 0;
}
void h_24(double *state, double *unused, double *out_9022773902133125623) {
   out_9022773902133125623[0] = state[4];
   out_9022773902133125623[1] = state[5];
}
void H_24(double *state, double *unused, double *out_4520271851034976936) {
   out_4520271851034976936[0] = 0;
   out_4520271851034976936[1] = 0;
   out_4520271851034976936[2] = 0;
   out_4520271851034976936[3] = 0;
   out_4520271851034976936[4] = 1;
   out_4520271851034976936[5] = 0;
   out_4520271851034976936[6] = 0;
   out_4520271851034976936[7] = 0;
   out_4520271851034976936[8] = 0;
   out_4520271851034976936[9] = 0;
   out_4520271851034976936[10] = 0;
   out_4520271851034976936[11] = 0;
   out_4520271851034976936[12] = 0;
   out_4520271851034976936[13] = 1;
   out_4520271851034976936[14] = 0;
   out_4520271851034976936[15] = 0;
}
void h_30(double *state, double *unused, double *out_5327908520537131135) {
   out_5327908520537131135[0] = state[4];
}
void H_30(double *state, double *unused, double *out_4975304948914773313) {
   out_4975304948914773313[0] = 0;
   out_4975304948914773313[1] = 0;
   out_4975304948914773313[2] = 0;
   out_4975304948914773313[3] = 0;
   out_4975304948914773313[4] = 1;
   out_4975304948914773313[5] = 0;
   out_4975304948914773313[6] = 0;
   out_4975304948914773313[7] = 0;
}
void h_26(double *state, double *unused, double *out_8117162285511574836) {
   out_8117162285511574836[0] = state[7];
}
void H_26(double *state, double *unused, double *out_5944829733989251112) {
   out_5944829733989251112[0] = 0;
   out_5944829733989251112[1] = 0;
   out_5944829733989251112[2] = 0;
   out_5944829733989251112[3] = 0;
   out_5944829733989251112[4] = 0;
   out_5944829733989251112[5] = 0;
   out_5944829733989251112[6] = 0;
   out_5944829733989251112[7] = 1;
}
void h_27(double *state, double *unused, double *out_1189140390810094755) {
   out_1189140390810094755[0] = state[3];
}
void H_27(double *state, double *unused, double *out_6262886936751398625) {
   out_6262886936751398625[0] = 0;
   out_6262886936751398625[1] = 0;
   out_6262886936751398625[2] = 0;
   out_6262886936751398625[3] = 1;
   out_6262886936751398625[4] = 0;
   out_6262886936751398625[5] = 0;
   out_6262886936751398625[6] = 0;
   out_6262886936751398625[7] = 0;
}
void h_29(double *state, double *unused, double *out_6271418373438793570) {
   out_6271418373438793570[0] = state[1];
}
void H_29(double *state, double *unused, double *out_1753074429408748441) {
   out_1753074429408748441[0] = 0;
   out_1753074429408748441[1] = 1;
   out_1753074429408748441[2] = 0;
   out_1753074429408748441[3] = 0;
   out_1753074429408748441[4] = 0;
   out_1753074429408748441[5] = 0;
   out_1753074429408748441[6] = 0;
   out_1753074429408748441[7] = 0;
}
void h_28(double *state, double *unused, double *out_3201592284613370746) {
   out_3201592284613370746[0] = state[5];
   out_3201592284613370746[1] = state[6];
}
void H_28(double *state, double *unused, double *out_6908639744438079814) {
   out_6908639744438079814[0] = 0;
   out_6908639744438079814[1] = 0;
   out_6908639744438079814[2] = 0;
   out_6908639744438079814[3] = 0;
   out_6908639744438079814[4] = 0;
   out_6908639744438079814[5] = 1;
   out_6908639744438079814[6] = 0;
   out_6908639744438079814[7] = 0;
   out_6908639744438079814[8] = 0;
   out_6908639744438079814[9] = 0;
   out_6908639744438079814[10] = 0;
   out_6908639744438079814[11] = 0;
   out_6908639744438079814[12] = 0;
   out_6908639744438079814[13] = 0;
   out_6908639744438079814[14] = 1;
   out_6908639744438079814[15] = 0;
}
}

extern "C"{
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



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
