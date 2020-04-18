// Copyright 2014 - 2015 The MathWorks, Inc.

#include "kalmanfilter.h"

static double x_est[6];
static double p_est[36];
static void kalmanfilter_init(void);
static void kalmanfilter_init(void)
{
  int i;
  for (i = 0; i < 6; i++) {
    x_est[i] = 0.0;
  }

  memset(&p_est[0], 0, 36U * sizeof(double));
}

void kalmanfilter(const double z[2], double y[2])
{
  signed char Q[36];
  int r2;
  double a[36];
  int k;
  double x_prd[6];
  static const signed char b_a[36] = { 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,
    1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1 };

  int i0;
  double p_prd[36];
  double a21;
  int r1;
  static const signed char b[36] = { 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1,
    0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 };

  double c_a[12];
  double S[4];
  static const signed char d_a[12] = { 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 };

  static const signed char b_b[12] = { 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 };

  static const short R[4] = { 1000, 0, 0, 1000 };

  double B[12];
  double a22;
  double b_z[2];
  double klm_gain[12];
  for (r2 = 0; r2 < 36; r2++) {
    Q[r2] = 0;
  }

  for (k = 0; k < 6; k++) {
    Q[k + 6 * k] = 1;
    x_prd[k] = 0.0;
    for (r2 = 0; r2 < 6; r2++) {
      x_prd[k] += (double)b_a[k + 6 * r2] * x_est[r2];
      a[k + 6 * r2] = 0.0;
      for (i0 = 0; i0 < 6; i0++) {
        a[k + 6 * r2] += (double)b_a[k + 6 * i0] * p_est[i0 + 6 * r2];
      }
    }
  }

  for (r2 = 0; r2 < 6; r2++) {
    for (i0 = 0; i0 < 6; i0++) {
      a21 = 0.0;
      for (r1 = 0; r1 < 6; r1++) {
        a21 += a[r2 + 6 * r1] * (double)b[r1 + 6 * i0];
      }

      p_prd[r2 + 6 * i0] = a21 + (double)Q[r2 + 6 * i0];
    }
  }

  for (r2 = 0; r2 < 2; r2++) {
    for (i0 = 0; i0 < 6; i0++) {
      c_a[r2 + (i0 << 1)] = 0.0;
      for (r1 = 0; r1 < 6; r1++) {
        c_a[r2 + (i0 << 1)] += (double)d_a[r2 + (r1 << 1)] * p_prd[i0 + 6 * r1];
      }
    }

    for (i0 = 0; i0 < 2; i0++) {
      a21 = 0.0;
      for (r1 = 0; r1 < 6; r1++) {
        a21 += c_a[r2 + (r1 << 1)] * (double)b_b[r1 + 6 * i0];
      }

      S[r2 + (i0 << 1)] = a21 + (double)R[r2 + (i0 << 1)];
    }

    for (i0 = 0; i0 < 6; i0++) {
      B[r2 + (i0 << 1)] = 0.0;
      for (r1 = 0; r1 < 6; r1++) {
        B[r2 + (i0 << 1)] += (double)d_a[r2 + (r1 << 1)] * p_prd[i0 + 6 * r1];
      }
    }
  }

  if (fabs(S[1]) > fabs(S[0])) {
    r1 = 1;
    r2 = 0;
  } else {
    r1 = 0;
    r2 = 1;
  }

  a21 = S[r2] / S[r1];
  a22 = S[2 + r2] - a21 * S[2 + r1];
  for (k = 0; k < 6; k++) {
    c_a[1 + (k << 1)] = (B[r2 + (k << 1)] - B[r1 + (k << 1)] * a21) / a22;
    c_a[k << 1] = (B[r1 + (k << 1)] - c_a[1 + (k << 1)] * S[2 + r1]) / S[r1];
  }

  for (r2 = 0; r2 < 2; r2++) {
    a21 = 0.0;
    for (i0 = 0; i0 < 6; i0++) {
      klm_gain[i0 + 6 * r2] = c_a[r2 + (i0 << 1)];
      a21 += (double)d_a[r2 + (i0 << 1)] * x_prd[i0];
    }

    b_z[r2] = z[r2] - a21;
  }

  for (r2 = 0; r2 < 6; r2++) {
    a21 = 0.0;
    for (i0 = 0; i0 < 2; i0++) {
      a21 += klm_gain[r2 + 6 * i0] * b_z[i0];
    }

    x_est[r2] = x_prd[r2] + a21;
    for (i0 = 0; i0 < 6; i0++) {
      a[r2 + 6 * i0] = 0.0;
      for (r1 = 0; r1 < 2; r1++) {
        a[r2 + 6 * i0] += klm_gain[r2 + 6 * r1] * (double)d_a[r1 + (i0 << 1)];
      }
    }

    for (i0 = 0; i0 < 6; i0++) {
      a21 = 0.0;
      for (r1 = 0; r1 < 6; r1++) {
        a21 += a[r2 + 6 * r1] * p_prd[r1 + 6 * i0];
      }

      p_est[r2 + 6 * i0] = p_prd[r2 + 6 * i0] - a21;
    }
  }

  for (r2 = 0; r2 < 2; r2++) {
    y[r2] = 0.0;
    for (i0 = 0; i0 < 6; i0++) {
      y[r2] += (double)d_a[r2 + (i0 << 1)] * x_est[i0];     // To inject error, multiple this value by 1.1
    }
  }
}

int sum_operation(int a, int b)
{
    return a + b;
}

void kalmanfilter_initialize(void)
{
  kalmanfilter_init();
}

void kalmanfilter_terminate(void)
{
}
