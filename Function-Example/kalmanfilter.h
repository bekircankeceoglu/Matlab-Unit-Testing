// Copyright 2014 - 2015 The MathWorks, Inc.

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

extern void kalmanfilter(const double z[2], double y[2]);
extern void kalmanfilter_initialize(void);
extern void kalmanfilter_terminate(void);
extern int sum_operation(int a, int b);