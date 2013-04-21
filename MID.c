/*
#------------------------------------------------------------------------------------------------------#
# MID (Mutual Information Dimension) for measuring statistical dependence between two random variables #
#------------------------------------------------------------------------------------------------------#
  Copyright (C) 2013 Mahito Sugiyama

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Contact: Mahito Sugiyama <mahito.sugiyama@tuebingen.mpg.de>
  
  Please refer the following article in your published research:
  Sugiyama, M., Borgwardt, K.M.: Measuring Statistical Dependence via the Mutual Information Dimension,  
  Proceedings of the 23rd International Joint Conference on Artificial Intelligence (IJCAI 2013), Beijing, China, Aug., 2013.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#define ZERO 0
#define ONE 1
#define BASE 2
#define BASE4 4
#define ROW_LENGTH 100
#define err(x) {printf("%s\n", (x));return;}

// Get the number of lines
long int getFile(FILE **fp, char *filename) {
  char *p, row[ROW_LENGTH];
  long int n = 1;

  if ((*fp=fopen(filename,"r")) == NULL) {
    printf("%s: cannot open the file\n", filename);
    return 0;
  }

  fgets(row, sizeof(row), *fp);
  p = strtok(row, ",\n");
  while (p != NULL) {
    p = strtok(NULL, ",\n");
  }

  while (fgets(row, sizeof(row), *fp)) n++;

  return n;
}

// Read a data file
void readFile(FILE *fp, double *x, double *y) {
  char *p, *ends, row[ROW_LENGTH];
  long int i = 0;

  rewind(fp);
  // input each line
  while(fgets(row, sizeof(row), fp) != NULL) {
    // divide at commas and LF
    p = strtok(row, ",\n");
    while (p != NULL) {
      if (i % BASE == 0) x[i / BASE] = strtod(p, &ends);
      else y[(i - 1) / BASE] = strtod(p, &ends);
      i++;
      p = strtok(NULL, ",\n");
    }
  }
}

// Min-max normalization to [0, 1] for each axis
void normalize(double *x, long int n) {
  long int i;
  double x_max, x_min;

  x_min = x[0];
  x_max = x[0];
  for (i = 0; i < n; i++) {
    x_max = x[i] > x_max ? x[i] : x_max;
    x_min = x[i] > x_min ? x_min : x[i];
  }
  if (x_min == x_max) {
    for (i = 0; i < n; i++) x[i] = 0.5;
  } else {
    for (i = 0; i < n; i++) {
      x[i] = (x[i] - x_min) / (x_max - x_min);
    }
  }
}

// Discretization
void discretize(double *x, int *codes, long int n, int k) {
  long int i;

  for (i = 0; i < n; i++) {
    codes[i] = floor(x[i] * k);
    if (codes[i] == k) {
      codes[i] = k - 1;
    }
  }
}

// Calculate entropy
double entropyEach(long int *B, long int n, int m) {
  int i;
  double p, result = 0;

  for (i = 0; i < m; i++) {
    if (B[i] != n && B[i] != 0) {
      p = B[i] / (double)n;
      result += -1 * p * (log(p) / log(BASE));
    }
  }

  return(result);
}

// Calculate joint entropy
double entropyCov(long int **B, long int n, int m) {
  int i, j;
  double p, result = 0;

  for (i = 0; i < m; i++) {
    for (j = 0; j < m; j++) {
      if (B[i][j] != n && B[i][j] != 0) {
        p = B[i][j] / (double)n;
        result += -1 * p * (log(p) / log(BASE));
      }
    }
  }
  return(result);
}

// Simple linear regression
void linearRegression(int n, double *x, double *y, double *a, double *b, double *rsq) {
  int i;
  double error = 0, tot = 0, xave = 0, yave = 0, xvar = 0, xyvar = 0;

  for (i = 0; i < n; i++) {
    xave += x[i];
    yave += y[i];
  }
  xave = xave / n;
  yave = yave / n;
  for (i = 0; i < n; i++) {
    xvar += pow(xave - x[i], 2);
    tot += pow(yave - y[i], 2);
    xyvar += (xave - x[i]) * (yave - y[i]);
  }
  xvar = xvar / n;
  xyvar = xyvar / n;

  *a = xyvar / xvar;
  *b = yave - *a * xave;

  for (i = 0; i < n; i++) {
    error += pow((*a * x[i] + *b) - y[i], 2);
  }

  *rsq = 1 - error / tot;
}

// Estimation of information dimension
double estimate(int xnum, double *yall, int width, _Bool cov, double minent) {
  _Bool flag;
  int i, j, i_end;
  double a, b, rsq, rsq_before = 0, coef = 0, *x, *y;

  i_end = xnum - width + 1;
  x = (double *)malloc(sizeof(double) * width);
  y = (double *)malloc(sizeof(double) * width);

  for (i = 1; i <= i_end; i++) {
    flag = 0;
    for (j = 0; j < width; j++) {
      x[j] = i + j;
      y[j] = yall[i + j - 1];
      if (j > 0 && y[j] == y[j - 1]) flag = 1;
    }
    if (flag == 1) {
      a = 0; rsq = 0;
    } else {
      linearRegression(width, x, y, &a, &b, &rsq);
    }
    if (cov == 0) {
      if (rsq > rsq_before) coef = a;
    } else {
      if (rsq > rsq_before && a > minent) coef = a;
    }
    rsq_before = rsq;
  }

  if (cov == 1 && coef == 0) coef = minent;
  return coef;
}

// Estimation of information dimension for two variables
double estimateCov(int xnum, double *yall, int width, _Bool cov, double minent) {
  double res;

  do {
    if (width > 1) {
      res = estimate(xnum, yall, width--, cov, minent);
    } else {
      res = 0;
      break;
    }
  } while (res <= 0);
  return(res);
}

double keepmax(double x, double y, double xy) {
  double tmp;

  tmp = x > y ? x : y;
  return(xy > tmp ? xy : tmp);
}

// Calculate information dimension for x, y, and xy
void idim(double *x, double *y, long int n, int level_max, int level_max_cov, double *idim_x, double *idim_y, double *idim_xy) {
  int level, k, m = 0, width, *codes_x, *codes_y;
  long int i, j, *B_x, *B_y, **B_xy;
  double *result_x, *result_y, *result_xy;

  codes_x = (int *)malloc(sizeof(int) * n);
  codes_y = (int *)malloc(sizeof(int) * n);

  result_x = (double *)malloc(sizeof(double) * level_max);
  result_y = (double *)malloc(sizeof(double) * level_max);
  result_xy = (double *)malloc(sizeof(double) * level_max_cov);

  // Normalization
  normalize(x, n);
  normalize(y, n);

  for (level=1; level <= level_max; level++) {
    // 2^level
    k = 1 << level;

    // Preparation
    B_x = (long int *)malloc(sizeof(long int) * k);
    B_y = (long int *)malloc(sizeof(long int) * k);
    memset(B_x, 0, sizeof(long int) * k);
    memset(B_y, 0, sizeof(long int) * k);
    if (level <= level_max_cov) {
      B_xy = (long int **)malloc(sizeof(long int *) * k);
      for (i = 0; i < k; i++) {
        B_xy[i] = (long int *)malloc(sizeof(long int) * k);
        memset(B_xy[i], 0, sizeof(long int) * k);
      }
    }

    // Discretization
    discretize(x, codes_x, n, k);
    discretize(y, codes_y, n, k);

    // Counting
    for (i = 0; i < n; i++) {
      (B_x[codes_x[i]])++;
      (B_y[codes_y[i]])++;
      if (level <= level_max_cov)
        (B_xy[codes_x[i]][codes_y[i]])++;
    }

    // Calculating
    result_x[m] = entropyEach(B_x, n, k);
    result_y[m] = entropyEach(B_y, n, k);
    if (level <= level_max_cov) {
      result_xy[m] = entropyCov(B_xy, n, k);
    }
    m++;

    // Free
    free(B_x); free(B_y);
    if (level <= level_max_cov) {
      for (i = 0; i < k; i++)
        free(B_xy[i]);
      free(B_xy);
    }
  }

  free(codes_x); free(codes_y);

  // Estimation of MID
  width = ceil(log(n) / log(BASE4));
  *idim_x = estimate(level_max, result_x, width, ZERO, ZERO);
  *idim_y = estimate(level_max, result_y, width, ZERO, ZERO);
  *idim_xy = estimateCov(level_max_cov, result_xy, width, ONE, *idim_x < *idim_y ? *idim_x : *idim_y);
  *idim_xy = keepmax(*idim_x, *idim_y, *idim_xy);

  free(result_x); free(result_y); free(result_xy);
}

// Main function
int main(int argc, char *argv[]) {
  int level_max, level_max_cov;
  long int i, n;
  double idim_x, idim_y, idim_xy, mid, *x, *y;
  FILE *fp;

  if (argc < 2) err("Error. Please input a filename");

  // Get the number of objects n
  n = getFile(&fp, argv[1]);
  if (n == 0) err("Error. There exist no data in the file.");

  // Prepare memory for a dataset
  x = (double *)malloc(sizeof(double) * n);
  y = (double *)malloc(sizeof(double) * n);

  // Read a data file
  readFile(fp, x, y);
  fclose(fp);

  // The maximum level
  level_max = floor(log(n) / log(BASE));
  level_max_cov = floor(log(n) / log(BASE4)) + 4;

  // MAIN PART
  idim(x, y, n, level_max, level_max_cov, &idim_x, &idim_y, &idim_xy);
  mid = idim_x + idim_y - idim_xy;
  if (mid > 1) mid = 1;
  if (mid < 0) mid = 0;

  printf("idim_x:  %f\n", idim_x);
  printf("idim_y:  %f\n", idim_y);
  printf("idim_xy: %f\n", idim_xy);
  printf("MID:     %f\n", mid);

  return 0;
}
