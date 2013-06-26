/*
 * doker_spe.c
 *
 * autor: Łukasz Bieniasz-Krzywiec
 *
 */


#include <spu_intrinsics.h>
#include <spu_mfcio.h>
#include <math.h>
#include <simdmath.h>

#ifdef TIMER
# include<spu_timer.h>
#endif

#include "../modules/err.h"
#include "../3dhit.h"

/*
 * PRZYDATNE MAKRA:
 */

#define waitag(t) mfc_write_tag_mask(1 << t); mfc_read_tag_status_all();

#undef TICS
#undef TICE
#undef TICP

#ifdef TIMER
  uint64_t timebase = 0, stics[16], etics[16], tics[16];
# define TICS(i) stics[(i)] = spu_clock_read()
# define TICE(i) do {\
    etics[(i)] = spu_clock_read();\
    tics[(i)] += etics[(i)] - stics[(i)];\
  } while(0)
# define TICP(i, s) do {\
    fprintf(stderr, "%s: %lld t\n", (s), tics[(i)]);\
    fprintf(stderr, "%s: %.4f s\n", (s), (float)tics[(i)] / (float)timebase);\
    fprintf(stderr, "%s: %.4f %%\n", (s), (float)tics[(i)] / (float)tics[0] * 100);\
  } while(0)
#else
# define TICS(i)
# define TICE(i)
# define TICP(i, s)
#endif


/*
 * DEKLARACJE STAŁYCH:
 */

#define Zero        0.0000001
#define Max_orth    0.001
#define Huge        (1.0e10)
#define MOVE_NO
#define TEST_MIRROR_NO
#define PRINT_MIRROR_NO
#define ENABLE_MIRROR_NO
#define TEST_ORTHOGONALITY
#define WARN_ORTHOGONALITY_ERROR_NO

#define AA          1           /* ali to ali */
#define AG          2           /* ali to gap */
#define VH          4           /* vertical */
#define ZE          8           /* zero */

#define PAIR
#define PAIRFULL

#define NEG_SCALE   10000.0
#define POS_SCALE   30000.0


/*
 * DEKLARACJE TYPÓW:
 */

/* Białko. */
typedef struct spe_protein_t_ {
  int     shift, a, b;
  short   num[SEGSIZE2];
  char    seq[SEGSIZE2];
  coor_t  coor[SEGSIZE2];
  coor_t  c[SEGSIZE2];
} spe_protein_t;

/* Krótkie białko. */
typedef struct protein_len_tmp_t_ {
  uchar num[ceil16(LEN) * sizeof(short) + 16];
  uchar seq[ceil16(LEN) * sizeof(char) + 16];
  uchar coor[ceil16(LEN) * sizeof(coor_t) + 16];
  uchar c[ceil16(LEN) * sizeof(coor_t) + 16];
} protein_len_tmp_t __attribute__ ((aligned(16)));


/*
 * DEKLARACJE ZMIENNYCH GLOBALNYCH:
 */

/* Tag IDs */
uint                    tag, tag_in, tag_out, tags_len[8];

/* Parametry wątku obliczeniowego. */
spe_worker_parameters_t params;

/* Wyjście wątku obliczeniowego. */
spe_worker_output_t     output;

/* Białka. */
spe_protein_t           prot1, prot2;

/* Tablica używane w wielu funkcjach. */
uchar                   memo[SEGSIZE2 * SEGSIZE2 / 2 + 1];

/* Pomocnicze przy ściąganiu: */
protein_len_tmp_t       plen1_tmp, plen2_tmp;
int                     off_num1, off_seq1, off_coor1, off_c1,
                        off_num2, off_seq2, off_coor2, off_c2;


/* Wyniki funkcji wmove(). */
coor_t                  yc;
coor_t                  xc;
double                  U[3][3];

/* Tablice używane w funkcji align2(). */
float                   al[SEGSIZE2] __attribute__ ((aligned(128)));
float                   gl[SEGSIZE2] __attribute__ ((aligned(128)));

/* Wyniki funkcji align2(). */
int                     pathn;
float                   pathv[SEGSIZE2] __attribute__ ((aligned(128)));

/* Tablice używane w jądrze obliczeniowym. */
short                   m[SEGSIZE2 * SEGSIZE2] __attribute__ ((aligned(128)));
coor_t                  xnew[SEGSIZE2] __attribute__ ((aligned(128)));
coor_t                  ynew[SEGSIZE2] __attribute__ ((aligned(128)));
coor_t                  cnew[SEGSIZE2] __attribute__ ((aligned(128)));
float                   movew[SEGSIZE2] __attribute__ ((aligned(128)));
int                     pathp[SEGSIZE2] __attribute__ ((aligned(128)));


/*
 * DEFINICJE PODPROGRAMÓW:
 */

extern int sys_nerr;

void syserr(const char *fmt, ...) {
  va_list fmt_args;

  fprintf(stderr, "ERROR: ");

  va_start(fmt_args, fmt);
  vfprintf(stderr, fmt, fmt_args);
  va_end (fmt_args);
  fprintf(stderr," (%d; %s)\n", errno, strerror(errno));
  exit(1);
}

void fatal(const char *fmt, ...) {
  va_list fmt_args;

  fprintf(stderr, "ERROR: ");

  va_start(fmt_args, fmt);
  vfprintf(stderr, fmt, fmt_args);
  va_end (fmt_args);

  fprintf(stderr,"\n");
  exit(1);
}

double wmove(const coor_t *y_old, coor_t *x_old, const float *movew, const int n) {
  int     i, j, k;
  double  R[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
  double  r_r[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
  double  eval[3] = {0.0,0.0,0.0};
  double  evec[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
  double  b[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
  double  rms = 0.0, wc = 0.0;

  static double Pi = 0.0;

  float * __restrict wmv_x0 = (float*)(memo);
  float * __restrict wmv_x1 = (float*)(memo + 1 * sizeof(float) * ceil16(SEGSIZE2));
  float * __restrict wmv_x2 = (float*)(memo + 2 * sizeof(float) * ceil16(SEGSIZE2));
  float * __restrict wmv_y0 = (float*)(memo + 3 * sizeof(float) * ceil16(SEGSIZE2));
  float * __restrict wmv_y1 = (float*)(memo + 4 * sizeof(float) * ceil16(SEGSIZE2));
  float * __restrict wmv_y2 = (float*)(memo + 5 * sizeof(float) * ceil16(SEGSIZE2));
  float * __restrict x_old0 = (float*)(memo + 6 * sizeof(float) * ceil16(SEGSIZE2));
  float * __restrict x_old1 = (float*)(memo + 7 * sizeof(float) * ceil16(SEGSIZE2));
  float * __restrict x_old2 = (float*)(memo + 8 * sizeof(float) * ceil16(SEGSIZE2));
  float * __restrict y_old0 = (float*)(memo + 9 * sizeof(float) * ceil16(SEGSIZE2));
  float * __restrict y_old1 = (float*)(memo + 10 * sizeof(float) * ceil16(SEGSIZE2));
  float * __restrict y_old2 = (float*)(memo + 11 * sizeof(float) * ceil16(SEGSIZE2));
  float * __restrict mvw    = (float*)(memo + 12 * sizeof(float) * ceil16(SEGSIZE2));
  float * __restrict tmp0   = (float*)(memo + 13 * sizeof(float) * ceil16(SEGSIZE2));
  float * __restrict tmp1   = (float*)(memo + 14 * sizeof(float) * ceil16(SEGSIZE2));
  float * __restrict tmp2   = (float*)(memo + 15 * sizeof(float) * ceil16(SEGSIZE2));
  float * __restrict tmp3   = (float*)(memo + 16 * sizeof(float) * ceil16(SEGSIZE2));
  float * __restrict tmp4   = (float*)(memo + 17 * sizeof(float) * ceil16(SEGSIZE2));
  float * __restrict tmp5   = (float*)(memo + 18 * sizeof(float) * ceil16(SEGSIZE2));
  float * __restrict tmp6   = (float*)(memo + 19 * sizeof(float) * ceil16(SEGSIZE2));
  float * __restrict tmp7   = (float*)(memo + 20 * sizeof(float) * ceil16(SEGSIZE2));
  float * __restrict tmp8   = (float*)(memo + 21 * sizeof(float) * ceil16(SEGSIZE2));
  float * __restrict _a     = (float*)(memo + 22 * sizeof(float) * ceil16(SEGSIZE2));

  __align_hint(wmv_x0, 16, 0);
  __align_hint(wmv_x1, 16, 0);
  __align_hint(wmv_x2, 16, 0);
  __align_hint(wmv_y0, 16, 0);
  __align_hint(wmv_y1, 16, 0);
  __align_hint(wmv_y2, 16, 0);
  __align_hint(x_old0, 16, 0);
  __align_hint(x_old1, 16, 0);
  __align_hint(x_old2, 16, 0);
  __align_hint(y_old0, 16, 0);
  __align_hint(y_old1, 16, 0);
  __align_hint(y_old2, 16, 0);
  __align_hint(mvw, 16, 0);
  __align_hint(tmp0, 16, 0);
  __align_hint(tmp1, 16, 0);
  __align_hint(tmp2, 16, 0);
  __align_hint(tmp3, 16, 0);
  __align_hint(tmp4, 16, 0);
  __align_hint(tmp5, 16, 0);
  __align_hint(tmp6, 16, 0);
  __align_hint(tmp7, 16, 0);
  __align_hint(tmp8, 16, 0);
  __align_hint(_a, 16, 0);

  vec_float4 *wmv_x0v = (vec_float4*)wmv_x0;
  vec_float4 *wmv_x1v = (vec_float4*)wmv_x1;
  vec_float4 *wmv_x2v = (vec_float4*)wmv_x2;
  vec_float4 *wmv_y0v = (vec_float4*)wmv_y0;
  vec_float4 *wmv_y1v = (vec_float4*)wmv_y1;
  vec_float4 *wmv_y2v = (vec_float4*)wmv_y2;
  vec_float4 *x_old0v = (vec_float4*)x_old0v;
  vec_float4 *x_old1v = (vec_float4*)x_old1v;
  vec_float4 *x_old2v = (vec_float4*)x_old2v;
  vec_float4 *y_old0v = (vec_float4*)y_old0v;
  vec_float4 *y_old1v = (vec_float4*)y_old1v;
  vec_float4 *y_old2v = (vec_float4*)y_old2v;
  vec_float4 *mvwv = (vec_float4*)mvw;
  vec_float4 *tmp0v = (vec_float4*)tmp0;
  vec_float4 *tmp1v = (vec_float4*)tmp1;
  vec_float4 *tmp2v = (vec_float4*)tmp2;
  vec_float4 *tmp3v = (vec_float4*)tmp3;
  vec_float4 *tmp4v = (vec_float4*)tmp4;
  vec_float4 *tmp5v = (vec_float4*)tmp5;
  vec_float4 *tmp6v = (vec_float4*)tmp6;
  vec_float4 *tmp7v = (vec_float4*)tmp7;
  vec_float4 *tmp8v = (vec_float4*)tmp8;
  vec_float4 *_av = (vec_float4*)_a;

  // note: not vectorized: unsupported unaligned store.
  for (i = 0; i < n; ++i) {
    y_old0[i] = y_old[i][0];
    y_old1[i] = y_old[i][1];
    y_old2[i] = y_old[i][2];
    x_old0[i] = x_old[i][0];
    x_old1[i] = x_old[i][1];
    x_old2[i] = x_old[i][2];
    mvw[i] = movew[i];
  }

  for (i = 0; i < 3; ++i) {
    yc[i] = 0.0;
    xc[i] = 0.0;
    for (j = 0; j < 3; ++j) {
      U[i][j] = 0.0;
    }
  }

  /* centers of Y,X */
  if (Pi == 0.0) {
    Pi = 2.0 * asin((double)1.0);
  }
  /* lbk
  for (i = 0; i < n; i++) {
    wc += movew[i];
    for (j = 0; j < 3; ++j) {
      yc[j] += y_old[i][j] * movew[i];
      xc[j] += x_old[i][j] * movew[i];
    }
  }
  */
  /* lbk
  // note: not vectorized: unsupported use in stmt.
  for (i = 0; i < n; ++i) {
    wc += mvw[i];
    yc[0] += y_old0[i] * mvw[i];
    yc[1] += y_old1[i] * mvw[i];
    yc[2] += y_old2[i] * mvw[i];
    xc[0] += x_old0[i] * mvw[i];
    xc[1] += x_old1[i] * mvw[i];
    xc[2] += x_old2[i] * mvw[i];
  }
  */
  // VECTORIZED!
  for (i = 0; i < n; ++i) {
    tmp0[i] = y_old0[i] * mvw[i];
    tmp1[i] = y_old1[i] * mvw[i];
    tmp2[i] = y_old2[i] * mvw[i];
    tmp3[i] = x_old0[i] * mvw[i];
    tmp4[i] = x_old1[i] * mvw[i];
    tmp5[i] = x_old2[i] * mvw[i];
  }
  // note: not vectorized: unsupported use in stmt.
  for (i = 0; i < n; ++i) {
    wc += mvw[i];
    yc[0] += tmp0[i];
    yc[1] += tmp1[i];
    yc[2] += tmp2[i];
    xc[0] += tmp3[i];
    xc[1] += tmp4[i];
    xc[2] += tmp5[i];
  }

  if (wc <= 0.0) {
    return -Huge;
  }
  for (j = 0; j < 3; ++j) {
    yc[j] /= wc;
    xc[j] /= wc;
  }
  /* lbk
  for (i = 0; i < n; ++i) {
    for (j = 0; j < 3; ++j) {
      wmv_y[i][j] = (y_old[i][j] - yc[j]);
      wmv_x[i][j] = (x_old[i][j] - xc[j]);
    }
  }
  */
  // VECTORIZED!
  for (i = 0; i < n; ++i) {
    wmv_y0[i] = (y_old0[i] - yc[0]);
    wmv_y1[i] = (y_old1[i] - yc[1]);
    wmv_y2[i] = (y_old2[i] - yc[2]);
    wmv_x0[i] = (x_old0[i] - xc[0]);
    wmv_x1[i] = (x_old1[i] - xc[1]);
    wmv_x2[i] = (x_old2[i] - xc[2]);
  }
  /* end of centers of sets */

  /* correlation matrix R */
  R[0][0] = R[0][1] = R[0][2] = R[1][0] = R[1][1] = R[1][2] = R[2][0] = R[2][1] = R[2][2] = 0.0;
  /* lbk
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      for (k = 0; k < n; ++k) {
        R[i][j] += wmv_y[k][i] * wmv_x[k][j] * movew[k];
      }
    }
  }
  */
  /* lbk
  // note: not vectorized: possible dependence between data-refs R[0][0] and R[0][0]
  for (k = 0; k < n; ++k) {
    R[0][0] += wmv_y0[k] * wmv_x0[k] * mvw[k];
    R[0][1] += wmv_y0[k] * wmv_x1[k] * mvw[k];
    R[0][2] += wmv_y0[k] * wmv_x2[k] * mvw[k];
    R[1][0] += wmv_y1[k] * wmv_x0[k] * mvw[k];
    R[1][1] += wmv_y1[k] * wmv_x1[k] * mvw[k];
    R[1][2] += wmv_y1[k] * wmv_x2[k] * mvw[k];
    R[2][0] += wmv_y2[k] * wmv_x0[k] * mvw[k];
    R[2][1] += wmv_y2[k] * wmv_x1[k] * mvw[k];
    R[2][2] += wmv_y2[k] * wmv_x2[k] * mvw[k];
  }
  */
  /* lbk
  // note: not vectorized: unsupported unaligned store. ???
  for (k = 0; k < n; ++k) {
    tmp0[k] = wmv_y0[k] * wmv_x0[k] * mvw[k];
    tmp1[k] = wmv_y0[k] * wmv_x1[k] * mvw[k];
    tmp2[k] = wmv_y0[k] * wmv_x2[k] * mvw[k];
    tmp3[k] = wmv_y1[k] * wmv_x0[k] * mvw[k];
    tmp4[k] = wmv_y1[k] * wmv_x1[k] * mvw[k];
    tmp5[k] = wmv_y1[k] * wmv_x2[k] * mvw[k];
    tmp6[k] = wmv_y2[k] * wmv_x0[k] * mvw[k];
    tmp7[k] = wmv_y2[k] * wmv_x1[k] * mvw[k];
    tmp8[k] = wmv_y2[k] * wmv_x2[k] * mvw[k];
  }
  */
  for (k = 0; k < n >> 2; ++k) {
    tmp0v[k] = spu_mul(spu_mul(wmv_y0v[k], wmv_x0v[k]), mvwv[k]);
    tmp1v[k] = spu_mul(spu_mul(wmv_y0v[k], wmv_x1v[k]), mvwv[k]);
    tmp2v[k] = spu_mul(spu_mul(wmv_y0v[k], wmv_x2v[k]), mvwv[k]);
    tmp3v[k] = spu_mul(spu_mul(wmv_y1v[k], wmv_x0v[k]), mvwv[k]);
    tmp4v[k] = spu_mul(spu_mul(wmv_y1v[k], wmv_x1v[k]), mvwv[k]);
    tmp5v[k] = spu_mul(spu_mul(wmv_y1v[k], wmv_x2v[k]), mvwv[k]);
    tmp6v[k] = spu_mul(spu_mul(wmv_y2v[k], wmv_x0v[k]), mvwv[k]);
    tmp7v[k] = spu_mul(spu_mul(wmv_y2v[k], wmv_x1v[k]), mvwv[k]);
    tmp8v[k] = spu_mul(spu_mul(wmv_y2v[k], wmv_x2v[k]), mvwv[k]);
  }
  for (k = n & n >> 2 << 2; k < n; ++k) {
    tmp0[k] = wmv_y0[k] * wmv_x0[k] * mvw[k];
    tmp1[k] = wmv_y0[k] * wmv_x1[k] * mvw[k];
    tmp2[k] = wmv_y0[k] * wmv_x2[k] * mvw[k];
    tmp3[k] = wmv_y1[k] * wmv_x0[k] * mvw[k];
    tmp4[k] = wmv_y1[k] * wmv_x1[k] * mvw[k];
    tmp5[k] = wmv_y1[k] * wmv_x2[k] * mvw[k];
    tmp6[k] = wmv_y2[k] * wmv_x0[k] * mvw[k];
    tmp7[k] = wmv_y2[k] * wmv_x1[k] * mvw[k];
    tmp8[k] = wmv_y2[k] * wmv_x2[k] * mvw[k];
  }
  for (k = 0; k < n; ++k) {
    R[0][0] += tmp0[k];
    R[0][1] += tmp1[k];
    R[0][2] += tmp2[k];
    R[1][0] += tmp3[k];
    R[1][1] += tmp4[k];
    R[1][2] += tmp5[k];
    R[2][0] += tmp6[k];
    R[2][1] += tmp7[k];
    R[2][2] += tmp8[k];
  }
  /* end correlation matrix */

  /* r_r = correlation matrix ** 2 */
  for (i = 2; i >= 0; --i) {
    for (j = i; j >= 0; --j) {
      r_r[i][j] = 0.0;
      for (k = 2; k >= 0; --k) {
        r_r[i][j] += R[k][i] * R[k][j];
      }
      r_r[j][i] = r_r[i][j];
    }
  }
  /* end r_r = correlation matrix ** 2 */

  /* eigenvalue , eigenvector */
  {
    /* det(r_r - eval * I) = 0 */
    double  _a = r_r[0][0];
    double  _b = r_r[0][1];
    double  _c = r_r[0][2];
    double  _d = r_r[2][2];
    float   _e = r_r[1][1];
    float   _f = r_r[1][2];
    float   _g, _h, PHI;
    /* A= -1; */
    double   B = -_a -_e -_d;
    double   C = - _f * _f - _b * _b - _c * _c + _a * _d + _e * _d + _a * _e;
    double   D = - _a * _e * _d - 2.0 * _b * _c * _f + _a * _f * _f + _e * _c * _c + _d * _b * _b;
    double   P = (C - B * B / 3.0) / 3.0;
    double   Q = (2.0 * B * B * B / 27.0 -( B * C / 3.0) + D) / 2.0;
    double _r = sqrt((double) fabs(P));

    if (Q < 0.0) {
      _r = -_r;
    }
    if (fabs(_r) > Zero && fabs(Q) > Zero) {
      double _x = Q / (_r *_r * _r);
      if (_x > 1.0) {
        _x = 1.0;
      }
      if (_x < -1.0) {
        _x = -1.0;
      }
      PHI = acos(_x) / 3.0;
    } else {
      PHI = Pi / 6.0;
    }
    eval[0] = -2.0 * _r * cos(PHI) - (B / 3.0);
    eval[1] = 2.0 * _r * cos(Pi / 3.0 - PHI) - B / 3.0;
    eval[2] = 2.0 * _r * cos(Pi / 3.0 + PHI) - B / 3.0;

    /* found eigenvalues using det = 0 */

    /* sorting evalues */
    if (eval[1] > eval[0] && eval[1] > eval[2]) {
      double _x = eval[0];
      eval[0] = eval[1];
      eval[1] = _x;
    } else {
      if (eval[2] > eval[0]) {
        double _x = eval[0];
        eval[0] = eval[2];
        eval[2] = _x;
      }
    }
    if (eval[2] > eval[1]) {
      double _x = eval[1];
      eval[1] = eval[2];
      eval[2] = _x;
    }
    /* end sorting eigenvalues */

    /* eigenvectors using GAUS */
    _b = r_r[0][1];
    _c = r_r[0][2];
    _d = r_r[1][0];
    _f = r_r[1][2];
    _g = r_r[2][0];
    _h = r_r[2][1];
    for (i = 0; i < 2; ++i) {
      double _x, _y, _z, _i;
      _a = r_r[0][0] - eval[i];
      _e = r_r[1][1] - eval[i];
      _i = r_r[2][2] - eval[i];
      _x = _a * _e - _d * _b;
      _y = _a * _i - _c * _g;
      _z = _e * _i - _f * _h;

      /* overwriting B, C, D */
      B = fabs(_x);
      C = fabs(_y);
      D = fabs(_z);
      if (B + C + D > Zero) {
        if (B >= C && B >= D) {
          evec[i][0] = (_e * _c - _f * _b) / _x;
          evec[i][1] = (_a * _f - _c * _d) / _x;
          evec[i][2] = -1.0;
        } else {
          if ( C > D || (C == D && i)) {
            evec[i][0] = (_b * _i - _c * _h) / _y;
            evec[i][2] = (_a * _h - _b * _g) / _y;
            evec[i][1] = -1.0;
          } else {
            evec[i][1] = (_d * _i - _f * _g) / _z;
            evec[i][2] = (_e * _g - _d * _h) / _z;
            evec[i][0] = -1.0;
          }
        }
      } else {
        if (!i) {
          evec[0][0] = 1.0;
          evec[0][1] = 0.0;
          evec[0][2] = 0.0;
        } else {
          /* assumed i = 0 or 1 */
          /* create evec[i] orthogonal to evec[0] */
          B = fabs(evec[0][0]);
          C = fabs(evec[0][1]);
          D = fabs(evec[0][2]);
          if (B > C && B > D) {
            evec[i][0] = evec[0][1];
            evec[i][1] = -evec[0][0];
            evec[i][2] = 0.0;
          } else {
            evec[i][0] = 0.0;
            if (C > D) {
              evec[i][1] = evec[0][2];
              evec[i][2] = -evec[0][1];
            } else {
              evec[i][1] = -evec[0][2];
              evec[i][2] = evec[0][1];
            }
          }
        }
      }

      /* normalizing eigenvectors */
      _x = sqrt(evec[i][0] * evec[i][0] + evec[i][1] * evec[i][1] + evec[i][2] * evec[i][2]);
      for (j = 0; j < 3; ++j) {
        evec[i][j] /= _x;
      }
    }
  }
  /* end eigenvalue , eigenvector */

  /* right handed system */
  evec[2][0] = evec[0][1] * evec[1][2] - evec[1][1] * evec[0][2];
  evec[2][1] = evec[0][2] * evec[1][0] - evec[0][0] * evec[1][2];
  evec[2][2] = evec[0][0] * evec[1][1] - evec[1][0] * evec[0][1];

  /* calculate b = R * evec/.. */
  for (i = 0; i < 2; ++i) {
    double _y = fabs(eval[i]);
    for (j = 0; j < 3; ++j) {
      double _x = 0.0;
      /* lbk
      for (k = 0; k < 3; ++k) {
        _x += R[j][k] * evec[i][k];
      }
      */
      _x = R[j][0] * evec[i][0] + R[j][1] * evec[i][1] + R[j][2] * evec[i][2];
      if (_y <= Zero) {
        if (!i) {
          b[0][0] = 1.0;
          b[0][1] = 0.0;
          b[0][2] = 0.0;
        } else {
          /* assumad i = 0 or 1 */
          /* create b[i] orthogonal to b[0] */
          double B = fabs(b[0][0]);
          double C = fabs(b[0][1]);
          double D = fabs(b[0][2]);
          if (B > C && B > D) {
            b[i][0] = b[0][1];
            b[i][1] = -b[0][0];
            b[i][2] = 0.0;
          } else {
            b[i][0] = 0.0;
            if (C > D) {
              b[i][1] = b[0][2];
              b[i][2] = -b[0][1];
            } else {
              b[i][1] = -b[0][2];
              b[i][2] = b[0][1];
            }
          }
          _x = sqrt(b[i][0] * b[i][0] + b[i][1] * b[i][1] + b[i][2] * b[i][2]);
          /* lbk
          for(j = 0; j < 3; ++j) {
            b[i][j] /= _x;
          }
          */
          b[i][0] /= _x;
          b[i][1] /= _x;
          b[i][2] /= _x;
        }
      } else {
        b[i][j] = _x / sqrt(_y);
      }
    }
  }
  b[2][0] = b[0][1] * b[1][2] - b[1][1] * b[0][2];
  b[2][1] = b[0][2] * b[1][0] - b[0][0] * b[1][2];
  b[2][2] = b[0][0] * b[1][1] - b[1][0] * b[0][1];
  /* end calculating b */

  /* calculatin U */
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      double _x = 0.0;
      for (k = 0; k < 3; ++k) {
        _x += b[k][i] * evec[k][j];
      }
      U[i][j] = _x;
    }
  }
  /* end calculation of U */

  /* moving x onto y */
  /* lbk
  float _x0, _x1, _x2;
  // note: not vectorized: complicated access pattern. ???
  // note: not vectorized: stmt not supported: D.5217_359 = (double) D.5216_358 ???
  for (i = 0; i < n; ++i) {
    _x0 = U[0][0] * wmv_x0[i] + U[0][1] * wmv_x1[i] + U[0][2] * wmv_x2[i];
    _x1 = U[1][0] * wmv_x0[i] + U[1][1] * wmv_x1[i] + U[1][2] * wmv_x2[i];
    _x2 = U[2][0] * wmv_x0[i] + U[2][1] * wmv_x1[i] + U[2][2] * wmv_x2[i];
    _a[i] = (_x0 - wmv_y0[i]) * (_x0 - wmv_y0[i]) * mvw[i] +
            (_x1 - wmv_y1[i]) * (_x1 - wmv_y1[i]) * mvw[i] +
            (_x2 - wmv_y2[i]) * (_x2 - wmv_y2[i]) * mvw[i];
  }
  */
  vec_float4 _x0v, _x1v, _x2v;
  for (i = 0; i < n >> 2; ++i) {
    _x0v = spu_madd(spu_splats((float)U[0][0]), wmv_x0v[i], spu_add(spu_mul(spu_splats((float)U[0][1]), wmv_x1v[i]), spu_mul(spu_splats((float)U[0][2]), wmv_x2v[i])));
    _x1v = spu_madd(spu_splats((float)U[1][0]), wmv_x0v[i], spu_add(spu_mul(spu_splats((float)U[1][1]), wmv_x1v[i]), spu_mul(spu_splats((float)U[1][2]), wmv_x2v[i])));
    _x2v = spu_madd(spu_splats((float)U[2][0]), wmv_x0v[i], spu_add(spu_mul(spu_splats((float)U[2][1]), wmv_x1v[i]), spu_mul(spu_splats((float)U[2][2]), wmv_x2v[i])));
    _av[i] = spu_madd((_x0v - wmv_y0v[i]), (_x0v - wmv_y0v[i]) * mvwv[i],
             (_x1v - wmv_y1v[i]) * (_x1v - wmv_y1v[i]) * mvwv[i] +
             (_x2v - wmv_y2v[i]) * (_x2v - wmv_y2v[i]) * mvwv[i]);
  }
  float _x0, _x1, _x2;
  for (i = n >> 2 << 2; i < n; ++i) {
    _x0 = U[0][0] * wmv_x0[i] + U[0][1] * wmv_x1[i] + U[0][2] * wmv_x2[i];
    _x1 = U[1][0] * wmv_x0[i] + U[1][1] * wmv_x1[i] + U[1][2] * wmv_x2[i];
    _x2 = U[2][0] * wmv_x0[i] + U[2][1] * wmv_x1[i] + U[2][2] * wmv_x2[i];
    _a[i] = (_x0 - wmv_y0[i]) * (_x0 - wmv_y0[i]) * mvw[i] +
            (_x1 - wmv_y1[i]) * (_x1 - wmv_y1[i]) * mvw[i] +
            (_x2 - wmv_y2[i]) * (_x2 - wmv_y2[i]) * mvw[i];
  }
  for (i = 0; i < n; ++i) {
    rms += _a[i];
  }

  rms = sqrt(rms / wc);

  return rms;

}

int matrix1w(
        short *m,
        const coor_t * p1, const coor_t * c1, const short * num1, const int lp1,
        const coor_t * p2, const coor_t * c2, const short * num2, const int lp2) {
  int   l1, l2, k = 0, n = 0;
  float c, d;

  for (l1 = 0; l1 < lp1; ++l1) {
    for (l2 = 0; l2 < lp2; ++l2) {
      m[k] = 0;
      if (num1[l1] > 0 && num2[l2] > 0) {
        c = c1[l1][0] * c2[l2][0] + c1[l1][1] * c2[l2][1] + c1[l1][2] * c2[l2][2];
        d = (p1[l1][0] - p2[l2][0]) * (p1[l1][0] - p2[l2][0]) +
            (p1[l1][1] - p2[l2][1]) * (p1[l1][1] - p2[l2][1]) +
            (p1[l1][2] - p2[l2][2]) * (p1[l1][2] - p2[l2][2]);

        if (c >= MINCOS && d <= SEGRMSD1 * SEGRMSD1) {
          m[k] = 1;
          n += 1;
        }
      }
      k += 1;
    }
  }

  return n;
}

float align1(short *m, const int lp1, const int lp2, int *pathp) {
  int l1, l2, k = 0;

  if (lp1 <= 0 || lp2 <= 0) return 0.0;

  /* forward */
  for (l2 = 1, ++k; l2 < lp2; ++l2, ++k) {
    if (m[k-1] > m[k]) {
      m[k] = m[k-1];
    }
  }
  for (l1 = 1; l1 < lp1; ++l1) {
    if (m[k-lp2] > m[k]) {
      m[k] = m[k-lp2];
    }
    for (l2 = 1, ++k; l2 < lp2; ++l2, ++k) {
      m[k] += m[k-lp2-1];
      if (m[k-1] > m[k]) {
        m[k] = m[k-1];
      }
      if (m[k-lp2] > m[k]) {
        m[k] = m[k-lp2];
      }
    }
  }

  /* backward */
  if (pathp != NULL) {
    for (--k; k >= lp2 && k % lp2; ) {
      if (m[k] == m[k-lp2]) {
        pathp[k/lp2] = -1;
        k -= lp2;
      } else if (m[k] == m[k-1]) {
        k -= 1;
      } else /* if m[k] == m[k-lp2-1] */ {
        pathp[k/lp2] = k % lp2;
        k -= lp2 + 1;
      }
    }
    if (m[k] > 0.0) {
      if (k < lp2) {
        for ( ; k > 0 && m[k] == m[k-1]; --k);
        pathp[0] = k;
      } else {
        for ( ; k > 0 && m[k] == m[k-lp2]; k -= lp2) {
          pathp[k/lp2] = -1;
        }
        pathp[k/lp2] = k % lp2;
        for(k -= lp2; k >= 0; k -= lp2) {
          pathp[k/lp2] = -1;
        }
      }
    } else {
      for (; k >= 0; k -= lp2) {
        pathp[k/lp2] = -1;
      }
    }
  }

  return m[lp1 * lp2 - 1];
}

int matrix2w(
        short *m,
        const coor_t *p1, const coor_t *c1, const short *num1, const int lp1,
        const coor_t *p2, const coor_t *c2, const short *num2, const int lp2) {
  int   l1, l2, k = 0, n = 0;
  float c, d, sigma = 2.0 * SEGSIGMA * SEGSIGMA;

  /*
   * Założenia: paraGap + paraExt <= 30000 / (2 * NEG_SCALE)
   */

  for (l1 = 0; l1 < lp1; ++l1) {
    for (l2 = 0; l2 < lp2; ++l2) {
      m[k] = (short) floor(-2.0 * (paraGap + paraExt) * NEG_SCALE);
      if (num1[l1] > 0 && num2[l2] > 0) {
        c = c1[l1][0] * c2[l2][0] +
            c1[l1][1] * c2[l2][1] +
            c1[l1][2] * c2[l2][2];
        d = (p1[l1][0] - p2[l2][0]) * (p1[l1][0] - p2[l2][0]) +
            (p1[l1][1] - p2[l2][1]) * (p1[l1][1] - p2[l2][1]) +
            (p1[l1][2] - p2[l2][2]) * (p1[l1][2] - p2[l2][2]);
        if (c >= MINCOS && d <= SEGRMSD2 * SEGRMSD2) {
#ifdef SEGRMSD0
          m[k] = (short)floor(exp(-d / sigma) * POS_SCALE);
#else
          m[k] = (short)(1.0 * POS_SCALE);
#endif
          n += 1;
        }
      }
      k += 1;
    }
  }

  return n;
}

#if (defined(ALIGLOB) || defined(ALIGLOLOC))
# define Loctest(a, k)
# define Glotest(a, g, k) do {\
    if ((a) > score) score = (a), be = (k);\
    if ((g) > score) score = (g), be = -(k);\
  } while (0)
#else
# define Glotest(a, g, k)
# define Loctest(a, k) do {\
    if ((a) <= 0.0) {\
      a = 0.0;\
      set_half_byte_array(memo, (k), get_half_byte_array(memo, (k)) | ZE);\
    } else {\
      if ((a) > score) {\
        score = (a);\
        be = (k);\
      }\
    }\
  } while(0)
#endif

#define scale_value(f) (((f) >= 0) ? (f) / POS_SCALE : (f) / NEG_SCALE)

/*
 * Operacje bitowe.
 *
 * ZAŁOŻENIA:
 *  value = 0,1,2,3,4,5,6,7,8
 */

uchar mask[2] = {240, 15};

void initialize_half_byte_array(uchar *hba) {
  int i;
  for (i = 0; i < SEGSIZE2 * SEGSIZE2 / 2 + 1; ++i) hba[i] = 0;
}

inline void set_half_byte_array(uchar *hba, int i, uchar value) {
  hba[i / 2] &= mask[!(i % 2)];
  hba[i / 2] |= (i % 2) ? value : value << 4;
}

inline uchar get_half_byte_array(uchar *hba, int i) {
  return (i % 2) ? hba[i/2] & mask[i % 2] : (hba[i/2] & mask[i % 2]) >> 4;
}

inline void or_half_byte_array(uchar *hba, int i, uchar value) {
  hba[i / 2] |= (i % 2) ? value : value << 4;
}

float align2(short *m, const int lp1, const int lp2, int* pathp) {
  int   segsize2 = 0;
  float score;
  int   be = 0;
  int   l1, l2, k, stat;

  if (lp1 > segsize2 || lp2 > segsize2) {
    if (segsize2 < SEGSIZE2) {
      segsize2 = SEGSIZE2;
    }
    if (lp1 > segsize2) {
      segsize2 = lp1 + 100;
    }
    if (lp2 > segsize2) {
      segsize2 = lp2 + 100;
    }
    initialize_half_byte_array(memo);
  }

#if (defined(ALIGLOB) || defined(ALIGLOLOC))
  score = -segsize2 * segsize2;
#else
  score = scale_value(m[0]);
#endif

  /* forward */
  for (l2 = k = 0; l2 < lp2; ++l2, ++k) {
    al[l2] = scale_value(m[k]);
    gl[l2] = -paraExt; /* changed from paraGap to paraExt */
#if (defined(ALIGLOB) || defined(ALIGLOLOC))
    al[l2] -= (float)l2 * paraExt;
    gl[l2] -= (float)l2 * paraExt;
#endif
    set_half_byte_array(memo, k, AA | VH);
    Loctest(al[l2], k);
  }
#ifdef ALIGLOLOC
  Glotest(al[l2-1], gl[l2-1], k-1);
#elif defined(ALIGLOB)
  Glotest(al[l2-1] - paraExt * (float)(lp1-1),
          gl[l2-1] - paraExt * (float)(lp1-1),
          k-1);
#endif
  for (l1 = 1; l1 < lp1; ++l1) {
    float ad = scale_value(m[k]);
    float gd = -paraExt; /* changed from paraGap to paraExt */
#ifdef ALIGLOB
    ad -= (float)l1 * paraExt;
    gd -= (float)l1 * paraExt;
#endif
    set_half_byte_array(memo, k, AA);
    Loctest(ad, k);
    for (l2 = 1, ++k; l2 < lp2; ++l2, ++k) {
      float an = scale_value(m[k]);
      /* position gapped */
      float gn = gd - paraExt;
      set_half_byte_array(memo, k, 0);
      if (gn < gl[l2] - paraExt) {
        gn = gl[l2] - paraExt;
        set_half_byte_array(memo, k, VH);
      }
      if (gn < ad - paraGap) {
        gn = ad - paraGap;
        set_half_byte_array(memo, k, AG);
      }
      if (gn < al[l2] - paraGap) {
        gn = al[l2] - paraGap;
        set_half_byte_array(memo, k, AG | VH);
      }
      /* position aligned */
      if (al[l2-1] > gl[l2 - 1]) {
        an += al[l2 - 1];
        or_half_byte_array(memo, k, AA);
      } else {
        an += gl[l2  -1];
      }
      /* moving up */
      al[l2 - 1] = ad;
      gl[l2 - 1] = gd;
      ad = an;
      gd = gn;
      Loctest(ad,k);
    }
    al[l2 - 1] = ad;
    gl[l2 - 1] = gd;
#ifdef ALIGLOLOC
    Glotest(ad, gd, k-1);
#elif defined(ALIGLOB)
    Glotest(ad - paraExt * (float)(lp1 - 1 - l1),
            gd - paraExt * (float)(lp1 - 1 - l1),
            k - 1)
#endif
    }
#if (defined(ALIGLOB) || defined(ALIGLOLOC))
  for (l2 = lp2 - 2; l2 >= 0; l2--) {
    Glotest(al[l2] - paraExt * (float)(lp2 - 1 - l2),
            gl[l2] - paraExt * (float)(lp2 - 1 - l2),
            k - (lp2 - l2))
  }
#endif

  /* backward */
  if (pathp != NULL) {
    stat = 1;
    pathn = 0;
    if (be < 0) {
      be = -be;
      stat = 0;
    }
    l1 = be / lp2;
    l2 = be % lp2;
    for (k = l1 + 1; k < lp1; ++k) {
      pathp[k] = -1;
    }
    if (be) {
      for ( ; ; ) {
        if (stat) {
#ifndef ALIGLOB
          if (get_half_byte_array(memo, be) & ZE) {
            /* TEST fprintf(stderr,"ZERO: %d,%d\n",be/lp2,be%lp2); */
            break;
          }
#endif
          pathp[l1] = l2;
          pathn += 1;
          pathv[l1] = scale_value(m[l1 * lp2 + l2]); /* save pathv */
          /* SaveScore(be,pathv[l1]); */
          stat = get_half_byte_array(memo, be) & AA;
          if (--l1 < 0) break;
          if (--l2 < 0) break;
          be -= lp2 + 1;
        } else {
          pathp[l1] = -1;
          stat = get_half_byte_array(memo, be) & AG;
          if (get_half_byte_array(memo, be) & VH) {
            if(--l1 < 0) break;
            be -= lp2;
          } else {
            if(--l2 < 0) break;
            be -= 1;
          }
        }
      }
    }
    for (; l1 >= 0; --l1) {
      pathp[l1] = -1;
    }
  }

  return (float)score;
}

inline void coor_assign(coor_t a, const coor_t b) {
  a[0] = b[0];
  a[1] = b[1];
  a[2] = b[2];
}

inline int coor_equal(const coor_t a, const coor_t b) {
  return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}

inline int shift(int p) {
  int a = p + (LEN - SEGSIZE2) / 2;
  if (a < 0) a = 0;
  return a;
}

inline int max(int a, int b) { return a > b ? a : b; }

inline void address(uintptr_t ss, uintptr_t se, uintptr_t *as, uintptr_t *ae, uint *s, int *o) {
  *as = ss - (ss % 16);
  *ae = (se + 16) - (se % 16);
  *o  = ss - *as;
  *s  = *ae - *as;

  ASSERT(MODULO(*as, 4) == 0);
  ASSERT(MODULO(*ae, 4) == 0);
  ASSERT(0 <= *o && *o < 16);
}

void get_num(short *d, short *s, uint n) {
  uintptr_t addr_start, addr_end;
  uint      size, i;
  int       off;
  char      *tmp = (char*)memo;

  address((uintptr_t)s, (uintptr_t)(s + n), &addr_start, &addr_end, &size, &off);

  mfc_get(tmp, addr_start, size, tag, 0, 0);
  waitag(tag);

  for (i = 0; i < n; ++i) {
    d[i] = ((short*)(tmp + off))[i];
  }
}

void get_seq(char *d, char *s, uint n) {
  uintptr_t addr_start, addr_end;
  uint      size, i;
  int       off;
  char      *tmp = (char*)memo;

  address((uintptr_t)s, (uintptr_t)(s + n), &addr_start, &addr_end, &size, &off);

  mfc_get(tmp, addr_start, size, tag, 0, 0);
  waitag(tag);

  for (i = 0; i < n; ++i) {
    d[i] = (tmp + off)[i];
  }
}

void get_coor(coor_t *d, coor_t *s, uint n) {
  uintptr_t addr_start, addr_end;
  uint      size, i;
  int       off;
  char      *tmp = (char*)memo;

  address((uintptr_t)s, (uintptr_t)(s + n), &addr_start, &addr_end, &size, &off);

  mfc_get(tmp, addr_start, size, tag, 0, 0);
  waitag(tag);

  for (i = 0; i < n; ++i) {
    coor_assign(d[i], ((coor_t*)(tmp + off))[i]);
  }
}

void get_protein(spe_protein_t *pd, prot_t ps, int p, int segsize) {
  int a, b, as, bs;

  a = (p + pd->shift) + (LEN - segsize) / 2;
  if (a < 0) a = 0;

  b = a + segsize;
  if (b > ps.len) b = ps.len;

  as = a - pd->shift;
  bs = b - pd->shift;

  if (as < pd->a) {
    get_num(pd->num + as, ps.num + a, pd->a - as);
    get_coor(pd->coor + as, ps.coor + a, pd->a - as);
    get_coor(pd->c + as, ps.c + a, pd->a - as);
    pd->a = as;
  }
  if (pd->b < bs) {
    get_num(pd->num + pd->b, ps.num + pd->b + pd->shift, b - pd->b - pd->shift);
    get_coor(pd->coor + pd->b, ps.coor + pd->b + pd->shift, b - pd->b - pd->shift);
    get_coor(pd->c + pd->b, ps.c + pd->b + pd->shift, b - pd->b - pd->shift);
    pd->b = bs;
  }
}

void get13(const int order) {
  int       p1 = order / (params.protein2.len - LEN + 1);
  int       p2 = order % (params.protein2.len - LEN + 1);
  int       a1 = p1;
  int       a2 = p2;
  uintptr_t as, ae;
  uint      size;

  address((uintptr_t)(params.protein1.num + a1),
          (uintptr_t)(params.protein1.num + a1 + LEN),
          &as, &ae, &size, &off_num1);
  mfc_get(plen1_tmp.num, as, size, tags_len[0], 0, 0);

  address((uintptr_t)(params.protein1.seq + a1),
          (uintptr_t)(params.protein1.seq + a1 + LEN),
          &as, &ae, &size, &off_seq1);
  mfc_get(plen1_tmp.seq, as, size, tags_len[1], 0, 0);

  address((uintptr_t)(params.protein1.coor + a1),
          (uintptr_t)(params.protein1.coor + a1 + LEN),
          &as, &ae, &size, &off_coor1);
  mfc_get(plen1_tmp.coor, as, size, tags_len[2], 0, 0);

  address((uintptr_t)(params.protein1.c + a1),
          (uintptr_t)(params.protein1.c + a1 + LEN),
          &as, &ae, &size, &off_c1);
  mfc_get(plen1_tmp.c, as, size, tags_len[3], 0, 0);

  address((uintptr_t)(params.protein2.num + a2),
          (uintptr_t)(params.protein2.num + a2 + LEN),
          &as, &ae, &size, &off_num2);
  mfc_get(plen2_tmp.num, as, size, tags_len[4], 0, 0);

  address((uintptr_t)(params.protein2.seq + a2),
          (uintptr_t)(params.protein2.seq + a2 + LEN),
          &as, &ae, &size, &off_seq2);
  mfc_get(plen2_tmp.seq, as, size, tags_len[5], 0, 0);

  address((uintptr_t)(params.protein2.coor + a2),
          (uintptr_t)(params.protein2.coor + a2 + LEN),
          &as, &ae, &size, &off_coor2);
  mfc_get(plen2_tmp.coor, as, size, tags_len[6], 0, 0);

  address((uintptr_t)(params.protein2.c + a2),
          (uintptr_t)(params.protein2.c + a2 + LEN),
          &as, &ae, &size, &off_c2);
  mfc_get(plen2_tmp.c, as, size, tags_len[7], 0, 0);
}

void copy13(spe_protein_t *pd, protein_len_tmp_t ps, int p,
            int off_num, int off_seq, int off_coor, int off_c) {
  int a = p, b = a + LEN, i;

  for (i = 0; i < LEN; ++i) {
    pd->num[p + i] = ((short*)(ps.num + off_num))[i];
  }
  for (i = 0; i < LEN; ++i) {
    pd->seq[p + i] = ((char*)(ps.seq + off_seq))[i];
  }
  for (i = 0; i < LEN; ++i) {
    coor_assign(pd->coor[p + i], ((coor_t*)(ps.coor + off_coor))[i]);
  }
  for (i = 0; i < LEN; ++i) {
    coor_assign(pd->c[p + i], ((coor_t*)(ps.c + off_c))[i]);
  }

  pd->a = a;
  pd->b = b;
}

uint f(const int order1, const int order2, spe_worker_output_t *output) {
  uint    mask = 0;
  int     p1 = order1 / (params.protein2.len - LEN + 1);
  int     p2 = order1 % (params.protein2.len - LEN + 1);

  int     a, a1, a2, b, b1, b2, aa, i, j, k, l, o, n;
  int     len, len1, len2, result = 0;

  int     mat1, matbest1 = 0, matbest2 = 0;
  float   ali1;

  coor_t  *x, *y;

  output->hitn = output->hitn1 = output->hitn2 = output->hitn3 = 0;
  output->alibest = -1.0;

  prot1.a = -1;
  prot1.b = -1;
  prot1.shift = shift(p1);
  prot2.a = -1;
  prot2.b = -1;
  prot2.shift = shift(p2);

  p1 = p1 - prot1.shift;
  p2 = p2 - prot2.shift;

  a1 = p1 + (LEN - SEGSIZE) / 2; /* start of segment of protein 1 */
  b1 = a1 + SEGSIZE;             /* end of segment of protein 1 */
  a2 = p2 + (LEN - SEGSIZE) / 2; /* start of segment of protein 2 */
  b2 = a2 + SEGSIZE;             /* end of segment of protein 2 */

  a = b = 0;

  if (a1 < -a) {
    a = -a1;
  }
  if (a2 < -a) {
    a = -a2;
  }

  if (b1 - (params.protein1.len - prot1.shift) > b) {
    b = b1 - (params.protein1.len - prot1.shift);
  }
  if (b2 - (params.protein2.len - prot2.shift) > b) {
    b = b2 - (params.protein2.len - prot2.shift);
  }

  a1 += a;        /* start in protin (query) */
  a2 += a;        /* end in protin (query) */

  b1 -= b;        /* start in template */
  b2 -= b;        /* end in template ; size = b1 - a1 = b2 - a2 */

  aa = p1 + LEN;  /* start of part to the right from hit */

  len = b1 - a1;
  len1 = p1 - a1;
  len2 = b1 - aa;

  ASSERT(0 <= len && len <= SEGSIZE2);
  ASSERT(0 <= len1 && len1 <= SEGSIZE2);
  ASSERT(0 <= len2 && len2 <= SEGSIZE2);

TICS(3);
  for (i = 0; i < 8; ++i) {
    mask |= (1 << tags_len[i]);
  }
  mfc_write_tag_mask(mask);
  mfc_read_tag_status_all();

  copy13(&prot1, plen1_tmp, p1, off_num1, off_seq1, off_coor1, off_c1);
  copy13(&prot2, plen2_tmp, p2, off_num2, off_seq2, off_coor2, off_c2);

  if (order2 >= 0) {
    get13(order2);
  }
TICE(3);

  x = prot1.coor + p1;
  y = prot2.coor + p2;

  ASSERT(0 <= p1 && p1 + LEN < SEGSIZE2);
  ASSERT(0 <= p2 && p2 + LEN < SEGSIZE2);

  for (l = 0, i = LEN - 1; i >= 0; --i) {
    ASSERT(0 <= i && i < LEN);
     if (prot1.num[p1 + i] && prot2.num[p2 + i]) {
      movew[i] = 1.0;
      l++;
    } else {
      movew[i] = 0.0;
    }
  }
  if (l < LEN / 2) {
    return 0;
  }

TICS(5);
  mat1 = (int)wmove((const coor_t*)y, x, movew, LEN);
TICE(5);

  if (prot1.seq[p1] == prot1.seq[p2]) {
    if (mat1 >= RMSDSEQ) {
      return 0;
    }
  } else {
#ifdef PAIRFULL
    if (mat1 >= RMSDSEQ) {
#else
    if (mat1 >= RMSD) {
#endif
      return 0;
    }
  }

  output->hitn += 1;

TICS(4);
  get_protein(&prot1, params.protein1, p1, SEGSIZE);
  get_protein(&prot2, params.protein2, p2, SEGSIZE);
TICE(4);

  ASSERT(0 <= a1 && a1 + len <= SEGSIZE2);
  ASSERT(0 <= a2 && a2 + len <= SEGSIZE2);

TICS(10);
  /* moving the segment */
  x = prot1.coor + a1;
  y = prot2.coor + a2;
  for (n = a1, i = len-1; i >= 0; --i, ++n) {
    ASSERT(0 <= i && i < len);
    ASSERT(a1 <= n && n < a1 + len);
    if (prot1.num[a1 + i]) {
      float _x[3] = {0.0, 0.0, 0.0};
      for (j = 2; j >= 0; --j) {
        xnew[i][j] = x[i][j] - xc[j];
      }
      for (j = 2; j >= 0; --j) {
        for(k = 0; k < 3; ++k) {
          _x[j] += U[j][k] * xnew[i][k];
        }
      }
      for (j = 2; j >= 0; --j) {
        xnew[i][j] = _x[j] + yc[j];
      }
      for (j = 0; j < 3; ++j) {
        cnew[i][j] = 0.f;
      }
      for (j = 0; j < 3; ++j) {
        for (k = 0; k < 3; ++k) {
          cnew[i][j] += prot1.c[n][k] * U[j][k];
        }
      }
    }
  }
TICE(10);

TICS(6);
  /* alignment score */
  /* calculate only parts of the (len * len) matrix */
  mat1 = l +
  matrix1w(m,
          (const coor_t*)xnew,           (const coor_t*)cnew,                     prot1.num + a1,           len1,
          (const coor_t*)y,              (const coor_t*)(prot2.c + a2),           prot2.num + a2,           len1)
           +
  matrix1w(m + SEGSIZE * (aa - a1) + aa - a1,
          (const coor_t*)xnew + aa - a1, (const coor_t*)(cnew + aa - a1),         prot1.num + aa,           len2,
          (const coor_t*)y    + aa - a1, (const coor_t*)(prot2.c + a2 + aa - a1), prot2.num + a2 + aa - a1, len2);
TICE(6);

  if (mat1 >= SEGHITS1) {
    output->hitn1 += 1;

TICS(7);
    ali1 = LEN +
           align1(m, len1, len1, pathp) +
           align1(m + SEGSIZE * (aa - a1) + aa - a1, len2, len2, pathp + aa - a1);
TICE(7);

    if ((int)ali1 > SEGHITS1) { /* change if changing gap parameters */
      float   ali2 = 0.0;
      int     mat2 = 0;
      int     pathm = 0;

      output->hitn2 += 1;

TICS(10);
      /* optimizing move, extract close atoms */
      for (k = 0, i = 0; i < len1; ++i) { /* part 1 */
        if (pathp[i] >= 0) {
          ASSERT(0 <= i && i < len);
          ASSERT(0 <= pathp[i] && pathp[i] < len);
          for (j = 2; j >= 0; --j) {
            movew[k] = 1.0;
            xnew[k][j] = x[i][j];
            ynew[k][j] = y[pathp[i]][j];
          }
          ++k;
        }
      }
      for(; i < aa - a1; ++i) { /* fixed alignment of hit segment */
        ASSERT(0 <= i && i < len);
        if (prot1.num[a1 + i] && prot2.num[a2 + i]){
          for(j = 2; j >= 0; --j){
            movew[k] = 1.0;
            xnew[k][j] = x[i][j];
            ynew[k][j] = y[i][j];
          }
          ++k;
        }
      }
      for (; i < len; ++i) { /* part 2 */
        if (pathp[i] >= 0) {
          ASSERT(0 <= i && i < len);
          ASSERT(0 <= pathp[i] + aa - a1 && pathp[i] + aa - a1 < len);
          for (j = 2; j >= 0; --j) {
            movew[k] = 1.0;
            xnew[k][j] = x[i][j];
            ynew[k][j] = y[pathp[i] + aa - a1][j];
          }
          ++k;
        }
      }
TICE(10);

      for (o = 1; ; ++o) {
        const int segsize = (int)(SEGSIZE + (SEGSIZE2 - SEGSIZE) * (float)o / OPTIMIZE);

        a1 = p1 + (LEN - segsize) / 2;  /* start of segment of protein 1 */
        b1 = a1 + segsize;              /* end of segment of protein 1 */

        a2 = p2 + (LEN - segsize) / 2;  /* start of segment of protein 2 */
        b2 = a2 + segsize;              /* end of segment of protein 2 */

        if (a1 < 0) {
          a1 = 0;
        }
        if (a2 < 0) {
          a2 = 0;
        }

        if (b1 > (params.protein1.len - prot1.shift)) {
          b1 = (params.protein1.len - prot1.shift);
        }
        if (b2 > (params.protein2.len - prot2.shift)) {
          b2 = (params.protein2.len - prot2.shift);
        }

        len1 = b1 - a1;
        len2 = b2 - a2;

        ASSERT(0 <= len1 && len1 <= segsize);
        ASSERT(0 <= len2 && len2 <= segsize);

        ASSERT(0 <= a1 && a1 + segsize <= SEGSIZE2);

TICS(4);
        get_protein(&prot1, params.protein1, p1, segsize);
        get_protein(&prot2, params.protein2, p2, segsize);
TICE(4);

        x = prot1.coor + a1;
        y = prot2.coor + a2;

TICS(5);
        wmove((const coor_t*)ynew, xnew, movew, k);
TICE(5);

TICS(10);
        /* move again */
        for (n = a1, i = 0; i < len1; ++i, ++n) {
          float _x[3] = {0.0, 0.0, 0.0};
          ASSERT(0 <= i && i < segsize);
          ASSERT(a1 <= n && n < a1 + segsize);
          for (j = 2; j >= 0; --j) {
            xnew[i][j] = x[i][j] - xc[j];
          }
          for (j = 2; j >= 0; --j) {
            for (k = 0; k < 3; ++k) {
              _x[j] += U[j][k] * xnew[i][k];
            }
          }
          for (j = 2; j >= 0; --j) {
            xnew[i][j] = _x[j] + yc[j];
          }
          for (j = 0; j < 3; ++j) {
            cnew[i][j] = 0.f;
          }
          for (j = 0; j < 3; ++j) {
            for (k = 0; k < 3; ++k) {
              cnew[i][j] += prot1.c[n][k] * U[j][k];
            }
          }
        }
TICE(10);

TICS(8);
        mat2 = matrix2w(m,
                       (const coor_t*)xnew, (const coor_t*)cnew,           prot1.num + a1, len1,
                       (const coor_t*)y   , (const coor_t*)(prot2.c + a2), prot2.num + a2, len2);
TICE(8);

TICS(9);
        ali2 = align2(m, len1, len2, pathp);
TICE(9);

        if (o >= OPTIMIZE || pathn < SEGHITS2) {
          break;
        }

        output->hitn3 += 1;

TICS(10);
        for (k = 0, i = 0; i < len1; ++i) {
          if (pathp[i] >= 0) {
            ASSERT(0 <= i && i < segsize);
            ASSERT(0 <= pathp[i] && pathp[i] < segsize);
            for (j = 2; j >= 0; --j) {
              movew[k] = pathv[i];
              xnew[k][j] = x[i][j];
              ynew[k][j] = y[pathp[i]][j];
            }
            ++k;
          }
        }
TICE(10);
      }

      result = 1;

      output->s = a1 + prot1.shift;
      output->t = (a1 + prot1.shift) + len1;

      for (i = len1 - 1; i >= 0; --i) {
        if (pathp[i] >= 0) {
          int j = (i + (a1 + prot1.shift)) * params.protein2.len + (pathp[i] + (a2 + prot2.shift));
          if (pathv[i] >= 0.0) {
            pathm++;
          }
          output->pfull[i] = j;
        } else {
          output->pfull[i] = -1;
        }
      }

      if (ali2 > output->alibest) {
        for (i = 0; i < len1; ++i) {
          if (pathp[i] < 0) {
            output->pbest[i] = -1;
          } else {
            output->pbest[i] = pathp[i] + (a2 + prot2.shift);
          }
        }

        for (i = 0; i < 3; ++i) {
          for(j = 0; j < 3; ++j) {
            output->Ubest[i][j] = (float)U[i][j];
          }
          output->xbest[i] = xc[i];
          output->ybest[i] = yc[i];
        }

        output->pathn2 = pathn;
        output->pathn3 = pathm;
        output->alibest = ali2;

        matbest1 = mat1;
        matbest2 = mat2;
      }
    }
  }

  output->result = result;

  return result;
}

inline int p(int msg) {
  int p1 = msg / (params.protein2.len - LEN + 1);
  int p2 = msg % (params.protein2.len - LEN + 1);
  return ((p1 + 43 * p2) % (params.protein1.len - LEN + 1)) * (params.protein2.len - LEN + 1) + p2;
}

uint32_t out_status __attribute__((aligned(128)));

/*
 * PROGRAM GŁÓWNY:
 */

int main(unsigned long long speid  __attribute__ ((unused)),
         unsigned long long argp   __attribute__ ((unused)),
         unsigned long long envp   __attribute__ ((unused)))
{
  uint  i, j;
  int   msg1, msg2;

#ifdef TIMER
  spu_slih_register(MFC_DECREMENTER_EVENT, spu_clock_slih);
  spu_clock_start();
  timebase = spu_timebase();
#endif

  TICS(0);

  /* Rezerwujemy tag ID. */
  tag = mfc_tag_reserve();
  if (tag == MFC_TAG_INVALID) {
    fatal("SPE: ERROR can't allocate tag ID.\n");
    return -1;
  }

  /* Rezerwujemy tag_in. */
  tag_in = mfc_tag_reserve();
  if (tag_in == MFC_TAG_INVALID) {
    fatal("SPE: ERROR can't allocate tag ID.\n");
    return -1;
  }

  /* Rezerwujemy tags_len. */
  for (i = 0; i < 8; ++i) {
    tags_len[i] = mfc_tag_reserve();
    if (tags_len[i] == MFC_TAG_INVALID) {
      fatal("SPE: ERROR can't allocate tag ID.\n");
      return -1;
    }
  }

  /* Rezerwujemy tag_out. */
  tag_out = mfc_tag_reserve();
  if (tag_out == MFC_TAG_INVALID) {
    fatal("SPE: ERROR can't allocate tag ID.\n");
    return -1;
  }

TICS(1);
  /* Pobieramy parametry */
  mfc_get(&params, argp, sizeof(spe_worker_parameters_t), tag_in, 0, 0);
  waitag(tag_in);
TICE(1);

  msg1 = params.from;
  get13(p(msg1));

  for (i = params.from, j = 0; i < params.to; ++i) {
    msg2 = (i < params.to) ? (int)(i + 1) : -1;

TICS(11);
    f(p(msg1), p(msg2) , &output);
TICE(11);

TICS(2);
    if (output.result) {
      mfc_put(&output,
              (uint32_t)(&params.output[j]),
              sizeof(spe_worker_output_t),
              tag_out, 0, 0);
      waitag(tag_out);
      j += 1;
    }
TICE(2);

    msg1 = msg2;
  }

  out_status = j + 1;

TICS(12);
  mfc_put(&out_status, (uint32_t)(params.output_status), sizeof(uint32_t), tag_out, 0, 0);
TICE(12);

  /* Zwalniamy tag ID. */
  mfc_tag_release(tag);

  /* Zwalniamy tag_in. */
  mfc_tag_release(tag_in);

  /* Zwalniamy tags_len. */
  for (i = 0; i < 8; ++i) {
    mfc_tag_release(tags_len[i]);
  }

  /* Zwalniamy tag_out. */
  mfc_tag_release(tag_out);

  TICE(0);

  TICP(0, "total_spe");
  TICP(1, "msg");
  TICP(2, "out");
  TICP(3, "in");
  TICP(4, "rin");
  TICP(5, "wmove");
  TICP(6, "matrix1w");
  TICP(7, "align1");
  TICP(8, "matrix2w");
  TICP(9, "align2");
  TICP(10, "moving");
  TICP(11, "f");
  TICP(12, "status");

  return 0;
}

