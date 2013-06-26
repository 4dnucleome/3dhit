/*
 * 3dhit.h
 *
 * autor: Łukasz Bieniasz-Krzywiec
 *
 */

#ifndef _3DHIT_H_
#define _3DHIT_H_

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "modules/err.h"

/*
 * DEKLARACJE STAŁYCH:
 */

#define LEN         13
#define RMSD        2.0   /* max rmsd */
#define ENDS        2.0   /* max distance of ends */
#define RMSDSEQ     3.0   /* 3.0 max rmsd */
#define ENDSSEQ     3.0   /* 3.0 max distance of ends */
#define SEGSIZE     99    /* superimposed segment */
#define SEGSIZE2    259   /* superimposed segment */
#define INCREMENT

#define DISTCUTOFF  8.0   /* cutoff for contact analysis, not used */
#define DISTSIGMA   3.0   /* distanse scaling factor */
#define MINSEQDIFF  6     /* min sequence separation */
#define MOVEMINSEG  4     /* initial segment to calculate superposition */
#define MOVEMINHIT  8.0   /* minimum superposition score to continue */
#define MOVESTEPS   4     /* minimum superposition score to continue */

#define MINCOS      0.707 /* acos(0.707)=45 => dAngle=90 */
#define SEGSIGMA    3.0
#define MVGSIGMA    9.0
#define MVGALIGN    10
#define MVGPRGAP    1.02
#define MVGPREXT    0.02

                          /* FAST - OK  */
#define SEGHITS1    18    /* 25     18  */
#define SEGRMSD1    3.0   /* 2.5    3.0 */
#define SEGHITS2    25    /* 30     25  */
#define SEGRMSD2    8.0   /* 8.0    8.0 */
#define SEGRMSD0    3.0   /* 2.5    3.0 */
#define MOVEW
#define TWOSTEP
#define OPTIMIZE    2     /* optimization steps */
#ifdef PAIR
#endif
#define ALIGNMENT
#define ALIMIN      0.0   /* min score for a position in final alignment */

#define ALIGLOB_NO
#define paraZero    0.0
#define paraGap     1.02
#define paraExt     0.02

#define PAIRSTEP    1

/* Maksymalna długość linii w pliku wejściowym. */
#define MAX_LINE_LENGTH 1000

/* Kody zleceń: */
#define MSG_STOP    0     /* zakończ pracę */


/*
 * DEKLARACJE TYPÓW:
 */

/* unsigned char */
typedef unsigned char uchar;

/* unsigned int */
typedef unsigned int uint;

/* unsigned long int */
typedef unsigned long int luint;

/* long long */
typedef long long LL;

/* unsigned long long */
typedef unsigned long long ULL;

/* Współrzędne. */
typedef float coor_t[3];

/* Białko. */
typedef struct prot_t_ {
  int     len;
  char    *name;
  short   *num;
  char    *seq;
  coor_t  *coor;
  coor_t  *c;
} prot_t;

/* Kontakty. */
typedef struct contacts_t_ {
  float   *map;
} contacts_t;

/* Wyjście wątku obliczeniowego. */
typedef struct spe_worker_output_t_ {
  int           result;
  int           hitn, hitn1, hitn2, hitn3;
  float         alibest;
  float         Ubest[3][3], xbest[3], ybest[3];
  int           pathn2, pathn3;
  int           s, t, pbest[SEGSIZE2], pfull[SEGSIZE2];
  uchar         pad[4 + 0 * 8];
} spe_worker_output_t __attribute__((aligned(128)));

/* Parametry wątku obliczeniowego. */
typedef struct spe_worker_parameters_t_ {
  prot_t              protein1;
  prot_t              protein2;
  spe_worker_output_t *output;
  volatile uint32_t   *output_status;
  uint                from;
  uint                to;
} spe_worker_parameters_t __attribute__((aligned(64)));


/*
 * DEKLARACJE ZMIENNYCH GLOBALNYCH:
 */


/*
 * DEKLARACJE PODPROGRAMÓW:
 */

/*
 * Funkcje z pliku modules/alignment.c
 */

float align2_f(float*, const int, const int, int*);
void calcvect(coor_t*, const coor_t*, int);
int calc_within(int*, float, int*, float, int*, int, coor_t*, coor_t*, float[3][3], float[3], float[3]);
float mvg_3Dscore(int*, float[3][3], coor_t*, int, float[3], coor_t*, int, float[3]);
void print_alin1(FILE*, int, int, char*, char*, int*);
void ascending_sort(float*, int*, int*, int, int);
void print_alin(FILE*, const int, const int, const char*, const char*, const int*, const int, const int, int);
void print_whole_pdb1(FILE*, char**, char**, float[3][3], float[3], float[3]);
void print_whole_pdb2(FILE*, char**, char**, float[3][3], float[3], float[3]);

/*
 * Funkcje z pliku modules/contacts.c
 */

void setcontacts(contacts_t*, coor_t*, int, short*);
void contactoverlap(contacts_t*, contacts_t*, int, float*, float*);

/*
 * Funkcje z pliku modules/input.c
 */

prot_t *readpdb(char*, char ***xpdb);

/*
 * Funkcje z pliku modules/memory.c
 */

void *mem(void*, size_t);


/*
 * PRZYDATNE MAKRA
 */

/* MODULO(a,w) = a mod 2^w */
#define MODULO(a, w) ((a) & ((1 << (w)) - 1))

/*
 * Macros for rounding input value to the next lower multiple of either
 * 16, 32 or 128.
 */
#define floor16(value)  ((value) - ((value) & ~15))
#define floor32(value)  ((value) - ((value) & ~31))
#define floor128(value) ((value) - ((value) & ~127))

/*
 * Macros for rounding input value to the next higher multiple of either
 * 16, 32 or 128.
 */
#define ceil16(value)   (((value) + 15) & ~15)
#define ceil32(value)   (((value) + 31) & ~31)
#define ceil128(value)  (((value) + 127) & ~127)


/*
 * MACROS TO HELP DEBUGGING:
 */

#undef PDEBUG
#undef ASSERT

#ifdef DEBUG
# ifdef __KERNEL__
    /* This one if debugging is on, and kernel space */
#   define PDEBUG(fmt, args...) printk(KERN_CRIT ": " fmt, ## args)
# else
    /* This one for user space */
#   define PDEBUG(fmt, args...) fprintf(stderr, fmt, ## args)
# endif

# define ASSERT(x) do {\
    if (!(x)) {\
      PDEBUG("Assertion failed at line %d in file %s: " #x "\n", __LINE__, __FILE__);\
      exit(1);\
    }\
  } while(0)
#else
# define PDEBUG(fmt, args...) /* not debugging: nothing */
# define ASSERT(x)
#endif

#endif /* _3DHIT_H_ */
