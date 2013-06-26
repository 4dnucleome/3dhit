#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <time.h>
#define LEN 13 
//#define SEQFILT
#define RMSD 2.0 /* max rmsd */
#define ENDS 2.0 /* max distance of ends */
#define RMSDSEQ 3.0 /* 3.0 max rmsd */
#define ENDSSEQ 3.0 /* 3.0 max distance of ends */
#define SEGSIZE 99 /* superimposed segment */
#define SEGSIZE2 259 /* superimposed segment */
#define INCREMENT

#define DISTCUTOFF 8.0 /* cutoff for contact analysis, not used */
#define DISTSIGMA  3.0 /* distanse scaling factor */ 
#define MINSEQDIFF 6 /* min sequence separation */ 
#define MOVEMINSEG 4 /* initial segment to calculate superposition */
#define MOVEMINHIT 8.0 /* minimum superposition score to continue */
#define MOVESTEPS 4 /* minimum superposition score to continue */

//added by mvg 
#define MINCOS   0.707  // acos(0.707)=45 => dAngle=90
#define SEGSIGMA 3.0
#define MVGSIGMA 9.0
#define MVGALIGN 10
#define MVGPRGAP 1.02
#define MVGPREXT 0.02
//end by mvg
                     /*FAST - OK*/
#define SEGHITS1 18  /* 25  18  */
#define SEGRMSD1 3.0 /* 2.5 3.0 */
#define SEGHITS2 25  /* 30  25  */
#define SEGRMSD2 8.0 /* 8.0 8.0 */
#define SEGRMSD0 3.0 /* 2.5 3.0 */
#define MOVEW
#define TWOSTEP
#define OPTIMIZE 2 /* optimization steps */
#ifdef PAIR
#endif
#define ALIGNMENT
#define ALIMIN   0.0 /* min score for a position in final alignment */

#define ALIGLOB_NO
#define paraZero 0.0
#define paraGap  1.02
#define paraExt  0.02

#define MAX_THREADS 16 /* by LBK */

// macros
#define DIST2(a,b) ((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))
typedef float coor_t[3];
typedef short pos_t[2];
typedef struct prot_s{
	char	*name;
	int	len;
	char	*seq;
	short	*num;
	coor_t	*coor;} prot_t; // in DATA/prot.c
typedef struct neib_s{
   float   d;
	short	prot1;
	short	prot2;
	short	res1;
	short	res2;} neib_t;
typedef struct cluster_s{
	char	c; // amino acid in the center
	float	d; // distance between ends
	int	size; // size of the cluster
	char	*name; // source of the fragment protein:position
	coor_t	coor[LEN];} cluster_t; // in DATA/clusters.c
typedef struct location_s{ // order of res and prot is important for sorting
	short	res;
	short	prot;} location_t;
typedef struct hash_s{
	int	size;
	location_t *loc;} hash_t;
typedef struct contacts_s{
	//short	*neib;
	//short	*neibn;
	float 	*map;} contacts_t;
        
/* by LBK: */
        
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

/* Wyjście wątku obliczeniowego. */
typedef struct spe_worker_output_t_ {
  int           result;
  int           hitn, hitn1, hitn2, hitn3;
  float         alibest;
  float         Ubest[3][3], xbest[3], ybest[3];
  int           pathn2, pathn3;
  int           s, t, pbest[SEGSIZE2], pfull[SEGSIZE2];
} spe_worker_output_t;

void* mem(void* p,long size);
float rmsd(coor_t *a,coor_t *b);
float dist(const float *x,const float *y);
double wmove(const coor_t* y_old,coor_t* x_old,const float* movew,const int n, float[3], float[3], double[3][3]);
//int matrix1(float *m,const coor_t *p1,const int lp1,const coor_t *p2,const int lp2);
//int matrix2(float *m,const coor_t *p1,const int lp1,const coor_t *p2,const int lp2);
int matrix1(float *m, //changed by mvg
		            const coor_t *p1,const coor_t *c1,const int lp1,
			                const coor_t *p2,const coor_t *c2,const int lp2);
int matrix2(float *m, //changed by mvg
		            const coor_t *p1,const coor_t *c1,const int lp1,
			                const coor_t *p2,const coor_t *c2,const int lp2);
float align1(float *m,const int lp1,const int lp2,int* pathp);
float align2(float *m,const int lp1,const int lp2,int* pathp, int *pathn, float **pathv, int *, char **, float **, float **);
void print_ali(FILE* file,const int qlen,const int len2,const char* seq1,const char* seq2,const int* pathp);
void print_alin(FILE* file,const int qlen,const int len2,const char* seq1,const char* seq2,const int* pathp,const int offset1,const int offset2);
void clusfromto(char aa,coor_t *co,int* c_from,int* c_to,float ends);
void setcontacts(contacts_t *cont,coor_t *coor,int l,short *seqnum);
void contactoverlap(contacts_t *cont1,contacts_t *cont2,int l,float *conta,float *contt);
//added by mvg
void profprof(float *pp, const float *p1, int n1, const float *p2, int n2);
void calcvect(coor_t * v, const coor_t *c, int n);
void print_pdb(FILE* out,  int *pfull, float U[3][3],
               coor_t* p1, int  n1,    float  xc[3],
                              coor_t* p2, int n2, float yc[3]);
float ligand_pair(FILE* out, char** pdb1, char** pdb2, char* lig);
void print_whole_pdb(FILE* out, char** pdb1, char** pdb2, 
		float U[3][3], float xc[3], float yc[3]);
float mvg_3Dscore(int *pfull, float U[3][3],
		              coor_t* p1, int n1, float xc[3],
	                      coor_t* p2, int n2, float yc[3]);
void ascending_sort(float* a, int* s, int* y, int l, int p);
//end by mvg
