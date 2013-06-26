#define PAIR
#define PAIRFULL

#ifndef PAIRSTEP
#define PAIRSTEP 1
#endif

#include "modul.h"

#include <omp.h>

#undef PROFS
#undef PROFE
#undef PROFP

/* by LBK, begin */
#ifdef PROF
# include <sys/time.h>
# define COUNTERS 16
  struct timeval  __timev1[COUNTERS], __timev2[COUNTERS];
  float           __time[COUNTERS];
# define PROFS(i) gettimeofday(&__timev1[(i)], NULL)
# define PROFE(i) do {\
    gettimeofday(&__timev2[(i)], NULL);\
    __time[(i)] += __timev2[(i)].tv_sec - __timev1[(i)].tv_sec;\
    __time[(i)] += 0.000001 * (__timev2[(i)].tv_usec - __timev1[(i)].tv_usec);\
  } while(0)
# define PROFP(i, s) do {\
    fprintf(stderr, "%s: %f s\n", (s), __time[(i)]);\
    fprintf(stderr, "%s: %f %%\n", (s), __time[(i)] / __time[0] * 100.0);\
  } while(0)
#else
# define PROFS(i)
# define PROFE(i)
# define PROFP(i, s)
#endif
/* by LBK, end */

typedef float float3[3];
typedef float float33[3][3];
typedef double double33[3][3];

short     max_len;
int       clus_a[21];
cluster_t cluster[1];

// functions
prot_t* readpdb(char* pdbfn,char*** xpdb);
int matrix1w(float *m, //changed by mvg, lbk float -> short
       const coor_t *p1,const coor_t *c1,const short *num1,const int lp1,
       const coor_t *p2,const coor_t *c2,const short *num2,const int lp2);
int matrix2w(float *m, //changed by mvg, lbk float -> short
       const coor_t *p1,const coor_t *c1,const short *num1,const int lp1,
       const coor_t *p2,const coor_t *c2,const short *num2,const int lp2);
int calc_within(int* xwin1, float win1, int* xwin2, float win2, int *pfull, int len1,
  coor_t* p1, coor_t* p2, float U[3][3], float xc[3], float yc[3]);
void print_alin1(FILE* file,int len1,int len2,char* seq1,char* seq2,int* pfull);
void print_alin2(FILE* file,int len1,int len2,char* seq1,char* seq2,int* pfull);
void print_whole_pdb1(FILE* out, char** pdb1, char** pdb2, float U[3][3], float xc[3], float yc[3]);
void print_whole_pdb2(FILE* out, char** pdb1, char** pdb2, float U[3][3], float xc[3], float yc[3]);

void f( const int p1, const int p2,
        const prot_t *protin1, const prot_t *protin2,
        const coor_t *c1, const coor_t *c2,
        int *pathn, float **pathv, int *segsize2, char **smpp, float **al, float **gl,
        spe_worker_output_t **output,
        uint *output_status)
{
  uint hitn = 0, hitn1 = 0, hitn2 = 0, hitn3 = 0;

  int l, matbest1=0, matbest2=0; /* by LBK:*/

  float3 yc;
  float3 xc;
  double33 U;
  coor_t *x,*y;
  coor_t xnew[SEGSIZE2];
  coor_t cnew[SEGSIZE2];
  int n;
  int a1=p1+(LEN-SEGSIZE)/2; // start of segment of protein 1
  int b1=a1+SEGSIZE; // end of segment of protein 1
  int a2=p2+(LEN-SEGSIZE)/2; // start of segment of protein 2
  int b2=a2+SEGSIZE; // end of segment of protein 2
  int aa,a=0,b=0,len,len1,len2,i,j,k,o;
  float m[SEGSIZE2*SEGSIZE2];
  float ali1;
  int mat1;
  float movew[SEGSIZE];

  if(a1<-a){ a=-a1; }
  if(a2<-a){ a=-a2; }
  if(b1-protin1->len>b){ b=b1-protin1->len; }
  if(b2-protin2->len>b){ b=b2-protin2->len; }
  a1+=a; // start in protin (query)
  a2+=a; // end in protin (query)
  b1-=b; // start in template
  b2-=b; // end in template ; size=b1-a1=b2-a2;
  aa=p1+LEN; // start of part to the right from hit
  len1=p1-a1;
  len2=b1-aa;
  len=b1-a1;
  /* superimposing LEN residues */
  x=protin1->coor+p1;
  y=protin2->coor+p2;
  /* marking missing residues */
  for(l=0,i=LEN-1;i>=0;i--){
    if(protin1->num[p1+i]&&protin2->num[p2+i]){
      movew[i]=1.0;
      l++;
    }else{
      movew[i]=0.0;
    }
  }
  if(l<LEN/2){ return;} 
  mat1=(int)(wmove((const coor_t*)y,x,movew,LEN, xc, yc, U));
  if(protin1->seq[p1]==protin2->seq[p2]){
    if(mat1>=RMSDSEQ){ return; }
  }else{
#ifdef PAIRFULL
    if(mat1>=RMSDSEQ){ return; }
#else
    if(mat1>=RMSD){ return; }
#endif
  }
  hitn += 1;
  /* moving the segment */
  x=protin1->coor+a1;
  y=protin2->coor+a2;
  // move((const coor_t*)y,x,len);// try this to rotate all len (99) atoms
  // now move all len (99) atoms
  for(n=a1,i=len-1;i>=0;i--,n++){ //changed by mvg
    if(protin1->num[a1+i]){
      float _x[3]={0.0,0.0,0.0};
      for(j=2;j>=0;j--){
        xnew[i][j]=x[i][j]-xc[j];
      }
      for(j=2;j>=0;j--){
        for(k=0;k<3;k++){
          _x[j]+=U[j][k]*xnew[i][k];
        }
      }
      for(j=2;j>=0;j--){
        xnew[i][j]=_x[j]+yc[j];
      }
      for(j=0;j<3;j++) cnew[i][j]=0.f;
      for(j=0;j<3;j++) for(k=0;k<3;k++) cnew[i][j]+=c1[n][k]*U[j][k];
    }
  }
  /* alignment score */
	// calculate only parts of the (len*len) matrix
	mat1=l+ //changed by mvg
	  matrix1w(m,                                                       // part 1 of the matrix
		   (const coor_t*)xnew      ,(const coor_t*) cnew        ,protin1->num+a1      ,len1,
		   (const coor_t*)y         ,(const coor_t*)(c2+a2)      ,protin2->num+a2      ,len1)+
	  matrix1w(m+SEGSIZE*(aa-a1)+aa-a1,                                 // part 2 (4)
		   (const coor_t*)xnew+aa-a1,(const coor_t*)(cnew +aa-a1),protin1->num+aa      ,len2,
		   (const coor_t*)y   +aa-a1,(const coor_t*)(c2+a2+aa-a1),protin2->num+a2+aa-a1,len2);

  if(mat1>=SEGHITS1){
    int pathp[SEGSIZE2];
    float movew[SEGSIZE2]; // use pathv
    hitn1 += 1;
    ali1=LEN+
      align1(m,len1,len1,pathp)+ // align part 1
      align1(m+SEGSIZE*(aa-a1)+aa-a1,len2,len2,pathp+aa-a1);//align part 2
    if((int)ali1>SEGHITS1){ // change if changing gap parameters
      coor_t ynew[SEGSIZE2];
      float ali2=0.0;
      int mat2=0;
      int pathm=0;
      hitn2 += 1;
      // optimizing move, extract close atoms
      for(k=0,i=0;i<len1;i++){ // part 1
        if(pathp[i]>=0){
          for(j=2;j>=0;j--){
            movew[k]=1.0;
            xnew[k][j]=x[i][j];
            ynew[k][j]=y[pathp[i]][j];
          }
          k++;
        }
      }
      for(;i<aa-a1;i++){ // fixed alignment of hit segment
        if(protin1->num[a1+i]&&protin2->num[a2+i]){
          for(j=2;j>=0;j--){
            movew[k]=1.0;
            xnew[k][j]=x[i][j];
            ynew[k][j]=y[i][j];
          }
          k++;
        }
      }
      for(;i<len;i++){ // part 2
        if(pathp[i]>=0){
          for(j=2;j>=0;j--){
            movew[k]=1.0;
            xnew[k][j]=x[i][j];
            ynew[k][j]=y[pathp[i]+aa-a1][j];
          }
          k++;
        }
      }
      for(o=1;;o++){
        const int segsize=(int)(SEGSIZE+(SEGSIZE2-SEGSIZE)*(float)o/OPTIMIZE);
        len1=protin1->len;
        len2=protin2->len;
        a1=p1+(LEN-segsize)/2; // start of segment of protein 1
        b1=a1+segsize; // end of segment of protein 1
        a2=p2+(LEN-segsize)/2; // start of segment of protein 2
        b2=a2+segsize; // end of segment of protein 2
        if(a1<0){a1=0;}
        if(a2<0){a2=0;}
        if(b1>len1){b1=len1;}
        if(b2>len2){b2=len2;}
        len1=b1-a1;
        len2=b2-a2;
        x=protin1->coor+a1;
        y=protin2->coor+a2;
        wmove((const coor_t*)ynew,xnew,movew,k,xc,yc,U);
        // move again
        for(n=a1,i=0;i<len1;i++,n++){ //changed by mvg
          float _x[3]={0.0,0.0,0.0};
          for(j=2;j>=0;j--){
            xnew[i][j]=x[i][j]-xc[j];
          }
          for(j=2;j>=0;j--){
            for(k=0;k<3;k++){
              _x[j]+=U[j][k]*xnew[i][k];
            }
          }
          for(j=2;j>=0;j--){
            xnew[i][j]=_x[j]+yc[j];
          }
          for(j=0;j<3;j++) cnew[i][j]=0.f;
          for(j=0;j<3;j++) for(k=0;k<3;k++) cnew[i][j]+=c1[n][k]*U[j][k];
	}
	mat2=matrix2w(m, //changed by mvg
            (const coor_t*)xnew,(const coor_t*) cnew  ,protin1->num+a1,len1,
            (const coor_t*)y   ,(const coor_t*)(c2+a2),protin2->num+a2,len2);
        ali2=align2(m,len1,len2,pathp, pathn, pathv, segsize2, smpp, al, gl);
        if(o>=OPTIMIZE||*pathn<SEGHITS2){break;}
        hitn3 += 1;
        for(k=0,i=0;i<len1;i++){ // part 1
          if(pathp[i]>=0){//&&(*pathv)[i]>0.0){ // use pathv
            for(j=2;j>=0;j--){
              movew[k]=(*pathv)[i];
              xnew[k][j]=x[i][j];
              ynew[k][j]=y[pathp[i]][j];
            }
            k++;
          }
        }
      }

      *output = (spe_worker_output_t*)mem(*output, sizeof(spe_worker_output_t) * (*output_status  + 1));

      (*output)[*output_status].hitn = hitn;
      (*output)[*output_status].hitn = hitn1;
      (*output)[*output_status].hitn = hitn2;
      (*output)[*output_status].hitn = hitn3;
      (*output)[*output_status].alibest = -1.0;

      (*output)[*output_status].s = a1;
      (*output)[*output_status].t = a1 + len1;

      for (i = len1 - 1; i >= 0; --i) {
        if (pathp[i] >= 0) {
          int j = (i + a1) * protin2->len + (pathp[i] + a2);
          if ((*pathv)[i] >= 0.0) {
            pathm++;
          }
          (*output)[*output_status].pfull[i] = j;
        } else {
          (*output)[*output_status].pfull[i] = -1;
        }
      }

      if (ali2 > (*output)[*output_status].alibest) {
        for (i = 0; i < len1; ++i) {
          if (pathp[i] < 0) {
            (*output)[*output_status].pbest[i] = -1;
          } else {
            (*output)[*output_status].pbest[i] = pathp[i] + a2;
          }
        }

        for (i = 0; i < 3; ++i) {
          for(j = 0; j < 3; ++j) {
            (*output)[*output_status].Ubest[i][j] = (float)U[i][j];
          }
          (*output)[*output_status].xbest[i] = xc[i];
          (*output)[*output_status].ybest[i] = yc[i];
        }

        (*output)[*output_status].pathn2 = *pathn;
        (*output)[*output_status].pathn3 = pathm;
        (*output)[*output_status].alibest = ali2;

        matbest1 = mat1;
        matbest2 = mat2;
      }

      *output_status += 1;

      /* LBK begin
      fprintf(stderr, "(%3d, %3d): hitn = %d, %d, %d, %d, ali = %6f\n",
        p1, p2,
        hitn, hitn1, hitn2, hitn3,
        (*output)[*output_status - 1].alibest);
      LBK end */
    }
  }
}

  int hitn = 0;
  int hitn1 = 0;
  int hitn2 = 0 ;
  int hitn3 = 0;
  int pathn2 = 0;
  int pathn3 = 0;
  float alibest2 = 0.0;

  int   *pfull;
  float *mfull;
  int   *pbest;

  float33 Ubest;
  float3  xbest, ybest;

  spe_worker_output_t **output;
  uint *output_status;

void copyProtein(prot_t *d, prot_t *s) {
  int i, j;
  d->name = s->name;
  d->len = s->len;
  d->seq = (char*)mem(NULL, sizeof(char) * s->len);
  d->num = (short*)mem(NULL, sizeof(short) * s->len);
  d->coor = (coor_t*)mem(NULL, sizeof(coor_t) * s->len);
  for (i = 0; i < s->len; ++i) {
    d->seq[i] = s->seq[i];
    d->num[i] = s->num[i];
    for (j = 0; j < 3; ++j) d->coor[i][j] = s->coor[i][j];
  }
}

void handle_output(int num_threads, int len) {
  int a, b, i, j;

  for (a = 0; a < num_threads; ++a) {
    for (b = 0; b < (int)output_status[a]; ++b) {

      hitn += output[a][b].hitn;
      hitn1 += output[a][b].hitn1;
      hitn2 += output[a][b].hitn2;
      hitn3 += output[a][b].hitn3;

      if (output[a][b].alibest != -1.0) {
        for (i = 0; i < output[a][b].t - output[a][b].s; ++i) {
          if (output[a][b].pfull[i] != -1 && mfull[output[a][b].pfull[i]] < output[a][b].alibest) {
            mfull[output[a][b].pfull[i]] = output[a][b].alibest;
          }
        }
      }

      if (alibest2 < output[a][b].alibest) {
        for (i = 0; i < output[a][b].s; ++i) {
          pbest[i] = -1;
        }
        for (i = output[a][b].s, j = 0; i < output[a][b].t; ++i, ++j) {
          pbest[i] = output[a][b].pbest[j];
        }
        for (i = output[a][b].t; i < len; ++i) {
          pbest[i] = -1;
        }

        for (i = 0; i < 3; ++i) {
          for(j = 0; j < 3; ++j) {
            Ubest[i][j] = output[a][b].Ubest[i][j];
          }
          xbest[i] = output[a][b].xbest[i];
          ybest[i] = output[a][b].ybest[i];
        }

        pathn2 = output[a][b].pathn2;
        pathn3 = output[a][b].pathn3;
        alibest2 = output[a][b].alibest;
      }
    }
  }
}

int main(int argc,char** argv) {

PROFS(0);

  int p1, p2, len1, len2, l;

  float conta=0.0;
  float contt=0.0;

  prot_t *protin1, *protin1t;
  prot_t *protin2, *protin2t;
  char** pdb1; // whole 1st PDB file
  char** pdb2; // whole 2nd PDB file
  coor_t *c1,*c2;

  int pathn[MAX_THREADS];
  float *pathv[MAX_THREADS];
  int segsize2[MAX_THREADS];
  char *smpp[MAX_THREADS];
  float *al[MAX_THREADS];
  float *gl[MAX_THREADS];

  if (argc < 3) {
    fprintf(stderr, "usage: %s pdb_file1 pdb_file2\n", argv[0]);
    exit(-1);
  }

  protin1 = readpdb(argv[1], &pdb1); // read protein
  protin2 = readpdb(argv[2], &pdb2); // read protein

  c1 = (coor_t*)mem(NULL, sizeof(coor_t) * protin1->len);
  c2 = (coor_t*)mem(NULL, sizeof(coor_t) * protin2->len);
  calcvect(c1, (const coor_t*)protin1->coor, protin1->len);
  calcvect(c2, (const coor_t*)protin2->coor, protin2->len);

  len1 = protin1->len;
  len2 = protin2->len;
  max_len = protin2->len;

  if (protin1->len < SEGHITS2 || protin2->len < SEGHITS2) return 0;

  mfull = (float*)mem(NULL, sizeof(float) * (protin1->len * protin2->len + 1));
  pfull = (int*)mem(NULL, sizeof(int) * (protin1->len + 1));
  pbest = (int*)mem(NULL, sizeof(int) * protin1->len);

  for(l = protin1->len * protin2->len; --l >= 0; ) mfull[l] = 0.0;

PROFS(1);

int num_threads = 0, thread_num, ii;

#pragma omp parallel private(thread_num, ii)
{
  #pragma omp single
  {
    num_threads = omp_get_num_threads();
    fprintf(stderr, "Liczba uruchomionych watkow: %d\n", num_threads);

    protin1t = (prot_t*)mem(NULL, sizeof(prot_t) * num_threads);
    protin2t = (prot_t*)mem(NULL, sizeof(prot_t) * num_threads);
    output = (spe_worker_output_t**)mem(NULL, sizeof(spe_worker_output_t*) * num_threads);
    output_status = (uint*)mem(NULL, sizeof(uint) * num_threads);
  }
  thread_num = omp_get_thread_num();

  copyProtein(&protin1t[thread_num], protin1);
  copyProtein(&protin2t[thread_num], protin2);

  pathn[thread_num] = 0;
  pathv[thread_num] = NULL;
  segsize2[thread_num] = 0;
  smpp[thread_num] = NULL;
  al[thread_num] = NULL;
  gl[thread_num] = NULL;

  output[thread_num] = NULL;
  output_status[thread_num] = 0;

  #pragma omp for schedule (dynamic)
  for (ii = 0; ii < (protin1->len - LEN + 1) * (protin2->len - LEN + 1); ++ii) {
    int p1 = ii / (protin2->len - LEN + 1);
    int p2 = ii % (protin2->len - LEN + 1);

    if (!protin1->num[p1]) continue;
#ifdef PAIRALI
    if(protin1->num[p1] != protin2->num[p2]) continue;
#elif defined(SEQFILT)
    if(protin1->seq[p1] != protin2->seq[p2]) continue;
#else
    if (!protin2->num[p2]) continue;
#endif

    thread_num = omp_get_thread_num();

    f(p1, p2,
      &protin1t[thread_num], &protin2t[thread_num],
      (const coor_t *)c1, (const coor_t *)c2,
      &pathn[thread_num], &pathv[thread_num], &segsize2[thread_num], &smpp[thread_num], &al[thread_num], &gl[thread_num],
      &output[thread_num],
      &output_status[thread_num]);
  }
}

PROFS(2);

handle_output(num_threads, protin1->len);

PROFE(2);

PROFE(1);

  // clean alignment matrix
  for(p1=0;p1<protin1->len;p1++){
    float top=ALIMIN;
    for(p2=0;p2<protin2->len;p2++){
      if(mfull[p1*protin2->len+p2]>top){
        top=mfull[p1*protin2->len+p2];
      }
    }
    for(p2=0;p2<protin2->len;p2++){
      if(mfull[p1*protin2->len+p2]<top){
        mfull[p1*protin2->len+p2]=2.0*(-paraGap-paraExt);
      }
      else{
        mfull[p1*protin2->len+p2] /= alibest2;
      }
    }
  }

  float tmp = align2(mfull,protin1->len,protin2->len,pfull, &pathn[0], &pathv[0], &segsize2[0], &smpp[0], &al[0], &gl[0]);

  // align
  if(tmp > 0.0){ // always true !
    // paste aligned coordinates
    coor_t *coor1=(coor_t*)mem(NULL,sizeof(coor_t)*(protin1->len));
    coor_t *coor2=(coor_t*)mem(NULL,sizeof(coor_t)*(protin1->len));
    short *seqnum1=(short int*)mem(NULL,sizeof(short)*(protin1->len));
    short *seqnum2=(short int*)mem(NULL,sizeof(short)*(protin1->len));
    contacts_t contacts1;
    contacts_t contacts2;
    contacts1.map=NULL;
    contacts2.map=NULL;
    for(l=0,p1=0;p1<protin1->len;p1++){
      if(pfull[p1]>=0){
        memcpy(coor1[l],protin1->coor[p1],sizeof(coor_t));
        memcpy(coor2[l],protin2->coor[pfull[p1]],sizeof(coor_t));
        seqnum1[l]=protin1->num[p1];
        seqnum2[l]=protin2->num[pfull[p1]];
        l++;
      }
    }
    setcontacts(&contacts1,coor1,l,seqnum1);
    setcontacts(&contacts2,coor2,l,seqnum2);
    contactoverlap(&contacts1,&contacts2,l,&conta,&contt);
  }

  //print  results
  if(argc==4) {
  FILE* out=fopen(argv[3],"wt");
    fprintf(out, "%s\t%d:%d:%d:%d\t%d\t%.1f\t%d\t%d\t%.1f\t%.1f\t",
            protin2->name,
            hitn,
            hitn1,
            hitn2,
            hitn3,
            protin2->len,
            alibest2,
            pathn3,
            pathn2,
            conta,
            contt);
    print_alin(out,protin1->len,protin2->len,protin1->seq,protin2->seq,pfull,protin1->num[0],protin2->num[0]);
  }

  if(argc==5) {
    FILE* out1=fopen(argv[3],"wt");
    FILE* out2=fopen(argv[4],"wt");
    print_whole_pdb1(out1,pdb1,pdb2,Ubest,xbest,ybest);
    print_whole_pdb1(out2,pdb1,pdb2,Ubest,xbest,ybest);
  } else {
    int i,j,win1,win2,ok,mypfull[protin1->len];
    int i1=strlen(protin1->name)-5; if(i1<0) i1=0;
    int i2=strlen(protin2->name)-5; if(i2<0) i2=0;
    for(i=0;i<protin1->len;i++) mypfull[i]=pfull[i];
    float mvg3Ds = mvg_3Dscore(mypfull,Ubest,protin1->coor,protin1->len,xbest,protin2->coor,protin2->len,ybest);
    calc_within(&win1,3.0,&win2,5.0,pfull,protin1->len,protin1->coor,protin2->coor,Ubest,xbest,ybest);
    for(ok=0;ok<protin1->len;ok++) if(pfull[ok]>=0) break;
    if(ok!=protin1->len) {
      printf("%s\t",&protin2->name[i2]);
      printf("%.1f\t%d\t%.1f\t%d\t%d\t%.1f\t%.1f\t",mvg3Ds,protin2->len,alibest2,win1,win2,conta,contt);
      print_alin1(stdout,protin1->len,protin2->len,protin1->seq,protin2->seq,pfull);
      for(i=0;i<3;i++) printf("\t%1.7e",xbest[i]);
      for(i=0;i<3;i++) for(j=0;j<3;j++) printf("\t%1.7e",Ubest[i][j]);
      for(i=0;i<3;i++) printf("\t%1.7e",ybest[i]);
      printf("\n");
    }
  }

PROFE(0);

PROFP(0, "total");
PROFP(1, "kernel");
PROFP(2, "handle_out");

  return 0;
}

prot_t* readpdb(char* pdbfn,char*** xpdb){
  int max_len=1000;
  char  *seq  =(char*)       mem(NULL,sizeof(char )*(max_len+2));
  short *num  =(short int*)  mem(NULL,sizeof(short)*(max_len+1));
  coor_t *coor=(float (*)[3])mem(NULL,sizeof(coor_t)*(max_len+1));
  prot_t *prot=(prot_t*)     mem(NULL,sizeof(prot_t));
  const char num2seq[21]={
    'G','A','S','C','V','T','I','P','M','D',
    'N','L','K','E','Q','R','H','F','Y','W','X'};
  const char num2name[21][4]={
    "GLY","ALA","SER","CYS","VAL","THR","ILE","PRO",
    "MET","ASP","ASN","LEU","LYS","GLU","GLN","ARG",
    "HIS","PHE","TYR","TRP","GAP"};
  int i,len=0;
  char line[1001];
  FILE *pdbfile=fopen(pdbfn,"rt");
  //added by mvg
  int npdb=0,mnpdb=1000;
  char** pdb=(char**)mem(NULL,sizeof(char*)*(mnpdb+1));
  line[1000]='\0';
  //end by mvg
  if(pdbfile==NULL) {
    fprintf(stderr,"There is no %s PDB file.\n\n",pdbfn);
    exit(-1);
  }
  while(fgets(line,1000,pdbfile)!=NULL){
    //added by mvg
    if(npdb>=mnpdb) {
	    mnpdb+=500;
	    pdb=(char**)mem(pdb,sizeof(char*)*(mnpdb+1));
    }
    i=strlen(line);
    pdb[npdb]=(char*)malloc((i+1)*sizeof(char));
    strncpy(pdb[npdb],line,i);
    pdb[npdb][i]='\0';
    npdb++;
    //end by mvg
    if(strncmp(line+13,"CA",2)||strncmp(line,"ATOM",4)) continue;
    if(len>=max_len) {
      max_len+=200;
      seq =(char*)       mem(seq,sizeof(char)*(max_len+1));
      num =(short int*)  mem(num,sizeof(short)  *max_len);
      coor=(float (*)[3])mem(coor,sizeof(coor_t)*max_len);
    }
    for(i=0;i<20;i++) if(!strncmp(line+17,num2name[i],3)) break;
    seq[len]=num2seq[i];
    num[len]=atoi(line+22);
    coor[len][0]=atof(line+30);
    coor[len][1]=atof(line+38);
    coor[len][2]=atof(line+46);
    line[13]='X';
    len++;
  }
  prot->len=len;
  if(len>0){
    prot->seq=seq;
    prot->seq[len]='\0';
  }
  prot->num=num;
  prot->coor=coor;
  prot->name=pdbfn;
  //added by mvg
  for(i=npdb;i<mnpdb;i++) pdb[i]=NULL;
  *xpdb=pdb;
  //end by mvg
  return prot;
}

int matrix1w(float *m,
	     const coor_t *p1,const coor_t *c1,const short *num1,const int lp1,
	     const coor_t *p2,const coor_t *c2,const short *num2,const int lp2) {
  float c,d;
  int l1,l2,k=0,n=0;
  for(l1=0;l1<lp1;l1++){
    for(l2=0;l2<lp2;l2++,k++){
      if(num1[l1]<=0||num2[l2]<=0) m[k]=0;
      else{
	c=c1[l1][0]*c2[l2][0]+c1[l1][1]*c2[l2][1]+c1[l1][2]*c2[l2][2];
       	if(c<MINCOS)           {m[k]=0;continue;}
        d=(p1[l1][0]-p2[l2][0])*(p1[l1][0]-p2[l2][0]);
        if(d>SEGRMSD1*SEGRMSD1){m[k]=0;continue;}
        d+=(p1[l1][1]-p2[l2][1])*(p1[l1][1]-p2[l2][1]);
        if(d>SEGRMSD1*SEGRMSD1){m[k]=0;continue;}
        d+=(p1[l1][2]-p2[l2][2])*(p1[l1][2]-p2[l2][2]);
        if(d>SEGRMSD1*SEGRMSD1){m[k]=0;continue;}
        m[k]=1;
        n++;
      }
    }
  }
  return n;
}

int matrix2w(float *m,
	     const coor_t *p1,const coor_t *c1,const short *num1,const int lp1,
	     const coor_t *p2,const coor_t *c2,const short *num2,const int lp2) {
  float c,d,sigma=2.0*SEGSIGMA*SEGSIGMA;
  int l1,l2,k=0,n=0;
  for(l1=0;l1<lp1;l1++){
    for(l2=0;l2<lp2;l2++,k++){
      if(num1[l1]<=0||num2[l2]<=0) m[k]=2.0*(-paraGap-paraExt);
      else{
	c=c1[l1][0]*c2[l2][0]+c1[l1][1]*c2[l2][1]+c1[l1][2]*c2[l2][2];
       	if(c<MINCOS)           {m[k]=2.0*(-paraGap-paraExt);continue;}
	d=(p1[l1][0]-p2[l2][0])*(p1[l1][0]-p2[l2][0]);
        if(d>SEGRMSD2*SEGRMSD2){m[k]=2.0*(-paraGap-paraExt);continue;}
        d+=(p1[l1][1]-p2[l2][1])*(p1[l1][1]-p2[l2][1]);
        if(d>SEGRMSD2*SEGRMSD2){m[k]=2.0*(-paraGap-paraExt);continue;}
        d+=(p1[l1][2]-p2[l2][2])*(p1[l1][2]-p2[l2][2]);
        if(d>SEGRMSD2*SEGRMSD2){m[k]=2.0*(-paraGap-paraExt);continue;}
#ifdef SEGRMSD0
	m[k]=exp(-d/sigma);
        //m[k]=(d<=SEGRMSD0*SEGRMSD0?1.0:(SEGRMSD2*SEGRMSD2-d)/(SEGRMSD2*SEGRMSD2-SEGRMSD0*SEGRMSD0));
	//printf("%f\t%f\t%f\t%f\n",sqrt(d),m[k],(d<=SEGRMSD0*SEGRMSD0?1.0:(SEGRMSD2*SEGRMSD2-d)/(SEGRMSD2*SEGRMSD2-SEGRMSD0*SEGRMSD0)),exp(-d/sigma));
#else
        m[k]=1.0;
#endif
        n++;
      }
    }
  }
  return n;
}
