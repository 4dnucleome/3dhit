/*
 * 3dhit.c
 *
 * autor: Łukasz Bieniasz-Krzywiec
 *
 */

#include <pthread.h>
#include <libspe2.h>

#include "3dhit.h"
#include "modules/err.h"

/*
 * PRZYDATNE MAKRA:
 */

#undef PROFS
#undef PROFE
#undef PROFP

#ifdef PROF
# include <sys/time.h>
# define COUNTERS 8
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

  
/*
 * DEKLARACJE STAŁYCH:
 */

/* Maksymalna liczba tworzonych wątków obliczeniowych. */
#define MAX_SPE_THREADS 16


/*
 * DEKLARACJE TYPÓW:
 */

/* Dane wątku obliczeniowego. */
typedef struct spe_worker_data_t_ {
  spe_context_ptr_t ctx;
  void              *argp;
  pthread_t         thread;
} spe_worker_data_t;


/*
 * DEKLARACJE ZMIENNYCH GLOBALNYCH:
 */

/* Program wątku obliczeniowego. */
extern spe_program_handle_t worker_spe;

/* Tablica z danymi wątków obliczeniowych. */
spe_worker_data_t           worker[MAX_SPE_THREADS];

/* Maksymalna liczba tworzonych wątków obliczeniowych. */
int                         max_spe_threads = MAX_SPE_THREADS;

/* Liczba wykorzystywanych SPE */
int                         liczba_spe = 0;

/* Parametry wątków obliczeniowych. */
spe_worker_parameters_t     params[MAX_SPE_THREADS];


/* Wyjścia wątków obliczeniowych. */
spe_worker_output_t         *output[MAX_SPE_THREADS];
volatile uint32_t           *output_status[MAX_SPE_THREADS];

/* Białka. */
prot_t                      *protin1, *protin2;

/* Rozmiary. */
uint                        size, msgs_size;


/*
 * DEFINICJE PODPROGRAMÓW:
 */

/* Funkcja zwraca losową liczbę całkowitą z przedziału [a,b]. */
inline int rand_int(int a, int b) {
  /* a <= b */
  return a + (int) ((b - a + 1) * (rand() / (RAND_MAX + 1.0)));
}

inline int min(int a, int b) { return (a < b) ? a : b; }

/*
 * Uruchamia wątek obliczeniowy.
 *
 * Argumenty:
 *  arg - wskaźnik na dane wątku obliczeniowego
 *
 * Wynik:
 *  NULL
 */
void *worker_function(void *arg) {
  spe_worker_data_t *worker = (spe_worker_data_t*)arg;
  uint              entry = SPE_DEFAULT_ENTRY;
  int               flags = 0, rc;

  PDEBUG("--> worker_function()\n");

  rc = spe_context_run(worker->ctx, &entry, flags, worker->argp, NULL, NULL);
  if (rc < 0) {
    syserr("spe_context_run()\n");
  }

  PDEBUG("<-- worker_function()\n");

  pthread_exit(NULL);
}

/*
 * Uruchamia wątki obliczeniowe.
 */
void create_spe_workers() {
  int i, rc;

  PDEBUG("--> create_spe_workers()\n");

  /* Pytamy o liczbę dostępnych SPE. */
  liczba_spe = spe_cpu_info_get(SPE_COUNT_USABLE_SPES, -1);

  if (liczba_spe < 0) {
    syserr("spe_cpu_info_get()\n");
  }
  if (liczba_spe > max_spe_threads) {
    liczba_spe = max_spe_threads;
  }

  /* Uruchamiamy wątki obliczeniowe. */
  for (i = 0; i < liczba_spe; ++i) {
    if ((worker[i].ctx = spe_context_create(0, NULL)) == NULL) {
      syserr("spe_context_create()\n");
    }

    rc = spe_program_load(worker[i].ctx, &worker_spe);
    if (rc != 0) {
      syserr("spe_program_load()\n");
    }

    params[i].protein1 = *protin1;
    params[i].protein2 = *protin2;
    params[i].output = output[i];
    params[i].output_status = output_status[i];
    params[i].from = i * size;
    params[i].to = min((i + 1) * size, msgs_size);

    *output_status[i] = 0;

    worker[i].argp = &params[i];

    rc = pthread_create(&worker[i].thread, NULL, &worker_function, &worker[i]);
    if (rc != 0) {
      syserr("pthread_create()\n");
    }
  }

  PDEBUG("liczba uruchomionych wątków obliczeniowych: %d\n", liczba_spe);

  PDEBUG("<-- create_spe_workers()\n");
}

/*
 * Zamyka wątki obliczeniowe.
 */
void destroy_spe_workers() {
  int i, rc;

  PDEBUG("--> destroy_spe_workers()\n");

  for (i = 0; i < liczba_spe; ++i) {
    rc = pthread_join(worker[i].thread, NULL);
    if (rc != 0) {
      syserr("pthread_join()\n");
    }

    rc = spe_context_destroy(worker[i].ctx);
    if (rc != 0) {
      syserr("spe_context_destroy()\n");
    }
  }

  PDEBUG("<-- destroy_spe_workers()\n");
}

int     len1;       /* długość 1 białka */
int     len2;       /* długość 2 białka */
int     max_len;

float   alibest2 = 0.0;
float   conta = 0.0;
float   contt = 0.0;
int     hitn = 0;
int     hitn1 = 0;
int     hitn2 = 0;
int     hitn3 = 0;
int     pathn2 = 0;
int     pathn3 = 0;

float   *mfull;     /* tablica rozmiaru len1*len2 + 1 */
int     *pfull;     /* tablica rozmiaru len1 + 1 */

char    **pdb1;     /* whole 1st PDB file */
char    **pdb2;     /* whole 2nd PDB file */

float   Ubest[3][3], xbest[3], ybest[3];
int     *pbest;

void create_io(int len1, int len2) {
  int i;

  msgs_size = (len1 - LEN + 1) * (len2 - LEN + 1);
  size = msgs_size / max_spe_threads + (msgs_size % max_spe_threads ? 1 : 0);

  for (i = 0; i < max_spe_threads; ++i) {
    posix_memalign((void*)(&output[i]), 128, sizeof(spe_worker_output_t) * size);
    posix_memalign((void*)(&output_status[i]), 128, sizeof(uint32_t));
  }
}

void handle_output() {
  int a, b, i, j;

  for (a = 0; a < liczba_spe; ++a) {
    for (b = 0; b < (int)*output_status[a]; ++b) {
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
        for (i = output[a][b].t; i < protin1->len; ++i) {
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


/*
 * Wypisuje sposób użycia.
 */
void usage(char **argv) {
  fprintf(stderr,"Usage:\n\t%s pdb_file1 pdb_file2\n", argv[0]);
  exit(1);
}


/*
 * PROGRAM GŁÓWNY:
 */

int main(int argc, char **argv) {
  PROFS(0);

  int p1, p2, l;

  if (argc < 3) {
    usage(argv);
  }

  protin1 = readpdb(argv[1], &pdb1); /* read protein 1 */
  protin2 = readpdb(argv[2], &pdb2); /* read protein 2 */

  len1 = protin1->len;
  len2 = protin2->len;
  max_len = protin2->len;

  mfull = (float*)mem(NULL, sizeof(float) * (len1 * len2 + 1));
  pfull = (int*)mem(NULL, sizeof(int) * (len1 + 1));
  pbest = (int*)mem(NULL, sizeof(int) * (protin1->len));

  calcvect(protin1->c, (const coor_t*)protin1->coor, len1);
  calcvect(protin2->c, (const coor_t*)protin2->coor, len2);

  if (protin1->len < SEGHITS2 || protin2->len < SEGHITS2) {
    return 0;
  }

  for (l = protin1->len * protin2->len; --l >= 0; ) {
    mfull[l] = 0.0;
  }

  /* AKCELERACJA */

PROFS(5);

PROFS(1);
  create_io(len1, len2);
PROFE(1);

PROFS(2);
  create_spe_workers();
PROFE(2);

PROFS(3);
  destroy_spe_workers();
PROFE(3);

PROFS(4);
  handle_output();
PROFE(4);

PROFE(5);

  /* AKCELERACJA */

  /* Reszta funkcji nie zmmodyfikowana. */

  /* clean alignment matrix */
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
        mfull[p1*protin2->len+p2]/=alibest2;
      }
    }
  }

  float tmp;
  tmp = align2_f(mfull, protin1->len, protin2->len, pfull);

  /* align */
  if (tmp > 0.0){ /* always true ! */
    /* paste aligned coordinates */
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

  /* print results */
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
    print_alin(out,protin1->len,protin2->len,protin1->seq,protin2->seq,pfull,protin1->num[0],protin2->num[0],max_len);
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

  PROFP(0, "total_ppe");
  PROFP(1, "create_io");
  PROFP(2, "create_spe");
  PROFP(3, "destroy_spe");
  PROFP(4, "handle_out");
  PROFP(5, "kernel");

  return 0;
}
