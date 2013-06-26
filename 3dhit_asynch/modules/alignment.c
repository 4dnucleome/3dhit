/*
 * alignment.c
 *
 * autor: Łukasz Bieniasz-Krzywiec
 *
 */

#include <ctype.h>
#include <math.h>

#include "../3dhit.h"

/* original begin */

#define AA 1 /* ali to ali */
#define AG 2 /* ali to gap */
#define VH 4 /* vertical */
#define ZE 8 /* zero */
#if (defined(ALIGLOB)||defined(ALIGLOLOC))
#define Loctest(a,k)
#define Glotest(a,g,k) \
  if(a>score){score=a;be= (k);} \
  if(g>score){score=g;be=-(k);}
#else
#define Glotest(a,g,k)
#define Loctest(a,k) \
  if(a<=0.0){a=0.0;smpp[k]|=ZE;} \
  else{if(a>score){score=a;be=k;}}
#endif

int   pathn;          /* number of aligned atoms */
float pathv[SEGSIZE2];  /* new */

float align2_f(float *m, const int lp1, const int lp2, int* pathp) {
  //char mp[1000*1000];
  //float al[1000];
  //float gl[1000];
  static int segsize2=0;
  static char *smpp=NULL;
  static float *al=NULL;
  static float *gl=NULL;
  float score;
  int be=0;
  if(lp1>segsize2||lp2>segsize2){
    if(segsize2<SEGSIZE2){
      segsize2=SEGSIZE2;}
    if(lp1>segsize2){
      segsize2=lp1+100;}
    if(lp2>segsize2){
      segsize2=lp2+100;}
    //pathv=(float*)mem(pathv,sizeof(float)*segsize2);
    smpp =(char* )mem(smpp, sizeof(char )*segsize2*segsize2);
    al   =(float*)mem(al,   sizeof(float)*segsize2);
    gl   =(float*)mem(gl,   sizeof(float)*segsize2);}
#if (defined(ALIGLOB)||defined(ALIGLOLOC))
  score=-segsize2*segsize2;
#else
  score=m[0];
#endif
  /* forward */
  { int l1,l2,k=0;
    for(l2=0;l2<lp2;l2++,k++){
      al[l2]=m[k];
      gl[l2]=-paraExt; // changed from paraGap to paraExt
#if (defined(ALIGLOB)||defined(ALIGLOLOC))
      al[l2]-=(float)l2*paraExt;
      gl[l2]-=(float)l2*paraExt;
#endif
      smpp[k]=AA|VH;
      Loctest(al[l2],k)}
#ifdef ALIGLOLOC
    Glotest(al[l2-1],gl[l2-1],k-1)
#elif defined(ALIGLOB)
    Glotest(al[l2-1]-paraExt*(float)(lp1-1),gl[l2-1]-paraExt*(float)(lp1-1),k-1)
#endif
    for(l1=1;l1<lp1;l1++){
      float ad=m[k];
      float gd=-paraExt; // changed from paraGap to paraExt
#ifdef ALIGLOB
      ad-=(float)l1*paraExt;
      gd-=(float)l1*paraExt;
#endif
      smpp[k]=AA;
      Loctest(ad,k)
      for(l2=1,k++;l2<lp2;l2++,k++){
        float an=m[k];
        /* position gapped */
        float gn=gd-paraExt;
          smpp[k]=0;
        if(gn<gl[l2]-paraExt){
          gn=gl[l2]-paraExt;
          smpp[k]=VH;}
        if(gn<ad-paraGap){
          gn=ad-paraGap;
          smpp[k]=AG;}
        if(gn<al[l2]-paraGap){
          gn=al[l2]-paraGap;
          smpp[k]=AG|VH;}
        /* position aligned */
        if(al[l2-1]>gl[l2-1]){
          an+=al[l2-1];
          smpp[k]|=AA;}
        else
          an+=gl[l2-1];
        /* moving up */
        al[l2-1]=ad;
        gl[l2-1]=gd;
        ad=an;
        gd=gn;
        Loctest(ad,k)}
      al[l2-1]=ad;
      gl[l2-1]=gd;
#ifdef ALIGLOLOC
      Glotest(ad,gd,k-1)
#elif defined(ALIGLOB)
      Glotest(ad-paraExt*(float)(lp1-1-l1),gd-paraExt*(float)(lp1-1-l1),k-1)
#endif
      }
#if (defined(ALIGLOB)||defined(ALIGLOLOC))
    for(l2=lp2-2;l2>=0;l2--){
      Glotest(al[l2]-paraExt*(float)(lp2-1-l2),
        gl[l2]-paraExt*(float)(lp2-1-l2),k-(lp2-l2))}
#endif
  }
  /* backward */
  if(pathp!=NULL)
  { int stat=1,l1,l2,k;
    pathn=0;
    if(be<0){
      be=-be;
      stat=0;}
    l1=be/lp2;
    l2=be%lp2;
    for(k=l1+1;k<lp1;k++){
      pathp[k]=-1;
      }
    if(be){
      for(;;){
        if(stat){
#ifndef ALIGLOB
          if(smpp[be]&ZE){
            /* TEST fprintf(stderr,"ZERO: %d,%d\n",be/lp2,be%lp2); */
            break;}
#endif
          pathp[l1]=l2;
          pathn++;
          pathv[l1]=m[l1*lp2+l2]; // save pathv
          //SaveScore(be,pathv[l1]);
          stat=smpp[be]&AA;
          if(--l1<0)
            break;
          if(--l2<0)
            break;
          be-=lp2+1;
          }
        else{
          pathp[l1]=-1;
          stat=smpp[be]&AG;
          if(smpp[be]&VH){
            if(--l1<0)
              break;
            be-=lp2;}
          else{
            if(--l2<0)
              break;
            be--;}}}}
    for(;l1>=0;l1--){
      pathp[l1]=-1;
      }
  }
  return((float)score);
}

#undef Loctest
#undef Glotest
#undef AA
#undef AG
#undef VH
#undef ZE

void calcvect(coor_t *v, const coor_t *c, int n) {
  //calculate normalized c1 vectors for protin
  int   i,j,k;
  float x,y,z,d;
  //calculate normalized c vectors
  for(i=0,j=1,k=2;k<n;k++,j++,i++) {
    x=c[k][0]+c[i][0]-2.f*c[j][0];
    y=c[k][1]+c[i][1]-2.f*c[j][1];
    z=c[k][2]+c[i][2]-2.f*c[j][2];
    d=1.0/sqrt(x*x+y*y+z*z);
    v[j][0]=d*x;
    v[j][1]=d*y;
    v[j][2]=d*z;
  }
  //calculation for the first residue in prot[p]
  x=2*v[1][0]-v[2][0];
  y=2*v[1][1]-v[2][1];
  z=2*v[1][2]-v[2][2];
  d=1.0/sqrt(x*x+y*y+z*z);
  v[0][0]=d*x;
  v[0][1]=d*y;
  v[0][2]=d*z;
  i=n-3; j=n-2; k=n-1;
  //calculation for the last  residue in prot[p]
  x=2*v[j][0]-v[i][0];
  y=2*v[j][1]-v[i][1];
  z=2*v[j][2]-v[i][2];
  d=1.0/sqrt(x*x+y*y+z*z);
  v[k][0]=d*x;
  v[k][1]=d*y;
  v[k][2]=d*z;
}

int calc_within(int* xwin1, float win1, int* xwin2, float win2, int *pfull, int len1,
  coor_t* p1, coor_t* p2, float U[3][3], float xc[3], float yc[3]) {
  int i1,i2,j,k,c1=0,c2=0;
  float d,p1rot[3],p1xyz[3];
  win1*=win1;
  win2*=win2;
        for(i1=0;i1<len1;i1++) {
    if(pfull[i1]<0) continue;
    i2=pfull[i1];
    for(j=0;j<3;j++) p1rot[j]=0.f;
    for(j=0;j<3;j++) p1xyz[j]=p1[i1][j]-xc[j];
    for(j=0;j<3;j++) for(k=0;k<3;k++) p1rot[j]+=U[j][k]*p1xyz[k];
    for(j=0;j<3;j++) p1xyz[j]=p1rot[j]+yc[j];
    for(d=0,j=0;j<3;j++) d+=(p1xyz[j]-p2[i2][j])*(p1xyz[j]-p2[i2][j]);
    if(d<win2) {
      if(d<win1) c1++;
      c2++;
    }
  }
  xwin1[0]=c1;
  xwin2[0]=c2;
  return 0;
}

float mvg_3Dscore(int *pfull, float U[3][3],
  coor_t* p1, int  n1,    float   xc[3],
  coor_t* p2, int  n2,    float   yc[3]) {
  int   i,j,k,m,n,p,q,r,m1,m2;
  int   a1[n1], a2[n2], x1[n1], xa[n1], x2[n2];
  float tmp, min, p1rot[3], p1xyz[3], dist[n1];
  float score=0.0, sigma=2.0*MVGSIGMA*MVGSIGMA;
  for(i=0;i<n1;i++) a1[i]=-1;
  for(i=0;i<n2;i++) a2[i]=-1;
  for(m1=0,i=0;i<n1;i++) { if(pfull[i]>=0) a2[pfull[i]]=i; else x1[m1++]=i; }
  for(m2=0,i=0;i<n2;i++) { if(   a2[i]>=0) a1[   a2[i]]=i; else x2[m2++]=i; }
  if(m2*m1) {
    for(i=0;i<m1;i++) {
      for(j=0;j<3;j++) p1rot[j]=0.f;
      for(j=0;j<3;j++) p1xyz[j]=p1[x1[i]][j]-xc[j];
      for(j=0;j<3;j++) for(k=0;k<3;k++) p1rot[j]+=U[j][k]*p1xyz[k];
      for(j=0;j<3;j++) p1xyz[j]=p1rot[j]+yc[j];
      for(m=x2[0],min=0.0,k=0;k<3;k++) min+=(p1xyz[k]-p2[x2[0]][k])*(p1xyz[k]-p2[x2[0]][k]);
      for(j=1;j<m2;j++) {
        for(tmp=0.0,k=0;k<3;k++) tmp+=(p1xyz[k]-p2[x2[j]][k])*(p1xyz[k]-p2[x2[j]][k]);
        if(tmp<min) { min=tmp; m=x2[j]; }
      }
      dist[i]=min;
      xa[i]=m;
    }
    ascending_sort(dist, x1, xa, 0, m1-1);
    for(i=0;i<m1;i++) {
      if(a2[xa[i]]<0) {
        a1[x1[i]]=xa[i];
        a2[xa[i]]=x1[i];
      }
    }
    for(i=0;i<n1;i++) {
      if(pfull[i]>=0) continue;
      for(j=i+1;j<n1;j++) if(pfull[j]>=0) break;
      for(k=i;k<j;k++) if(a1[k]>=0) break;
      if(k==j) { i=j; continue; }
      for(n=1,m=k,p=k+1;p<j;p++) {
        if(a1[p]<0) continue;
        if(a1[m]<a1[p]) { n++; m=p; continue; }
        for(q=p;q<j;q++) if(a1[q]>=0) break;
        if(q==j) break;
        if(n<MVGALIGN) {
          if(i==0 || j==n1) { for(r=k;r<p;r++) a1[r]=-1; }
          else { for(r=k;r<p;r++) if(a1[r]<pfull[i-1] || pfull[j]<a1[r]) a1[r]=-1; }
        }
        n=1; k=q; m=q; p=q;
      }
      if(n<MVGALIGN) {
        if(i==0 || j==n1) { for(r=k;r<p;r++) a1[r]=-1; }
        else { for(r=k;r<p;r++) if(a1[r]<pfull[i-1] || pfull[j]<a1[r]) a1[r]=-1; }
      }
      i=j;
    }

  }
  for(j=0,i=0;i<n1;i++) if(pfull[i]>=0) j++;
  for(k=0,i=0;i<n1;i++) if(   a1[i]>=0) k++;
  if(j<k) {
    for(i=0;i<n1;i++) pfull[i]=a1[i];
    return mvg_3Dscore(pfull,U,p1,n1,xc,p2,n2,yc);
  }
  for(i=0;i<n1;i++) {
    if(pfull[i]<0) continue;
    for(j=0;j<3;j++) p1rot[j]=0.f;
    for(j=0;j<3;j++) p1xyz[j]=p1[i][j]-xc[j];
    for(j=0;j<3;j++) for(k=0;k<3;k++) p1rot[j]+=U[j][k]*p1xyz[k];
    for(j=0;j<3;j++) p1xyz[j]=p1rot[j]+yc[j];
    for(tmp=0.0,j=0;j<3;j++) tmp+=(p1xyz[j]-p2[pfull[i]][j])*(p1xyz[j]-p2[pfull[i]][j]);
    score+=exp(-tmp/sigma);
  }
  for(i=0;i<n1;i++) {
    if(pfull[i]<0) {
      score-=(MVGPRGAP+MVGPREXT);
      for(i++;i<n1;i++) {
        if(pfull[i]<0) score-=MVGPREXT;
        else break;
      }
    }
  }
  return score;
}

void ascending_sort(float* a, int* s, int* y, int l, int p) {
  float w,x=a[(l+p)/2];
  int i=l,j=p;
  int t,z;
  while(i<j) {
    while(a[i]<x) i++;
    while(x<a[j]) j--;
    if(i<=j) {
      w=a[i]; a[i]=a[j]; a[j]=w;
      t=s[i]; s[i]=s[j]; s[j]=t;
      z=y[i]; y[i]=y[j]; y[j]=z;
      i++;
      j--;
    }
  }
  if(l<j) ascending_sort(a,s,y,l,j);
  if(i<p) ascending_sort(a,s,y,i,p);
}

void print_alin1(FILE* file,int len1,int len2,char* seq1,char* seq2,int* pfull) {
  char alin[len1+len2];
  int i,j,a1,a2,b1=0,b2=0,lst,x=0,id=0,tt=0; /* lbk */
  for(a1=0;a1<len1;a1++) if(pfull[a1]>=0) {
    b1=pfull[a1];
    break;
  }
  if(a1==len1) return;
  for(a2=len1-1;a2>=0;a2--) if(pfull[a2]>=0) {
    b2=pfull[a2];
    break;
  }
  for(lst=b1-1,i=a1;i<=a2;i++) {
    if((j=pfull[i])>=0) {
      for(lst++;lst<j;lst++) alin[x++]=tolower(seq2[lst]);
      if(seq1[i]==seq2[j]) id++;
      alin[x++]=seq2[j];
      tt++;
    }
    else alin[x++]='-';
  }
  alin[x]='\0';
  fprintf(file,"%.0f%%\t%i\t[%i:%i-%i:%i]\t%s",100.*id/tt,tt,a1+1,a2+1,b1+1,b2+1,alin);
}

void print_alin(FILE* file,const int qlen,const int len2,const char* seq1,const char* seq2,const int* pathp,const int offset1,const int offset2, int max_len)
{ static char *ali=NULL;
  static int alilen=0;
  int lp=-1,a,a1,a2=0,b,b1,b2=0,ll=0;
  int id=0,match=0;

  if(qlen+len2+1>alilen){
    if(alilen==0){
      alilen=qlen+max_len+1;
      ali=(char*)mem(NULL,sizeof(char)*alilen);}
    else{
      free(ali);
      alilen=qlen+len2+100;
      ali=(char*)mem(NULL,sizeof(char)*alilen);}}
  for(a=0;pathp[a]<0;){
    if(++a>=qlen){
      fprintf(file,"\n");
      return;}}
  b=pathp[a];
  a1=a;
  b1=b;
  for(;a<qlen;){
    if(b==pathp[a]){
      ali[++lp]=toupper(seq2[b]);
      if(seq1[a]==seq2[b]){
        id++;}
      match++;
      a2=a;
      b2=b;
      a++;
      b++;
      ll=lp;
      continue;}
    for(;pathp[a]<0;){
      ali[++lp]='-';
      if(++a>=qlen){
        goto DONE;}}
    for(;b<pathp[a];){
      ali[++lp]=tolower(seq2[b]);
      if(++b>=len2){
        goto DONE;}}}
DONE:
  ali[ll+1]='\0';
  a=a2-a1;
  b=b2-b1;
  if(a<b){
    lp=a+1;}
  else{
    lp=b+1;}
  fprintf(file,"%d%%\t%d\t[%d:%d-%d:%d]\t%s\n",
    id*100/match,match,a1+offset1,a2+offset1,b1+offset2,b2+offset2,ali);
}

void print_whole_pdb1(FILE* out, char** pdb1, char** pdb2,
                float U[3][3], float xc[3], float yc[3]) {
  int n,j,k;
  float p1rot[3];
  float p1xyz[3];
  char str30[31];
  str30[30]='\0';
  for(n=0;pdb1[n]!=NULL;n++) {
    if( strncmp(pdb1[n],"ATOM",4) && strncmp(pdb1[n],"HETATM",6) ) {
      fprintf(out,"%s",pdb1[n]);
      continue;
    }
    strncpy(str30,pdb1[n],30);
    p1xyz[0]=atof(pdb1[n]+30)-xc[0];
    p1xyz[1]=atof(pdb1[n]+38)-xc[1];
    p1xyz[2]=atof(pdb1[n]+46)-xc[2];
    for(j=0;j<3;j++) p1rot[j]=0.f;
    for(j=0;j<3;j++) for(k=0;k<3;k++) p1rot[j]+=U[j][k]*p1xyz[k];
    for(j=0;j<3;j++) p1xyz[j]=p1rot[j]+yc[j];
    fprintf(out,"%s%8.3f%8.3f%8.3f%s",str30,p1xyz[0],p1xyz[1],p1xyz[2],pdb1[n]+54);
  }
  for(n=0;pdb2[n]!=NULL;n++) fprintf(out,"%s",pdb2[n]);
}

void print_whole_pdb2(FILE* out, char** pdb1, char** pdb2,
                float U[3][3], float xc[3], float yc[3]) {
  int n,j,k;
  float p2rot[3];
  float p2xyz[3];
  char str30[31];
  str30[30]='\0';
  for(n=0;pdb1[n]!=NULL;n++) fprintf(out,"%s",pdb1[n]);
  for(n=0;pdb2[n]!=NULL;n++) {
    if( strncmp(pdb2[n],"ATOM",4) && strncmp(pdb2[n],"HETATM",6) ) {
      fprintf(out,"%s",pdb2[n]);
      continue;
    }
    strncpy(str30,pdb2[n],30);
    p2xyz[0]=atof(pdb2[n]+30)-yc[0];
    p2xyz[1]=atof(pdb2[n]+38)-yc[1];
    p2xyz[2]=atof(pdb2[n]+46)-yc[2];
    for(j=0;j<3;j++) p2rot[j]=0.f;
    for(j=0;j<3;j++) for(k=0;k<3;k++) p2rot[j]+=U[k][j]*p2xyz[k];
    for(j=0;j<3;j++) p2xyz[j]=p2rot[j]+xc[j];
    fprintf(out,"%s%8.3f%8.3f%8.3f%s",str30,p2xyz[0],p2xyz[1],p2xyz[2],pdb2[n]+54);
  }
}

/* original end */
