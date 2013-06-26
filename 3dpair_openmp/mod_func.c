#include "modul.h"

void* mem(void* p,long size)
{
  p=realloc(p,size);
  if(p==NULL)
  { fprintf(stderr,"no memory for %ld bytes, fatal\n",size);
    exit(-1);
  }
  return(p);
}
float dist(const float *x,const float *y)
{ return(sqrtf(
    (x[0]-y[0])*(x[0]-y[0])+
    (x[1]-y[1])*(x[1]-y[1])+
    (x[2]-y[2])*(x[2]-y[2])
    ));
}
float dist2(const float *x,const float *y)
{ return(
    (x[0]-y[0])*(x[0]-y[0])+
    (x[1]-y[1])*(x[1]-y[1])+
    (x[2]-y[2])*(x[2]-y[2]));
}
int min_int(const int* i,int l)
{ int n=(--l)-1;
  for(;n>=0;n--){
    if(i[n]<i[l]){
      l=n;}}
  return(l);
}
int min_float(const float* i,int l)
{ int n=(--l)-1;
  for(;n>=0;n--){
    if(i[n]<i[l]){
      l=n;}}
  return(l);
}
int max_float(const float* i,int l)
{ int n=(--l)-1;
  for(;n>=0;n--){
    if(i[n]>i[l]){
      l=n;}}
  return(l);
}
void clusfromto(char aa,coor_t *co,int* c_from,int* c_to,float ends)
{  
  extern const int clus_a[21];
  extern const cluster_t cluster[];
  float d=dist(co[0],co[LEN-1]),da=d-ends,db=d+ends;
  const char num2seq[21]={
    'G','A','S','C','V','T','I','P','M','D',
    'N','L','K','E','Q','R','H','F','Y','W','.'};

  int a,ca=0,cb=0,c1,c2,c;
  for(a=0;a<20;a++){ // select cluster range for given amino acid
    if(num2seq[a]==aa){
      ca=clus_a[a];
      cb=clus_a[a+1]-1;
      break;}}
  if(a>=20){
    *c_from=0;
    *c_to=0;
    return;}
  if(db<cluster[ca].d||da>cluster[cb].d){
    *c_from=0;
    *c_to=0;
    return;}
  c1=ca;
  c2=cb;
  if(da>cluster[c1].d){
    while(cluster[c1].d<cluster[c2].d){
      c=(int)(c1+(c2-c1)*(da-cluster[c1].d)/(cluster[c2].d-cluster[c1].d));
      if(cluster[c].d<=da){
        c1=c+1;
        if(cluster[c1].d>=da){
          c1=c2;
          break;}}
      else{
        c2=c-1;
        if(cluster[c2].d<=da){
          c1=c2;
          break;}}}}
  *c_from=c1;
  c2=cb;
  if(db<cluster[c2].d){
    while(cluster[c1].d<cluster[c2].d){
      c=(int)(c1+(c2-c1)*(db-cluster[c1].d)/(cluster[c2].d-cluster[c1].d));
      if(cluster[c].d<db){
        c1=c+1;
        if(cluster[c1].d>=db){
          c1=c2;
          break;}}
      else{
        c2=c-1;
        if(cluster[c2].d<=db){
          c1=c2;
          break;}}}}
  *c_to=c2+1;
  return;
}
