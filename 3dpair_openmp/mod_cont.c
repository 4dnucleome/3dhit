#include "modul.h"

   // distsigma=2.0*DISTSIGMA*DISTSIGMA (LSerrano)
float distsigma=1.442695*DISTSIGMA*DISTSIGMA;
int   minseqdiff=MINSEQDIFF;

void setcontacts(contacts_t *cont,coor_t *coor,int l,short *seqnum)
{ int i,j,k;
  if(cont->map==NULL&&l>0){
    cont->map=(float*)mem(NULL,sizeof(float)*(l-1)*l/2);
    for(k=0,i=0;i<l;i++){
      for(j=0;j<i;j++,k++){
        if(seqnum[i]>0&&seqnum[j]>0&&abs(seqnum[i]-seqnum[j])>=minseqdiff){
          //float dist2=DIST2(coor[i],coor[j]);
          //if(dist2<DISTCUTOFF*DISTCUTOFF){
          //  cont->map[k]=exp(-dist2/distsigma);}
          //else{
          //  cont->map[k]=0.0;}}
          cont->map[k]=exp(- DIST2(coor[i],coor[j])/distsigma);}
        else{
          cont->map[k]=0.0;}}}}
}
void contactoverlap(contacts_t *cont1,contacts_t *cont2,int l,float *conta,float *contt) { 
  int i,j,k;
  float *cc=(float*)mem(NULL,sizeof(float)*l*3);
  float *c1=cc+l,t1=0.0;
  float *c2=c1+l,t2=0.0;
  float a=0.0,t=0.0;
  memset(cc,0,sizeof(float)*l*3);
  for(k=0,i=0;i<l;i++){
    for(j=0;j<i;j++,k++){
      float c=(cont1->map[k]<cont2->map[k]?cont1->map[k]:cont2->map[k]);
      cc[i]+=c;
      c1[i]+=cont1->map[k];
      c2[i]+=cont2->map[k];
      cc[j]+=c;
      c1[j]+=cont1->map[k];
      c2[j]+=cont2->map[k];
    }
  }
  for(k=0,i=0;i<l;i++){
    if(c1[i]>0.0&&c2[i]>0.0){
      k++;
      t+=cc[i];
      a+=2.0*cc[i]/(c1[i]+c2[i]);
      t1+=c1[i];
      t2+=c2[i];
    }
  }
  *conta=a; /* overlap by atom */
  if(t>0.0) t*=(float)k*2.0/(t1+t2);
  *contt=t; /* total overlap times atoms */
  free(cc);
}
