#include "modul.h"

#define Zero 0.0000001
#define Max_orth 0.001
#define Huge (1.0e10)
#define MOVE_NO
#define TEST_MIRROR_NO
#define PRINT_MIRROR_NO
#define ENABLE_MIRROR_NO
#define TEST_ORTHOGONALITY
#define WARN_ORTHOGONALITY_ERROR_NO

double wmove(const coor_t* y_old,coor_t* x_old,const float* movew,const int n, float xc[3], float yc[3], double U[3][3])
{ int i,j,k;
  double R[3][3]={{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
  double r_r[3][3]={{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
  double eval[3]={0.0,0.0,0.0};
  double evec[3][3]={{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
  double b[3][3]={{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
  double rms=0.0,wc=0.0;
#ifdef TEST_MIRROR
  double mirror=0.0;
#endif
  double Pi=0.0;
  coor_t *y=(float (*)[3])mem(NULL,sizeof(coor_t)*n);
  coor_t *x=(float (*)[3])mem(NULL,sizeof(coor_t)*n);

  for(i=0;i<3;i++){
    yc[i]=0.0;
    xc[i]=0.0;
    for(j=0;j<3;j++){
      U[i][j]=0.0;}} 

  /* centers of Y,X */
  if(Pi==0.0)
    Pi=2.0*asin((double)1.0);
  for(i=0;i<n;i++){
    wc+=movew[i];
    for(j=0;j<3;j++){
      yc[j]+=y_old[i][j]*movew[i];
      xc[j]+=x_old[i][j]*movew[i];}}
  if(wc<=0.0){
    return(-Huge);}
  for(j=0;j<3;j++){
    yc[j]/=wc;
    xc[j]/=wc;}
  for(i=0;i<n;i++){
    for(j=0;j<3;j++){
      y[i][j]=(y_old[i][j]-yc[j]);   //*movew[i];
      x[i][j]=(x_old[i][j]-xc[j]);}} //*movew[i];}}

#ifdef MOVETEST
  fprintf(stderr,"\nnew coordinates of x,y:\n");
  for(i=0;i<n;i++){
    for(j=0;j<3;j++)
      fprintf(stderr,"\t%.2f",x[i][j]);
    putc('\t',stderr);
    for(j=0;j<3;j++)
      fprintf(stderr,"\t%.2f",y[i][j]);
    putc('\n',stderr);}
  fprintf(stderr,"centres of x,y:\n");
  for(j=0;j<3;j++)
    fprintf(stderr,"\t%.2f",xc[j]);
  putc('\t',stderr);
  for(j=0;j<3;j++)
    fprintf(stderr,"\t%.2f",yc[j]);
  putc('\n',stderr);
  fprintf(stderr,"\ncorrelations matrix:\n");
#endif
  /* end of centers of sets */
  /* correlation matrix R */
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      double _x=0.0;
      for(k=0;k<n;k++)
        _x+=y[k][i]*x[k][j]*movew[k]; // added !
      R[i][j]=_x;
#ifdef MOVETEST
      fprintf(stderr,"\t%.2f",R[i][j]);
#endif
    }
#ifdef MOVETEST
    fprintf(stderr,"\n");
#endif
  }
  /* end correlation matrix */
  /* r_r = correlation matrix ** 2 */
  for(i=2;i>=0;i--)
    for(j=i;j>=0;j--){
      r_r[i][j]=0.0;
      for(k=2;k>=0;k--)
       	r_r[i][j]+=R[k][i]*R[k][j];
      r_r[j][i]=r_r[i][j];}
#ifdef MOVETEST
  fprintf(stderr,"correlation matrix ** 2\n");
  for(i=0;i<3;i++){
    for(j=0;j<3;j++)
      fprintf(stderr,"\t%.2f",r_r[i][j]);
    fprintf(stderr,"\n");}
#endif
  /* end r_r = correlation matrix ** 2 */
  /* eigenvalue , eigenvector */
  { /* det(r_r-eval*I)=0 */
    double _a=r_r[0][0];
    double _b=r_r[0][1];
    double _c=r_r[0][2];
    double _d=r_r[2][2];
    double _e=r_r[1][1];
    double _f=r_r[1][2];
    double _g,_h,PHI;
    /* A= -1; */
    double B=-_a-_e-_d;
    double C=-_f*_f-_b*_b-_c*_c+_a*_d+_e*_d+_a*_e;
    double D=-_a*_e*_d-2.0*_b*_c*_f+_a*_f*_f+_e*_c*_c+_d*_b*_b;
    double P=(C-B*B/3.0)/3.0;
    double Q=(2.0*B*B*B/27.0-(B*C/3.0)+D)/2.0;
    double _r=sqrt((double)fabs(P));
#ifdef MOVETEST
    fprintf(stderr,"det=x**3+(%-.2f)*x**2+(%-.2f)*x+(%-.2f)\n",B,C,D);
    fprintf(stderr,"P=%.7f\nQ=%.7f\n0>D=%.2f\n",P,Q,Q*Q+P*P*P);
#endif

    if(Q<0.0)
      _r=-_r;
    if(fabs(_r)>Zero&&fabs(Q)>Zero){
      double _x=Q/(_r*_r*_r);
#ifdef MOVETEST
      fprintf(stderr,"r=%.7f\nQ/r/r/r=%.2f\n",_r,_x);
#endif
      if(_x>1.0) _x=1.0;
      if(_x<-1.0) _x= -1.0;
      PHI=acos(_x)/3.0;}
    else
      PHI=Pi/6.0;
    eval[0]=-2.0*_r*cos(PHI)-(B/3.0);
    eval[1]= 2.0*_r*cos(Pi/3.0-PHI)-B/3.0;
    eval[2]= 2.0*_r*cos(Pi/3.0+PHI)-B/3.0;
#ifdef MOVETEST
    fprintf(stderr,"r=%.2f\nPHI=%.2f\n",_r,PHI*3.0);
    fprintf(stderr,"PHI/3=%.2f\n",PHI);
    fprintf(stderr,"PHI/3+1/3=%.2f\n",PHI+1.0/3.0);
    fprintf(stderr,"PHI/3-1/3=%.2f\n",PHI-1.0/3.0);
    fprintf(stderr,"(Q/r**3)=%.2f = %.2f=cos(PHI)\n",Q/_r/_r/_r,cos(PHI*3.0));
    fprintf(stderr,
      "cos(PHI/3)=%.2f\ncos(Pi/3-PHI/3)=%.2f\ncos(Pi/3+PHI/3)=%.2f\n",
      cos(PHI),cos(-PHI+Pi/3.0),cos(PHI+Pi/3.0));
    fprintf(stderr,"%.2f*%.2f*%.2f-(%.2f)=(%.2f)\n",-2.0,_r,
      cos(PHI),B/3.0,-2.0*_r*cos(PHI)-(B/3.0));
    fprintf(stderr,"eigenvalues :\n");
    for(i=0;i<3;i++){
      fprintf(stderr,"\t%.2f\n",eval[i]);
      fprintf(stderr,"det=x**3+(%-.2f)*x**2+(%-.2f)*x+(%-.2f)=%f\n",B,C,D,
	eval[i]*eval[i]*eval[i]+B*eval[i]*eval[i]+C*eval[i]+D);}
    fprintf(stderr,"\n");
#endif
    /* found eigenvalues using det=0 */
    /* sorting evalues */
    if(eval[1]>eval[0]&&eval[1]>eval[2]){
      double _x=eval[0];
      eval[0]=eval[1];
      eval[1]=_x;}
    else
      if(eval[2]>eval[0]){
        double _x=eval[0];
	eval[0]=eval[2];
	eval[2]=_x;}
    if(eval[2]>eval[1]){
      double _x=eval[1];
      eval[1]=eval[2];
      eval[2]=_x;}
    /* end sorting eigenvalues */
    /* eigenvectors using GAUS */
    /* might be it is not nessecery */
    _b=r_r[0][1];
    _c=r_r[0][2];
    _d=r_r[1][0];
    _f=r_r[1][2];
    _g=r_r[2][0];
    _h=r_r[2][1];
#ifdef MOVETEST
    for(i=0;i<3;i++)
#else
    for(i=0;i<2;i++)
#endif
    { double _x,_y,_z,_i;
      _a=r_r[0][0]-eval[i];
      _e=r_r[1][1]-eval[i];
      _i=r_r[2][2]-eval[i];
      _x=_a*_e-_d*_b;
      _y=_a*_i-_c*_g;
      _z=_e*_i-_f*_h;
      /* overwriting B,C,D */
      B=fabs(_x);
      C=fabs(_y);
      D=fabs(_z);
#ifdef MOVETEST /* check if linear */
      if(B>Zero)
	fprintf(stderr,"x=%f,[%f],[%f],[%f]\n",_x,
	  (_e*_c-_f*_b)/_x,(_a*_f-_c*_d)/_x,-1.0);
      if(C>Zero)
	fprintf(stderr,"y=%f,[%f],[%f],[%f]\n",_y,
	  (_b*_i-_c*_h)/_y,-1.0,(_a*_h-_b*_g)/_y);
      if(D>Zero)
	fprintf(stderr,"z=%f,[%f],[%f],[%f]\n",_z,
	  -1.0,(_d*_i-_f*_g)/_z,(_e*_g-_d*_h)/_z);
#endif
      if(B+C+D>Zero){
        if(B>=C && B>=D){
	  evec[i][0]=(_e*_c-_f*_b)/_x;
	  evec[i][1]=(_a*_f-_c*_d)/_x;
	  evec[i][2]= -1.0;}
	else{
	  if(C>D||(C==D&&i)){
	    evec[i][0]=(_b*_i-_c*_h)/_y;
	    evec[i][2]=(_a*_h-_b*_g)/_y;
	    evec[i][1]= -1.0;}
	  else{
	    evec[i][1]=(_d*_i-_f*_g)/_z;
	    evec[i][2]=(_e*_g-_d*_h)/_z;
	    evec[i][0]= -1.0;}}}
      else{
        if(!i){
	  evec[0][0]=1.0;
	  evec[0][1]=0.0;
	  evec[0][2]=0.0;}
	/* assumed i = 0 or 1 */
	/* create evec[i] orthogonal to evec[0] */
	else{
	  B=fabs(evec[0][0]);
	  C=fabs(evec[0][1]);
	  D=fabs(evec[0][2]);
	  if(B>C&&B>D){
	    evec[i][0]=  evec[0][1];
	    evec[i][1]= -evec[0][0];
	    evec[i][2]=0.0;}
	  else{
	    evec[i][0]=0.0;
	    if(C>D){
	      evec[i][1]=  evec[0][2];
	      evec[i][2]= -evec[0][1];}
	    else{
	      evec[i][1]= -evec[0][2];
	      evec[i][2]=  evec[0][1];}}}}
      /* normalizing eigenvectors */
      _x=
        sqrt(evec[i][0]*evec[i][0]+evec[i][1]*evec[i][1]+evec[i][2]*evec[i][2]);
      for(j=0;j<3;j++)
        evec[i][j]/=_x;
#ifdef MOVETEST
      for(j=0;j<3;j++){
        r_r[j][j]-=eval[i];
	fprintf(stderr,"0=%.2f*%.2f+%.2f*%.2f+%.2f*%.2f=%.2f\n",
	  r_r[j][0],evec[i][0],r_r[j][1],evec[i][1],r_r[j][2],evec[i][2],
	  r_r[j][0]*evec[i][0]/_x+r_r[j][1]*evec[i][1]/_x+
	  r_r[j][2]*evec[i][2]/_x);
	r_r[j][j]+=eval[i];}
#endif
    }
  }
  /* right handed system */
  evec[2][0]=evec[0][1]*evec[1][2]-evec[1][1]*evec[0][2];
  evec[2][1]=evec[0][2]*evec[1][0]-evec[0][0]*evec[1][2];
  evec[2][2]=evec[0][0]*evec[1][1]-evec[1][0]*evec[0][1];
#ifdef MOVETEST
  fprintf(stderr,"evec[2] again\n");
  for(j=0;j<3;j++){
    r_r[j][j]-=eval[2];
    fprintf(stderr,"0=%.2f*%.2f+%.2f*%.2f+%.2f*%.2f=%.2f\n",
    r_r[j][0],evec[2][0],r_r[j][1],evec[2][1],r_r[j][2],evec[2][2],
    r_r[j][0]*evec[2][0]+r_r[j][1]*evec[2][1]+r_r[j][2]*evec[2][2]);
    r_r[j][j]+=eval[2];}
  fprintf(stderr,"eigenvalues :\n");
  for(i=0;i<3;i++)
    fprintf(stderr,"\t%.2f",eval[i]);
  fprintf(stderr,"\neigenvectors :\n");
  for(i=0;i<3;i++){
    for(j=0;j<3;j++)
      fprintf(stderr,"\t%.2f",evec[i][j]);
    fprintf(stderr,"\n");}
#endif
  /* end eigenvalue , eigenvector */
  /* calculate b=R*evec/.. */
  for(i=0;i<2;i++){
    double _y=fabs(eval[i]);
    for(j=0;j<3;j++){
      double _x=0.0;
      for(k=0;k<3;k++)
        _x+=R[j][k]*evec[i][k];
      if(_y<=Zero){
        if(!i){
	  b[0][0]=1.0;
	  b[0][1]=0.0;
	  b[0][2]=0.0;}
	/* assumad i = 0 or 1 */
	/* create b[i] orthogonal to b[0] */
	else{
	  double B=fabs(b[0][0]);
	  double C=fabs(b[0][1]);
	  double D=fabs(b[0][2]);
	  if(B>C&&B>D){
	    b[i][0]= b[0][1];
	    b[i][1]=-b[0][0];
	    b[i][2]=0.0;}
	  else{
	    b[i][0]=0.0;
	    if(C>D){
	      b[i][1]= b[0][2];
	      b[i][2]=-b[0][1];}
	    else{
	      b[i][1]=-b[0][2];
	      b[i][2]= b[0][1];}}
	  _x=sqrt(b[i][0]*b[i][0]+b[i][1]*b[i][1]+b[i][2]*b[i][2]);
	  for(j=0;j<3;j++)
	    b[i][j]/=_x;}}
      else
	b[i][j]=_x/sqrt(_y);}}
  b[2][0]=b[0][1]*b[1][2]-b[1][1]*b[0][2];
  b[2][1]=b[0][2]*b[1][0]-b[0][0]*b[1][2];
  b[2][2]=b[0][0]*b[1][1]-b[1][0]*b[0][1];
#ifdef TEST_MIRROR
  for(j=0;j<3;j++){
    double _x=0.0;
    for(k=0;k<3;k++)
      _x+=R[j][k]*evec[2][k];
    mirror+=b[2][j]*_x;}
  if(mirror<-Zero){
#ifdef PRINT_MIRROR
    fprintf(stderr,"\nmirror transformation !(%f)\n",mirror);
#endif
#ifdef ENABLE_MIRROR
    for(k=0;k<3;k++)
      b[2][k]=-b[2][k];
#endif
  }
#endif
#ifdef MOVETEST
  fprintf(stderr,"R*eigenvectors/sqrt((double)eigenvalues) :\n");
  for(i=0;i<3;i++){
    for(j=0;j<3;j++)
      fprintf(stderr,"\t%.2f",b[i][j]);
    fprintf(stderr,"\n");}
#endif
  /* end calculating b */
  /* calculatin U */
  for(i=0;i<3;i++)
    for(j=0;j<3;j++){
      double _x=0.0;
      for(k=0;k<3;k++)
        _x+=b[k][i]*evec[k][j];
      U[i][j]=_x;}
#ifdef MOVETEST
  fprintf(stderr,"U :\n");
  for(i=0;i<3;i++){
    for(j=0;j<3;j++)
      fprintf(stderr,"\t%.2f",U[i][j]);
    fprintf(stderr,"\n");}
  fprintf(stderr,"U*U :\n");
#endif
#ifdef TEST_ORTHOGONALITY
  { double orth=0.0;
    for(i=0;i<3;i++){
      for(j=0;j<3;j++){
        double _x=0.0;
        for(k=0;k<3;k++)
          _x+=U[i][k]*U[j][k];
#ifdef MOVETEST
        fprintf(stderr,"\t%.2f",_x);
#endif
        if(i!=j)
          orth+=fabs(_x);
        else
          orth+=fabs(_x-1.0);}
#ifdef MOVETEST
      fprintf(stderr,"\n");
#endif
    }
#ifdef MOVETEST
    fprintf(stderr,"orthogonality error %f\n",orth);
#endif
    if(orth>Max_orth){
#ifdef WARN_ORTHOGONALITY_ERROR
      fprintf(stderr,"\northogonality error !(%f)\n",orth);
#endif
      /*U[0][0]=U[1][1]=U[2][2]=1.0;
      U[0][1]=U[1][0]=U[2][0]=U[0][2]=U[1][2]=U[2][1]=0.0;*/
      return(-Huge);}
  }
#endif
  /* end calculation of U */
  /* moving x onto y */
#ifdef MOVETEST
  fprintf(stderr,"xnew and rms:\n");
#endif
  for(i=0;i<n;i++){
    double _a=0.0;
    for(j=0;j<3;j++){
      double _x=0.0;
      for(k=0;k<3;k++)
        _x+=U[j][k]*x[i][k];
      /* moving x_old */
#ifdef WMOVE
      x_old[i][j]=_x+yc[j];
#endif
      _a+=(_x-y[i][j])*(_x-y[i][j])*movew[i];}
    rms+=_a;
#ifdef MOVETEST
    _a=sqrt(_a);
    for(j=0;j<3;j++)
      fprintf(stderr,"\t%.2f",x_old[i][j]);
    fprintf(stderr,"\trms %.2f",_a);
    for(j=0;j<3;j++)
      fprintf(stderr,"\t%.2f",y_old[i][j]);
    fprintf(stderr,"\n");
#endif
  }
  free(x);
  free(y);
  rms=sqrt(rms/wc);
#ifdef ENABLE_MIRROR
  if(mirror<-Zero)
    rms=-rms;
#endif
#ifdef MOVETEST
  fprintf(stderr,"\n -> rms %f\n",rms);
#endif
  return(rms);
}

