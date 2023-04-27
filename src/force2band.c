#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <unistd.h>
#include <string.h>

typedef double mat33_t[3][3];   /* type of 3x3 matrix */
typedef double vec3_t[3];

int usage(char *argv[]){
  printf("%s natom density_file", argv[0]);
  return 0;
}

int main(int argc, char *argv[]){
  const double hbarc = 0.19732e4;
  const double matssunit = 9.315e8;
  const double kb=8.614e-5;
  
  int i,j,k,na,mlen;
  double x,y,z;
  gsl_complex elem, pfactor;
  vec3_t *pos;
  mat33_t basis, kbasis;
  gsl_matrix_complex *force;
  
  FILE *fid = fopen(argv[1], "r");
  int linelen=100;
  char line[linelen];
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      fscanf(fid, "%lf", &(basis[i][j]));
  fscanf(fid, "%d", &na);

  pos = (vec3_t *) malloc(sizeof(vec3_t)*na);
  mlen = na*3;
  force = gsl_matrix_complex_alloc(mlen, mlen);
  fgets(line, linelen, fid);
  for (i=0;i<na;i++){
    fscanf(fid, "%lf %lf %lf", &x, &y, &z);
    pos[i][0] = x*basis[0][0]+y*basis[1][0]+z*basis[2][0];
    pos[i][1] = x*basis[0][1]+y*basis[1][1]+z*basis[2][1];
    pos[i][2] = x*basis[0][2]+y*basis[1][2]+z*basis[2][2];
    fgets(line, linelen, fid);
  }
  fclose(fid);

  fid = fopen("recip_vectors", "r");
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      fscanf(fid, "%lf", &(kbasis[i][j]));
  fclose(fid);
  
  fid = fopen("force_matrix.dat", "r");
  for (i=0;i<mlen;i++)
    for (j=0;j<mlen;j++){
      fscanf(fid, "%lf", &x);
      GSL_REAL(elem)=x;
      GSL_IMAG(elem)=0;
      gsl_matrix_complex_set(force, i,j, elem);
    }
  fclose(fid);

  double T;
  fid = fopen("Trun", "r");
  fscanf(fid, "%lf", &T);
  fclose(fid);
  double beta=1/(T*kb), mass;
  int ntype, atom_num;
  fid = fopen("Mass", "r");
  fscanf(fid, "%d", &ntype);
  fscanf(fid, "%d", &atom_num);
  fscanf(fid, "%lf", &mass);
  fclose(fid);
  
  //  for (i=0;i<na;i++)
  //    printf("%lf %lf %lf\n", pos[i][0], pos[i][1], pos[i][2]);
  
  vec3_t kpoint;
  int s;
  gsl_matrix_complex *dyn, *dyn_cell,*evec;
  dyn = gsl_matrix_complex_alloc(3, 3);
  dyn_cell = gsl_matrix_complex_alloc(3, 3);
  evec = gsl_matrix_complex_alloc(3, 3);
  gsl_eigen_hermv_workspace *w;
  w = gsl_eigen_hermv_alloc(3);
  gsl_vector *eval;
  gsl_vector_complex *vec1,*vec2,*vec3;
  eval = gsl_vector_alloc(3);
  vec1 = gsl_vector_complex_alloc(3);
  vec2 = gsl_vector_complex_alloc(3);
  vec3 = gsl_vector_complex_alloc(3);

  int u,v;
  int npath, ninterv=100;
  double dstart =0, interv;
  double cost,sint,dx, dy, dz, p1,p2,p3, v1, v2, v3;
  double xreal, ximag, x1real, x2real, x3real, x1imag, x2imag, x3imag;
  vec3_t *kpath, vec;
  fid = fopen("kpath", "r");
  fscanf(fid, "%d", &npath);
  kpath = (vec3_t *) malloc(sizeof(vec3_t)*npath);
  for (i=0;i<npath;i++)
    fscanf(fid, "%lf %lf %lf", &(kpath[i][0]), &(kpath[i][1]), &(kpath[i][2]));
  fclose(fid);
  for (i=1;i<npath;i++){
    dx = (kpath[i][0]-kpath[i-1][0]);
    dy = (kpath[i][1]-kpath[i-1][1]);
    dz = (kpath[i][2]-kpath[i-1][2]);
    interv = sqrt(dx*dx+dy*dy+dz*dz);
    dx/=ninterv; dy/=ninterv; dz/=ninterv;
    interv/=ninterv;
    for (j=0;j<=ninterv;j++){
      kpoint[0] = kpath[i-1][0]+dx*j;
      kpoint[1] = kpath[i-1][1]+dy*j;
      kpoint[2] = kpath[i-1][2]+dz*j;
      gsl_matrix_complex_set_zero(dyn);
      for (s=0;s<na;s++){
	for (u=0;u<3;u++)
	  for (v=0;v<3;v++)
	    gsl_matrix_complex_set(dyn_cell, u,v, gsl_matrix_complex_get(force, u, 3*s+v));
	x = pos[s][0]*kpoint[0]*kbasis[0][0]+pos[s][1]*kpoint[1]*kbasis[1][0]+pos[s][2]*kpoint[2]*kbasis[2][0];
	y = pos[s][0]*kpoint[0]*kbasis[0][1]+pos[s][1]*kpoint[1]*kbasis[1][1]+pos[s][2]*kpoint[2]*kbasis[2][1];
	z = pos[s][0]*kpoint[0]*kbasis[0][2]+pos[s][1]*kpoint[1]*kbasis[1][2]+pos[s][2]*kpoint[2]*kbasis[2][2];
	pfactor = gsl_complex_polar(1.0, (x+y+z));
	gsl_matrix_complex_scale(dyn_cell, pfactor);
	gsl_matrix_complex_add(dyn, dyn_cell);
      }
      gsl_eigen_hermv(dyn, eval, evec, w);
      //printf("%lf %lf %lf ", kpoint[0], kpoint[1], kpoint[2]);

      x = kpoint[0]*kbasis[0][0]+kpoint[1]*kbasis[1][0]+kpoint[2]*kbasis[2][0];
      y = kpoint[0]*kbasis[0][1]+kpoint[1]*kbasis[1][1]+kpoint[2]*kbasis[2][1];
      z = kpoint[0]*kbasis[0][2]+kpoint[1]*kbasis[1][2]+kpoint[2]*kbasis[2][2];
      xreal = sqrt(x*x+y*y+z*z);
      if (xreal > 0) {x/=xreal; y/=xreal; z/=xreal;}

      x1real = GSL_REAL(gsl_matrix_complex_get(evec, 0, 0));
      x1imag = GSL_IMAG(gsl_matrix_complex_get(evec, 0, 0));
      x2real = GSL_REAL(gsl_matrix_complex_get(evec, 1, 0));
      x2imag = GSL_IMAG(gsl_matrix_complex_get(evec, 1, 0));
      //xreal = GSL_REAL(gsl_matrix_complex_get(evec, 0, 1));
      //ximag = GSL_IMAG(gsl_matrix_complex_get(evec, 0, 1));
      x3real = GSL_REAL(gsl_matrix_complex_get(evec, 2, 0));
      x3imag = GSL_IMAG(gsl_matrix_complex_get(evec, 2, 0));
      //xreal = GSL_REAL(gsl_matrix_complex_get(evec, 0, 2));
      //ximag = GSL_IMAG(gsl_matrix_complex_get(evec, 0, 2));
      xreal = x1real+x2real+x3real;
      ximag = x1imag+x2imag+x3imag;
      cost = xreal/sqrt(xreal*xreal+ximag*ximag);
      sint = ximag/sqrt(xreal*xreal+ximag*ximag);
      vec[0] = x1real*cost+x1imag*sint;
      vec[1] = x2real*cost+x2imag*sint;
      vec[2] = x3real*cost+x3imag*sint;
      xreal = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
      vec[0] /= xreal;
      vec[1] /= xreal;
      vec[2] /= xreal;
      vec[0] = x1real;
      vec[1] = x2real;
      vec[2] = x3real;
      p1 = fabs(x*vec[0]+y*vec[1]+z*vec[2]);

      x1real = GSL_REAL(gsl_matrix_complex_get(evec, 0, 1));
      x1imag = GSL_IMAG(gsl_matrix_complex_get(evec, 0, 1));
      //xreal = GSL_REAL(gsl_matrix_complex_get(evec, 1, 0));
      //ximag = GSL_IMAG(gsl_matrix_complex_get(evec, 1, 0));
      x2real = GSL_REAL(gsl_matrix_complex_get(evec, 1, 1));
      x2imag = GSL_IMAG(gsl_matrix_complex_get(evec, 1, 1));
      x3real = GSL_REAL(gsl_matrix_complex_get(evec, 2, 1));
      x3imag = GSL_IMAG(gsl_matrix_complex_get(evec, 2, 1));
      //xreal = GSL_REAL(gsl_matrix_complex_get(evec, 1, 2));
      //ximag = GSL_IMAG(gsl_matrix_complex_get(evec, 1, 2));
      xreal = x1real+x2real+x3real;
      ximag = x1imag+x2imag+x3imag;
      cost = xreal/sqrt(xreal*xreal+ximag*ximag);
      sint = ximag/sqrt(xreal*xreal+ximag*ximag);
      vec[0] = x1real*cost+x1imag*sint;
      vec[1] = x2real*cost+x2imag*sint;
      vec[2] = x3real*cost+x3imag*sint;
      xreal = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
      vec[0] /= xreal;
      vec[1] /= xreal;
      vec[2] /= xreal;
      vec[0] = x1real;
      vec[1] = x2real;
      vec[2] = x3real;
      p2 = fabs(x*vec[0]+y*vec[1]+z*vec[2]);

      x1real = GSL_REAL(gsl_matrix_complex_get(evec, 0, 2));
      x1imag = GSL_IMAG(gsl_matrix_complex_get(evec, 0, 2));
      //xreal = GSL_REAL(gsl_matrix_complex_get(evec, 2, 0));
      //ximag = GSL_IMAG(gsl_matrix_complex_get(evec, 2, 0));
      x2real = GSL_REAL(gsl_matrix_complex_get(evec, 1, 2));
      x2imag = GSL_IMAG(gsl_matrix_complex_get(evec, 1, 2));
      //xreal = GSL_REAL(gsl_matrix_complex_get(evec, 2, 1));
      //ximag = GSL_IMAG(gsl_matrix_complex_get(evec, 2, 1));
      x3real = GSL_REAL(gsl_matrix_complex_get(evec, 2, 2));
      x3imag = GSL_IMAG(gsl_matrix_complex_get(evec, 2, 2));
      xreal = x1real+x2real+x3real;
      ximag = x1imag+x2imag+x3imag;
      cost = xreal/sqrt(xreal*xreal+ximag*ximag);
      sint = ximag/sqrt(xreal*xreal+ximag*ximag);
      vec[0] = x1real*cost+x1imag*sint;
      vec[1] = x2real*cost+x2imag*sint;
      vec[2] = x3real*cost+x3imag*sint;
      xreal = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
      vec[0] /= xreal;
      vec[1] /= xreal;
      vec[2] /= xreal;
      vec[0] = x1real;
      vec[1] = x2real;
      vec[2] = x3real;
      p3 = fabs(x*vec[0]+y*vec[1]+z*vec[2]);

      //p1 = fabs(GSL_REAL(gsl_matrix_complex_get(evec,0,0))+GSL_REAL(gsl_matrix_complex_get(evec,1,0)));
      //p2 = fabs(GSL_REAL(gsl_matrix_complex_get(evec,0,1))+GSL_REAL(gsl_matrix_complex_get(evec,1,1)));
      //p3 = fabs(GSL_REAL(gsl_matrix_complex_get(evec,0,2))+GSL_REAL(gsl_matrix_complex_get(evec,1,2)));
      if ((p3>p2) && (p3>p1)) {
	v3=gsl_vector_get(eval, 2);
	gsl_matrix_complex_get_col(vec3, evec,2);
	if ( p2>p1 ){
	  v2 = gsl_vector_get(eval, 1);
	  v1 = gsl_vector_get(eval, 0);
	  gsl_matrix_complex_get_col(vec2, evec,1);
	  gsl_matrix_complex_get_col(vec1, evec,0);
	}
	else{
	  v2 = gsl_vector_get(eval, 0);
	  v1 = gsl_vector_get(eval, 1);
	  gsl_matrix_complex_get_col(vec2, evec,0);
	  gsl_matrix_complex_get_col(vec1, evec,1);
	}
      }
      if ((p2>p3) && (p2>p1)) {
	v3=gsl_vector_get(eval, 1);
	gsl_matrix_complex_get_col(vec3, evec,1);
	if (p3>p1){
	  v2 = gsl_vector_get(eval, 2);
	  v1 = gsl_vector_get(eval, 0);
	  gsl_matrix_complex_get_col(vec2, evec,2);
	  gsl_matrix_complex_get_col(vec1, evec,0);
	}
	else{
	  v2 = gsl_vector_get(eval, 0);
	  v1 = gsl_vector_get(eval, 2);
	  gsl_matrix_complex_get_col(vec2, evec,0);
	  gsl_matrix_complex_get_col(vec1, evec,2);
	}
      }
      if ((p1>p3) && (p1>p2)){
	v3=gsl_vector_get(eval, 0);
	gsl_matrix_complex_get_col(vec3, evec,0);
	if (p3>p2){
	  v2 = gsl_vector_get(eval, 2);
	  v1 = gsl_vector_get(eval, 1);
	  gsl_matrix_complex_get_col(vec2, evec,2);
	  gsl_matrix_complex_get_col(vec1, evec,1);
	}
	else{
	  v2 = gsl_vector_get(eval, 1);
	  v1 = gsl_vector_get(eval, 2);
	  gsl_matrix_complex_get_col(vec2, evec,1);
	  gsl_matrix_complex_get_col(vec1, evec,2);
	}
      }
      gsl_vector_set(eval, 0, v3);
      gsl_vector_set(eval, 1, v2);
      gsl_vector_set(eval, 2, v1);
      gsl_matrix_complex_set_col(evec, 0, vec3);
      gsl_matrix_complex_set_col(evec, 1, vec2);
      gsl_matrix_complex_set_col(evec, 2, vec1);
      //*
      printf("kvector: %lf %lf %lf\n", x, y, z);
      printf("--------------------------------------------------------------------------------\n");
      for (u=0;u<3;u++){
	for (v=0;v<3;v++){
	  xreal = GSL_REAL(gsl_matrix_complex_get(evec, u, v));
	  ximag = GSL_IMAG(gsl_matrix_complex_get(evec, u, v));
	  printf("%lf%+lfi ", xreal, ximag);
	}
	printf("\n");
      }
      printf("--------------------------------------------------------------------------------\n");
      //*/
      //printf("%lf %lf %lf ", p1, p2, p3);
      printf("dispersion: %lf ", dstart+j*interv);
      for (s=0;s<3;s++){
	x = gsl_vector_get(eval, s);
	if ( x < 0)
	  printf("%.8e ", -hbarc*sqrt(-x)/sqrt(matssunit*mass)*1000);
	else
	  printf("%.8e ", hbarc*sqrt(x)/sqrt(matssunit*mass)*1000);
      }
      printf("\n");
    }
    dstart+=interv*ninterv;
  }
  gsl_matrix_complex_free(force);
  gsl_matrix_complex_free(dyn);
}
