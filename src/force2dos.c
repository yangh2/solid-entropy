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
  gsl_matrix_complex *force, *force2;
  gsl_vector *eigen;
  
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
  force2 = gsl_matrix_complex_alloc(mlen, mlen);
  eigen = gsl_vector_alloc(mlen);
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
      gsl_matrix_complex_set(force2, i,j, elem);
    }
  fclose(fid);
  double hmax=0;
  

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
  
  int kdim=atoi(argv[2]);
  vec3_t kpoint;
  int s;
  gsl_matrix_complex *dyn, *dyn_cell;
  gsl_eigen_herm_workspace *w;
  gsl_vector *eval;

  w = gsl_eigen_herm_alloc(mlen);
  gsl_eigen_herm(force2, eigen, w);
  for (i=0;i<mlen;i++) {
    x = gsl_vector_get(eigen, i);
    hmax = x>hmax?x:hmax;
  }
  gsl_eigen_herm_free(w);

  dyn = gsl_matrix_complex_alloc(3, 3);
  dyn_cell = gsl_matrix_complex_alloc(3, 3);
  w = gsl_eigen_herm_alloc(3);
  eval = gsl_vector_alloc(3);
  
  int nbins=400, count=0;
  double hist[nbins]; 		/* histogram: [-0.5*hmax,1.5*hmax] */
  for (i=0;i<nbins;i++) hist[i]=0;
  hmax=hbarc*sqrt(hmax)/sqrt(matssunit*mass)*1000;
  double hwidth=2*hmax/nbins, omega;
  
  int u,v;
  for (i=0;i<kdim;i++)
    for (j=0;j<kdim;j++)
      for (k=0;k<kdim;k++){
	kpoint[0] = (double) i/kdim;
	kpoint[1] = (double) j/kdim;
	kpoint[2] = (double) k/kdim;
	gsl_matrix_complex_set_zero(dyn);
	for (s=0;s<na;s++){
	  for (u=0;u<3;u++)
	    for (v=0;v<3;v++)
	      gsl_matrix_complex_set(dyn_cell, u,v, gsl_matrix_complex_get(force, u, 3*s+v));
	  x = pos[s][0]*kpoint[0]*kbasis[0][0]+pos[s][1]*kpoint[1]*kbasis[1][0]+pos[s][2]*kpoint[2]*kbasis[2][0];
	  y = pos[s][0]*kpoint[0]*kbasis[0][1]+pos[s][1]*kpoint[1]*kbasis[1][1]+pos[s][2]*kpoint[2]*kbasis[2][1];
	  z = pos[s][0]*kpoint[0]*kbasis[0][2]+pos[s][1]*kpoint[1]*kbasis[1][2]+pos[s][2]*kpoint[2]*kbasis[2][2];
	  pfactor = gsl_complex_polar(1.0, x+y+z);
	  gsl_matrix_complex_scale(dyn_cell, pfactor);
	  gsl_matrix_complex_add(dyn, dyn_cell);
	}
	gsl_eigen_herm(dyn, eval, w);
	for (s=0;s<3;s++){
	  x = gsl_vector_get(eval, s);
	  if ( x < 0)
	    omega=-hbarc*sqrt(-x)/sqrt(matssunit*mass)*1000;
	  else
	    omega=hbarc*sqrt(x)/sqrt(matssunit*mass)*1000;
	  hist[(int)((omega+0.5*hmax)/hwidth)]++;
	  count ++;
	}
      }
  for (i=0;i<nbins;i++) hist[i]/=(count*2*hmax/nbins/3);
  for (i=0;i<nbins;i++)
    printf("%lf %lf\n", hwidth*i-0.5*hmax, hist[i]);
  gsl_matrix_complex_free(force);
  gsl_matrix_complex_free(dyn);
}
