#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#include <unistd.h>
#include <string.h>

typedef double mat33_t[3][3]; 	/* type of 3x3 matrix */
typedef double vec3_t[3];
typedef struct pair_type{
  int atom1, atom2;		/* atom1, atom2: dx=atom2-atom1; */
  vec3_t dv;
  int sym;			/* symmetry operation index */
  int weight;
  struct pair_type *next;
}pair_t;

const double epsilon=1e-6;

int Usage(char* cmd){
  printf("Usage %s: -f xyzfile XDATCARs [-s symmetryfile] [-r rcut] [-m sample_max] [-i nskip]\n", cmd);
  printf("-f xyzfile: xyz like file\n" );
  printf("-s symmetryfile: symmetry operation. Use 'aflow-sym' to generate.");
  printf("-m sample_max: sample maximum.");
  printf("-i nskip: skip n sample runs.");
  printf("-r rcut: set rcut and quantum correction fails.");
  return 0;
}

double hypot (const double x, const double y) {
  return sqrt(x*x+y*y);
}

int dir2cart(mat33_t basis, vec3_t dv, vec3_t cv){
  cv[0] = basis[0][0]*dv[0]+basis[1][0]*dv[1]+basis[2][0]*dv[2];
  cv[1] = basis[0][1]*dv[0]+basis[1][1]*dv[1]+basis[2][1]*dv[2];
  cv[2] = basis[0][2]*dv[0]+basis[1][2]*dv[1]+basis[2][2]*dv[2];
  return 0;
}

int symoper(mat33_t basis, vec3_t dv, vec3_t cv){
  cv[0] = basis[0][0]*dv[0]+basis[0][1]*dv[1]+basis[0][2]*dv[2];
  cv[1] = basis[1][0]*dv[0]+basis[1][1]*dv[1]+basis[1][2]*dv[2];
  cv[2] = basis[2][0]*dv[0]+basis[2][1]*dv[1]+basis[2][2]*dv[2];
  return 0;
}

int matrix_transform(mat33_t T, mat33_t A, mat33_t B){ /* B=T^-1 A T; */
  mat33_t Ti, stack;
  int i,j,k;
  for (i=0;i<3;i++)
    for(j=0;j<3;j++){
      Ti[i][j]=T[j][i];
      stack[i][j]=0;
      B[i][j]=0;
    }
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      for (k=0;k<3;k++){
	stack[i][j]+=A[i][k]*T[k][j];
      }
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      for (k=0;k<3;k++){
	B[i][j]+=Ti[i][k]*stack[k][j];
      }
  return 0;  
}

int compare_vector(vec3_t v1, vec3_t v2, double tol){
  int i,j,k;
  for (i=-1;i<2;i++)
    for (j=-1;j<2;j++)
      for (k=-1;k<2;k++)
	if ( (fabs(v1[0]-v2[0]+i)<tol) && (fabs(v1[1]-v2[1]+j)<tol) && (fabs(v1[2]-v2[2]+k)<tol) ) return 1;
  //if ( (fabs(v1[0]-v2[0])<tol) && (fabs(v1[1]-v2[1])<tol) && (fabs(v1[2]-v2[2])<tol) ) return 1;
  //if ( (fabs(v1[0]+v2[0])<tol) && (fabs(v1[1]+v2[1])<tol) && (fabs(v1[2]+v2[2])<tol) ) return 1;
  return 0;
}

int main(int argc, char * argv[]){
  const double pi=3.1425926;
  const double e=2.7183;
  const double kb=8.614e-5;	    /* eV/T */
  const double hbarc=0.197326e4; /* eV*A */
  const double massunit=9.314941e8; /* eV */

  int i,j,k,t,h;
  FILE *fid;
  const int linelen=1000;
  char line[linelen];

  int natom, count;
  double *dxcenter, *dycenter, *dzcenter,  dx, dy, dz, dxsym, dysym, dzsym;		/* centers of atoms */
  double *dxreg, *dyreg, *dzreg;
  double *dxtraj, *dytraj, *dztraj;
  double dxa, dxb, dxc, dya, dyb, dyc;
  int idRx, idRy, idRz;

  int c, sym=0;
  double rc=100000;
  char fxyz[100], fxyzscell[100], symfile[100];
  double tol=1e-4;		/* minimum distance of a pair; unit: Ang */
  int max_sample = 100000, interv=1, residual=0, printmatrix=0;
  while ((c = getopt (argc, argv, "f:r:s:t:m:i:z:p")) != -1)
    switch (c)
      {
      case 'f':
	strcpy(fxyz, optarg);
	break;
      case 'l':
	strcpy(fxyzscell, optarg);
	break;
      case 'r':
        rc = atof(optarg);
        break;
      case 's':
	sym=1;
	strcpy(symfile, optarg);
	break;
      case 't':
	tol=atof(optarg);
	break;
      case 'm':
	max_sample=atoi(optarg);
	break;
      case 'i':
	interv=atoi(optarg);
	break;
      case 'z':
	residual=atoi(optarg);
	break;
      case 'p':
	printmatrix=1;
	break;
      default:
	Usage(argv[0]);
        abort ();
      }

  //++++++++++++++++++++ INPUT SYMMETRY OPERATIONS ++++++++++++++++++++
  int symnum=1;
  mat33_t *symmatrix;
  if (sym){
    fid = fopen(symfile, "r");
    fscanf(fid, "%d", &symnum);
    symmatrix = (mat33_t *)malloc(sizeof(mat33_t)*symnum);
    for (i=0;i<symnum;i++){
      for (j=0;j<3;j++)
	for (k=0;k<3;k++)
	  fscanf(fid, "%lf", &(symmatrix[i][j][k]));
    }
  }
  else {
    symmatrix = (mat33_t *)malloc(sizeof(mat33_t)*symnum);
    for (j=0;j<3;j++)
      for (k=0;k<3;k++)
	if (j==k)
	  symmatrix[0][j][k]=1;
	else
	  symmatrix[0][j][k]=0;
  }

  //-------------------- INPUT SYMMETRY OPERATIONS --------------------

  mat33_t basis, basis_inv, basis_tmp;
  vec3_t *sitelist, dsite, ssite;
  int npair=0; pair_t *pairlist, *pair;
  int *typelist;
  
  fid = fopen(fxyz, "r");
  //for (i=0;i<2;i++) fgets(line, linelen, fid);
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      fscanf(fid, "%lf", &(basis[i][j]));
  fscanf(fid, "%d", &natom);

  sitelist = (vec3_t *) malloc(sizeof(vec3_t)*natom);
  typelist = (int *) malloc(sizeof(int)*natom);
  pairlist = (pair_t *) malloc(sizeof(pair_t)*natom*natom);
  for (i=0;i<natom*natom;i++) {pairlist[i].next=NULL;pairlist[i].weight=1;}
  
  fgets(line, linelen, fid);
  for (i=0;i<natom;i++){
    fscanf(fid, "%lf %lf %lf %d", &(dsite[0]), &(dsite[1]), &(dsite[2]), &(typelist[i]));
    //printf("%lf %lf %lf %d\n", (dx), (dy), (dz), (typelist[i]));
    //dir2cart(basis, dsite, sitelist[i]);
    sitelist[i][0]=dsite[0];
    sitelist[i][1]=dsite[1];
    sitelist[i][2]=dsite[2];
    fgets(line, linelen, fid);
  }

  int idx;
  count =0;
  for (i=0;i<natom;i++)
    for (j=0;j<natom;j++){
      dsite[0]=sitelist[i][0]-sitelist[j][0];
      dsite[1]=sitelist[i][1]-sitelist[j][1];
      dsite[2]=sitelist[i][2]-sitelist[j][2];
      if (dsite[0]<-0.5-tol) dsite[0] ++; if (dsite[0]>0.5+tol) dsite[0]--;
      if (dsite[1]<-0.5-tol) dsite[1] ++; if (dsite[1]>0.5+tol) dsite[1]--;
      if (dsite[2]<-0.5-tol) dsite[2] ++; if (dsite[2]>0.5+tol) dsite[2]--;
      idx=-1;
      if (npair == 0){
	pairlist[npair].dv[0] = dsite[0];
	pairlist[npair].dv[1] = dsite[1];
	pairlist[npair].dv[2] = dsite[2];
	if ( fabs(fabs(dsite[0])-0.5)<tol) pairlist[npair].weight *=2;
	if ( fabs(fabs(dsite[1])-0.5)<tol) pairlist[npair].weight *=2;
	if ( fabs(fabs(dsite[2])-0.5)<tol) pairlist[npair].weight *=2;
	idx = npair;
	npair ++;
      }
      else{
	for (t=0;t<npair;t++){
	  for (k=0;k<symnum;k++){
	    symoper(symmatrix[k], dsite, ssite);
	    if (compare_vector(ssite, pairlist[t].dv, tol)){
	      idx = t;
	      break;
	    }
	    if ( idx >=0) break;
	  }
	}
      }

      if (idx < 0){
	pairlist[npair].dv[0] = dsite[0];
	pairlist[npair].dv[1] = dsite[1];
	pairlist[npair].dv[2] = dsite[2];
	if ( fabs(fabs(dsite[0])-0.5)<tol) pairlist[npair].weight *=2;
	if ( fabs(fabs(dsite[1])-0.5)<tol) pairlist[npair].weight *=2;
	if ( fabs(fabs(dsite[2])-0.5)<tol) pairlist[npair].weight *=2;
	idx = npair;
	npair++;
      }

      for (k=0;k<symnum;k++){
	symoper(symmatrix[k], dsite, ssite);
	if (compare_vector(ssite, pairlist[idx].dv, tol)){
	  pair = (pair_t *) malloc(sizeof(pair_t));
	  pair->atom1=i; pair->atom2=j;
	  pair->dv[0]= dsite[0]; pair->dv[1] = dsite[1]; pair->dv[2] = dsite[2];
	  pair->sym = k;
	  pair->next = pairlist[idx].next;
	  pairlist[idx].next = pair;
	  count ++;
	}
      }
    }

  pair_t * ptr=NULL;
  for (i=0;i<npair;i++){
    printf("pair %d %lf %lf %lf\n", i, pairlist[i].dv[0], pairlist[i].dv[1], pairlist[i].dv[2]);
    count=0;
    ptr = pairlist[i].next;
    while (ptr!=NULL){
      count++;
      ptr=ptr->next;
    }
    printf("%d, count=%d\n", i, count*pairlist[i].weight);
  }
  
  mat33_t *cov;
  int *countlist;
  cov = (mat33_t *) malloc(sizeof(mat33_t)*npair);
  for (i=0;i<npair;i++)
    for (j=0;j<3;j++)
      for (k=0;k<3;k++)
	cov[i][j][k] = 0;
  countlist = (int *) malloc(sizeof(int)*npair);
  //++++++++++++++++++++ CALCULATE AVERAGE POSITION ++++++++++++++++++++
  double temperature, beta;
  int u,v;
  fid = fopen(argv[optind], "r");
  for (i=0;i<2;i++) fgets(line, linelen, fid);
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      fscanf(fid, "%lf", &(basis[i][j]));

  int signum;
  gsl_permutation * p = gsl_permutation_alloc (3);
  gsl_matrix *m1,*m2;
  m1= gsl_matrix_alloc(3,3);
  m2= gsl_matrix_alloc(3,3);
  for (u=0;u<3;u++)
    for(v=0;v<3;v++)
      gsl_matrix_set(m1,u,v,basis[u][v]);
  gsl_linalg_LU_decomp(m1, p, &signum);
  gsl_linalg_LU_invert(m1,p, m2);
  for (u=0;u<3;u++)
    for(v=0;v<3;v++)
      basis_inv[u][v]=gsl_matrix_get(m2,u,v);
  printf("----------------------------------------\n");
  for (u=0;u<3;u++){
  	for (v=0;v<3;v++)
  	  printf("%lf ", (basis_inv[u][v]));
  	printf("\n");
  }
  printf("----------------------------------------\n");
  
  for (i=0;i<2;i++) fgets(line, linelen, fid);
  fscanf(fid, "%d", &natom);
  fclose(fid);

  fid = fopen("Trun", "r");
  fscanf(fid, "%lf", &temperature);
  fclose(fid);
  beta = 1/(temperature * kb);

  int ntype;
  int types[100];
  double masses[100];
  int nsamples=0, ncount=0;
  fid = fopen("Mass", "r");
  fscanf(fid, "%d", &ntype);
  for (i=0;i<ntype;i++) fscanf(fid, "%d", &(types[i]));
  for (i=0;i<ntype;i++) fscanf(fid, "%lf", &(masses[i]));
  fclose(fid);
  
  dxcenter = (double *) malloc(natom*sizeof(double));
  dycenter = (double *) malloc(natom*sizeof(double));
  dzcenter = (double *) malloc(natom*sizeof(double));
  for (i=0;i<natom;i++) {dxcenter[i] =0; dycenter[i]=0; dzcenter[i]=0;}

  
  dxreg = (double *) malloc(natom*sizeof(double));
  dyreg = (double *) malloc(natom*sizeof(double));
  dzreg = (double *) malloc(natom*sizeof(double));
  for (i=0;i<natom;i++) {dxreg[i] =0; dyreg[i]=0; dzreg[i]=0;}

  dxtraj = (double *) malloc(natom*sizeof(double));
  dytraj = (double *) malloc(natom*sizeof(double));
  dztraj = (double *) malloc(natom*sizeof(double));
  for (i=0;i<natom;i++) {dxtraj[i] =0; dytraj[i]=0; dztraj[i]=0;}
  count=0;
  printf("Collecting position list...\n");
  int iargc;
  for (iargc=optind;iargc<argc;iargc++){
    printf("%s\n", argv[iargc]);
    fid = fopen(argv[iargc], "r");
    //for (i=0;i<7;i++) fgets(line, linelen, fid);
    while (! feof(fid)){
      fgets(line, linelen, fid);
      if (feof(fid))break;
      fgets(line, linelen, fid);
      for (u=0;u<3;u++)
	for (v=0;v<3;v++)
	  fscanf(fid, "%lf", &(basis_tmp[u][v]));
      /* printf("----------------------------------------\n"); */
      /* for (u=0;u<3;u++){ */
      /* 	for (v=0;v<3;v++) */
      /* 	  printf("%lf ", (basis_tmp[u][v])); */
      /* 	printf("\n"); */
      /* } */
      /* printf("----------------------------------------\n"); */
      fgets(line, linelen, fid);
      fgets(line, linelen, fid);
      fgets(line, linelen, fid);
      fgets(line, linelen, fid);

      for (i=0;i<natom;i++){
	fscanf(fid, "%lf %lf %lf", &dx, &dy, &dz);
	fgets(line, linelen, fid);
	if ( count >= 1 ){	/* solve for periodic boundary */
	  if ( fabs(dx+1-dxtraj[i])<0.5 ) dx=dx+1;
	  if ( fabs(dx-1-dxtraj[i])<0.5 ) dx=dx-1;
	  if ( fabs(dy+1-dytraj[i])<0.5 ) dy=dy+1;
	  if ( fabs(dy-1-dytraj[i])<0.5 ) dy=dy-1;
	  if ( fabs(dz+1-dztraj[i])<0.5 ) dz=dz+1;
	  if ( fabs(dz-1-dztraj[i])<0.5 ) dz=dz-1;
	}
	dxtraj[i]=dx; dytraj[i]=dy; dztraj[i]=dz;
      }

      ncount ++;
    
      if (ncount % interv == residual){ /* valid sample */
	count++;
	nsamples++;
	for (i=0;i<natom;i++){
	  dxcenter[i]+=dxtraj[i]; dycenter[i]+=dytraj[i]; dzcenter[i]+=dztraj[i];
	}
      }
      if (nsamples > max_sample) break;
    }
    fclose(fid);
    if (nsamples > max_sample) break;
  }
  for (i=0;i<natom;i++){ dxcenter[i]=dxcenter[i]/count; dycenter[i]=dycenter[i]/count; dzcenter[i]=dzcenter[i]/count;}
  
  //-------------------- CALCULATE AVERAGE POSITION --------------------

  nsamples=0; ncount=0;
  mat33_t cov_tmp;
  mat33_t stack;

  int flag1=0;
  for (iargc=optind;iargc<argc;iargc++){
    printf("%s\n", argv[iargc]);
    fid = fopen(argv[iargc], "r");
    //for (i=0;i<7;i++) fgets(line, linelen, fid);
    //fgets(line, linelen, fid);
    flag1=0;
    while (! feof(fid)){
      //if (flag1>0) fgets(line, linelen, fid);
      flag1++;
      fgets(line, linelen, fid);
      if (feof(fid))break;
      fgets(line, linelen, fid);
      for (u=0;u<3;u++)
	for (v=0;v<3;v++)
	  fscanf(fid, "%lf", &(basis_tmp[u][v]));
      /* printf("----------------------------------------\n"); */
      /* for (u=0;u<3;u++){ */
      /* 	for (v=0;v<3;v++) */
      /* 	  printf("%lf ", (basis_tmp[u][v])); */
      /* 	printf("\n"); */
      /* } */
      /* printf("----------------------------------------\n"); */
      fgets(line, linelen, fid);
      fgets(line, linelen, fid);
      fgets(line, linelen, fid);
      fgets(line, linelen, fid);
      
      for (i=0;i<natom;i++){
	fscanf(fid, "%lf %lf %lf", &dx, &dy, &dz);
	//printf("%lf %lf %lf\n", dx, dy, dz);
	fgets(line, linelen, fid);
	dxtraj[i]=dx-dxcenter[i]; dytraj[i]=dy-dycenter[i]; dztraj[i]=dz-dzcenter[i];

	/* solve periodic boundary condition */
	if (fabs(dxtraj[i]+1)<0.5) dxtraj[i] ++; if (fabs(dxtraj[i]-1)<0.5) dxtraj[i]--;
	if (fabs(dytraj[i]+1)<0.5) dytraj[i] ++; if (fabs(dytraj[i]-1)<0.5) dytraj[i]--;
	if (fabs(dztraj[i]+1)<0.5) dztraj[i] ++; if (fabs(dztraj[i]-1)<0.5) dztraj[i]--;
	//if ( i==0) printf("%lf\n", dztraj[2]);  
      }

      //++++++++++++++++++++++ reduce drift ++++++++++++++++++
      dx=0; dy=0; dz=0;
      for (i=0;i<natom;i++){
	dx += dxtraj[i]; dy += dytraj[i]; dz += dztraj[i];
      }
      dx /= natom; dy /= natom; dz /= natom;
      //printf("drift: %.4e %.4e %.4e\n", dx, dy, dz);
      for (i=0;i<natom;i++){
	dxtraj[i] -= dx; dytraj[i] -= dy; dztraj[i] -= dz;
      }
      dxtraj[0] = 0; dytraj[0] = 0; dztraj[0] = 0;
      for (i=1;i<natom;i++){
	dxtraj[0] -= dxtraj[i]; dytraj[0] -= dytraj[i]; dztraj[0] -= dztraj[i];
      }
      dx=0; dy=0; dz=0;
      for (i=0;i<natom;i++){
	dx += dxtraj[i]; dy += dytraj[i]; dz += dztraj[i];
      }
      dx /= natom; dy /= natom; dz /= natom;
      //printf("drift after: %.4e %.4e %.4e\n", dx, dy, dz);
      //---------------------- reduce drift ------------------
      ncount ++;
    
      if (ncount % interv == residual){ /* valid samples */
	count++;
	nsamples++;
	dx=0;
	for (t=0;t<npair;t++){
	  ptr = pairlist[t].next;
	  while (ptr != NULL ){
	    i = ptr->atom1;
	    j = ptr->atom2;
	    k = ptr->sym;
	    dxa = dxtraj[i]*symmatrix[k][0][0] +dytraj[i]*symmatrix[k][0][1] +dztraj[i]*symmatrix[k][0][2];
	    dxb = dxtraj[i]*symmatrix[k][1][0] +dytraj[i]*symmatrix[k][1][1] +dztraj[i]*symmatrix[k][1][2];
	    dxc = dxtraj[i]*symmatrix[k][2][0] +dytraj[i]*symmatrix[k][2][1] +dztraj[i]*symmatrix[k][2][2];
	    dya = dxtraj[j]*symmatrix[k][0][0] +dytraj[j]*symmatrix[k][0][1] +dztraj[j]*symmatrix[k][0][2];
	    dyb = dxtraj[j]*symmatrix[k][1][0] +dytraj[j]*symmatrix[k][1][1] +dztraj[j]*symmatrix[k][1][2];
	    dyc = dxtraj[j]*symmatrix[k][2][0] +dytraj[j]*symmatrix[k][2][1] +dztraj[j]*symmatrix[k][2][2];

	    stack[0][0]=dxa*dya;
	    stack[0][1]=dxa*dyb;
	    stack[0][2]=dxa*dyc;
	    stack[1][0]=dxb*dya;
	    stack[1][1]=dxb*dyb;
	    stack[1][2]=dxb*dyc;
	    stack[2][0]=dxc*dya;
	    stack[2][1]=dxc*dyb;
	    stack[2][2]=dxc*dyc;
	    
	    
	    for ( u=0;u<3;u++)
	      for ( h=0;h<3;h++)
		cov_tmp[u][h]=0;
	    for ( u=0;u<3;u++)
	      for ( h=0;h<3;h++)
		for ( v=0;v<3;v++)
		  cov_tmp[u][h] += stack[u][v]*basis_tmp[v][h];
	    for ( u=0;u<3;u++)
	      for ( h=0;h<3;h++)
		stack[u][h]=0;
	    for ( u=0;u<3;u++)
	      for ( h=0;h<3;h++)
		for ( v=0;v<3;v++)
		  stack[u][h] += basis_tmp[v][u]*cov_tmp[v][h];

	    for ( u=0;u<3;u++)
	      for ( h=0;h<3;h++)
		cov_tmp[u][h]=0;
	    for ( u=0;u<3;u++)
	      for ( h=0;h<3;h++)
		for ( v=0;v<3;v++)
		  cov_tmp[u][h] += stack[u][v]*basis_inv[h][v];
	    for ( u=0;u<3;u++)
	      for ( h=0;h<3;h++)
		stack[u][h]=0;
	    for ( u=0;u<3;u++)
	      for ( h=0;h<3;h++)
		for ( v=0;v<3;v++)
		  stack[u][h] += basis_inv[u][v]*cov_tmp[v][h];


	    cov[t][0][0]+=stack[0][0];
	    cov[t][0][1]+=stack[0][1];
	    cov[t][0][2]+=stack[0][2];
	    cov[t][1][0]+=stack[1][0];
	    cov[t][1][1]+=stack[1][1];
	    cov[t][1][2]+=stack[1][2];
	    cov[t][2][0]+=stack[2][0];
	    cov[t][2][1]+=stack[2][1];
	    cov[t][2][2]+=stack[2][2];
	    dx+=dxa*dyc;
	    countlist[t]++;
	    ptr=ptr->next;
	  }
	}
	//printf("test dx=%.8e\n", dx);
      }
      if (nsamples > max_sample) break;
    }
    fclose(fid);
    if (nsamples > max_sample) break;
  }

  for (t=0;t<npair;t++)
    for (i=0;i<3;i++)
      for(j=0;j<3;j++){
	cov[t][i][j] /= countlist[t];
	//cov[t][i][j] /= countlist[t]*pairlist[t].weight;
      }
  dx=0;
  for (t=0;t<npair;t++)
    dx += cov[t][0][0];
  //printf("test dx=%.4e\n", dx);
  
  int matrixlen=natom*3;
  gsl_matrix *mx=gsl_matrix_alloc(matrixlen,matrixlen);

  for (i=0;i<matrixlen;i++)
    for (j=0;j<matrixlen;j++)
      gsl_matrix_set(mx,i,j,0);

  fid = fopen("covariance_pair.dat", "w");
  for (t=0;t<npair;t++){
    ptr = pairlist[t].next;
    i = ptr->atom1;
    j = ptr->atom2;
    k = ptr->sym;
    matrix_transform(symmatrix[k], cov[t], stack);
    for ( u=0;u<3;u++)
      for ( h=0;h<3;h++)
	cov_tmp[u][h]=0;
    for ( u=0;u<3;u++)
      for ( h=0;h<3;h++)
	for ( v=0;v<3;v++)
	  cov_tmp[u][h] += stack[u][v]*basis[v][h];
    for ( u=0;u<3;u++)
      for ( h=0;h<3;h++)
	stack[u][h]=0;
    for ( u=0;u<3;u++)
      for ( h=0;h<3;h++)
	for ( v=0;v<3;v++)
	  stack[u][h] += basis[v][u]*cov_tmp[v][h];
    fprintf(fid, "No. %d\n", t);
    dir2cart(basis, pairlist[t].dv, ssite);
    fprintf(fid, "r2-r1=%lf %lf %lf #directional space\n", ssite[0], ssite[1], ssite[2]);
    dy = sqrt(ssite[0]*ssite[0] +ssite[1]*ssite[1] +ssite[2]*ssite[2]);
    fprintf(fid, "|r2-r1|=%lf\n", dy);
    dx=0;
    for (u=0;u<3;u++)
      for (v=0;v<3;v++) dx += fabs(stack[u][v]);
    fprintf(fid, "cov matrix M: |M|= %lf\n", dx);
    fprintf(fid, "%.6e %.6e %.6e\n", stack[0][0], stack[0][1], stack[0][2]);
    fprintf(fid, "%.6e %.6e %.6e\n", stack[1][0], stack[1][1], stack[1][2]);
    fprintf(fid, "%.6e %.6e %.6e\n", stack[2][0], stack[2][1], stack[2][2]);

    if ( dy > rc ) continue;
    while (ptr != NULL ){
      i = ptr->atom1;
      j = ptr->atom2;
      k = ptr->sym;
      matrix_transform(symmatrix[k], cov[t], stack);

      for ( u=0;u<3;u++)
	for ( h=0;h<3;h++)
	  cov_tmp[u][h]=0;
      for ( u=0;u<3;u++)
	for ( h=0;h<3;h++)
	  for ( v=0;v<3;v++)
	    cov_tmp[u][h] += stack[u][v]*basis[v][h];
      for ( u=0;u<3;u++)
	for ( h=0;h<3;h++)
	  stack[u][h]=0;
      for ( u=0;u<3;u++)
	for ( h=0;h<3;h++)
	  for ( v=0;v<3;v++)
	    stack[u][h] += basis[v][u]*cov_tmp[v][h];

      for (u=0;u<3;u++)
	for (v=0;v<3;v++){
	  gsl_matrix_set(mx, i*3+u, j*3+v, stack[u][v]);
	}
      ptr=ptr->next;
    }    
  }
  for (i=0;i<matrixlen;i++)
    for (j=0;j<matrixlen;j++){
      dx = (gsl_matrix_get(mx, i, j)+gsl_matrix_get(mx, j, i))/2;
      gsl_matrix_set(mx, i, j, dx);
      gsl_matrix_set(mx, j, i, dx);
    }

  printf("NPAIR=%d\n", npair);
  fid=fopen("density_matrix_row1.dat", "w");
  for (j=0;j<matrixlen;j++)
    fprintf(fid, "%.8e\n", gsl_matrix_get(mx, 0,j));
  fclose(fid);

  fid=fopen("density_matrix.dat", "w");
  for (i=0;i<matrixlen;i++){
    for (j=0;j<matrixlen;j++)
      fprintf(fid, "%.8e ", gsl_matrix_get(mx, i,j));
    fprintf(fid, "\n");
  }
  fclose(fid);

  /* for (i=0;i<matrixlen;i++){ */
  /*   dx=0; */
  /*   dy=0; */
  /*   for (j=1;j<matrixlen;j=j+3){ */
  /*     dx+=gsl_matrix_get(mx, i,j); */
  /*     dy+=gsl_matrix_get(mx, j,i); */
  /*   } */
  /*   printf("%d sum=%.8e %.8e\n",i,dx, dy); */
  /* } */

  gsl_vector *diag;
  gsl_matrix *eval, *exec, *force;
  diag = gsl_vector_alloc(matrixlen);
  exec = gsl_matrix_alloc(matrixlen, matrixlen);
  eval = gsl_matrix_alloc(matrixlen, matrixlen);
  force = gsl_matrix_alloc(matrixlen, matrixlen);
  
  gsl_eigen_symmv_workspace *w;
  w = gsl_eigen_symmv_alloc(matrixlen);
  //gsl_eigen_symm(mx, diag, w);
  gsl_eigen_symmv(mx, diag, exec, w);
  gsl_eigen_symmv_free(w);

  gsl_matrix_set_zero(eval);
  for (i=0;i<matrixlen;i++){
    dx = gsl_vector_get(diag,i);
    if (fabs(dx) < 1e-10)
      gsl_matrix_set(eval, i,i, 0);
    else
      gsl_matrix_set(eval, i,i, 1/beta/gsl_vector_get(diag, i));
  }

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, exec, eval, 0, force);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, force, exec, 0, eval);
  //gsl_linalg_ldlt_decomp (mx);

  fid=fopen("force_matrix.dat", "w");
  for (i=0;i<matrixlen;i++){
    for (j=0;j<matrixlen;j++)
      fprintf(fid, "%.15e ", gsl_matrix_get(eval, i,j));
    fprintf(fid, "\n");
  }
  fclose(fid);

  double mass1, mass2;
  for (i=0;i<natom;i++)
    for (j=0;j<natom;j++){
      for ( k=0;k<ntype;k++ )
	if ( types[k] == typelist[i] ){
	  mass1 = masses[k];
	  break;
	}
      for ( k=0;k<ntype;k++ )
	if ( types[k] == typelist[j] ){
	  mass2 = masses[k];
	  break;
	}
    }

  double lambda=sqrt(hbarc*hbarc*2*pi/mass1/massunit*beta);
  double detALL=0, sigma2, omegahbar;
  double Squantum=0;
  printf("Lambda=%lf\n", lambda);
  for (i=0;i<matrixlen;i++){
    //sigma2 = gsl_matrix_get(mx, i, i);
    sigma2 = gsl_vector_get(diag,i);
    for ( k=0;k<ntype;k++ )
      if ( types[k] == typelist[i] ){
	mass1 = masses[k];
	break;
      }
    //detALL += log(fabs((2*pi*e)*omega2))+log(2*pi*e*mass1*massunit*beta/(hbarc*hbarc*(2*pi*2*pi)));
    //if ( omega2 <= 0) continue;
    if (sigma2 <=1e-10) continue;
    detALL += log(fabs(sigma2));
    omegahbar = hbarc/sqrt(sigma2*beta*mass1*massunit);
    Squantum += -log(1-exp(-beta*omegahbar))+beta*omegahbar/(exp(beta*omegahbar)-1);
  }
  //*/
  //------------------------------ CLASSIC ------------------------------
	 
  lambda=sqrt(hbarc*hbarc*2*pi/mass1/massunit*beta);

  fid=fopen("Eigenvalues.dat", "w");
  for (i=0;i<matrixlen;i++){
    //sigma2 = gsl_matrix_get(mx, i, i);
    sigma2 = gsl_vector_get(diag,i);
    omegahbar = hbarc/sqrt(sigma2*beta*mass1*massunit);
    fprintf(fid, "%.8e %.8e\n", sigma2, omegahbar);
  }
  fclose(fid);
  /* detALL=gsl_linalg_LU_det(m, s); */

  //printf("rc=%lf\tS=%.4e\n", rc, 0.5*detALL/natom_sc);
  printf("T=%.4lf\tln(DET)=%.4lf\tS_classic=%.4lf\tS_quantum=%.4lf\n", temperature, detALL/natom, 0.5*(detALL/natom-2*3*log(lambda)+3*log(2*pi))+3-(double)3/natom, Squantum/natom);
  gsl_matrix_free(mx);
  free(dxcenter); free(dycenter); free(dzcenter);
  free(dxreg); free(dyreg); free(dzreg);
  free(dxtraj); free(dytraj); free(dztraj);

  return 0;
}
