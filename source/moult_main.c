/*
c------------------------------------------------------------------------------
c***************** MAIN PROGRAM - MOULT/ MIGRATION VERSION ********************
c------------------------------------------------------------------------------
c
c AUTHOR:  BOB WELHAM, UNIVERSITY OF BRISTOL			OCTOBER 1999
c STRICTLY CONFIDENTIAL - ALL RIGHTS RESERVED

rewritten in C by me
*/


#include"moult.h"


#define MAXN 100

double ucritcl[MT][ME+1];
double unocare[MR+1][MF-MO+1][MF-MO+1][ME+1];
double vnocare[MR+1][MF-MO+1][MF-MO+1][ME+1];

double repval[MT+1][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
double states[MT+1][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
double youngs[MT+1][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
double uvalue[MT][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
int    action[MT][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
int    moult1[MT][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
int    moult2[MT][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];

extern double G[MT+1][ME+1];
extern double G_o[MT+1][ME+1];
extern double theta[ME+1];
extern double FOOD[MT+1];

type_par P;

long idum=(-13);

int Sz;

FILE *OutPut3;

extern void maxuv(int, int, int, int, int, int, int,
	   double, double, double *, double *);
extern void initialise(int, double, int);
extern void savevalues_back(void);
extern void savevalues_for(double);
extern void normalise(double *, double *, int, int);
extern void stepforwards(int, int, int, double);
extern void iniFP(void);
extern void inigain(double);
extern void Ind(char *, int);
extern void ReadIn(char *);

void maxV(int, int, int, int, int, int, double *, double *, int *);
void Backward(void);
void Forward(double *, double *);
void old_Forward(void);
void OutPut(char *);
void Cohort(char *);
void inityoung(double, int);
void PolicyOut(char *);
void BinWrite(char *);
void hiba(char *);
void Treatment(int, t_ind *);
void Food_Treatment(int, int);
void dump_repval(void);
void dump_states(void);
double back_lambda(double);
void Update_food(void);

/************************* main *************************/
int main(int argc, char *argv[])
{
  char filename[30];
  char filen[30];
  double lambda_back, diff_back;
  double lambda_for = 1.0, diff_for;
  double prev3lam;
  int I;

  idum = (long) time(NULL);
  if (idum > 0) idum = -idum;
  
  if (argc == 1) strcpy(filename, "moult");
  else if (argc == 2) {
    if (strcmp(argv[argc-1],"2bin")!=0) strcpy(filename, argv[1]);
    else strcpy(filename, "moult");
  } else {
    strcpy(filename, argv[1]);
  }
  printf("filename: %s\n", filename);
  strcpy(filen,""); strcat(filen,filename); strcat(filen,".sum");
  if ((OutPut3 = fopen(filen,"wt")) == NULL) {
    printf("ERROR: Can't open the output file\n");
    exit (1);
  }
  fprintf(OutPut3,"Seed: %ld\n", idum);
  ReadIn(filename);
  iniFP();
  initialise(0, 1.0, BACKWARD);
  initialise(0, 1.0, FORWARD);
  lambda_back = 100.0;
  Sz = 0.0;  
  inigain(1.0);
  Backward(); 
  normalise(&lambda_back,&diff_back,MAX,NORM);
  printf("backpass done, sz= %d lambda= %f diff = %f n= %f\n", 
	 Sz, lambda_back, diff_back, P.n); 
  do {
    /*    do {
      Sz++;
      prev2lam = lambda_back;*/
/*      Update_food();*/
    Backward(); 
    normalise(&lambda_back,&diff_back,MAX,NORM);
    printf("backpass done, sz= %d lambda= %f diff = %f n= %f\n", 
	   Sz, lambda_back, diff_back, P.n); 
      /*if(Sz > 1000) {
	fprintf(OutPut3, "ERROR: no convergence, Sz= %d\n", Sz);
	break;
      }
    } while (diff_back > P.tol && (fabs(lambda_back-1.0)>P.tol ? 
    fabs(lambda_back-prev2lam) > limit*P.tol : 1));*/
    I = 0;
    do {
      Sz++;
      I++;
      if ((I % 17) == 0) prev3lam = diff_for;
      Forward(&lambda_back,&diff_for);
      normalise(&lambda_for, &diff_for, SUM,NONORM);
      printf("forwardpass done, lambda=%f, diff=%f, n=%f\n", 
	     lambda_for, diff_for, P.n);
      if(fabs(prev3lam-diff_for)< 0.1*P.tol) {
	fprintf(OutPut3, "ERROR: cycling, Sz= %d\n", Sz);
	break;
      }
    } while (diff_for > P.tol);
    if((Sz > 1000) || ((lambda_for < P.tol) && (lambda_back < 1.0))) {
      fprintf(OutPut3, "ERROR: no convergence, Sz= %d\n", Sz);
      break;
    }
  } while (/*((diff_back+diff_for) > 2.0*P.tol) && */
     (fabs(1.0-lambda_back) > P.tol));
  old_Forward();
  PolicyOut(filename);
  if (strcmp(argv[argc-1], "2bin")==0) BinWrite(filename);
  Ind(filename, -1);
  Cohort(filename);
  OutPut(filename);
  fclose(OutPut3);
  return 0;
}
/************************* end of main ******************/

#define simit 0.1

/************************* Update_food *************************/
void Update_food(void)
{
  int t, e;

  for(t=0; t<=MT; t++)
    for(e=0; e<=ME; e++) {
      G_o[t][e] = G[t][e] = G[t][e]*simit + G_o[t][e]*(1.0-simit);
    }
}
/************************* end of Update_food ******************/



/************************* Backward *************************/
void Backward(void)
{
  int t, e, a, r, f1, f2;
  double vmax, umax, bv, bu;
  int bm1, bm2, ba, acmax, m1, m2;
  int Kapcsolo = 0;
  int ex;
  
  for(t = 0; t < MT; t++)
    for(ex = 0; ex <= ME; ex++){
	ucritcl[t][ex] = P.Cost.keep/G[t][ex];
    }
  savevalues_back();
  for(t = MT-1; t >= 0; t--)
  for(a = 0; a <= MA; a++)
  for(f2 = MO; f2 <= MF; f2++)
  for(f1 = MO; f1 <= MF; f1++)
  for(r = 0; r <= MR; r++)
  for(e = 0; e <= ME; e++)	{
    bv = bu = 0.0;
    ba = DELAY;
    bm1 = bm2 = NOMOULT;
    m2 = NOMOULT-1;
    do {
      if (f2>=0) m2++;
      else m2 = NOMOULT;
      m1 = NOMOULT-1;
      do {
	if (f1 >= 0) m1++;
	else m1 = NOMOULT;
	if (m1==MOULT1 && m2==MOULT1) {
	  vmax = repval[t][r][MO-MO][MO-MO][e][a];
	  umax = uvalue[t][r][MO-MO][MO-MO][e][a];
	  acmax = action[t][r][MO-MO][MO-MO][e][a];
	}
	else if (m1==MOULT1) {
	  vmax = repval[t][r][MO-MO][f2-MO][e][a];
	  umax = uvalue[t][r][MO-MO][f2-MO][e][a];
	  acmax = action[t][r][MO-MO][f2-MO][e][a];
	}
	else if (m2==MOULT1) {
	  vmax = repval[t][r][f1-MO][MO-MO][e][a];
	  umax = uvalue[t][r][f1-MO][MO-MO][e][a];
	  acmax = action[t][r][f1-MO][MO-MO][e][a];
	}
	else maxV(t,r,f1,f2,e,a,&vmax,&umax,&acmax);
	if (vmax > bv) {
	  /* if vmax == 0; acmax is not defined very accurately here !!!*/
	  bv = vmax;
	  bu = umax;
	  ba = acmax;
	  bm1 = m1;
	  bm2 = m2;
	}
      } while (m1 < MOULT1 && f1>=0);
    } while (m2 < MOULT1 && f2>=0);
    repval[t][r][f1-MO][f2-MO][e][a] = bv;
    uvalue[t][r][f1-MO][f2-MO][e][a] = bu;
    action[t][r][f1-MO][f2-MO][e][a] = ba;
    moult1[t][r][f1-MO][f2-MO][e][a] = bm1;
    moult2[t][r][f1-MO][f2-MO][e][a] = bm2;
  }
  Kapcsolo = 0;
  if (Kapcsolo == 1) {dump_repval();}
}
/************************* end of Backward *****************/


/************************* Forward *************************/
void Forward(double *lambda, double *diff)
{
  int t;
  double l, l3;

  l = *lambda;
  l3 = (l-1.0) * (l-1.0) * (l-1.0);
  l = l + 10.0*l3;
  savevalues_for(l);
  for(t = 0; t < MT; t++) stepforwards(t, ALL, INIT, 1.0);
}
/************************* end of Forward ******************/

/************************* old_Forward *************************/
void old_Forward(void)
{
  int Kapcsolo, t;
  double lambda, diff;

  /*  initialise(0,1.0,FORWARD);

  do {
    savevalues_for(1.0);
    for(t = 0; t < MT; t++) stepforwards(t, ALL, INIT, 1.0);
    normalise(&lambda, &diff, SUM,NORM);
  } while (diff > P.tol);
  */
  Kapcsolo = 0;
  if (Kapcsolo == 1) dump_repval();
  inityoung(0.0, INIT);
  savevalues_for(1.0);
  for(t = 0; t < MT; t++) stepforwards(t, ALL, INIT, 1.0);
  normalise(&lambda, &diff, SUM,NONORM);
  fprintf(OutPut3,"forwards convergence, lambda=%f, diff=%f, n=%f\n", 
	  lambda, diff, P.n);
  
}
/************************* end of old_Forward ******************/

/************************* inityoung *************************/
void inityoung(double initval, int code)
{
  int t, r, f1, f2, e, a;
  
  for(t = 0; t <= MT; t++)
  for(r = 0; r <= MR; r++)
  for(f1 = MO; f1 <= MF; f1++)
  for(f2 = MO; f2 <= MF; f2++)
  for(e = 0; e <= ME; e++)
  for(a = 0; a <= MA; a++){
    if (code == INIT) {
      if (r == 0) youngs[t][r][f1-MO][f2-MO][e][a] = 0.0;
      else youngs[t][r][f1-MO][f2-MO][e][a] = initval;
    }
    else {
      states[t][r][f1-MO][f2-MO][e][a] = 
	youngs[t][r][f1-MO][f2-MO][e][a];
    }
  }
}
/************************* end of inityoung ******************/

/************************* Cohort *************************/
void Cohort(char *filename)
{
  char filen[30];
  int year, t, r, f1, f2, e, a, act, i, mol1, mol2;
  double s, total, absr, absf1, absf2, Nstart, Naban, Nm1, Nm2;
  double Ntot;
  double alive, breed;
  double moult_1, moult_2;
  double start, aban;
  double res[2], fe1[2], fe2[2], uv[2], ex[2];
  FILE *OutPut1, *OutPut2;
  

  strcpy(filen,""); strcat(filen,filename); strcat(filen,".C.pop");
  if ((OutPut1 = fopen(filen,"wt")) == NULL) {
    printf("ERROR: Can't open the output file: %s\n", filen);
    exit (1);
  }
  strcpy(filen,""); strcat(filen,filename); strcat(filen,".C.sta");
  if ((OutPut2 = fopen(filen,"wt")) == NULL) {
    printf("ERROR: Can't open the output file: %s\n", filen);
    exit (1);
  }
  inityoung(0.0, COPY);
  Ntot = Nstart = Naban = Nm1 = Nm2 = 0.0;
  for(t = 0; t < MT; t++){
    total = 0.0;
    alive = breed = 0.0;
      moult_1 = moult_2 = start = aban = 0.0;
      for(i = 0; i < 2; i++){
	res[i] = fe1[i] = fe2[i] = 0.0;
	uv[i] = ex[i] = 0.0;
      }
    for(a = 0; a <= MA; a++)
    for(r = 0; r <= MR; r++){
      absr = ((double) r)/TOPR;
    for(f1 = MO; f1 <= MF; f1++){
      absf1 = ((double) f1)/TOPF;
    for(f2 = MO; f2 <= MF; f2++){
      absf2 = ((double) f2)/TOPF;
    for(e = 0; e <= ME; e++){
      total += (s = states[t][r][f1-MO][f2-MO][e][a]);
      act = action[t][r][f1-MO][f2-MO][e][a];
      mol1 = moult1[t][r][f1-MO][f2-MO][e][a];
      mol2 = moult2[t][r][f1-MO][f2-MO][e][a];
      alive += s;
      if (mol1 == MOULT1) moult_1 += s;
      else if (f1<0) moult_1 += s;
      if (mol2 == MOULT1) moult_2 += s;
      else if (f2<0) moult_2 += s;
      if (a != 0) breed += s;
      if (act == START) start += s;
      else if (act == ABANDON) aban += s;
      if (a == 0){
	res[0] += absr*s;
	fe1[0] += absf1*s;
	fe2[0] += absf2*s;
	uv[0] += uvalue[t][r][f1-MO][f2-MO][e][a]*s;
	ex[0] += e*s;
      }
      else {
	res[1] += absr*s;
	fe1[1] += absf1*s;
	fe2[1] += absf2*s;
	uv[1] += uvalue[t][r][f1-MO][f2-MO][e][a]*s;
	ex[1] += e*s;
      }
    }
    }
    }
    }
    fprintf(OutPut1, "%d ", t);
    fprintf(OutPut2, "%d ", t);
    Ntot += alive;
    if (alive > 0.0) {
      Nstart += start/alive;
      Naban += aban/alive;
      Nm1 += moult_1/alive;
      Nm2 += moult_2/alive;
    }
    fprintf(OutPut1, "%f %f %f %f %f %f ", 
	    alive, breed, start, aban, moult_1, moult_2);
    fprintf(OutPut2, "%f %f ", alive, breed);
    for(i = 0; i < 2; i++){
      fprintf(OutPut2, "%f %f %f %f %f ", res[i], fe1[i], fe2[i], 
	      uv[i], ex[i]);
    }
    
    fprintf(OutPut1, "\n");
    fprintf(OutPut2, "\n");
  }
  fprintf(OutPut1, "\n\n");
  fprintf(OutPut2, "\n\n");
  fprintf(OutPut3, "%f %f %f %f %f\n",
	  Ntot, Nstart, Naban, Nm1, Nm2);
  for(year = 0; year <= 5; year++){
    Ntot = Nstart = Naban = Nm1 = Nm2 = 0.0;
    for(t = 0; t < MT-1; t++){
      if (year == 0) stepforwards(t, COPY, COPY, 1.0);
      else stepforwards(t, COPY, INIT, 1.0);
      /*      savevalues_for(1.0);*/
    }
    for(t = 0; t < MT; t++){
      total = 0.0;
      alive = breed = 0.0;
      moult_1 = moult_2 = start = aban = 0.0;
      for(i = 0; i < 2; i++){
	res[i] = fe1[i] = fe2[i] = 0.0;
	uv[i] = ex[i] = 0.0;
      }
      for(a = 0; a <= MA; a++)
      for(r = 0; r <= MR; r++){
	absr = ((double) r)/TOPR;
      for(f1 = MO; f1 <= MF; f1++){
	absf1 = ((double) f1)/TOPF;
      for(f2 = MO; f2 <= MF; f2++){
	absf2 = ((double) f2)/TOPF;
      for(e = 0; e <= ME; e++){
	total += (s = states[t][r][f1-MO][f2-MO][e][a]);
	act = action[t][r][f1-MO][f2-MO][e][a];
	mol1 = moult1[t][r][f1-MO][f2-MO][e][a];
	mol2 = moult2[t][r][f1-MO][f2-MO][e][a];
	alive += s;
	if (mol1 == MOULT1) moult_1 += s;
	else if (f1 < 0) moult_1 += s;
	if (mol2 == MOULT1) moult_1 += s;
	else if (f2 < 0) moult_2 += s;
	if (a != 0) breed += s;
	if (act == START) start += s;
	else if (act == ABANDON) aban += s;
	if (a == 0){
	  res[0] += absr*s;
	  fe1[0] += absf1*s;
	  fe2[0] += absf2*s;
	  uv[0] += uvalue[t][r][f1-MO][f2-MO][e][a]*s;
	  ex[0] += e*s;
	}
	else {
	  res[1] += absr*s;
	  fe1[1] += absf1*s;
	  fe2[1] += absf2*s;
	  uv[1] += uvalue[t][r][f1-MO][f2-MO][e][a]*s;
	  ex[1] += e*s;
	}
      }
      }
      }
      }
      fprintf(OutPut1, "%d ", t);
      fprintf(OutPut2, "%d ", t);
      Ntot += alive;
      if (alive > 0.0){
	Nstart += start/alive;
	Naban += aban/alive;
	Nm1 += moult_1/alive;
	Nm2 += moult_2/alive;
      }
      fprintf(OutPut1, "%f %f %f %f %f %f ", 
	      alive, breed, start, aban, moult_1,moult_2);
      fprintf(OutPut2, "%f %f ", alive, breed);
      for(i = 0; i < 2; i++){
	fprintf(OutPut2, "%f %f %f %f %f ", res[i], fe1[i], fe2[i], 
		uv[i], ex[i]);
      }
      fprintf(OutPut1, "\n");
      fprintf(OutPut2, "\n");
    }
    stepforwards(MT-1, COPY, INIT, 1.0);
    fprintf(OutPut1, "\n\n");
    fprintf(OutPut2, "\n\n");
    fprintf(OutPut3, "%f %f %f %f %f\n",
	   Ntot, Nstart, Naban, Nm1, Nm2);
  }
  fclose(OutPut1);
  fclose(OutPut2);
}
/************************* end of Cohort ******************/


/************************* OutPut *************************/
void OutPut(char *filename)
{
  char filen[30];
  int t, r, f1, f2, e, a, act, mol1, mol2;
  double s, total;
  double alive[ME+1], breed[ME+1], north[ME+1];
  double south[ME+1], moult_1[ME+1], moult_2[ME+1];
  double res[ME+1], fe1[ME+1], fe2[ME+1];
  FILE *OutPut1, *OutPut2;
  

  strcpy(filen,""); strcat(filen,filename); strcat(filen,".pop");
  if ((OutPut1 = fopen(filen,"wt")) == NULL) {
    printf("ERROR: Can't open the output file: %s\n", filen);
    exit (1);
  }
  strcpy(filen,""); strcat(filen,filename); strcat(filen,".sta");
  if ((OutPut2 = fopen(filen,"wt")) == NULL) {
    printf("ERROR: Can't open the output file: %s\n", filen);
    exit (1);
  }

  for(t = 0; t < MT; t++){
    total = 0.0;
    for(e = 0; e <= ME; e++){
      alive[e] = breed[e] = north[e] = south[e] = 0.0;
      moult_1[e] = moult_2[e] = 0.0;
      res[e] = fe1[e] = fe2[e] = 0.0;
    }
    for(a = 0; a <= MA; a++)
    for(r = 0; r <= MR; r++)
    for(f1 = MO; f1 <= MF; f1++)
    for(f2 = MO; f2 <= MF; f2++)
    for(e = 0; e <= ME; e++){
      total += (s = states[t][r][f1-MO][f2-MO][e][a]);
      act = action[t][r][f1-MO][f2-MO][e][a];
      mol1 = moult1[t][r][f1-MO][f2-MO][e][a];
      mol2 = moult2[t][r][f1-MO][f2-MO][e][a];
      alive[e] += s;
      if (f1<0) moult_1[e] += s;
      else if (f1 < 0) moult_1[e] += s;
      if (mol2 == MOULT1) moult_2[e] += s;
      else if (f2 < 0) moult_2[e] += s;
      if (act == KEEP || act == START) breed[e] += s;
      res[e] += r*s;
      if (f1 > -1) fe1[e] += f1*s;
      if (f2 > -1) fe2[e] += f2*s;
    }
    if (total > 0.0) {
      fprintf(OutPut1, "%d %f ", t, total);
      fprintf(OutPut2, "%d %f ", t, total);
      for(e = 0; e <= ME; e++){
	fprintf(OutPut1, "%f %f %f %f ", alive[e], breed[e], 
		moult_1[e], moult_2[e]);
	fprintf(OutPut2, "%f %f %f %f ", alive[e], res[e], 
		fe1[e], fe2[e]);
      }
      fprintf(OutPut1, "\n");
      fprintf(OutPut2, "\n");

    }
  }
  fclose(OutPut1);
  fclose(OutPut2);
}
/************************* end of OutPut ******************/

/************************* PolicyOut *************************/
void PolicyOut(char *filename)
{
  int t, r, f1, f2, e, a;
  char filen[30];
  FILE *OutPut1;
  FILE *OutPut2;
  FILE *OutPut3;
  double rf, pr;
  int rlo, rhi;
  double N_for[ME+1], n_fora;

  strcpy(filen,""); strcat(filen,filename); strcat(filen,".pol");
  if ((OutPut1 = fopen(filen,"wt")) == NULL) {
    printf("ERROR: Can't open the output file: %s\n", filen);
    exit (1);
  }
  strcpy(filen,""); strcat(filen,filename); strcat(filen,".mol");
  if ((OutPut2 = fopen(filen,"wt")) == NULL) {
    printf("ERROR: Can't open the output file: %s\n", filen);
    exit (1);
  }
  strcpy(filen,""); strcat(filen,filename); strcat(filen,".yrv");
  if ((OutPut3 = fopen(filen,"wt")) == NULL) {
    printf("ERROR: Can't open the output file: %s\n", filen);
    exit (1);
  }
  fprintf(OutPut1,"%d %d %d %d %d\n",MT, MR+1, MF-MO+1,MF-MO+1,MA+1);
  fprintf(OutPut2,"%d %d %d %d %d %d\n",MT,MR+1,MF-MO+1,MF-MO+1,MA+1, NFEATHER);
  
  rf = 0.5*TOPR;
  rlo = (int) rf;
  rhi = rlo + 1;
  pr = rf - ((double) rlo);

  for(t = 0; t < MT; t++){
    rf = repval[t][rlo][MF-MO][MF-MO][0][0]*(1.0-pr)+
      repval[t][rhi][MF-MO][MF-MO][0][0]*pr;
    fprintf(OutPut3,"%f\t%f\t",rf,FOOD[t]);
  
    n_fora = 0.0;
    for(e = 0; e <= ME; e++){
      N_for[e] = 0.0;
      for(r = 0; r <= MR; r++)
      for(f1 = MO; f1 <= MF; f1++)
      for(f2 = MO; f2 <= MF; f2++)
      for(a = 0; a <= MA; a++){
	N_for[e] += states[t][r][f1-MO][f2-MO][e][a];
      }
      n_fora += N_for[e]*theta[e];
    }
    fprintf(OutPut3,"%f\t", n_fora);
    for(e = 0; e <= ME; e++) {
      G[t][e] = (FOOD[t]/(n_fora+P.N0)) * theta[e];
      fprintf(OutPut3,"%f\t%f\t%f\t",G[t][e],theta[e],N_for[e]);
    }
    fprintf(OutPut3,"\n");
    for(r = 0; r <= MR; r++)
    for(f1 = MO; f1 <= MF; f1++)
    for(f2 = MO; f2 <= MF; f2++)
    for(a = 0; a <= MA; a++){
    fprintf(OutPut1,"%d ", action[t][r][f1-MO][f2-MO][ME][a]);
    fprintf(OutPut2,"%d %d ", moult1[t][r][f1-MO][f2-MO][ME][a],
	    moult2[t][r][f1-MO][f2-MO][ME][a]);
  }
}
  fclose(OutPut1);
  fclose(OutPut2);
  fclose(OutPut3);
}
/************************* PolicyOut ******************/


/************************* BinWrite *************************/
void BinWrite(char *filename)
{
  char filen[30];
  FILE *BinOut;
  int n_y, n_a;

  strcpy(filen,""); strcat(filen,filename); strcat(filen,".bin");
  if ((BinOut = fopen(filen,"wb")) == NULL) {
    printf("ERROR: Can't open the output file: %s\n", filen);
    exit (1);
  }
  n_y = (MT+1)*(MR+1)*(MF-MO+1)*(MF-MO+1)*(ME+1)*(MA+1);
  n_a = MT*(MR+1)*(MF-MO+1)*(MF-MO+1)*(ME+1)*(MA+1);
  if (fwrite(&n_y,sizeof(int),1,BinOut) != 1) hiba("BinWrite 1");
  if (fwrite(repval, sizeof(double),n_y,BinOut) != n_y) hiba("BinWrite 2a");
  if (fwrite(youngs, sizeof(double),n_y,BinOut) != n_y) hiba("BinWrite 2b");
  if (fwrite(&n_a,sizeof(int),1,BinOut) != 1) hiba("BinWrite 3");
  if (fwrite(uvalue, sizeof(double),n_a,BinOut) != n_a) hiba("BinWrite 4");
  if (fwrite(action, sizeof(int),n_a,BinOut) != n_a) hiba("BinWrite 5");
  if (fwrite(moult1, sizeof(int),n_a,BinOut) != n_a) hiba("BinWrite 6");
  if (fwrite(moult2, sizeof(int),n_a,BinOut) != n_a) hiba("BinWrite 7");
  fclose(BinOut);
}
/************************* BinWrite *************************/

/************************* hiba *************************/
void hiba(char *leiras)
{
  fprintf(stderr, "ERROR: %s\n", leiras);
  exit(1);
}
/************************* hiba *************************/

/************************* Treatment *************************/
void Treatment(int t, t_ind *Bird)
{
  /* dump function for simulated treatments, 
     working version in `moult_treat.c' */
}
/************************* Treatment *************************/

/************************* Treatment *************************/
void Food_Treatment(int ido, int t)
{
  /* dump function for simulated treatments, 
     working version in `moult_treat.c' */
}
/************************* Treatment *************************/

/************************* dump_repval *************************/
void dump_repval(void)
{
  
  char *filen="repval_dump.txt";
  FILE *OutPut;
  int t, r, f1, f2, a, e;
  int tm, rm, f1m, f2m, am, em;
  int tm0, rm0, f1m0, f2m0, am0, em0;
  int tmT, rmT, f1mT, f2mT, amT, emT;
  double rvm, rvm0, rvmT;

  if ((OutPut = fopen(filen,"wt")) == NULL) {
    printf("ERROR: Can't open the dump file: %s\n", filen);
    exit (1);
  }
  rvm = rvm0 = rvmT = 0.0;
  for(t = 0; t <= MT; t++)
  for(f1 = MO; f1 <= MF; f1++)
  for(f2 = MO; f2 <= MF; f2++)
  for(a = 0; a <= MA; a++)
  for(r = 0; r <= MR; r++)
  for(e = 0; e <= ME; e++){
    if (t==0) {
      if (rvm0<repval[t][r][f1-MO][f2-MO][e][a]) {
	rvm0 = repval[t][r][f1-MO][f2-MO][e][a];
	tm0 = t;
	rm0 = r;
	f1m0 = f1;
	f2m0 = f2;
	am0 = a;
	em0 = e;
      }
    }
    if (t==MT) {
      if (rvmT<repval[t][r][f1-MO][f2-MO][e][a]) {
	rvmT = repval[t][r][f1-MO][f2-MO][e][a];
	tmT = t;
	rmT = r;
	f1mT = f1;
	f2mT = f2;
	amT = a;
	emT = e;
      }
    }
    if (rvm<repval[t][r][f1-MO][f2-MO][e][a]) {
      rvm = repval[t][r][f1-MO][f2-MO][e][a];
      tm = t;
      rm = r;
      f1m = f1;
      f2m = f2;
      am = a;
      em = e;
      }
    if (t==MT) fprintf(OutPut, "%d %d %d %d %d %d %.15f %s %s %s %s\n", 
		       t,r,f1,f2,e,a,repval[t][r][f1-MO][f2-MO][e][a],
		       "NA","NA","NA","NA");
    else fprintf(OutPut, "%d %d %d %d %d %d %.15f %d %d %d %f\n", 
		 t,r,f1,f2,e,a,repval[t][r][f1-MO][f2-MO][e][a],
		 action[t][r][f1-MO][f2-MO][e][a], 
		 moult1[t][r][f1-MO][f2-MO][e][a], 
		 moult2[t][r][f1-MO][f2-MO][e][a], 
		 uvalue[t][r][f1-MO][f2-MO][e][a]);

  }  
  fprintf(OutPut,"\n");
  fclose(OutPut);
  printf("t==0; t:\t%d, r:\t%d, f1:\t%d, f2:\t%d, e:\t%d, a:\t%d, rv:\t%f\n",
	 tm0, rm0, f1m0, f2m0, em0, am0, rvm0);
  printf("t!=0; t:\t%d, r:\t%d, f1:\t%d, f2:\t%d, e:\t%d, a:\t%d, rv:\t%f\n",
	 tm, rm, f1m, f2m, em, am, rvm);
  printf("t==T; t:\t%d, r:\t%d, f1:\t%d, f2:\t%d, e:\t%d, a:\t%d, rv:\t%f\n",
	 tmT, rmT, f1mT, f2mT, emT, amT, rvmT);
}
/************************* dump_repval *************************/

/************************* dump_states *************************/
void dump_states(void)
{
  
  char *filen="states_dump.txt";
  FILE *OutPut;
  int t, r, f1, f2, a, e;
  int tm, rm, f1m, f2m, am, em;
  int tm0, rm0, f1m0, f2m0, am0, em0;
  int tmT, rmT, f1mT, f2mT, amT, emT;
  double rvm, rvm0, rvmT;

  if ((OutPut = fopen(filen,"wt")) == NULL) {
    printf("ERROR: Can't open the dump file: %s\n", filen);
    exit (1);
  }
  rvm = rvm0 = rvmT = 0.0;
  for(t = 0; t <= MT; t++)
  for(f1 = MO; f1 <= MF; f1++)
  for(f2 = MO; f2 <= MF; f2++)
  for(a = 0; a <= MA; a++)
  for(r = 0; r <= MR; r++)
  for(e = 0; e <= ME; e++){
    if (t==0) {
      if (rvm0<states[t][r][f1-MO][f2-MO][e][a]) {
	rvm0 = states[t][r][f1-MO][f2-MO][e][a];
	tm0 = t;
	rm0 = r;
	f1m0 = f1;
	f2m0 = f2;
	am0 = a;
	em0 = e;
      }
    }
    if (t==MT) {
      if (rvmT<states[t][r][f1-MO][f2-MO][e][a]) {
	rvmT = states[t][r][f1-MO][f2-MO][e][a];
	tmT = t;
	rmT = r;
	f1mT = f1;
	f2mT = f2;
	amT = a;
	emT = e;
      }
    }
    if (rvm<states[t][r][f1-MO][f2-MO][e][a]) {
      rvm = states[t][r][f1-MO][f2-MO][e][a];
      tm = t;
      rm = r;
      f1m = f1;
      f2m = f2;
      am = a;
      em = e;
      }
    if (t==MT) fprintf(OutPut, "%d %d %d %d %d %d %.15f %s %s %s %s\n", 
		       t,r,f1,f2,e,a,states[t][r][f1-MO][f2-MO][e][a],
		       "NA","NA","NA","NA");
    else fprintf(OutPut, "%d %d %d %d %d %d %.15f %d %d %d %f\n", 
		 t,r,f1,f2,e,a,states[t][r][f1-MO][f2-MO][e][a],
		 action[t][r][f1-MO][f2-MO][e][a], 
		 moult1[t][r][f1-MO][f2-MO][e][a], 
		 moult2[t][r][f1-MO][f2-MO][e][a], 
		 uvalue[t][r][f1-MO][f2-MO][e][a]);

  }  
  fprintf(OutPut,"\n");
  fclose(OutPut);
  printf("t==0; t:\t%d, r:\t%d, f1:\t%d, f2:\t%d, e:\t%d, a:\t%d, rv:\t%f\n",
	 tm0, rm0, f1m0, f2m0, em0, am0, rvm0);
  printf("t!=0; t:\t%d, r:\t%d, f1:\t%d, f2:\t%d, e:\t%d, a:\t%d, rv:\t%f\n",
	 tm, rm, f1m, f2m, em, am, rvm);
  printf("t==T; t:\t%d, r:\t%d, f1:\t%d, f2:\t%d, e:\t%d, a:\t%d, rv:\t%f\n",
	 tmT, rmT, f1mT, f2mT, emT, amT, rvmT);
}
/************************* dump_states *************************/

/************************* maxV *************************/
void maxV(int t, int r, int f1, int f2, int e, int a,
	  double *vmax, double *umax, int *acmax)
{
  double u1, u2, v1, v2;
  /*  double u3 u4, u5, u6, v3, v4, v5, v6;*/
  
  *umax = 0.0;
  *acmax = DELAY;
  *vmax = 0.0;
  if (r == 0) return;
  if (a == 0) {
    maxuv(t,r,f1,f2,e,a,DELAY,0.0,1.0,&u1,&v1);
    maxuv(t,r,f1,f2,e,a,START,0.0,1.0,&u2,&v2);
    unocare[r][f1-MO][f2-MO][e] = u1;
    vnocare[r][f1-MO][f2-MO][e] = v1;
    
    if (v2>v1){
      *umax = u2;
      *acmax = START;
      *vmax = v2;
    }
    else {
      *umax = u1;
      *acmax = DELAY;
      *vmax = v1;
    }
  }
  else if (a>0 && a<MA) {
    if (ucritcl[t][e] > 1.0){
      *umax = unocare[r][f1-MO][f2-MO][e];
      *acmax = ABORT;
      *vmax = vnocare[r][f1-MO][f2-MO][e];
    }
    else {
      maxuv(t,r,f1,f2,e,a,KEEP,ucritcl[t][e],1.0,&u1,&v1);
      if (v1>vnocare[r][f1-MO][f2-MO][e]) {
	*umax = u1;
	*acmax = KEEP;
	*vmax = v1;
      }
      else {
	*umax = unocare[r][f1-MO][f2-MO][e];
	*acmax = ABORT;
	*vmax = vnocare[r][f1-MO][f2-MO][e];
      }
    }
  }
  else if (a == MA) {
    maxuv(t,r,f1,f2,e,a,ABANDON,0.0,1.0,&u1,&v1);
    *umax = u1;
    *acmax = ABANDON;
    *vmax = v1;
  }
  return;
}
/************************* maxV *************************/

