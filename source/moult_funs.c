

#include"moult.h"


extern type_par P;

extern double states[MT+1][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];

extern FILE *OutPut3;

double G[MT+1][ME+1];
double G_o[MT+1][ME+1];
double FP[NFEATHER][MF-MO+1];
double FOOD[MT+1];


/* these two arrays below introduced because of speed optimization */
double pred_dfe[MF-MO+1][MF-MO+1];
double fora_dfe[MF-MO+1][MF-MO+1];
double theta[ME+1];


double newres(double *, int, int, int, int, int, double, int);
double newfq1(int, int, int, double, int, int);
void iniFP(void);
void inigain(double);
void ReadIn(char *);


/************************* newfq1 *************************/
/* f1 contains the feather state under consideration; f2 is the other
state; feather gives the 'name' of the feather under consideration */

double newfq1(int r, int f1, int f2, double u, int action, int feather)
{
  double nfq1;

  /*  if (MF == 0) return 0.0;*/
  if (f1 == -1) return 1.0;				/* finishing moult */
  if (f1 < -1) {
    nfq1 = (((double) f1)+P.F[feather].NU)/TOPF;	/* growing feathers */
    return nfq1;
  }
  if (P.F[feather].fixed) nfq1 = ((double) f1)/TOPF - P.F[feather].fora;
  else {
    nfq1 = ((double) f1)/TOPF - P.F[feather].fora*u;
    if (f2 < 0) nfq1 -= P.F[feather].moult;
  }
  nfq1 = (nfq1 > 1.0 ? 1.0 : (nfq1 < 0.0 ? 0.0: nfq1));
  return nfq1;
}
/************************* end of newfq1 ******************/

/************************* newres *************************/
double newres(double *Pred, int t, int r, int f1, int f2, int e, double u, 
	      int action)
{

  double nr;
  double c;
  unsigned i1, i2;
  int fl1, fl2;
  double absr, a2, u2;

  absr = ((double) r)/TOPR;
  a2 = absr*absr;
  i1 = i2 = 0;
  fl1 = f1-MO;
  fl2 = f2-MO;
  if (f1 < 0) i1 = 1;
  if (f2 < 0) i2 = 1;
  c = P.Cost.base + P.Cost.mass * a2;
  u2 = u*u;
  if ((*Pred = P.Pred.base+u2*(0.1*a2+1.0)*pred_dfe[fl1][fl2]) > 1.0) 
    *Pred = 1.0;
  c += u2*(0.1*a2+1.0)*fora_dfe[fl1][fl2]; 
  c += i1*P.Cost.moult1+i2*P.Cost.moult2+(i1||i2)*P.Cost.moint;
  if (action==START) c += P.Cost.start;
  else if (action==KEEP) c += P.Cost.keep;

  nr = absr + u*G[t][e] - c;
  nr = (nr > 1.0 ? 1.0 : (nr < 0.0 ? 0.0: nr));
  return nr;
}
/************************* end of newres ******************/

#define MOULT_TIME 1

/************************* iniFP *************************/
void iniFP(void)
{
  int i, f, f1, f2;
  double absf, dfe1, dfe2, dfe;

  for(f=0; f<NFEATHER; f++)
    for(i=MO; i <=MF; i++){
      if (i < 0) {
	absf = ((i-MO)*(TOPF/((double) abs(MO))))/TOPF;
      } else {
	absf = ((double) i)/TOPF;
      }
      FP[f][i-MO] = P.F[f].worn+pow(absf,P.F[f].flight)*(1.0-P.F[f].worn);
    }

  for(f1 = MO; f1 <= MF; f1++) {
    dfe1 = FP[0][f1-MO];
    for(f2 = MO; f2 <= MF; f2++){
      dfe2 = FP[1][f2-MO];
      dfe = (P.F[0].def_coeff*dfe1 + P.F[1].def_coeff*dfe2 + 
	     P.Cost.dfeint*dfe1*dfe2)/P.Cost.sum_ability;
      pred_dfe[f1-MO][f2-MO] = P.Pred.mass/dfe;
      fora_dfe[f1-MO][f2-MO] = P.Cost.propNF*P.Cost.fora + 
	(P.Cost.propF*P.Cost.def)/dfe;
    }
  }

}
/************************* end of iniFP ******************/

/************************* inigain *************************/
void inigain(double lambda)
{
  double n_fora;
  int t, r, f1, f2, a, e;
  double N_for[MT+1][ME+1];
  static int KEZD = 0;

  for(t = 0; t <= MT; t++)
    for(e = 0; e <= ME; e++) {
      N_for[t][e] = 0.0;
    }
  if (KEZD != 0) {
    for(t = 0; t <= MT; t++)
    for(e = 0; e <= ME; e++) 
    for(a = 0; a <= MA; a++)
    for(f2 = MO; f2 <= MF; f2++)
    for(f1 = MO; f1 <= MF; f1++)
    for(r = 0; r <= MR; r++) {
      N_for[t][e] += states[t][r][f1-MO][f2-MO][e][a];
    }    
  }
  for(t=0; t<=MT; t++){
    n_fora = 0.0;
    for(e=0; e<=ME; e++){
      if(KEZD == 0) n_fora = 1.0;
      else n_fora += N_for[t][e]*theta[e];
    }
    /*    fprintf(OutPut3, "AAA %d ", t);*/
    for(e=0; e<=ME; e++){
      G_o[t][e] = G[t][e] = (FOOD[t]/(lambda*n_fora+P.N0)) * theta[e];
      /*      fprintf(OutPut3, "%f %f ", G[t][e], N_for[t][e]);*/
    }
    /*    fprintf(OutPut3, "\n");*/
  }
  if (KEZD == 0) KEZD = 1;
}
/************************* end of inigain ******************/

/************************************* ReadIn ******************************/
void ReadIn(char *filename)
{
  FILE *InPut;
  char filen[30];
  char *ext = ".ini";
  int i, f;
  char junk[50];
  int errors = 0;
  int t, r, fe, a, e, fo;

  strcpy(filen,""); strcat(filen,filename); strcat(filen,ext);
  if ((InPut = fopen(filen,"rt")) == NULL) {
    printf("ERROR: Can't open the input file: %s\n", filen);
    exit (1);
  }

  for(f = 0; f < NFEATHER; f++){ 
    if (2 != fscanf(InPut, "%s%d", junk, &(P.F[f].fixed))) errors++;
    if (2 != fscanf(InPut, "%s%lf", junk, &(P.F[f].fora))) errors++;
    if (2 != fscanf(InPut, "%s%lf", junk, &(P.F[f].moult))) errors++;
    if (2 != fscanf(InPut, "%s%lf", junk, &(P.F[f].migra))) errors++;
    if (2 != fscanf(InPut, "%s%lf", junk, &(P.F[f].flight))) errors++;
    if (2 != fscanf(InPut, "%s%lf", junk, &(P.F[f].worn))) errors++;
    if (2 != fscanf(InPut, "%s%lf", junk, &(P.F[f].moulted))) errors++;
    if (2 != fscanf(InPut, "%s%lf", junk, &(P.F[f].def_coeff))) errors++;
    if (2 != fscanf(InPut, "%s%lf", junk, &(P.F[f].NU))) errors++;
  }
				/* metabolisms */
  if (2 != fscanf(InPut, "%s%lf", junk, &(P.Cost.fora))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P.Cost.propF))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P.Cost.def))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P.Cost.start))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P.Cost.keep))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P.Cost.moult1))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P.Cost.moult2))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P.Cost.moint))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P.Cost.dfeint))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P.Cost.migra))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P.Cost.mass))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P.Cost.base))) errors++;
				/* predation */
  if (2 != fscanf(InPut, "%s%lf", junk, &(P.Pred.base))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P.Pred.mass))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P.THETA))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P.MAXEPS))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P.FOOD))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P.n))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P.N0))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P.pe))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P.pa))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P.tol))) errors++;
  if (2 != fscanf(InPut, "%s%lf", junk, &(P.alpha))) errors++;
  if (2 != fscanf(InPut, "%s%d", junk, &r)) errors++;
  if (2 != fscanf(InPut, "%s%d", junk, &t)) errors++;
  if (2 != fscanf(InPut, "%s%d", junk, &fe)) errors++;
  if (2 != fscanf(InPut, "%s%d", junk, &fo)) errors++;
  if (2 != fscanf(InPut, "%s%d", junk, &a)) errors++;
  if (2 != fscanf(InPut, "%s%d", junk, &e)) errors++;
  if (errors > 0) {
    printf("ERROR: in reading input parameters from file %s at %d\n", filen, 
	   errors);
    exit(1);
  }
  if ((MR!=r) | (MF!=fe) | (MT!=t) | (ME!=e) | (MA!=a) | (MO!=fo)) {
    printf("ERROR: the compiled grid not the same as the assumed one!\n");
    printf("\tCompiled:\tMT=%d\tMR=%d\tMF=%d\tME=%d\tMA=%d\tMO=%d\n",
	   MT, MR, MF, ME, MA, MO);
    printf("\tAssumed:\tMT=%d\tMR=%d\tMF=%d\tME=%d\tMA=%d\tMO=%d\n",
	   t, r, fe, e, a, fo);
    exit(1);
  }
  P.Cost.propNF = 1.0 - P.Cost.propF;
  if ((P.Cost.sum_ability = P.F[0].def_coeff + P.F[1].def_coeff+P.Cost.dfeint) 
      == 0) {
    printf("ERROR: P.Cost.sum_ability is zero: exiting\n");
    exit(1);
  }
  P.n_alpha = 1.0-3.0*P.alpha;
  P.alpha2 = 2.0*P.alpha;
  for (i=0; i<NFEATHER; i++){
    P.F[i].moulted = P.F[i].worn - P.F[i].moulted;
    P.F[i].NU = 1.0-P.F[i].NU;
  }
  P.FOOD *= P.N0 + 1;
  P.MAXEPS *= P.N0 +1;
  for(e=0; e<=ME; e++) 
    theta[e] = pow(((double) P.THETA), ((double) ME-e));
  for(t=0; t<=MT; t++){
    FOOD[t] = P.FOOD+P.MAXEPS*sin(2.0*PI*((double) t)/((double) MT) - PI/2.0);
  }
  fclose(InPut);
}
/************************************* ReadIn ***********************/


