
#include<string.h>
#include"moult.h"

#define NBIRDS 10000
#define NYEARS 5
#define FNBIRD ((float) NBIRDS)

extern double youngs[MT+1][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
extern double uvalue[MT][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
extern int    action[MT][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
extern int    moult1[MT][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
extern int    moult2[MT][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];

extern type_par P;

extern FILE *OutPut3;

t_ind Bird[NBIRDS];

extern float randf(float);
extern double newfq1(int, int, int, double, int, int);
extern double newres(double *, int, int, int, int, int, double, int);
extern void Treatment(int, t_ind *);
extern void Food_Treatment(int, int);

void ynormalise(double *, int );
void Ind(char *, int);
void MoultStr(int, char *);
void ActionStr(int, char *);
void WriteOut(int, int, int, char *, char *, char *, char *, t_ind, double, FILE *);
int smooth_r(double);

/************************* Ind *************************/
void Ind(char *filename, int tt)
{
  char filen[30], S1[30], S2[30], S3[30];
  int t, r, f1, f2, e, a, act, ido;
  double lambda, u, pmort;
  int year, i, ef1, ef2, ml1, ml2, al;
  FILE *OutPut1;
  double rf, ff1, ff2, pf1, pf2;
  int f1lo, f1hi, f2lo, f2hi, alive;
  int WRITE;

  strcpy(filen,""); strcat(filen,filename); strcat(filen,".dat");
  if ((OutPut1 = fopen(filen,"wt")) == NULL) {
    printf("ERROR: Can't open the output file: %s\n", filen);
    exit (1);
  }
  fprintf(OutPut1, "year\tweek\tbird\tevent\tmoult1\tmoult2\taction\tres\tfea1\tml1\tfea2\tml2\tbrood\texper\tuval\n");

  ynormalise(&lambda, SUM);

  for(year = 0; year < NYEARS; year++){
    alive = 0;
    WRITE = year >= 0;
    for(t = 0; t < MT; t++){
      ido = year*52 + t;
      if (tt != -1) Food_Treatment(ido, t);
      for(i = 0; i < NBIRDS; i++){
	if (year == 0 && Bird[i].birth == t) Bird[i].birth=-1;
	if (Bird[i].alive && Bird[i].birth == -1){
	  if (t==0) alive++;
	  if (tt != -1) Treatment(ido, &(Bird[i]));
	  ef1 = Bird[i].ef1;
	  ml1 = Bird[i].ml1;
	  ef2 = Bird[i].ef2;
	  ml2 = Bird[i].ml2;
	  al = Bird[i].al;
	  r = Bird[i].r;
	  f1 = Bird[i].f1;
	  f2 = Bird[i].f2;
	  e = Bird[i].e;
	  a = Bird[i].a;
	  u = uvalue[t][r][f1-MO][f2-MO][e][a];
	  act = action[t][r][f1-MO][f2-MO][e][a];
	  if (moult1[t][r][f1-MO][f2-MO][e][a]==MOULT1) {
	    f1 = MO;
	    strcpy(S1, "START");
	  }
	  else MoultStr(f1, S1);
	  if (moult2[t][r][f1-MO][f2-MO][e][a]==MOULT1) {
	    f2 = MO;
	    strcpy(S2, "START");
	  }
	  else MoultStr(f2, S2);
	  ActionStr(act,S3);
	  rf = TOPR*newres(&pmort,t,r,f1,f2,e,u,act);
	  if (Bird[i].r == 0.0) {
	    if (WRITE) WriteOut(year,t,i,"STARV",S1,S2,S3,Bird[i],u,OutPut1);
	    Bird[i].alive = 0;
	  }
	  else if (randf(1.0) < pmort) {
	    if (WRITE) WriteOut(year,t,i,"PREDA",S1,S2,S3,Bird[i],u,OutPut1);
	    Bird[i].alive = 0;
	  }
	  else {
	    if (WRITE) WriteOut(year,t,i,"ALIVE",S1,S2,S3,Bird[i],u,OutPut1);
	    if (f1<0) ml1++;
	    else ml1=0;
	    if (f2<0) ml2++;
	    else ml2=0;
	    if (a > 0) al++;
	    else al=0;
	    ff1 = TOPF*newfq1(r,f1,f2,u,act,0);
	    ff2 = TOPF*newfq1(r,f2,f1,u,act,1);
	    if (ff1 <= -0.9) {
	      f1lo = (int) ff1;
	      f1hi = f1lo-1;
	      pf1 = ((double) f1lo)-ff1;
	    }
	    else if (ff1 <= 0.0) {
	      f1lo = 0;
	      f1hi = 0;
	      pf1 = 1.0;
	    }
	    else if (ff1 >= TOPF) {
	      f1lo = MF;
	      f1hi = MF;
	      pf1 = 0.0;
	    }
	    else {
	      f1lo = (int) ff1;
	      f1hi = f1lo + 1;
	      pf1 = ff1 - ((double) f1lo);
	    }
	    if (ff2 <= -0.9) {
	      f2lo = (int) ff2;
	      f2hi = f2lo-1;
	      pf2 = ((double) f2lo)-ff2;
	    }
	    else if (ff2 <= 0.0) {
	      f2lo = 0;
	      f2hi = 0;
	      pf2 = 1.0;
	    }
	    else if (ff2 >= TOPF) {
	      f2lo = MF;
	      f2hi = MF;
	      pf2 = 0.0;
	    }
	    else {
	      f2lo = (int) ff2;
	      f2hi = f2lo + 1;
	      pf2 = ff2 - ((double) f2lo);
	    }
	    Bird[i].r = smooth_r(rf);
	    if (randf(1.0)<pf1) Bird[i].f1 = f1hi;
	    else Bird[i].f1 = f1lo;
	    if (randf(1.0)<pf2) Bird[i].f2 = f2hi;
	    else Bird[i].f2 = f2lo;
	    if (randf(1.0)<P.pe) Bird[i].e = (e+1 < ME ? e+1 : ME);
	    Bird[i].a = 0;
	    if (act == START) Bird[i].a = 1;
	    if (act == KEEP) {
	      if (randf(1.0) < P.pa) Bird[i].a = (a+1 > MA ? MA : a+1);
	      else Bird[i].a = a;
	    }
	    Bird[i].ef1 = f1;
	    Bird[i].ml1 = ml1;
	    Bird[i].ef2 = f2;
	    Bird[i].ml2 = ml2;
	    Bird[i].al = al;
	  }
	}
      }
    }
    fprintf(OutPut3,"%d, %d\n", year, alive);
  }
  fclose(OutPut1);
}
/************************* end of Ind ******************/

#define SF 4

/************************* get_r *************************/
/* mimic the smoothing of r */

int smooth_r(double rf) 
{
  double p[SF], coin, pr, sum;
  int r[SF], i;

  if (rf >= TOPR) {
    r[1] = r[2] = r[3] = MR;
    r[0] = MR - 1;
    p[2] = p[3] = 0.0;
    p[1] = 1.0 - P.alpha;
    p[0] = P.alpha;
  }
  else {
    r[1] = (int) rf;
    r[2] = r[1] + 1;
    r[3] = ((r[3]=r[2]+1) > MR ? MR : r[3]);
    r[0] = ((r[0]=r[1]-1) > 0 ? r[0] : 0);
    pr = rf - ((double) r[1]);
    p[0] = P.alpha*(1.0-pr);
    p[1] = 1.0 - pr*P.n_alpha - P.alpha2;
    p[2] = pr*P.n_alpha + P.alpha;
    p[3] = pr*P.alpha;
  }

  coin = randf(1.0);
  sum = 0.0;
  for(i = 0; i < SF; i++) {
    sum += p[i];
    if (coin <= sum) return r[i];
  }
  return r[3];
  /*  return r;
  p = randf(1.0);
  if (p < P.alpha) return ((--r) < 0 ? 0 : r);
  else if (p > (1.0-P.alpha)) return ((++r) > MR ? MR : r);
  else return r;*/
}
/************************* end of get_r ******************/


/************************* ynormalise *************************/
void ynormalise(double *lambda, int how)
{
  int t, r, f1, f2, e, a;
  double tmp, L1, L2;
  int i;

  *lambda = 0.0;
  for(t = 0; t <= MT; t++)
  for(r = 1; r <= MR; r++)
  for(f1 = MO; f1 <= MF; f1++)
  for(f2 = MO; f2 <= MF; f2++)
  for(e = 0; e <= ME; e++)
  for(a = 0; a <= MA; a++){
    tmp = youngs[t][r][f1-MO][f2-MO][e][a];
    if ((how == MAX) && (*lambda < tmp)) *lambda = tmp;
    if (how == SUM) *lambda += tmp;
  }
  L1 = L2 = 0.0;
  i = 0;
  if (*lambda > 0.0){
    for(t = 0; t <= MT; t++)
    for(r = 1; r <= MR; r++)
    for(f1 = MO; f1 <= MF; f1++)
    for(f2 = MO; f2 <= MF; f2++)
    for(e = 0; e <= ME; e++)
    for(a = 0; a <= MA; a++){
      L2 += (youngs[t][r][f1-MO][f2-MO][e][a] /= *lambda);
      while(((((float) i)/FNBIRD) < L2) && (i<NBIRDS)){
	Bird[i].alive = 1;
	Bird[i].r = r;
	Bird[i].f1 = f1;
	Bird[i].f2 = f2;
	Bird[i].e = e;
	Bird[i].a = a;
	Bird[i].ef1 = 0;
	Bird[i].ml1 = 0;
	Bird[i].ef2 = 0;
	Bird[i].ml2 = 0;
	Bird[i].al = 0;
	Bird[i].birth= t;
	i++;
      }
    }
  }
  else {
    fprintf(OutPut3,"ERROR: Something is wrong with lambda in ynormalise\n");
    exit(1);
  }
}
/************************* end of ynormalise ******************/


/************************* MoultStr *************************/
void MoultStr(int f, char *S)
{

  if(f<-1) strcpy(S, "MOULT");
  else if(f==-1) strcpy(S, "END");
  else strcpy(S, "NO");
  return;
}
/************************* end of MoultStr ******************/

/************************* MoultStr *************************/
void ActionStr(int act, char *S)
{

  switch (act) {
  case DELAY: strcpy(S, "DELAY"); break;
  case START: strcpy(S, "START"); break;
  case KEEP: strcpy(S, "KEEP"); break;
  case ABORT: strcpy(S, "ABORT"); break;
  case ABANDON: strcpy(S, "ABAND"); break;
  default: sprintf(S,"%d", act); break;
  }
}
/************************* end of MoultStr ******************/


/************************* WriteOut *************************/
void WriteOut(int year, int t, int i, char *S0, char *S1, char *S2, char *S3,
	      t_ind Bird, double u, FILE *OutPut)
{
  fprintf(OutPut, "%d\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n",
	  year, t, i, S0, S1, S2, S3, Bird.r, Bird.f1,
	  Bird.ml1, Bird.f2, Bird.ml2,Bird.al,Bird.e,u);

}
  /************************* WriteOut *************************/

