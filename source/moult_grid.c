/*
c-----------------------------------------------------------------------------
c ********************** GRID ROUTINES - MOULT VERSION ***********************
c-----------------------------------------------------------------------------
c  AUTHOR:  BOB WELHAM, UNIVERSITY OF BRISTOL			OCTOBER 1999
c  STRICTLY CONFIDENTIAL - ALL RIGHTS RESERVED
c
c-----------------------------------------------------------------------------

rewritten in C by me
*/

#include"moult.h"

extern double repval[MT+1][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
extern double states[MT+1][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
extern double youngs[MT+1][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];

extern FILE *OutPut3;

extern type_par P;

void initialise(int, double, int);
void savevalues_back(void);
void savevalues_for(double);
void setup_n_for(double);
double gridread(int, double, double, double, int, int);
double gridget(int, int, int, int, int, int, int, int);
void normalise(double *, double *, int, int);
void gridwrite(int, double, double, double, int, int, double, int);
void gridinc(int, int, int, int, int, int, double, int);

/************************* gridinc *************************/
void gridinc(int t, int r, int f1, int f2, int e, int a, 
	     double inc, int code)
{

  if (r == 0) return;
  if (code == YOUNGS) {
    youngs[t][r][f1-MO][f2-MO][e][a] += inc;
  }
  else {
    states[t][r][f1-MO][f2-MO][e][a] += inc;
  }
}
/************************* end of gridinc ******************/


/************************* gridwrite *************************/
void gridwrite(int t, double r, double f1, double f2, int e, 
	       int a, double inc, int code)
{
  double pr, pf1, pf2, qr, qf1, qf2;
  int f1lo, f2lo, f1hi, f2hi;
  double p1, p2, p3, p4;
  int r1, r2, r3, r4;

  if (inc == 0.0 || r <= 0.0) return;
  else if (r >= TOPR) {
    r2 = r3 = r4 = MR;
    r1 = MR - 1;
    p3 = p4 = 0.0;
    p2 = 1.0 - P.alpha;
    p1 = P.alpha;
  }
  else {
    r2 = (int) r;
    r3 = r2 + 1;
    r4 = ((r4=r3+1) > MR ? MR : r4);
    r1 = ((r1=r2-1) > 0 ? r1 : 0);
    pr = r - ((double) r2);
    p1 = P.alpha*(1.0-pr);
    p2 = 1.0 - pr*P.n_alpha - P.alpha2;
    p3 = pr*P.n_alpha + P.alpha;
    p4 = pr*P.alpha;
  }
  if (f1 <= -0.9) {
    f1lo = (int) f1;
    f1hi = f1lo-1;
    pf1 = ((double) f1lo)-f1;
  }
  else if (f1 <= 0.0) {
    f1lo = 0;
    f1hi = 0;
    pf1 = 1.0;
  }
  else if (f1 >= TOPF) {
    f1lo = MF;
    f1hi = MF;
    pf1 = 0.0;
  }
  else {
    f1lo = (int) f1;
    f1hi = f1lo + 1;
    pf1 = f1 - ((double) f1lo);
  }
  if (f2 <= -0.9) {
    f2lo = (int) f2;
    f2hi = f2lo-1;
    pf2 = ((double) f2lo)-f2;
  }
  else if (f2 <= 0.0) {
    f2lo = 0;
    f2hi = 0;
    pf2 = 1.0;
  }
  else if (f2 >= TOPF) {
    f2lo = MF;
    f2hi = MF;
    pf2 = 0.0;
  }
  else {
    f2lo = (int) f2;
    f2hi = f2lo + 1;
    pf2 = f2 - ((double) f2lo);
  }
  qr = 1.0 - pr;
  qf1 = 1.0 - pf1;
  qf2 = 1.0 - pf2;

  if(r2!=0) {
    gridinc(t,r1,f1lo,f2lo,e,a,p1*qf1*qf2*inc,code);
    gridinc(t,r1,f1hi,f2lo,e,a,p1*pf1*qf2*inc,code);
    gridinc(t,r1,f1lo,f2hi,e,a,p1*qf1*pf2*inc,code);
    gridinc(t,r1,f1hi,f2hi,e,a,p1*pf1*pf2*inc,code);
    gridinc(t,r2,f1lo,f2lo,e,a,p2*qf1*qf2*inc,code);
    gridinc(t,r2,f1hi,f2lo,e,a,p2*pf1*qf2*inc,code);
    gridinc(t,r2,f1lo,f2hi,e,a,p2*qf1*pf2*inc,code);
    gridinc(t,r2,f1hi,f2hi,e,a,p2*pf1*pf2*inc,code);
  }

  gridinc(t,r3,f1lo,f2lo,e,a,p3*qf1*qf2*inc,code);
  gridinc(t,r3,f1hi,f2lo,e,a,p3*pf1*qf2*inc,code);
  gridinc(t,r3,f1lo,f2hi,e,a,p3*qf1*pf2*inc,code);
  gridinc(t,r3,f1hi,f2hi,e,a,p3*pf1*pf2*inc,code);
  gridinc(t,r4,f1lo,f2lo,e,a,p4*qf1*qf2*inc,code);
  gridinc(t,r4,f1hi,f2lo,e,a,p4*pf1*qf2*inc,code);
  gridinc(t,r4,f1lo,f2hi,e,a,p4*qf1*pf2*inc,code);
  gridinc(t,r4,f1hi,f2hi,e,a,p4*pf1*pf2*inc,code);
}
/************************* end of gridwrite ******************/


/************************* normalise *************************/
void normalise(double *lambda, double *diff, int how, int act)
{
  int r, f1, f2, e, a;
  double tmp, tmp2;

  *lambda = 0.0;
  *diff = 0.0;
  for(r = 0; r <= MR; r++)
  for(f1 = MO; f1 <= MF; f1++)
  for(f2 = MO; f2 <= MF; f2++)
  for(e = 0; e <= ME; e++)
  for(a = 0; a <= MA; a++){
    if (how == MAX) {
      tmp = repval[0][r][f1-MO][f2-MO][e][a];
      if (*lambda < tmp) *lambda = tmp;
    }
    if (how == SUM) {
      *lambda += states[0][r][f1-MO][f2-MO][e][a];
    }
  }
  if (*lambda > 0.0){
    tmp = 0.0;
    for(r = 0; r <= MR; r++)
    for(f1 = MO; f1 <= MF; f1++)
    for(f2 = MO; f2 <= MF; f2++)
    for(e = 0; e <= ME; e++)
    for(a = 0; a <= MA; a++){
      if (how == MAX) {
	if (act==NORM) repval[0][r][f1-MO][f2-MO][e][a] /= *lambda;
	*diff += fabs(repval[0][r][f1-MO][f2-MO][e][a] -
		      repval[MT][r][f1-MO][f2-MO][e][a]);
      } else {
	if (act==NORM) states[0][r][f1-MO][f2-MO][e][a] /= *lambda;
	*diff += (tmp2 = fabs(states[0][r][f1-MO][f2-MO][e][a] -
		      states[MT][r][f1-MO][f2-MO][e][a]));
	if (tmp < tmp2) tmp = tmp2;
      }
    }
    /*    if (how == SUM) *diff /= tmp;*/
  }
  else {
    fprintf(OutPut3,
	    "ERROR: Something is wrong with lambda in normalise: %f\n", 
	    *lambda);
    /*    exit(1);*/
  }
}
/************************* end of normalise ******************/


/************************* gridread *************************/
double gridread(int t, double r, double f1, double f2, int e, int a)
{
  double pr, pf1, pf2, qr, qf1, qf2;
  double gr = 0;
  double p1, p2, p3, p4;
  /*  int rlo, rhi, rlu, rld, rhu, rhd;*/
  int r1, r2, r3, r4;
  int f1lo, f2lo, f1hi, f2hi;

  if (r <= 0.0) return 0.0;
  else if (r >= TOPR) {
    r2 = r3 = r4 = MR;
    r1 = MR - 1;
    p3 = p4 = 0.0;
    p2 = 1.0 - P.alpha;
    p1 = P.alpha;
  }
  else {
    r2 = (int) r;
    r3 = r2 + 1;
    r4 = ((r4=r3+1) > MR ? MR : r4);
    r1 = ((r1=r2-1) > 0 ? r1 : 0);
    pr = r - ((double) r2);
    p1 = P.alpha*(1.0-pr);
    p2 = 1.0 - pr*P.n_alpha - P.alpha2;
    p3 = pr*P.n_alpha + P.alpha;
    p4 = pr*P.alpha;
  }
  if (f1 <= -0.9) {
    f1lo = (int) f1;
    f1hi = f1lo-1;
    pf1 = ((double) f1lo)-f1;
  }
  else if (f1 <= 0.0) {
    f1lo = 0;
    f1hi = 0;
    pf1 = 1.0;
  }
  else if (f1 >= TOPF) {
    f1lo = MF;
    f1hi = MF;
    pf1 = 0.0;
  }
  else {
   f1lo = (int) f1;
   f1hi = f1lo + 1;
   pf1 = f1 - ((double) f1lo);
  }
  if (f2 <= -0.9) {
    f2lo = (int) f2;
    f2hi = f2lo-1;
    pf2 = ((double) f2lo)-f2;
  }
  else if (f2 <= 0.0) {
    f2lo = 0;
    f2hi = 0;
    pf2 = 1.0;
  }
  else if (f2 >= TOPF) {
    f2lo = MF;
    f2hi = MF;
    pf2 = 0.0;
  }
  else {
   f2lo = (int) f2;
   f2hi = f2lo + 1;
   pf2 = f2 - ((double) f2lo);
  }
  qr = 1.0 - pr;
  qf1 = 1.0 - pf1;
  qf2 = 1.0 - pf2;
  f1lo -= MO;
  f1hi -= MO;
  f2lo -= MO;
  f2hi -= MO;
  

  if(r2!=0) {
    gr = qf1*qf2*(p1*repval[t][r1][f1lo][f2lo][e][a]  +
		  p2*repval[t][r2][f1lo][f2lo][e][a]) +
         pf1*qf2*(p1*repval[t][r1][f1hi][f2lo][e][a]  +
		  p2*repval[t][r2][f1hi][f2lo][e][a]) + 
         qf1*pf2*(p1*repval[t][r1][f1lo][f2hi][e][a]  + 
		  p2*repval[t][r2][f1lo][f2hi][e][a]) + 
         pf1*pf2*(p1*repval[t][r1][f1hi][f2hi][e][a]  + 
		  p2*repval[t][r2][f1hi][f2hi][e][a]);
  }

  gr += qf1*qf2*(p3*repval[t][r3][f1lo][f2lo][e][a]  +
		 p4*repval[t][r4][f1lo][f2lo][e][a]) +
        pf1*qf2*(p3*repval[t][r3][f1hi][f2lo][e][a]  +
		 p4*repval[t][r4][f1hi][f2lo][e][a]) + 
        qf1*pf2*(p3*repval[t][r3][f1lo][f2hi][e][a]  + 
		 p4*repval[t][r4][f1lo][f2hi][e][a]) + 
        pf1*pf2*(p3*repval[t][r3][f1hi][f2hi][e][a]  + 
		 p4*repval[t][r4][f1hi][f2hi][e][a]);
  return gr;
}
/************************* end of gridread ******************/

/************************* savevalues_back *************************/
void savevalues_back(void)
{
  int r, f1, f2, e, a;

  for(r = 0; r <= MR; r++)
  for(f1 = MO; f1 <= MF; f1++)
  for(f2 = MO; f2 <= MF; f2++)
  for(e = 0; e <= ME; e++)
  for(a = 0; a <= MA; a++){
    repval[MT][r][f1-MO][f2-MO][e][a] = 
      repval[0][r][f1-MO][f2-MO][e][a];
  }
}
/************************* end of savevalues_back ******************/


/************************* savevalues_for *************************/
void savevalues_for(double lambda)
{
  int r, f1, f2, e, a;

  for(r = 0; r <= MR; r++)
  for(f1 = MO; f1 <= MF; f1++)
  for(f2 = MO; f2 <= MF; f2++)
  for(e = 0; e <= ME; e++)
  for(a = 0; a <= MA; a++){
    states[MT][r][f1-MO][f2-MO][e][a] = 
      states[0][r][f1-MO][f2-MO][e][a];
  }
}
/************************* end of savevalues_for ******************/

/************************* setup_n_for *************************/
void setup_n_for(double lambda)
{
  int r, f1, f2, e, a;
  double sum_n;

  sum_n = 0.0;
  for(r = 0; r <= MR; r++)
  for(f1 = MO; f1 <= MF; f1++)
  for(f2 = MO; f2 <= MF; f2++)
  for(e = 0; e <= ME; e++)
  for(a = 0; a <= MA; a++){
    sum_n += states[0][r][f1-MO][f2-MO][e][a];
  }
  lambda /= sum_n;
  for(r = 0; r <= MR; r++)
  for(f1 = MO; f1 <= MF; f1++)
  for(f2 = MO; f2 <= MF; f2++)
  for(e = 0; e <= ME; e++)
  for(a = 0; a <= MA; a++){
    states[0][r][f1-MO][f2-MO][e][a] *= lambda;
    states[MT][r][f1-MO][f2-MO][e][a] = 
      states[0][r][f1-MO][f2-MO][e][a];
  }
}
/************************* setup_n_for ******************/



/************************* initialise *************************/
void initialise(int t, double initval, int what)
{
  int r, f1, f2, e, a;

  for(e = 0; e <= ME; e++)
  for(r = 0; r <= MR; r++)
  for(f1 = MO; f1 <= MF; f1++)
  for(f2 = MO; f2 <= MF; f2++)
  for(a = 0; a <= MA; a++){
    if (what == BACKWARD) {
      if (r == 0) repval[t][r][f1-MO][f2-MO][e][a] = 0.0;
      else repval[t][r][f1-MO][f2-MO][e][a] = initval;
    } else {
      states[t][r][f1-MO][f2-MO][e][a] = 0.0;
      /*      if (r == 0) states[t][r][f1-MO][f2-MO][e][a] = 0.0;
	      else states[t][r][f1-MO][f2-MO][e][a] = initval;*/
    }
  }
  states[t][MR][MF][MF][ME][0] = initval;
}
/************************* end of initialise ******************/


