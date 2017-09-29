/*
c------------------------------------------------------------------------------
c ***** BACKWARDS AND FORWARDS VALUE PROPAGATION - MOULT VERSION **************
c------------------------------------------------------------------------------
c    AUTHOR:  BOB WELHAM, UNIVERSITY OF BRISTOL			OCTOBER 1999
c    STRICTLY CONFIDENTIAL - ALL RIGHTS RESERVED
c
c------------------------------------------------------------------------------

rewritten in C by me
*/

#include"moult.h"

extern double repval[MT+1][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
extern double states[MT+1][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
extern double uvalue[MT][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
extern int action[MT][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
extern int moult1[MT][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
extern int moult2[MT][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
extern double theta[ME+1];
extern double FOOD[MT+1];
extern double G[MT+1][ME+1];

extern type_par P;

extern void initialise(int, double, int);
extern double newres(double *, int, int, int, int, int, double, int);
extern double newfq1(int, int, int, double, int, int);
extern double gridread(int, double, double, double, int, int);
extern void gridwrite(int, double, double, double, int, int, double, int);

double backvalue(int, int, int, int, int, int, double, int);
void stepforwards(int, int, int, double);

/************************* stepforwards ***********************/
void stepforwards(int t, int how, int code, double lambda)
{
  double newr, newf1, newf2, pl, u, w, pmort;
  int r, f1, f2, e, a, actn, t1, e1, a1, ai;
  int m1, m2;
  double N_for[ME+1], n_fora;
  
  n_fora = 0.0;
  for(e = 0; e <= ME; e++){
    N_for[e] = 0.0;
    for(r = 0; r <= MR; r++)
    for(f1 = MO; f1 <= MF; f1++)
    for(f2 = MO; f2 <= MF; f2++)
    for(a = 0; a <= MA; a++){
      N_for[e] += states[t][r][f1-MO][f2-MO][e][a];
    }
    n_fora += N_for[e] * theta[e];
  }
  for(e = 0; e <= ME; e++) {
    G[t][e] = (FOOD[t]/(lambda*n_fora+P.N0)) * theta[e];
  }
  t1 = (t+1) % MT;
  if (code == INIT) initialise(t1, 0.0, FORWARD);
  for(a = MA; a >= 0; a--)
  for(f1 = MF; f1 >= MO; f1--)
  for(f2 = MF; f2 >= MO; f2--)
  for(r = 0; r <= MR; r++)
  for(e = 0; e <= ME; e++){
    w = states[t][r][f1-MO][f2-MO][e][a];
    if (w > 0.0) {
      m1 = (moult1[t][r][f1-MO][f2-MO][e][a]==MOULT1);
      m2 = (moult2[t][r][f1-MO][f2-MO][e][a]==MOULT1);
      if (m1 && m2) {
	states[t][r][MO-MO][MO-MO][e][a] += w;
	states[t][r][f1-MO][f2-MO][e][a] = 0.0;
      }
      else if (m1) {
	states[t][r][MO-MO][f2-MO][e][a] += w;
	states[t][r][f1-MO][f2-MO][e][a] = 0.0;
      }
      else if (m2) {
	states[t][r][f1-MO][MO-MO][e][a] += w;
	states[t][r][f1-MO][f2-MO][e][a] = 0.0;
      }
      else {
	u = uvalue[t][r][f1-MO][f2-MO][e][a];
	actn = action[t][r][f1-MO][f2-MO][e][a];
	newr = newres(&pmort,t,r,f1,f2,e,u,actn)*TOPR;
	newf1 = newfq1(r,f1,f2,u,actn,0)*TOPF;
	newf2 = newfq1(r,f2,f1,u,actn,1)*TOPF;
	pl = 1.0 - pmort;
	e1 = (e+1 < ME ? e+1 : ME);
	a1 = 0;
	ai = a;
	if (actn == START) {
	    a1 = 1;
	    ai = 1;
	  }
	else if (actn == KEEP) a1 = ((a+1)<MA ? (a+1) : MA);
	else if (actn == ABANDON || actn == ABORT) a1 = ai = 0;
	gridwrite(t1,newr,newf1,newf2,e, a1,P.pa*(1.0-P.pe)*pl*w,NOYOUNGS);
	gridwrite(t1,newr,newf1,newf2,e1,a1,P.pa*P.pe*      pl*w,NOYOUNGS);
	gridwrite(t1,newr,newf1,newf2,e,ai,(1.0-P.pa)*(1.0-P.pe)*pl*w,NOYOUNGS);
	gridwrite(t1,newr,newf1,newf2,e1,ai,(1.0-P.pa)*P.pe*pl*w,NOYOUNGS);
	if (actn == ABANDON && how == ALL) {
	  gridwrite(t,0.5*TOPR,TOPF,TOPF,0,0,P.n*w, NOYOUNGS);
	  gridwrite(t,0.5*TOPR,TOPF,TOPF,0,0,P.n*w, YOUNGS);
	}
      }
    }
  }
}
/************************* end of stepforwards ******************/


/************************* backvalue *************************/
double backvalue(int t, int r, int f1, int f2, int e, int a, 
		 double u, int actn)
{
  double newr, newf1, newf2, pl, pmort;
  int t1, e1, a1;
  double bv;

  t1 = (t+1) % MT;
  newr = newres(&pmort,t,r,f1,f2,e,u,actn)*TOPR;
  if (newr <= 0.0) {
    return 0.0;
  }
  newf1 = newfq1(r,f1,f2,u,actn,0)*TOPF;
  newf2 = newfq1(r,f2,f1,u,actn,1)*TOPF;
  pl = 1.0-pmort;
  if (e < ME) {
    e1 = (ME > e+1 ? e+1 : ME);
    if (actn == DELAY || actn == ABANDON) {
      bv = (1.0-P.pe)*gridread(t1, newr, newf1, newf2, e,  0) +
	P.pe      *gridread(t1, newr, newf1, newf2, e1, 0); 
      bv *= pl;
      if (actn == ABANDON) 
	bv += P.n*gridread(t,0.5*TOPR, TOPF, TOPF, 0, 0);
    } 
    else if (actn == START) {
      bv = (1.0-P.pe)*gridread(t1, newr, newf1, newf2, e,  1) +
	P.pe      *gridread(t1, newr, newf1, newf2, e1, 1); 
      bv *= pl;
    }
    else if (actn == KEEP) {
      a1 = (MA > a+1 ? a+1 : MA);
      bv = (1.0-P.pe)*P.pa*      gridread(t1, newr, newf1, newf2, e,  a1) +
	P.pe      *P.pa*      gridread(t1, newr, newf1, newf2, e1, a1) +
	(1.0-P.pe)*(1.0-P.pa)*gridread(t1, newr, newf1, newf2, e,  a) +
	P.pe      *(1.0-P.pa)*gridread(t1, newr, newf1, newf2, e1, a);
      bv *= pl;
    }
  } else {
    if (actn == DELAY || actn == ABANDON) {
      bv = pl*gridread(t1, newr, newf1, newf2, e,  0);
      if (actn == ABANDON) 
	bv += P.n*gridread(t,0.5*TOPR, TOPF, TOPF, 0, 0);
    } 
    else if (actn == START) {
      bv = pl*gridread(t1, newr, newf1, newf2, e,  1);
    }
    else if (actn == KEEP) {
      a1 = (MA > a+1 ? a+1 : MA);
      bv = P.pa*      gridread(t1, newr, newf1, newf2, e,  a1) +
	(1.0-P.pa)*gridread(t1, newr, newf1, newf2, e,  a);
      bv *= pl;
    }
  }
  return bv;
}
/************************* end of backvalue ******************/
