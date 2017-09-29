
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<fpu_control.h>
#include<time.h>

/*states*/
#define MR 12
#define MF 10
#define MT 52
#define ME 2
#define MA 7
#define MO -6

#define TOPR ((double) MR)
#define TOPF ((double) MF)

#define NFEATHER 2
#define UGRAIN 10
#define UTOL 1.0e-3

#define PI 3.1415926535898

#define MAX 900
#define SUM 901
#define ALL 902
#define BACKWARD 903
#define FORWARD 904
#define YOUNGS 905
#define NOYOUNGS 906
#define INIT 907
#define COPY 908
#define NORM 909
#define NONORM 910

/*actions*/
#define DELAY 1
#define START 2
#define KEEP 3
#define ABANDON 4
#define ABORT 5
/*#define FLYNORTH 6
  #define FLYSOUTH 7*/
#define NOMOULT 8
#define MOULT1 9

typedef struct {
  int fixed;			/* abrasion fixed? 1: yes */
  double fora;			/* effect of foraging intensity on feather abrasion */
  double moult;			/* effect of moulting the other feather on abrasion */
  double migra;			/* effect of migration */
  double flight;		/* effect of feather quaility on the flight ability */
  double worn;			/* the flight ability of birds with feathers in very poor condition */
  double moulted;		/* the flight ability of birds with moulting feathers */
  double def_coeff;		/* the effect of feather abrasion on flight efficiency */
  double NU;			/* stochasticity of moult length  */
} t_f;

typedef struct {
  double fora;			/* cost of foraging */
  double propF;			/* the proportion of flight used during foraging */
  double propNF;		/* 1-propF; just to save computation time */
  double sum_ability;		/* def_coeff1+def_coeff2+def; just to save computation time */
  double def;			/* the cost of decreased flight eff. during foraging */
  double start;			/* cost of starting a brood */
  double keep;			/* amount of food needed by the brood */
  double moult1;		/* cost of moulting feather 1 alone */
  double moult2;		/* cost of moulting feather 2 alone */
  double moint;			/* additional cost when moulting both feathers */
  double dfeint;		/* interaction coeff for decrease in flight efficiency */
  double migra;			/* migration cost given feathers are in top quality */
  double mass;			/* mass dependent cost */
  double base;			/* basic cost */
} t_c;

typedef struct {
  double base;
  double mass;
} t_p;

typedef struct {
  t_f F[NFEATHER];
  t_c Cost;
  t_p Pred;
  double THETA;
  double FOOD;
  double MAXEPS;
  double N0;
  double pe;
  double pa;
  double tol;
  double n;
  double alpha;
  double alpha2;
  double n_alpha;
  double nlo;
  double nhi;
  double lambdalo;
  double lambdahi;
} type_par;


typedef struct /* structure used in moult_sim.c */
{
  int alive;
  int r;
  int f1;
  int ef1;
  int ml1;
  int f2;
  int ef2;
  int ml2;
  int e;
  int a;
  int birth;
  int al;
} t_ind;

