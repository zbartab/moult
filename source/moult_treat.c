

#include"moult.h"


extern double G[MT+1][ME+1];
double G_ori[MT+1][ME+1];


double youngs[MT+1][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
double repval[MT+1][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
double states[MT+1][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
double uvalue[MT][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
int    action[MT][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
int    moult1[MT][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];
int    moult2[MT][MR+1][MF-MO+1][MF-MO+1][ME+1][MA+1];

type_par P;

#define N_PAR 8
#define MAX_TREAT 20

double pa = 1.0; /* factor affecting the length of a breeding attempt */

long idum=(-13);

FILE *OutPut3;

double Treats[MAX_TREAT][N_PAR];
int n_tr;

extern void Ind(char *, int);
extern void iniFP(void);

void ReadIn(char *);
void BinRead(char *);
void hiba(char *);
void Treatment(int, t_ind *);
void Food_Treatment(int, int);
void Read_gain(char *) ;
void get_treat(int, char **, char *, char *, int *); 

/************************* main *************************/
int main(int argc, char *argv[])
{
  char filename[128];
  char filen[128];
  char ext[128];
  int i, j;
  float differ, pdiff;
  float In, D;

  idum = (long) time(NULL);
  if (idum > 0) idum = -idum;

  get_treat(argc, argv, filename, ext, &n_tr);
  strcpy(filen,""); strcat(filen,filename); strcat(filen,"-");
  strcat(filen,ext);strcat(filen,".sim");
  if ((OutPut3 = fopen(filen,"wt")) == NULL) {
    printf("Can't open the output file\n");
    exit (1);
  }
  fprintf(OutPut3,"Seed: %ld\n", idum);
  for(i = 0; i <= n_tr; i++) {
    for(j = 0; j < 5; j++) fprintf(OutPut3, "%lf ", Treats[i][j]);
    fprintf(OutPut3, "\n");
  }
  ReadIn(filename);
  Read_gain(filename);
  iniFP();
  BinRead(filename);
  strcpy(filen,""); strcat(filen,filename); strcat(filen,"-");
  strcat(filen,ext);
  Ind(filen, Treats[0][0]);
  fclose(OutPut3);
}
/************************* end of main ******************/

#define YES 1
#define NO 0

/************************* get_treat *************************/
void get_treat(int argc, char *argv[], char *filename, char *ext, int *n) 
     /* args format -t ido1 ido2 food f1 f2 filename extension */
{
  int i = 1, j = 0, l, p = -1;
  int READ = NO, name_read = NO;

  strcpy(filename,"baseline");
  strcpy(ext,"sim");
  while (i < argc) {
    if (READ == YES) {
      Treats[p][j] = atof(argv[i]);
      j++;
      if (j >= N_PAR) READ = NO;
    } else {
      if (strcmp(argv[i],"-t") == 0) {
	READ = YES;
	p++;
	j = 0;
      } else {
	if (name_read == NO) {
	  strcpy(filename,argv[i]);
	  name_read = YES;
	} else {
	  strcpy(ext,argv[i]);
	  break;
	}
      }
    }
    i++;
  }
  *n = p;
}
/************************* end of get_treat ******************/


/************************* BinRead *************************/
void BinRead(char *filename)
{
  int t, r, f1, e, a, l;
  char filen[30];
  FILE *BinOut;
  int n_y, n_a;
  int wrote;

  strcpy(filen,""); strcat(filen,filename); strcat(filen,".bin");
  if ((BinOut = fopen(filen,"rb")) == NULL) {
    printf("Can't open the input file\n");
    exit (1);
  }
  if (fread(&n_y,sizeof(int),1,BinOut) != 1) hiba("BinWrite 1");
  if (n_y != (MT+1)*(MR+1)*(MF-MO+1)*(MF-MO+1)*(ME+1)*(MA+1)) 
    hiba("Size of youngs is not correct");
  if (fread(repval, sizeof(double),n_y,BinOut) != n_y) hiba("BinWrite 2b");
  if (fread(youngs, sizeof(double),n_y,BinOut) != n_y) hiba("BinWrite 2a");
  if (fread(&n_a,sizeof(int),1,BinOut) != 1) hiba("BinWrite 3");
  if (n_a != MT*(MR+1)*(MF-MO+1)*(MF-MO+1)*(ME+1)*(MA+1))
    hiba("Size of action is not correct");
  if (fread(uvalue, sizeof(double),n_a,BinOut) != n_a) hiba("BinWrite 4");
  if (fread(action, sizeof(int),n_a,BinOut) != n_a) hiba("BinWrite 5");
  if (fread(moult1, sizeof(int),n_a,BinOut) != n_a) hiba("BinWrite 6");
  if (fread(moult2, sizeof(int),n_a,BinOut) != n_a) hiba("BinWrite 7");
  fclose(BinOut);
}
/************************* BinRead *************************/

/************************* hiba *************************/
void hiba(char *leiras)
{
  fprintf(stderr, "%s\n", leiras);
  exit(1);
}
/************************* hiba *************************/


/************************* Treatment *************************/
void Treatment(int t, t_ind *Bird)
/*void Treatment(int t, double *rf, double *ff1, double *ff2, double *u)*/
{
  double s;
  int i,ido;

  for(i = 0; i <= n_tr; i++) {
    if (t >= Treats[i][0] && t<=Treats[i][1]) {
      s = (*Bird).r + Treats[i][2];
      (*Bird).r = (s > MR ? MR : (s < 0 ? 0 : s));
      if ((*Bird).f1>=0) 
	(*Bird).f1 = ((s=(*Bird).f1+Treats[i][3])<0 ? 0 : (s>MF ? MF: s));
      if ((*Bird).f2>=0) 
	(*Bird).f2 = ((s=(*Bird).f2+Treats[i][4])<0 ? 0 : (s>MF ? MF: s));
      if ((*Bird).a == Treats[i][5]) (*Bird).a = Treats[i][6];
    }
  }
}
/************************* Treatment *************************/


/************************* Food_Treatment *************************/
void Food_Treatment(int ido, int t)
{
  double s;
  int i,e;

  for(e = 0; e <= ME; e++) {
    G[t][e] = G_ori[t][e];
  }
  for(i = 0; i <= n_tr; i++) {
    if (ido >= Treats[i][0] && ido<=Treats[i][1]) {
      for(e = 0; e <= ME; e++) {
	G[t][e] *= Treats[i][N_PAR-1];
      }
    }
  }
}
/************************* Treatment *************************/


/************************* Read_gain *************************/
void Read_gain(char *filename) 
{
  char filen[30];
  FILE *GAIN_FILE;
  int t, e;
  double j1, j2, j3;

  strcpy(filen,""); strcat(filen,filename); strcat(filen,".yrv");
  if ((GAIN_FILE = fopen(filen,"rt")) == NULL) {
    printf("Can't open the input file\n");
    exit (1);
  }
  for(t = 0; t < MT; t++) {
    if (fscanf(GAIN_FILE,"%lf%lf%lf",&j1,&j2,&j3) != 3) 
      hiba("Read in error in Read_gain:1");
    for(e = 0; e <= ME; e++) {
      if (fscanf(GAIN_FILE,"%lf%lf%lf",&(G_ori[t][e]),&j1,&j2) !=3)
      hiba("Read in error in Read_gain:2");
    }
  }
  fclose(GAIN_FILE);
}     
/************************* end of Read_gain ******************/

