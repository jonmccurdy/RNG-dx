#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include <R_ext/Random.h>
using namespace Rcpp;

#define K_MAX	3500
#define PP 2147483647  //2^31-1
#define IPP 4.6566129e-10
#define ROTR(x,n) ((((unsigned long) x)>>(n))|(((unsigned long) x)<<(32-(n))))
#define f1(x) (ROTR(x, 25))
#define f2(x) (ROTR(x, 23))
#define f3(x) (ROTR(x, 12))
#define f4(x) (ROTR(x, 10))
#define TSIZE 256 		// actual table size used
#define maxT 1024 		//maximum table size
#define maxB 15
#define B_X 32896  		
#define B_Y 32776
#define T_MOD TSIZE-1
#define MOD2_32 4294967295  
static unsigned long  K[8], IV[8];   /*KEY and IV vectors*/
static int TABLE_X[maxT], TABLE_Y[maxT];
static int I_Y, BURN_X, BURN_Y;

static long long B_X1 = 536869888; // will change based on K, S- can be done more efficiently
static long long B_X2 = 65011712;
static long long B_X3 = 67633152;
static long long B_X4 = 67108736;

static int I_X;           /* running index */

unsigned long MODP(unsigned long z) {return ((((z)&PP)+((z)>>31)) &PP);}

static int K_X = 47;
static int K_SX = 4;
static int K_SY = 5;
static int S_X = 1;
static int nseed = 5;
static int seed[K_MAX]={K_X,S_X,K_SX,K_SY};
static int XX[K_MAX], YY[K_MAX];
static double res;
int set=0;

int i, K12, K13, K23;

void dx_1(){
  K_X = seed[0];
  S_X = seed[1];
  int II0 = I_X;
  if(++I_X >= K_X)  I_X = 0;     /*wrap around running index */
  XX[I_X] = MODP(B_X1 * XX[I_X] + XX[II0]);
  res = (double) XX[I_X] * IPP;
}

void dx_2(){
  K_X = seed[0];
  S_X = seed[1];
  int II0 = I_X;
  if(++I_X >= K_X)  I_X = 0;     /*wrap around running index */
  XX[I_X] = MODP(B_X2 * (XX[I_X] + XX[II0]));
  res = (double) XX[I_X] * IPP;
}

void dx_3(){
  K_X = seed[0];
  S_X = seed[1];
  int II0 = I_X;
  if(++I_X >= K_X)  I_X = 0;     /*wrap around running index */
  if(++K12 >= K_X) K12 = 0;    /*wrap around K12*/
  XX[I_X] = MODP(B_X3 * (XX[I_X] + XX[K12] + XX[II0]));
  res= (double) XX[I_X] * IPP;
}

void dx_4(){
  K_X = seed[0];
  S_X = seed[1];
  int II0 = I_X;
  if(++I_X >= K_X)   I_X = 0;    /*wrap around running index */
  if(++K13 >= K_X) K13 = 0;    /*wrap around K13*/
  if(++K23 >= K_X) K23 = 0;    /*wrap around K23*/
  XX[I_X] = MODP(B_X4 * (XX[I_X]+XX[K13]+XX[K23]+XX[II0]));
  res= (double) XX[I_X] * IPP;
}

/*  SAFE */

void X_RNG_INIT(){
  int i;
  srand((seed[1]+K[0]) & PP);
  for (i=0; i<8; i++)
    XX[i] = (rand()+K[i]+IV[8-i]) & PP; 
  for (i=8; i<K_MAX; i++)
    XX[i] = rand() & PP; 
  I_X=K_SX-1;
}

int RNG_X(void){
  int II0 = I_X;
  if(++I_X >= K_SX)  I_X = 0;     /*wrap around running index */
  XX[I_X] = MODP(B_X * XX[I_X] + XX[II0]);
return XX[I_X];
}

void Y_RNG_INIT(){
  int i;
  srand((seed[1]+K[1]) & PP);
  for (i=0; i<8; i++)
    YY[i] = (rand()+K[8-i]+IV[i]) & PP;
  for (i=8; i<K_MAX; i++)
    YY[i] = rand() & PP; 
  I_Y=K_SY-1;
}

int RNG_Y(void){
  int II0 = I_Y;
  if(++I_Y >= K_SY)  I_Y = 0;     /*wrap around running index */
  YY[I_Y] = MODP(B_Y * YY[I_Y] + YY[II0]);
  return YY[I_Y];
}

void init_SAFE_X_Y(void){
  int i;
  
  srand(seed[4]);
  for (i = 0; i < 8; i++){
    K[i] = rand();
    IV[i] = rand();
  }
  
  X_RNG_INIT();
  Y_RNG_INIT();
  
  BURN_X = (K[0]+ K[3] + IV[1] + IV[2]) & maxB;
  BURN_Y = (K[1]+ K[2] + IV[0] + IV[3]) & maxB;
  for(i=0; i<BURN_X; i++) RNG_X(); 
  for(i=0; i<BURN_X; i++) RNG_Y(); 
  
  for(i=0; i<TSIZE; i++) 
    TABLE_X[i] = RNG_X(); 
  for(i=0; i<TSIZE; i++) 
    TABLE_Y[i] = RNG_Y(); 
}

void SAFE_X_Y(void){
  unsigned long x0, y0, vv, ww;
  int X0, Y0;
  x0 = RNG_X();
  y0 = RNG_Y();
  Y0=y0&T_MOD;
  X0=x0&T_MOD;
  
  vv = TABLE_X[Y0]; TABLE_X[Y0]=x0;
  ww = TABLE_Y[X0]; TABLE_Y[X0]=y0;
  
  res=(double)((f1(x0)+f2(y0)+f4(ww)) & MOD2_32)* IPP;
}
/*  SAFE */

void dx_init(){
  srand(seed[4]);
  for (i=0; i< K_MAX; i++){
    XX[i] = rand() & PP;
  }
}

void (*dx_gen)();

void generator_type(){
  switch(S_X){
  case 1: dx_gen=&dx_1;
    dx_init();
    init_SAFE_X_Y();
    break;
  case 2: dx_gen=&dx_2;
    dx_init();
    init_SAFE_X_Y();
    break;
  case 3: dx_gen=&dx_3;
    dx_init();
    init_SAFE_X_Y();
    break;
  case 4: dx_gen=&dx_4;
    dx_init();
    init_SAFE_X_Y();
    break;
  case 13: dx_gen=&SAFE_X_Y;
    dx_init();
    init_SAFE_X_Y();
    break;
  default: dx_gen=&dx_1;
    dx_init();
    init_SAFE_X_Y();
    S_X=1;
    seed[1]=S_X;
    Rcerr << "The value of 'S' that was chosen is not compatible with this package. \nBy default, a DX-" << K_X<< "-1 generator was chosen.\n";
    break;
  }
}

double * user_unif_rand(){
  dx_gen();
  return &res;
}

void user_unif_init(Int32 seed_in) {
  K_X = seed[0];
  S_X = seed[1];
  K_SX = seed[2];
  K_SY = seed[3];
  
  srand(seed_in);
  seed[4] = rand() & PP;
  
  generator_type();

  I_X = seed[0]-1;
  K12 = seed[0]/2-1;
  K13 = 2*seed[0]/3-1;   
  K23 = seed[0]/3-1;
}

int * user_unif_nseed() {return &nseed; }
int * user_unif_seedloc() { return (int *) &seed; }
