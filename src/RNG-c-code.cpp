#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include <R_ext/Random.h>
using namespace Rcpp;

#define K_MAX	3500
#define PP 2147483647  //2^31-1
#define IPP 4.6566129e-10
static long long B_X1 = 536869888; // will change based on K, S- can be done more efficiently
static long long B_X2 = 65011712;
static long long B_X3 = 67633152;
static long long B_X4 = 67108736;

static int I_X;           /* running index */

unsigned long MODP(unsigned long z) {return ((((z)&PP)+((z)>>31)) &PP);}

static int K_X = 47;
static int S_X = 1;
static int nseed = 3;
static Int32 seed[K_MAX];
static int XX[K_MAX];
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

#define B20(x) ( ((x)>>11) + ((x)<<20)&PP ) 
#define B9(x) ( ((x)>>22) + ((x)<<9) &PP )
unsigned long S ;

void dx_2_fast(){
  K_X = seed[0];
  S_X = seed[1];
  int II0 = I_X;
  if(++I_X >= K_X)  I_X = 0;     /*wrap around running index */
  S = MODP(XX[I_X] + XX[II0]);
  XX[I_X] = MODP(B20(S) + B9(S));
  res = (double) XX[I_X] * IPP;
}

void lcg_basic(){
  seed[2] = seed[2] + 1;
  res = seed[2] * 2.32830643653869e-10;
}

void (*dx_gen)();

void generator_type(){
  switch(S_X){
  case 1: dx_gen=&dx_1;
    break;
  case 2: dx_gen=&dx_2;
    break;
  case 3: dx_gen=&dx_3;
    break;
  case 4: dx_gen=&dx_4;
    break;
  case 5: dx_gen=&dx_2_fast;
    break;
  case 6: dx_gen=&lcg_basic;
    break;
  default: dx_gen=&dx_1;
  Rcerr << "The value of 'S' that was chosen is not compatible with this package. \nBy default, a DX-" << K_X<< "-1 generator was chosen.\n";
  S_X=1;
  break;
  }
}

double * user_unif_rand(){
  dx_gen();
  return &res;
}

void user_unif_init(Int32 seed_in) {
  generator_type();
  seed[0] = K_X;
  seed[1] = S_X;
  srand(seed_in);
  seed[2] = rand() & PP;
  for (i=0; i< seed[0]; i++){
    XX[i] = rand() & PP;
  } 
  I_X = seed[0]-1;
  K12 = seed[0]/2-1;
  K13 = 2*seed[0]/3-1;   
  K23 = seed[0]/3-1;
}

int * user_unif_nseed() {return &nseed; }
int * user_unif_seedloc() { return (int *) &seed; }
