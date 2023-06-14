#include <stdio.h>
#include <stdlib.h>
#include "cliquer.h"

#define HSP_DEBUGGING 1

void hsp_log2(int a, int* a1, int* l2)
// a = a1 * 2^l2;
{
  *a1=a;
  *l2=0;
  while(*a1 %2 == 0){
    *a1 >>= 1;
    *l2 +=1;
  }
}

int hsp_jacobi(int a, int n)
//Computes the Jacobi symbol (a|n) 
{
  int symbol = 1;
  int a0=a;
  int a1;
  int l2;
  
  while(1)
    {
      a0 = a0%n;
      if(a0 == 0){
	symbol = n == 1 ?  symbol :  0;
	return symbol;
      }
      hsp_log2(a0,&a1,&l2);
      if(l2%2 == 1 && (n%8 != 1 && n%8 != 7)) symbol=-symbol;
      if(a1%4 != 1 && n%4 !=1) symbol=-symbol;
      a0=n;
      n=a1;
    }
}

void hsp_paley(int p, int*** Q)
 // Q is a matrix of size pxp
{ 
  for(int i=0; i<p; i++){
    for(int j=0; j<p; j++){
      (*Q)[i][j]=(1-hsp_jacobi(p+(j-i),p))/2; // displayed as 0-1 matrix
    }
  }
}

void hsp_allocmat(int*** M, int m, int n)
// allocates a matrix of dimensions mxn
{
  (*M) = malloc(sizeof(int*)*m);
  for(int i=0; i<m; i++){
    (*M)[i]=malloc(sizeof(int)*n);
  }
}

void hsp_freemat(int*** M, int m, int n)
// frees a matrix of dimensions mxn
{
  for(int i=0; i<m; i++) free((*M)[i]);
  free(*M);
}

void hsp_dispmat(int*** M, int m, int n)
// display matrix
{
  for(int i=0; i<m; i++){
    for(int j=0; j<n; j++) printf("%d",(*M)[i][j]);
    printf("\n");
  }
}

int hsp_inArray(int i, int l, int* v)
// returns true if i is in the array v of length l
{
  for(int j=0; j<l; j++){
    if(i==v[j]) return 1;
  }
  return 0;
}

void hsp_paley_ix(int p, int h, int* ix, int*** Q)
// Q is at least a pxh array, and len(ix)=h
// assuming: ix[h-1]<p
// returns rows Q[i,ix] with i not in ix
{
  int nrows = 0; //number of rows added 
  for(int i=0; i<p; i++){
    if(!hsp_inArray(i,h,ix)){
      for(int j=0; j<h; j++){
	(*Q)[nrows][j]=(1-hsp_jacobi(p+(ix[j]-i),p))/2;
      }
      nrows++;
    }
  }
  if(HSP_DEBUGGING) printf("rows added = %d\n", nrows);
}

int hsp_areOrthogonal(int i, int j, int len, int**Q)
// returns 1 if rows i and j of Q are orthogonal, 0 otherwise
{
  int matches = 0;
  for(int k = 0; k<len; k++){
    if (Q[i][k] == Q[j][k]){
      matches++;
    }
    else{
      matches--;
    }
  }
  return matches==0;
}


int main(void){
  int** M;
  int ix[4] = {1,3,5,7};
  
  hsp_allocmat(&M,13,4);
  hsp_paley_ix(13,4,ix,&M);
  hsp_dispmat(&M,13-4,4);
  printf("Creating graph..\n");
  graph_t* g = graph_new(9);
  for(int i=0; i<9; i++){
    for(int j=i+1; j<9; j++){
      if(hsp_areOrthogonal(i,j,4,M)){
	GRAPH_ADD_EDGE(g,i,j);
      }
    }
  }
  ASSERT(graph_test(g,NULL));
  set_t s;
  s=clique_find_single(g,3,3,FALSE,NULL);
  set_print(s);
  
  hsp_freemat(&M,13,4);
  return 0;
}
