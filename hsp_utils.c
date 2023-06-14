#include "hsp_utils.h"

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
    for(int j=0; j<n; j++) printf("%d ",(*M)[i][j]);
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

void hsp_paley_ix(int p, int h, int* ix, int* rx, int*** Q)
// Q is at least a pxh array, and len(ix)=h
// assuming: ix[h-1]<p
// returns rows Q[i,ix] with i not in ix
// rx is set to be the indices of rows 
{
  int nrows = 0; //number of rows added 
  for(int i=0; i<p; i++){
    if(!hsp_inArray(i,h,ix)){
      for(int j=0; j<h; j++){
	(*Q)[nrows][j]=(1-hsp_jacobi(p+(ix[j]-i),p))/2;
      }
      rx[nrows]=i;
      nrows++;
    }
  }
}

int hsp_areOrthogonal(int i, int j, int len, int***Q)
// returns 1 if rows i and j of Q are orthogonal, 0 otherwise
{
  int matches = 0;
  for(int k = 0; k<len; k++){
    if ((*Q)[i][k] == (*Q)[j][k]){
      matches++;
    }
    else{
      matches--;
    }
  }
  return matches==0;
}

int hsp_areOrthogonal_ix(int i, int j, int len, int* ix, int***Q)
// returns 1 if rows i and j, restricted to ix, are orthogonal, 0 o.w.
// do not use if i or j are in ix
{
  int matches = 0;
  for(int k=0; k<len; k++){
    if((*Q)[i][ix[k]]==(*Q)[j][ix[k]]){
      matches++;
    }
    else{
      matches--;
    }
  }
  return matches==0;
}

int hsp_areOrthogonal_rc(int i, int j, int h, int* row_indices, int* col_indices, int*** Q){
  // ret 1 if row[i] and row[j] restricted to col_indices are orthogonal, 0 o.w
  int matches = 0;
  int row_i, row_j;
  row_i = row_indices[i];
  row_j = row_indices[j];
  
  int col_k;
  for(int k=0; k<h; k++){
    col_k = col_indices[k];
    if((*Q)[row_i][col_k]==(*Q)[row_j][col_k]){
      matches++;
    }
    else{
      matches --;
    }
  }
  return matches == 0;
}

void hsp_set2array(set_t s, int* v){
  // assuming v has at least set_size allocated
  int i=0;
  int item=-1;
  while((item=set_return_next(s,item)) >= 0){
    v[i]=item;
    i++;
  }
}

void hsp_disparray(int** v, int len){
  printf("[");
  for(int i=0; i<len-1; i++) printf("%d, ",(*v)[i]);
  printf("%d]\n",(*v)[len-1]);
}

void hsp_residuemat(int la, int lb, int p, int* a, int* b, int***M){
  for(int i=0; i<la; i++){
    for(int j=0; j<lb; j++){
      (*M)[i][j]=(p+(a[i]-b[j]))%p;
    }
  }
}

void hsp_randarray(int a, int b, int len, int* v)
// v is an array of length len, it is filled
// randomly with elements in the range [a..b]
{
  int drawn_number;
  int ndrawn=0;
  while(ndrawn<len){
    drawn_number=(rand() % (b-a+1))+a;
    if(! hsp_inArray(drawn_number,ndrawn,v)){
      v[ndrawn]=drawn_number;
      ndrawn++;
    }
  }
}
