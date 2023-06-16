#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "cliquer.h"
#include "hsp_utils.h"


#define HSP_DEBUGGING 0
#define BUFF_SIZE 20

void hsp_data2json(FILE* fdat, int p, int h, int* row, int* col){
  fprintf(fdat, "{\n");
  fprintf(fdat, "  \"prime\": %d,\n", p);
  fprintf(fdat, "  \"hadamard_order\": %d,\n", h);
  
  fprintf(fdat, "  \"row_indices\": [");
  for (int i = 0; i < h; i++) {
    fprintf(fdat, "%d", row[i]);
    if (i < h - 1) {
      fprintf(fdat, ", ");
    }
  }
  fprintf(fdat, "],\n");
  
  fprintf(fdat, "  \"col_indices\": [");
  for (int i = 0; i < h; i++) {
    fprintf(fdat, "%d", col[i]);
    if (i < h - 1) {
      fprintf(fdat, ", ");
    }
  }
  fprintf(fdat, "]\n");
  fprintf(fdat, "},\n");
}

void hsp_create_row_indices(int* ix, int* row_indices, int p, int h){
  int row_index=0;
  for(int i=0; i<p; i++){
    if (!hsp_inArray(i,h,ix)){
      row_indices[row_index]=i;
      row_index++;
    }
  }
}

void hsp_rand_row_indices(int* ix, int* row_indices, int N, int p, int h){
  int row_index=0;
  int drawn_number;
  
  while(row_index<N){
    drawn_number=rand()%p;
    if(!hsp_inArray(drawn_number,row_index,row_indices)
       && !hsp_inArray(drawn_number,h,ix))
      {
	row_indices[row_index]=drawn_number;
	row_index++;
      }
  }
}
    

int hsp_rand_check_prime(int p, int h, int NRand, FILE* fpass, FILE* ffail, FILE* fdat)
// randomised check for an hxh Hadamard submatrix in Paley(p)
// if the test passes p is written in fpass, the corresponding data in fdat
// if the test fails, p is written in ffail
{
  int** M;
  int* cols = malloc(sizeof(int)*h);
  int* row_indices = malloc(sizeof(int)*(p-h));

  // init M = Paley(Q)
  hsp_allocmat(&M,p,p);
  hsp_paley(p, &M);

  // generate random columns until success or tries == NRand

  int success = 0;
  int ntries = 0;

  graph_t* g=graph_new(p-h);
  ASSERT(graph_test(g,NULL));
  set_t s;
  int clique_size;
  int* rows = malloc(sizeof(int)*h);
  int* clique_vertices = malloc(sizeof(int)*h);

  while( !success && ntries < NRand ){
    hsp_randarray(0,p-1,h,cols);
    hsp_create_row_indices(cols, row_indices, p, h);
    if(HSP_DEBUGGING){
      hsp_disparray(&cols,h);
    }
    // fill in graph
    for(int i=0; i<p-h; i++){
      for(int j=i; j<p-h; j++){
	if(hsp_areOrthogonal_rc(i,j,h,row_indices,cols,&M)){
	  GRAPH_ADD_EDGE(g,i,j);
	}
      }
    }
    // call to cliquer
    // I redirect stdout to /dev/null and then restore
    // to avoid printing cliquer's output at every call
    freopen("/dev/null", "w", stdout); 
    s = clique_find_single(g,h,h,FALSE,NULL);
    freopen("/dev/tty", "w", stdout);  
    
    ntries++;
    if(s!=NULL){
      if(HSP_DEBUGGING) set_print(s);
      clique_size = set_size(s);
      if(clique_size==h){
	success=1;
	hsp_set2array(s,clique_vertices);
	hsp_disparray(&clique_vertices,h);
      }
    }
    g=graph_new(p-h);
  }
  
  if(success){
    fprintf(fpass,"%d\n",p);
    for(int i=0; i<h; i++){
      rows[i]=row_indices[clique_vertices[i]];
    }
    hsp_data2json(fdat, p, h, rows, cols);
  }
  else{
    fprintf(ffail,"%d : %d\n", ntries, p);
  }

  free(clique_vertices);
  free(rows);
  graph_free(g);
  hsp_freemat(&M,p,p);
  free(row_indices);
  free(cols);
  return success;
}

int hsp_rand_rows_check_prime(int p, int h, int Nrows, int NRand, FILE* fpass, FILE* ffail, FILE* fdat)
// randomised check for an hxh Hadamard submatrix in Paley(p)
// here we choose Nrows random rows and h columns
// if the test passes p is written in fpass, the corresponding data in fdat
// if the test fails, p is written in ffail
{
  int** M;
  int* cols = malloc(sizeof(int)*h);
  int* row_indices = malloc(sizeof(int)*(Nrows));

  // init matrix M
  hsp_allocmat(&M,Nrows,h);

  // generate random columns until success or tries == NRand

  int success = 0;
  int ntries = 0;

  graph_t* g=graph_new(Nrows);
  ASSERT(graph_test(g,NULL));
  set_t s;
  int clique_size;
  int* rows = malloc(sizeof(int)*h);
  int* clique_vertices = malloc(sizeof(int)*h);

  while( !success && ntries < NRand ){
    hsp_randarray(0,p-1,h,cols);
    hsp_rand_row_indices(cols, row_indices, Nrows, p, h);
    // fill in the matrix
    hsp_paley_small(p,Nrows,h,row_indices,cols,&M);
    
    if(HSP_DEBUGGING){
      hsp_disparray(&cols,h);
    }
    // fill in graph
    for(int i=0; i<Nrows; i++){
      for(int j=i; j<Nrows; j++){
	if(HSP_DEBUGGING){
	  printf("Attempting i=%d, j=%d, row[i]=%d, col[j]=%d\n",i,j,
		 row_indices[i],cols[j]);
	}
	if(hsp_areOrthogonal(i,j,h,&M)){
	  GRAPH_ADD_EDGE(g,i,j);
	}
      }
    }
    // call to cliquer
    // I redirect stdout to /dev/null and then restore
    // to avoid printing cliquer's output at every call
    freopen("/dev/null", "w", stdout); 
    s = clique_find_single(g,h,h,FALSE,NULL);
    freopen("/dev/tty", "w", stdout);  
    
    ntries++;
    if(s!=NULL){
      if(HSP_DEBUGGING){
	printf("Found clique!\n");
	set_print(s);
      }
      clique_size = set_size(s);
      if(clique_size==h){
	success=1;
	hsp_set2array(s,clique_vertices);
	hsp_disparray(&clique_vertices,h);
      }
    }
    g=graph_new(Nrows);
  }
  
  if(success){
    fprintf(fpass,"%d\n",p);
    for(int i=0; i<h; i++){
      rows[i]=row_indices[clique_vertices[i]];
    }
    hsp_data2json(fdat, p, h, rows, cols);
  }
  else{
    fprintf(ffail,"%d : %d\n", ntries, p);
  }

  free(clique_vertices);
  free(rows);
  graph_free(g);
  hsp_freemat(&M,Nrows,h);
  free(row_indices);
  free(cols);
  return success;
}

 
int main(int argc, char* argv[]){

  if (argc < 4) {
    printf("Usage: %s <hadamard_order> <number_of_tries> <primes_file>\n", argv[0]);
    printf("or\n");
    printf("Usage: %s <hadamard_order> <number_of_tries> <primes_file> <n_random_rows>\n", argv[0]);
    return 1; 
  }
  
  int h = atoi(argv[1]);
  int NRand = atoi(argv[2]);
  char* fin_name = argv[3];
  int Nrows;
  if(argc>4){
    Nrows=atoi(argv[4]);
  }
  else{
    Nrows=-1;
  }
  char fpass_name[BUFF_SIZE];
  char ffail_name[BUFF_SIZE];
  char fdata_name[BUFF_SIZE];

  sprintf(fpass_name, "data/pass_%d",h);
  sprintf(ffail_name, "data/fail_%d",h);
  sprintf(fdata_name, "data/data_%d.json",h);
  
  FILE* fpass = fopen(fpass_name, "w");
  FILE* ffail = fopen(ffail_name, "w");
  FILE* fdata = fopen(fdata_name, "w");
  FILE* fin = fopen(fin_name, "r");

  if(fin==NULL){
    printf("Error reading input file\n");
    return 1;
  }

  if(fpass==NULL || ffail==NULL || fdata==NULL){
    printf("Error opening output files.\n");
    return 1;
  }
  else{
    int seed = time(NULL);
    srand(seed);
    
    int p;
    int success;
    fprintf(fdata,"[");
    while (fscanf(fin, "%d", &p) == 1) {
	printf("Checking %d..\n",p);
	if(Nrows==-1){
	  success=hsp_rand_check_prime(p,h,NRand,fpass,ffail,fdata);
	}
	else{
	  success=hsp_rand_rows_check_prime(p,h,Nrows,NRand,fpass,ffail,fdata);
	}
	printf("Done: ");
	if(success) printf("SUCCESS!\n");
	else printf("FAILED!\n");
    }
    fseek(fdata,-2,SEEK_CUR);
    fputc(']',fdata);
    
    fclose(fin);
    fclose(fdata);
    fclose(ffail);
    fclose(fpass);
  }

  return 0;
}
