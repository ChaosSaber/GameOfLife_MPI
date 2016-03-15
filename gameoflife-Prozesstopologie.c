#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <endian.h>
#include <omp.h>
#include "mpi.h"

#define calcIndex(width, x,y)  ((y)*(width) + (x))

void show(unsigned* currentfield, int w, int h) {
  printf("\033[H");
  for (int y = 0; y < h; y++) {
    for (int x = 0; x < w; x++) printf(currentfield[calcIndex(w, x,y)] ? "\033[07m  \033[m" : "  ");
    printf("\033[E");
  }
  fflush(stdout);
}


float convert2BigEndian( const float inFloat )
{
   float retVal;
   char *floatToConvert = ( char* ) & inFloat;
   char *returnFloat    = ( char* ) & retVal;

   // swap the bytes into a temporary buffer
   returnFloat[0] = floatToConvert[3];
   returnFloat[1] = floatToConvert[2];
   returnFloat[2] = floatToConvert[1];
   returnFloat[3] = floatToConvert[0];

   return retVal;
}

void writeVTK(unsigned* currentfield, int w, int h, int t, char* prefix) {
  char name[1024] = "\0";
  sprintf(name, "%s_%d.vtk", prefix, t);
  FILE* outfile = fopen(name, "w");

  /*Write vtk header */                                                           
  fprintf(outfile,"# vtk DataFile Version 3.0\n");       
  fprintf(outfile,"frame %d\n", t);     
  fprintf(outfile,"BINARY\n");     
  fprintf(outfile,"DATASET STRUCTURED_POINTS\n");     
  fprintf(outfile,"DIMENSIONS %d %d %d \n", w, h, 1);        
  fprintf(outfile,"SPACING 1.0 1.0 1.0\n");//or ASPECT_RATIO                            
  fprintf(outfile,"ORIGIN 0 0 0\n");                                              
  fprintf(outfile,"POINT_DATA %d\n", h*w);                                    
  fprintf(outfile,"SCALARS data float 1\n");                              
  fprintf(outfile,"LOOKUP_TABLE default\n");         
 
  for (int y = 0; y < h; y++) {
    for (int x = 0; x < w; x++) {
      float value = currentfield[calcIndex(w, x,y)]; // != 0.0 ? 1.0:0.0;
      value = convert2BigEndian(value);
      fwrite(&value, 1, sizeof(float), outfile);
    }
  }
  fclose(outfile);
}
int coutLifingsPeriodic(unsigned* currentfield, int x , int y, int w, int h) {
   int n = 0;
   for (int y1 = y - 1; y1 <= y + 1; y1++) 
     for (int x1 = x - 1; x1 <= x + 1; x1++) 
       if (currentfield[calcIndex(w, (x1 + w) % w, (y1 + h) % h)]) 
         n++;
   return n;
}
 
int evolve(unsigned* currentfield, unsigned* newfield, int w, int h) {
int changes = 0;
  for (int y = 0; y < h; y++) {
    for (int x = 0; x < w; x++) {
      int n = coutLifingsPeriodic(currentfield, x , y, w, h);
      if (currentfield[calcIndex(w, x,y)]) 
        n--;
      newfield[calcIndex(w, x,y)] = (n == 3 || (n == 2 && currentfield[calcIndex(w, x,y)]));   
      if(newfield[calcIndex(w, x,y)] != currentfield[calcIndex(w, x,y)])
        changes++;
    }
  } 
  return changes;
}
 
void filling(unsigned* currentfield, int w, int h) {
  for (int i = 0; i < h*w; i++) {
    currentfield[i] = (rand() < RAND_MAX / 10) ? 1 : 0; ///< init domain randomly
  }
}
 
void game(int w, int h, int timesteps) {
  unsigned *currentfield = calloc(w*h, sizeof(unsigned));
  unsigned *newfield     = calloc(w*h, sizeof(unsigned));
  
  filling(currentfield, w, h);
  for (int t = 0; t < timesteps; t++) {
    show(currentfield, w, h);
    writeVTK(currentfield, w, h, t, "out/output");
    int changes = evolve(currentfield, newfield, w, h);
    if (changes == 0) {
    	sleep(3);
    	break;
    }
    
    usleep(200000);

    //SWAP
    unsigned *temp = currentfield;
    currentfield = newfield;
    newfield = temp;
  }
  
  free(currentfield);
  free(newfield);
}
 
int main(int c, char **v) {
     
  int w, h, timesteps = 10;
///< read height
  
   w = 30; ///< default width
   h = 30; ///< default height
  MPI_Init(&c, &v);
  MPI_Comm cart_comm;
  int dims[1];
  int periodics[1];
  dims[0]=3;
  periodics[0]=0;
  MPI_Cart_create( MPI_COMM_WORLD, 1, dims, periodics, 0, &cart_comm);
  int myrank;
  MPI_Comm_rank(cart_comm, &myrank);
  //printf("%d\n", myrank);
  int links,rechts;
  MPI_Cart_shift(cart_comm, 0, 1, &links, &rechts);
  if(myrank == 0)
     printf("Process %d: keinen linken Nachbar, rechts: %d\n",myrank, rechts);
   else if(myrank == (dims[0] - 1))
      printf("Process %d: links: %d, keinen rechten Nachbar\n",myrank, links);
    else
      printf("Process %d: links: %d,rechts: %d\n",myrank, links, rechts);
 // game(w, h, timesteps);
  MPI_Finalize();
}
