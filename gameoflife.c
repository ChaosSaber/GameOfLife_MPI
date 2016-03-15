#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <endian.h>
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

void writeVTK(unsigned* currentfield, int w, int h, int t, char* prefix, int myrank) 
{
    char name[1024] = "\0";
    sprintf(name, "%s_%d_%d.vtk", prefix, myrank, t);
    FILE* outfile = fopen(name, "w");

    /*Write vtk header */                                                           
    fprintf(outfile,"# vtk DataFile Version 3.0\n");       
    fprintf(outfile,"FRAME %d\n", t);     
    fprintf(outfile,"BINARY\n");     
    fprintf(outfile,"DATASET STRUCTURED_POINTS\n");     
    fprintf(outfile,"DIMENSIONS %d %d %d \n",  w, h, 1);        
    fprintf(outfile,"SPACING 1.0 1.0 1.0\n");//or ASPECT_RATIO                            
    fprintf(outfile,"ORIGIN %d %d 0\n", 0, (h - 1) * myrank );                                              
    fprintf(outfile,"POINT_DATA %d\n", w*h);
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

int coutLifingsPeriodic(unsigned* currentfield, unsigned* upperGhostLayer, unsigned* lowerGhostLayer, int x , int y, int w, int h) {
   int n = 0;
   for (int y1 = y - 1; y1 <= y + 1; y1++) 
     for (int x1 = x - 1; x1 <= x + 1; x1++) 
       if (y1 < 0 && upperGhostLayer[x1]) 
       {
         n++;
       }
       else if (y1 >= h && lowerGhostLayer[x1])
       {
        n++;
       }
       else if (currentfield[calcIndex(w, (x1 + w) % w, (y1 + h) % h)])
       {
        n++;
       }

   return n;
}
 
int evolve(unsigned* currentfield, unsigned* upperGhostLayer, unsigned* lowerGhostLayer, unsigned* newfield, int w, int h) 
{
  int changes = 0;
  for (int y = 0; y < h; y++) {
    for (int x = 0; x < w; x++) {
      int n = coutLifingsPeriodic(currentfield, upperGhostLayer, lowerGhostLayer, x , y, w, h);
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
 
void game(int w, int h, int timesteps, int myrank, unsigned* initialField, MPI_Comm cart_comm) 
{
  unsigned *currentfield = initialField;
  unsigned *newfield     = calloc(w * h, sizeof(unsigned));
  unsigned *upperGhostLayer = calloc(w, sizeof(unsigned));
  unsigned *lowerGhostLayer = calloc(w, sizeof(unsigned));

  int links,rechts;
  MPI_Cart_shift(cart_comm, 0, 1, &links, &rechts);

  MPI_Status *status = calloc(1, sizeof(MPI_Status));
  for (int t = 0; t < timesteps; t++) {
    //show(currentfield, w, h);
    writeVTK(currentfield, w, h, t, "out/output", myrank);


    if (myrank % 2== 0)
    {

   //  kommuniziere mit deinem Rechten 
      MPI_Sendrecv(currentfield + calcIndex(w, 0, h - 1), w,  MPI_UNSIGNED, rechts, 1, lowerGhostLayer, w, MPI_UNSIGNED, rechts, 1, cart_comm, status);
     // kommuniziere mit deinem linken
      MPI_Sendrecv(currentfield + calcIndex(w, 0, 0), w,  MPI_UNSIGNED, links, 1, upperGhostLayer, w, MPI_UNSIGNED, links, 1, cart_comm, status);
    }
    if(myrank & 2 == 1)
    {
   //   kommuniziere mit deinem linken
      MPI_Sendrecv(currentfield + calcIndex(w, 0, 0), w,  MPI_UNSIGNED, links, 1, upperGhostLayer, w, MPI_UNSIGNED, links, 1, cart_comm, status);
     // kommuniziere mit deinem rechten
        MPI_Sendrecv(currentfield + calcIndex(w, 0, h - 1), w,  MPI_UNSIGNED, rechts, 1, lowerGhostLayer, w, MPI_UNSIGNED, rechts, 1, cart_comm, status);
    }

    int changes = evolve(currentfield, upperGhostLayer, lowerGhostLayer, newfield, w, h);
    if (changes == 0) {
      sleep(3);
      break;
    }
    
    //usleep(200000);

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
  
  w = 30; ///< default width
  h = 30; ///< default height

  unsigned* initialField = calloc(w*h, sizeof(unsigned));

  filling(initialField, w, h);

  MPI_Init(&c, &v);
  MPI_Group MPI_GROUP_WORLD, not_group_world;
  MPI_Comm cart_comm;
  MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);

  int process[]={0,1,2,3};
  int dims[1];
  int periodics[1];
  dims[0]=3;
  periodics[0]=1;
  MPI_Cart_create( MPI_COMM_WORLD, 1, dims, periodics, 0, &cart_comm);

  int myrank;
  MPI_Comm_rank(cart_comm, &myrank);

  int links,rechts;
  MPI_Cart_shift(cart_comm, 0, 1, &links, &rechts);

  /*if(myrank == 0)
     printf("Process %d: keinen linken Nachbar, rechts: %d, Gebietsgröße: %d\n",myrank, rechts,  h*w/dims[0]);
   else if(myrank == (dims[0] -1))
      printf("Process %d: links: %d, keinen rechten Nachbar, Gebietsgröße: %d\n",myrank, links,  h*w/dims[0]);
    else*/
    //  printf("Process %d: links: %d,rechts: %d, Gebietsgröße: %d\n",myrank, links, rechts,  h*w/dims[0]);

  int blockHeight = h/dims[0];
  game(w, h/dims[0], timesteps, myrank, &initialField[myrank * blockHeight * w], cart_comm);
  MPI_Finalize();
}