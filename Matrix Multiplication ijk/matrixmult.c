/*
 * @file matrixmult.c
 * @author Long Nguyen <lnguyen34@mail.csuchico.edu>
 * @date 4/8/2018
 *
 * @section DESCRIPTION
 * This program is the fouth programming assignment for CSCI 551
 * (Numerical Methods and Parallel Programming) at CSU Chico
 *
 * The purpose of the program is to use MPI to run Matrix Multiplication
 * usng 3 forms: ijk, ikj and kij
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <mpi.h>

/*
* @name ijk
*
* @section DESCRIPTION
* this function will take matrix A and B calculate it
* using ijk and push it in to matrix c
*
* @param  n size of matrix
* @param  x matrix X
* @param  sendCount number of Count sending
*
* @return none
*/
void ijk(int n, int *A, int *B, int *c, int sendCount)
{
 for(int i = 0; i<sendCount/n; i++)
    for(int j = 0; j<n; j++)
        for(int k = 0; k < n; k++)
           c[i * n + j] += A[i * n + k] * B[k * n + j];
}

/*
* @name ikj
*
* @section DESCRIPTION
* this function will take matrix A and B calculate it
* using ikj and push it in to matrix c
*
* @param  n size of matrix
* @param  x matrix X
* @param  sendCount number of Count sending
*
* @return none
*/
void ikj(int n, int *A, int *B, int *c, int sendCount)
{
 for(int i = 0; i<sendCount/n; i++)
    for(int k = 0; k<n; k++)
        for(int j = 0; j < n; j++)
           c[i * n + j] += A[i * n + k] * B[k * n + j];
}

/*
* @name kij
*
* @section DESCRIPTION
* this function will take matrix A and B calculate it
* using kij and push it in to matrix c
*
* @param  n size of matrix
* @param  x matrix X
* @param  sendCount number of Count sending
*
* @return none
*/
void kij(int n, int *A, int *B, int *c, int sendCount)
{
 for(int k = 0; k<n; k++)
    for(int i = 0; i<sendCount/n; i++)
        for(int j = 0; j < n; j++)
           c[i * n + j] += A[i * n + k] * B[k * n + j];
}

/*
* @name print_matrix
*
* @section DESCRIPTION
* this fucntion will print out the matrix
*
* @param  n size of matrix
* @param  m matrix M
*
* @return none
*/
void print_matrix(int n, int *m)
{
  for(int i = 0; i< n; i++)
  {
    for(int j = 0; j<n; j++)
       printf("%i ", m[i * n + j]);
    printf("\n");
  }
}

/*
* @name create_matrix
*
* @section DESCRIPTION
* this fucntion will create the matrix either random or input
*
* @param  flag random or input
* @param  n    size of matrix
* @param  x    matrix X
*
* @return none
*/

void create_matrix(char flag, int n, int *A,  int *B)
{
  int input;
  //if R create a random matrix
  if (flag == 'R')
  {
      srand(time(NULL));
      for(int i = 0; i< n; i++)
      {
          for(int j = 0; j<n; j++)
          {
              A[i * n + j] = rand() % 101;
              B[i * n + j] = rand() % 101;
          }
      }
  }
  //if I user input matrix
  else if (flag == 'I')
  {
      for(int i = 0; i< n; i++)
      {
         for(int j = 0; j<n; j++)
         {
           scanf ("%i",&input);
           A[i * n + j] = input;
         }
      }
      for(int i = 0; i< n; i++)
      {
          for(int j = 0; j<n; j++)
          {
              scanf ("%i",&input);
              B[i * n + j] = input;
          }
      }
  }
}

/*
* @name main
*
* @section DESCRIPTION
* run the matrix Multiplication
*
* @param  none
*
* @return 0
*/
int main(void)
{
  int        comm_sz;               /* Number of processes    */
  int        my_rank;               /* My process rank        */

  int n;
  char flag, form[4];

  int *A, *B, *c, *sendCount, *displs;
  long double time_start, time_end;

 // Start MPI
 MPI_Init(NULL, NULL);

 // Get the number of proccesses running
 MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

 // Get the current process's rank
 MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

 /*Input from the user */
 if (my_rank == 0)
 {
   /* Input form */
   scanf (" %s",form);
   while ((strcmp(form,"ijk") != 0) && (strcmp(form,"ikj") != 0) && (strcmp(form,"kij") != 0))
   {
     printf("%s\n", "Not valid form. Enter either ijk, ikj or kij" );
     scanf (" %s",form);
   }

   /* Input flag*/
   scanf (" %c",&flag);
   while (flag != 'R' && flag != 'I')
   {
     printf("%s\n", "Not valid flag. Enter either R = random or I = Input" );
     scanf (" %c",&flag);
   }

   /* Input size*/
    scanf ("%i",&n);
    A = (int*)calloc(n*n, sizeof(int));
    B = (int*)calloc(n*n, sizeof(int));
    c = (int*)calloc(n*n, sizeof(int));

    create_matrix(flag,n,A,B);

  }

  /* MPI Time function*/
  MPI_Barrier(MPI_COMM_WORLD);
  /* time start*/
  if (my_rank == 0)
    time_start = MPI_Wtime();

  MPI_Bcast(&form, 4, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* Note: you never send matrix A out
  so u dont need to create matrix A when rank is not 0*/
  if (my_rank != 0)
  {
    B = (int*)calloc(n*n, sizeof(int)),
    c = (int*)calloc(n*n, sizeof(int));
  }

  sendCount = (int*)calloc(comm_sz, sizeof(int)),
  displs = (int*)calloc(comm_sz, sizeof(int));
  int mod = n%comm_sz;

  for (int i = 0; i < comm_sz; i++)
  {
    sendCount[i] = n / comm_sz * n;
    if (mod > 0)
    {
      sendCount[i] = sendCount[i]+ n;
      mod--;
    }
  }

  for (int i = 1; i < comm_sz; i++)
    displs[i] = displs[i - 1] + sendCount[i - 1];

  int *localA = (int*)calloc(sendCount[my_rank], sizeof(int));
  int *localc = (int*)calloc(sendCount[my_rank], sizeof(int));

  MPI_Bcast(B, n * n, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatterv(A, sendCount, displs, MPI_INT, localA, sendCount[my_rank], MPI_INT, 0, MPI_COMM_WORLD);

   if ((strcmp(form,"ijk") == 0))
     ijk(n, localA, B, localc, sendCount[my_rank]);
   else if ((strcmp(form,"ikj") == 0))
     ikj(n, localA, B, localc, sendCount[my_rank]);
   else if ((strcmp(form,"kij") == 0))
     kij(n, localA, B, localc, sendCount[my_rank]);

  MPI_Gatherv(localc, sendCount[my_rank], MPI_INT, c, sendCount, displs, MPI_INT, 0, MPI_COMM_WORLD);

    /* print all the information*/
   if (my_rank == 0)
   {
     /* Time end */
     time_end = MPI_Wtime();
     printf("running on %d processors\n", comm_sz);
     printf("elapsed time = %.6Le seconds\n", time_end - time_start);
     if (flag == 'I')
       print_matrix(n, c);
   }

  /* free all the data*/
  free(B);  free(c);
  free(sendCount);  free(displs);
  free(localA);  free(localc);

  /*close MPI*/
  MPI_Finalize();
  return 0;
}
