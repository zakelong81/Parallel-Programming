/*
 * @file mm.c
 * @author Long Nguyen <lnguyen34@mail.csuchico.edu>
 * @date 4/29/2018
 *
 * @section DESCRIPTION
 * This program is the 1st programming assignment for CSCI 551
 * (Numerical Methods and Parallel Programming) at CSU Chico
 *
 * The purpose of the program is to vectorlize matrix multiplication
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <sys/resource.h>

/*
 * @name print_function
 *
 * @section DESCRIPTION
 * this function will take matrix A and B and push it in to matrix c
 *
 * @param int n - The size of the matrix and vector
 * @param matrix[n][n] - The matrix A
 * @param matrix[n][n] - The matrix B
 * @param array[n] - The vector c
 *
 * @return None
 */

void matrix_multiplication(int n, float A[n][n], float B[n][n], float c[n][n])
{
 float tmp = 0.0;
 for(int i = 0; i<n; i++)
 {
    for(int j = 0; j<n; j++)
    {
        tmp = 0;
        for(int k = 0; k < n; k++)
        {
           tmp += A[i][k] * B[k][j];
        }
        //reduce dependency if place it here instead of inside for k loop
        c[i][j] = tmp;
    }
 }
}

/*
 * @name print_function
 *
 * @section DESCRIPTION
 * Prints a given matrix
 *
 * @param int n - The size of the matrix and vector
 * @param matrix[n][n] - The matrix to be printed
 *
 * @return None
 */
void print_matrix(int n, float m[n][n])
{
  for(int i = 0; i< n; i++)
  {
     for(int j = 0; j<n; j++)
     {
       printf("%.0f ", m[i][j]);
     }
     printf("\n");
  }
}

/*
 * @name print_function
 *
 * @section DESCRIPTION
 * this fucntion will print out usage
 *
 * @param None
 *
 * @return None
 */
void print_usage()
{
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);

  printf("%s ","user time: " );
  printf("%f ", (double)usage.ru_utime.tv_sec+((double)usage.ru_utime.tv_usec/(double)1000000));
  printf("%s\n","sec" );

  printf("%s ","system time: " );
  printf("%f ", (double)usage.ru_stime.tv_sec+((double)usage.ru_stime.tv_usec/(double)1000000));
  printf("%s\n","sec" );

  printf("%s ","sum of user and system time: " );
  printf("%f ",
  ((double)usage.ru_utime.tv_sec+((double)usage.ru_utime.tv_usec/(double)1000000))+
  ((double)usage.ru_stime.tv_sec+((double)usage.ru_stime.tv_usec/(double)1000000)));
  printf("%s\n","sec" );

  printf("%s ","max resident set: " );
  printf("%f ", (double)usage.ru_maxrss);
  printf("%s\n","kB" );
}

/*
 * @name main
 *
 * @section DESCRIPTION
 * create 2 matrix and then do matrix multiplication
 *
 * @param None
 *
 * @return 0
 */
int main()
{
    int n;
    char flag;
    float input;

    scanf ("%c",&flag);
    scanf ("%i",&n);

    // data aligning
    float A[n][n], B[n][n], c[n][n] __attribute__((aligned(32)));

    if (flag != 'R' & flag != 'I')
    {
      printf("%s\n", "Not valid flag R = random, I = Input" );
      exit(1);
    }


    //if R create a random matrix
    if (flag == 'R')
    {
        srand(time(NULL));
        for(int i = 0; i< n; i++)
        {
            for(int j = 0; j<n; j++)
            {
                A[i][j] = rand() % 101 + (-50);
                B[i][j] = rand() % 101 + (-50);
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
             scanf ("%f",&input);
             A[i][j] = input;
           }
        }
        for(int i = 0; i< n; i++)
        {
            for(int j = 0; j<n; j++)
            {
                scanf ("%f",&input);
                B[i][j] = input;
            }
        }
     }

     if (flag == 'R')
       matrix_multiplication(n, A, B, c);
     else if (flag == 'I')
     {
       matrix_multiplication(n, A, B, c);
       print_matrix(n, c);
     }

     print_usage();
     return 0;
}


// float *A[n];
// for (int i=0; i<n; i++)
//     A[i] = (float *)malloc(n * sizeof(int));
//
// float *B[n];
// for (int i=0; i<n; i++)
//     B[i] = (float *)malloc(n * sizeof(int));
//
// float *c[n];
// for (int i=0; i<n; i++)
//     c[i] = (float *)malloc(n * sizeof(int));
