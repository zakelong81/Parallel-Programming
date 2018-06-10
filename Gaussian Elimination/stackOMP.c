/*
 * @file sgauss.c
 * @author Long Nguyen <lnguyen34@mail.csuchico.edu>
 * @date 4/22/2018
 *
 * @section DESCRIPTION
 * This program is the 4th programming assignment for CSCI 551
 * (Numerical Methods and Parallel Programming) at CSU Chico
 *
 * The purpose of the program is to write Serial Gaussian Elimination
 */


#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include <sys/resource.h>
/*
 * @name print_function
 *
 * @section DESCRIPTION
 * Prints a given matrix and vector next to each other
 *
 * @param int n - The size of the matrix and vector
 * @param matrix[n][n] - The matrix to be printed
 * @param array[n] - The vector to be printed
 *
 * @return None
 */
void print_function(int n, double A[n][n], double b[n])
{
  int i,j;
  printf("\n");
  for(i = 0; i < n; i++)
  {
     for(j = 0; j <n; j++)
     {
       printf("%.10e ", A[i][j]);
     }
     printf("%.10e", b[i]);
     printf("\n");
   }
}

/*
 * @name print_function_vector
 *
 * @section DESCRIPTION
 * Prints a last interation x vector
 *
 * @param int n - The size of the vector
 * @param array[n] - The vector to be printed
 *
 * @return None
 */
void print_function_vector(int n, double x[n])
{
  int i;
  for(i = 0; i < n; i++)
  {
     printf("%.10e ", x[i]);
  }
  printf("\n");
}

/*
 * @name forward_elimination
 *
 * @section DESCRIPTION
 * do forward elimination
 *
 * @param int n - The size of the matrix and vector
 * @param matrix[n][n] - The matrix input
 * @param array[n] - The vector input
 * @param array[n] - The vector that store answer
 *
 * @return None
 */
void forward_elimination(int n, double A[n][n], double b[n])
{
    int swap, swaprow;
    int i,j,k;
    double absvalue;
    double tempA, tempb;

    for (i = 0; i < n-1; i++)
    {
      absvalue = fabs(A[i][i]);
      swap = 1;

      for (j = i+1; j < n; j++)
      {
        if (fabs(A[j][i]) > absvalue)
        {
          swap = 0; // set swap flag to true
          swaprow = j;
          absvalue = fabs(A[j][i]);
        }
      }

      // swap row A and row b when flag swap is true
      if (swap == 0)
      {
        for (j = i; j < n; j++)
        {
          tempA = A[i][j];
          A[i][j] = A[swaprow][j];
          A[swaprow][j] = tempA;
        }
          tempb = b[i];
          b[i] = b[swaprow];
          b[swaprow] = tempb;
      }

      #pragma omp parallel for private(j, k) //shared(A, b, i, n)
      for ( j = i+1; j < n; j++)
      {
        double mul = A[j][i] / A[i][i];
        for (k = 0; k < n+1; k++)
        {
          A[j][k] = A[j][k] - mul * A[i][k];
        }
        b[j] = b[j] - b[i] * mul;
      }

    }
}

/*
 * @name back_substitution
 *
 * @section DESCRIPTION
 * do back substitution
 *
 * @param int n - The size of the matrix and vector
 * @param matrix[n][n] - The matrix input
 * @param array[n] - The vector input
 * @param array[n] - The vector that store answer
 *
 * @return None
 */
void back_substitution(int n, double A[n][n], double b[n], double x[n])
{
  int i,j;
  double s;
    for (i = n-1; i >= 0; i--)
    {
        s = b[i];
        #pragma omp parallel for private(j) reduction(-:s) //shared(A, x, i, n)
        for (j=i+1; j<n; j++)
        {
            //#pragma omp critical
            s = s - A[i][j]*x[j];
        }
        //#pragma omp critical
        //#pragma omp barrier
        x[i] = s/A[i][i];
    }

}
/*
 * @name copy_original_matrix
 *
 * @section DESCRIPTION
 * make a copy of original matrix
 *
 * @param int n - The size of the matrix and vector
 * @param matrix[n][n] - The matrix you want to copy
 * @param array[n] - The vector you want to copy
 * @param matrix[n][n] - New copied matrix
 * @param array[n] - New copied vector
 *
 * @return None
 */

void copy_original_matrix(int n, double A[n][n], double b[n], double Ao[n][n], double bo[n])
{
  int i,j;
  for(i = 0; i< n; i++)
  {
      for(j = 0; j<n; j++)
      {
          Ao[i][j] = A[i][j];
      }
      bo[i] = b[i];
  }
}

/*
 * @name euclidean_norm
 *
 * @section DESCRIPTION
 * calculate euclidean norm
 *
 * @param int n - The size of the matrix and vector
 * @param matrix[n][n] - The matrix input original
 * @param array[n] - The vector input original
 * @param array[n] - The vector that store answer for back_substitution
 *
 * @return sqrt of norm
 */
double euclidean_norm(int n, double Ao[n][n], double bo[n], double x[n])
{
  double norm = 0;
  double r[n];
  int i,j;
  // r = Ax - b
  for (i = 0; i < n; i++)
  {
    r[i] = 0;
    for (j = 0; j < n; j++)
    {
      r[i] = r[i]+ Ao[i][j] * x[j];
    }
    r[i] = r[i] - bo[i];
    norm =  pow(r[i],2);
  }
  return sqrt(norm);
}

/*
 * @name print_usage
 *
 * @section DESCRIPTION
 * gather use and print all rusage function
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

  printf("%s ","max resident set: " );
  printf("%f ", (double)usage.ru_maxrss);
  printf("%s\n","kB" );

  printf("minor page faults: %ld\n", usage.ru_minflt);
  printf("major page faults: %ld\n", usage.ru_majflt);

}

/*
 * @name main
 *
 * @section DESCRIPTION
 * create matrix and array then apply back and forward elimination to find x
 *
 * @param None
 *
 * @return 0
 */
int main(int argc, char* argv[])
{
  if (argc != 2)
  {
    return 1;
  }

  int n;
  int p;
  //scanf ("%i",&n);
  n = atoi(argv[1]);

  double A[n][n];
  double b[n];
  double Ao[n][n];
  double bo[n];
  double x[n];
  double input;
  double norm;
  //int procCount;
  double starttime;

  #pragma omp parallel
  {
    p = omp_get_num_threads();
  }

  //procCount = omp_get_num_procs();
  //#pragma omp parallel

  printf("Total Number of Thread: %d\n", p);
  //printf("%d ", procCount);

  // n is less than 11 input the matrix and array
  int i,j;
  if (n <11)
  {
    for(i = 0; i< n; i++)
    {
        for(j = 0; j<n; j++)
        {
            scanf ("%lf",&input);
            A[i][j] = input;
            //Ao[i][j] = A[i][j];
        }
        scanf ("%lf",&input);
        b[i] = input;
        //bo[i] = b[i];
    }
  }
  //n is bigger than 10 random matrix and array
  else if (n > 10)
  {
    srand48(time(NULL));
    for(i = 0; i< n; i++)
    {
      for(j = 0; j<n; j++)
      {
          A[i][j] = drand48() * 1.0e6 - 1.0e6;
          //Ao[i][j] = A[i][j] ;
      }
      b[i] = drand48() * 1.0e6 - 1.0e6;
      //bo[i] =   b[i];
    }
  }
  starttime = omp_get_wtime();
  copy_original_matrix(n,A,b,Ao,bo);
  forward_elimination(n,A,b);
  back_substitution(n,A,b,x);
  norm = euclidean_norm(n,Ao,bo,x);

  print_usage();
  if (n < 11)
  {
    print_function(n,Ao,bo);
    print_function_vector(n,x);
  }
  printf("norm of the residual: %.10e\n", norm );

  double finaltime = omp_get_wtime() - starttime;
  printf("%lf\n", finaltime );
  return 0;
}
