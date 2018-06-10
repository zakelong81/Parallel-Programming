/*
 * @file p2.c
 * @author Long Nguyen <lnguyen34@mail.csuchico.edu>
 * @date 3/4/2018
 *
 * @section DESCRIPTION
 * This program is the second programming assignment for CSCI 551
 * (Numerical Methods and Parallel Programming) at CSU Chico
 *
 * The purpose of the program is to write Jacobi Iteration and Gauss-Seidel Iteration
 * to compare the speed and number of interation
 */


#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

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
void print_function(int n, float A[n][n], float b[n])
{
  printf("\n");
  for(int i = 0; i < n; i++)
  {
     for(int j = 0; j <n; j++)
       printf("%.10e ", A[i][j]);
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
void print_function_vector(int n, float x[n])
{
  for(int i = 0; i < n; i++)
     printf("%.10e ", x[i]);
  printf("\n");
}

/*
 * @name check_diagonal
 *
 * @section DESCRIPTION
 * check if a matrix is diagonal
 *
 * @param int n - The size of the matrix
 * @param matrix[n][n] - The matrix to be checked
 *
 * @return None
 */
void check_diagonal(int n, float A[n][n])
{
  float check;
  float jsum;
  int flag = 0;
  for(int i = 0; i< n; i++)
  {
    if (flag == 0)
    {
        check = A[i][i];
        jsum = 0;
        for(int j = 0; j<n; j++)
        {
          if (j != i)
            jsum = jsum + A[i][j];
          if(check < jsum)
              flag = 1;
        }
    }
    else
      break;
  }

  if (flag ==  0)
  {
    printf("\n");
    printf("%s\n", "the system is diagonally dominant");
  }
  else
  {
    printf("\n");
    printf("%s\n", "the system is not diagonally dominant");
  }
}


/*
 * @name jacobi
 *
 * @section DESCRIPTION
 * implementation of Jacobi Iteration
 *
 * @param int n - The size of the matrix and vector
 * @param matrix[n][n] - The matrix input
 * @param array[n] - The vector input
 * @param array[n] - The vector that store last Iteration
 *
 * @return None
 */
void jacobi(int n, float A[n][n], float b[n], float x[n])
{
  printf("\n");
  printf("%s\n", "Jacobi Iteration algorithm is running" );
  float previousX[n];
  float presentX[n];
  float arae;
  int flag = 1;

  //set all previous to 0 and present to 1
  for ( int i= 0; i < n; i++ )
  {
    previousX[i]=0;
    presentX[i]= 1;
  }

  //interation 100 times
  for ( int iteration = 1; iteration < 101; iteration++ )
  {
      //set current = previous
      for ( int i= 0; i < n; i++ )
        previousX[i] = presentX[i];

      //Jacobi implementation
      for ( int i= 0; i < n; i++ )
      {
        presentX[i]  = b[i];
        for (int j= 0; j < n; j++ )
        {
          if ( i != j)
            presentX[i] = presentX[i] - A[i][j]*previousX[j];
        }
        presentX[i] = presentX[i] / A[i][i];
      }

      //caculate the flag
      for ( int i =0; i < n; i++ )
          arae = fabs((presentX[i]-previousX[i])/ presentX[i]);

      //check for flag if ARAR < 1.0e-3
      for ( int i =0; i < n; i++ )
      {
        if (arae < 1.0e-3)
          flag = 0;
      }

      //if flag is set terminate the loop and print out # of iteration
      if (flag == 0)
      {
        //record the final iteration
        for (int i=0; i< n; i++)
            x[i] = presentX[i];

         printf("%s\n","Convergence" );
         printf("%s", "Number of Iteration: ");
         printf("%i\n",  iteration );
         iteration=1000;
      }
  }

 // if the flag is not set print out not Convergent
  if (flag == 1)
    printf("%s\n","Not Convergence" );

}

/*
 * @name gauss_seidel
 *
 * @section DESCRIPTION
 * implementation of Gauss Seidel Iteration
 *
 * @param int n - The size of the matrix and vector
 * @param matrix[n][n] - The matrix input
 * @param array[n] - The vector input
 * @param array[n] - The vector that store last Iteration
 *
 * @return None
 */
void gauss_seidel(int n, float A[n][n], float b[n], float x[n])
{
  printf("\n");
  printf("%s\n", "Gauss Seidel Iteration algorithm is running" );
  float previousX[n];
  float presentX[n];
  float arae;
  int flag =1;

  for ( int i= 0; i < n; i++ )
  {
    previousX[i]=0;
    presentX[i]= 1;
  }

  for ( int iteration = 1; iteration <= 100; iteration++ )
  {
      //set current = previous
      for ( int i= 0; i < n; i++ )
          previousX[i] = presentX[i];

      //Gauss Seidel implementation
      for ( int i= 0; i < n; i++ )
      {
        presentX[i] = b[i];
        for (int j= 0; j < n; j++ )
        {
          if ( i != j)
            presentX[i] = presentX[i] - A[i][j]*presentX[j];
        }
        presentX[i] = presentX[i] / A[i][i];
      }

      //caculate the flag
      for ( int i =0; i < n; i++ )
          arae = fabs((presentX[i]-previousX[i])/ presentX[i]);

      //check for flag if ARAR < 1.0e-3
      for ( int i =0; i < n; i++ )
      {
        if (arae < 1.0e-3)
          flag = 0;
      }

      //if flag is set terminate the loop and print out # of iteration
      if (flag == 0)
      {
        //record the final iteration
        for (int i=0; i< n; i++)
            x[i] = presentX[i];

         printf("%s\n","Convergence" );
         printf("%s", "Number of Iteration: ");
         printf("%i\n",  iteration );
         iteration=1000;
      }
    }

    // if the flag is not set print out not Convergence
    if (flag == 1)
      printf("%s\n","Not Convergence" );
}

/*
 * @name main
 *
 * @section DESCRIPTION
 * create matrix and array call Jacobi and Gauss Seidel then print them out
 *
 * @param None
 *
 * @return 0
 */
int main()
{
  int n;
  //seed srand48
  srand48(time(NULL));
  scanf ("%i",&n);

  float A[n][n];
  float b[n];
  float x1[n];
  float x2[n];
  float input;

  // n is less than 11 input the matrix and array
  if (n <11)
  {
    for(int i = 0; i< n; i++)
    {
        for(int j = 0; j<n; j++)
        {
            scanf ("%f",&input);
            A[i][j] = input;
        }
        scanf ("%f",&input);
        b[i] = input;
    }
    print_function(n,A,b);
  }
  //n is bigger than 10 random matrix and array
  else if (n > 10)
  {
    for(int i = 0; i< n; i++)
    {
      for(int j = 0; j<n; j++)
      {
          if (i==j)
            A[i][j] = drand48() * 1.0e6 - 1.0e6;
          else
            A[i][j] = drand48() * 1.0e3 - 1.0e3;
      }
      b[i] = drand48();
    }
  }

  check_diagonal(n,A);
  jacobi(n, A, b, x1);
  if (n < 11)
    print_function_vector(n,x1);

  gauss_seidel(n, A, b, x2);
  if (n < 11)
    print_function_vector(n,x2);

  return 0;
}
