/*
 * @file tmin.c
 * @author Long Nguyen <lnguyen34@mail.csuchico.edu>
 * @date 3/25/2018
 *
 * @section DESCRIPTION
 * This program is the third programming assignment for CSCI 551
 * (Numerical Methods and Parallel Programming) at CSU Chico
 *
 * The purpose of the program is to use MPI to run parallel trapezoids rule
 */
#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <math.h>

/*
* @name funtion
* @section DESCRIPTION
* take a current x and calculate
* @param  x
* @return function
*/
long double funtion(double x)
{
 return ( (-3)*cosl(x/6) + sqrtl( powl( ((5*sinl(x/9))+(8*sinl(x/2))),4) ) + 1 );
}

/*
* @name Trap
* @section DESCRIPTION
* Find area under the curve with upper and lower bound
* @param  a lower bound
* @param  b upper bound
* @param  num_trap # of trap
* @param  h lens
* @return trap under the curve value estimate
*/
long double Trap(long double a, long double b, long double num_trap)
{
 long double x, estimate,h;
  h = (b - a) /  num_trap;
 estimate = (funtion(a) + funtion(b)) / 2.0;
 for (int i = 1; i < num_trap; i++)
 {
   x = a + i * h;
   estimate = estimate + funtion(x);
 }
 estimate = estimate * h;
 return estimate;
}

/*
* @name main
* @section DESCRIPTION
* using MPI to parallel trapezoids rule

* @return 0
*/
int main(void)
{
   int my_rank, comm_sz, n, local_n;
   long double a, b, local_a, local_b;
   long double local_integrate, total_integrate;
   double time_start, time_end, total_time;
   long double arte;
   long double true_v = 4754.019228858818136649863;
   long double criteria = 0.5*pow(10, -14);

   /* Let the system do what it needs to start up MPI */
   MPI_Init(NULL, NULL);
   /* Get my process rank */
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   /* Find out how many processes are being used */
   MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);


   if (my_rank == 0)
   {
      printf("Enter a, b, and n\n");
      scanf("%Lf %Lf %d", &a, &b, &n);
      printf("\nRunning on %d cores.\n", comm_sz);
   }

   //giving out information
   MPI_Bcast(&a, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Bcast(&b, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

   //start timer
   if (my_rank == 0)
   {
    time_start = MPI_Wtime();
   }

   local_n = n/comm_sz;  /* So is the number of trapezoids  */

   local_a = a + my_rank*local_n*((b-a)/n);
   local_b = local_a + local_n*((b-a)/n);

   local_integrate = Trap(local_a, local_b, local_n);

   /* Add up the integrals calculated by each process */
   MPI_Reduce(&local_integrate, &total_integrate, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

   //stop the time
   if (my_rank == 0)
   {
    time_end = MPI_Wtime();
    total_time = time_end - time_start;
   }

   /* Print the result */
   if (my_rank == 0)
   {
     printf("Elapsed time = %.6e seconds\n", total_time);
     printf("With n = %d trapezoids, our estimate\n", n);
     printf("of the integral from %Lf to %Lf = %.13Le\n", a, b, total_integrate);
     printf("true value = %.19Le\n" ,true_v);
     arte = fabs((true_v - total_integrate) / true_v);
     printf("absolute relative true error = %.19Le\n" , arte);
     if (arte < criteria)
       printf(" is less than criteria = %.19Le\n",criteria);
     else
       printf(" is NOT less than criteria = %.19Le\n",criteria);
   }

   /* Shut down MPI */
   MPI_Finalize();

   return 0;
}
