/*
 * @file tmin.c
 * @author Long Nguyen <lnguyen34@mail.csuchico.edu>
 * @date 3/25/2018
 *
 * @section DESCRIPTION
 * This program is the third programming assignment for CSCI 551
 * (Numerical Methods and Parallel Programming) at CSU Chico
 *
 * The purpose of the program is to find the min # of trap
 */

#include <stdio.h>
#include <stdlib.h>
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
	return ((-3)*cosl(x/6) + sqrtl( powl( ((5*sinl(x/9))+(8*sinl(x/2))),4) ) + 1 );
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
long double Trap(double a, double b, double num_trap)
{
	long double x, h, estimate;
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
* Find the minimum number of trap needed so with the error less than 0.5e-14
* @param  None
* @return 0
*/
int main()
{
	//Given value
	long double abs_err = 0.5*pow(10,-14);
	long double true_v = 4754.019228858818136649863;
	//a = uperbound , b = underbound, num_trap = number of trap, prev_trap = previous number of trap
  long double a, b, num_trap, prev_trap;
	int adder = 1000000;
  long double arte, curr_v;

	printf("Insert a, b and t:\n");
  scanf("%Lf%Lf%Lf", &a, &b, &num_trap);
  prev_trap = num_trap;

	int check = 0;
  do
	{
    if (check == 0)
		{
      prev_trap = num_trap;
      num_trap = num_trap + adder;
    }
		else
      check = 1;

    curr_v = Trap(a, b, num_trap);
    arte  = fabs( (true_v - curr_v) / true_v);
    // printf("Number of Trap: %.0Lf Integrate: %.13Le ARTE: %0.19Le\n", num_trap, curr_v, arte );
    if (arte < abs_err)
		{
			//printf("Reduce adder by 10\n");
			if (adder > 1)
				arte = 1; // reset arte so it will continute
			num_trap = prev_trap;
			adder = adder / 10;
    }
  } while (arte > abs_err);

	printf("**Final result**\n");
	printf("Number of Trapazoids: %.0Lf\n", num_trap);
	printf("Integrate: %.13Le\n", curr_v);
	printf("ARTE: %0.19Le\n", arte );

  return 0;
}
