#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;

/*Creating two void functions for calculating the right side of the system and
another one for the Runge - Kuta method*/

void dervis(double x, double* y, double* dydx, double* m)
{


	double* force, * Ri, * Rj;
	int i, j, n;
	force = (double*)malloc(3 * sizeof(m));
	Ri = (double*)malloc(3 * sizeof(double));
	Rj = (double*)malloc(3 * sizeof(double));

	double gm = 0.0172020985;

	// Get the number of bodies from mass array size //

	n = sizeof(m);

	for (i = 0; i < n; i++)
	{
		force[i * 3] = 0.0;
		force[i * 3 + 1] = 0.0;
		force[i * 3 + 2] = 0.0;

		Ri[0] = y[i * 6];
		Ri[1] = y[i * 6 + 1];
		Ri[2] = y[i * 6 + 2];

		for (j = 0; j < n; j++)
		{
			if (i != j)
			{
				Rj[0] = y[j * 6];
				Rj[1] = y[j * 6 + 1];
				Rj[2] = y[j * 6 + 2];

				force[i * 3] = force[i * 3] + gm * gm * m[i] * m[j] * (Rj[0] - Ri[0]) / (pow((Ri[0] - Rj[0]) * (Ri[0] - Rj[0]) + (Ri[1] - Ri[1]) * (Ri[1] - Rj[1]) + (Ri[2] - Rj[2]) * (Ri[2] - Rj[2]), 1.5));
				force[i*3+1]=force[i*3+1] + gm * gm * m[i] * m[j] * (Rj[1] - Ri[1]) / (pow((Ri[0] - Rj[0]) * (Ri[0] - Rj[0]) + (Ri[1] - Ri[1]) * (Ri[1] - Rj[1]) + (Ri[2] - Rj[2]) * (Ri[2] - Rj[2]), 1.5));
				force[i*3+2]=force[i*3+2] + gm * gm * m[i] * m[j] * (Rj[2] - Ri[2]) / (pow((Ri[0] - Rj[0]) * (Ri[0] - Rj[0]) + (Ri[1] - Ri[1]) * (Ri[1] - Rj[1]) + (Ri[2] - Rj[2]) * (Ri[2] - Rj[2]), 1.5));
			}
		}
	}

	//Calculate the derivatives//

	for (i = 0; i < n; i++)
	{
		dydx[i * 6] = y[i * 6 + 3];
		dydx[i * 6 + 1] = y[i * 6 + 4];
		dydx[i * 6 + 2] = y[i * 6 + 5];

		dydx[i * 6 + 3] = force[i * 3]/m[i];
		dydx[i * 6 + 4] = force[i * 3 + 1]/m[i];
		dydx[i * 6 + 5] = force[i * 3 + 2]/m[i];
	}

}

void rk4(double* y, double* dydx, double x, double h, double* yout, void*, double* m)
{
	void dervis(double x, double* y, double* dydx, double* m);

	int ndum;
	double h6, hh, xh;

	double* dym, * dyt, * yt;

	dym = (double*)malloc(3 * sizeof(y));
	dyt = (double*)malloc(3 * sizeof(y));
	yt = (double*)malloc(3 * sizeof(y));

	ndum = sizeof(y);

	hh = h * 0.5;
	h6 = h / 6.0;
	xh = x + hh;

	int i, n;

	for (i = 0; i < n; i++)
	{
		yt[i] = y[i] + hh * dydx[i];
	}

	dervis(xh, yt, dyt, m);
	for (i = 0; i < n; i++)
	{
		yt[i] = y[i] + hh * dyt[i];
	}

	dervis(xh, yt, dym, m);
	for (i = 0; i < n; i++)
	{
		yt[i] = y[i] + h * dym[i];
		dym[i] = dyt[i] + dym[i];
	}

	dervis(x + h, yt, dyt, m);

	for (i = 0; i < n; i++)
	{
		yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2 * dym[i]);
	}
}

//Main program integrationg the equation of motion of n celestial obejects//

int main()
{
	double h, x, period, wrtout, maxit;
	double* y, * ynew, * dydx, * m;
	int i, n, j;

	char frmt (20);
	

	//Read the total integration period//

	cout << "Enter period:";
	cin >> period;

	//Read the integration step//

	cout << "Enter integration step:";
	cin >> h;
	if (h <= 0)
	{
		cout << "wrong h" << endl;
	}

	cin >> wrtout;

	//Read the number of objects//
	cout << "Enter number of objects";
	cin >> n;

	//Allocating required memory for working arrays based on - n//

	y = (double*)malloc(6 * n * sizeof(double));
	ynew = (double*)malloc(6 * n * sizeof(double));
	dydx = (double*)malloc(6 * n * sizeof(double));
	m = (double*)malloc(6 * n * sizeof(double));

	//Read inital positions, velocity and mass for each object//

	for (i = 0; i < n; i++)
	{
		cin >> y[i * 6];
		cin >> y[i * 6 + 1];
		cin >> y[i * 6 + 2];
		cin >> y[i * 6 + 3];
		cin >> y[i * 6 + 4];
		cin >> y[i * 6 + 5];
		cin >> m[i];
	}

	//Calculating iteration step//

	maxit = period / h;

	for (i = 0; i < maxit; i++)
	{
		x = x + h;
		dervis(x, y, dydx, m);
		rk4(y, dydx, x, h, ynew, dervis, m);

		for (i = 0; i < n; i++)
		{
			y[i] = ynew[i];
		}

		double c = i * h;
		if (fmod(c, wrtout) == 0)
		{
			cout << "x = " << x;
			for (j = 0; j < n; j++){
			cout << " y = " << y[j] << endl;
			}
		}
	}

	//Dealocating the arrays//

	free(y);
	free(ynew);
	free(dydx);
	free(m);

	return 0;
}