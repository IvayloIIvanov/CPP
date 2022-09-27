#include<iostream>
#include<cstdlib>
using namespace std;


double eps = 1.e-5;
double hx, hy;
double** Unew, ** Uold, ** f, ** Uexact;
int N = 64, M = 64;
double xmin = 0., xmax = 1., ymin = 0., ymax = 1.;

int main() {

	hx = (xmax - xmin) / (N + 1);
	hy = (ymax - ymin) / (M + 1);

	Unew = (double**)malloc((N + 2) * sizeof(double*));
	Uold = (double**)malloc((N + 2) * sizeof(double*));
	Uexact = (double**)malloc((N + 2) * sizeof(double*));
	f = (double**)malloc((N + 2) * sizeof(double*));

	for (int ii = 0; ii < N + 2; ii++)
	{
		Unew[ii] = (double*)malloc((M + 2) * sizeof(double));
		Uold[ii] = (double*)malloc((M + 2) * sizeof(double));
		Uexact[ii] = (double*)malloc((M + 2) * sizeof(double));
		f[ii] = (double*)malloc((M + 2) * sizeof(double));
	}

	for (int i = 0; i < (N + 2); i++)
	{
		double x = xmin + i * hx;
		for (int j = 0; j < (M + 2); j++)
		{
			double y = ymin + j * hy;
			f[i][j] = 2 * x * x * x - 6 * x * y * (1. - y);
			Uexact[i][j] = y * (1. - y) * x * x * x;
		}
	}

	for (int i = 0; i < (N + 1); i++)
	{
		for (int j = 0; j < (M + 1); j++)
		{
			Uold[i][j] = f[i][j];
		}
	}

	for (int i = 0, j = 0; j < (M + 2); j++)
	{
		Uold[i][j] = Uexact[i][j];
		Unew[i][j] = Uexact[i][j];
	}

	for (int i = N + 1, j = 0; j < (M + 2); j++)
	{
		Uold[i][j] = Uexact[i][j];
		Unew[i][j] = Uexact[i][j];
	}

	for (int i = 0, j = 0; i < (N + 2); i++)
	{
		Uold[i][j] = Uexact[i][j];
		Unew[i][j] = Uexact[i][j];
	}

	for (int i = 0, j = M + 1; i < (N + 2); i++)
	{
		Uold[i][j] = Uexact[i][j];
		Unew[i][j] = Uexact[i][j];
	}


	double diff = 0.;

	do {
		for (int i = 1; i < (N + 1); i++) {
			for (int j = 1; j < (M + 1); j++) {
				Unew[i][j] = (hx * hx * hy * hy) / (2 * (hx * hx + hy * hy)) * (f[i][j] + ((Uold[i - 1][j] + Uold[i + 1][j]) / hx * hx) + ((Uold[i][j + 1] + Uold[i][j - 1]) / hy * hy));
			}
		}
		diff = 0.;
		for (int i = 1; i < (N + 1); i++)
		{
			for (int j = 1; j < (M + 1); j++)
			{
				diff += (Unew[i][j] - Uold[i][j]) * (Unew[i][j] - Uold[i][j]);
			}
		}

		for (int i = 1; i < (N + 1); i++)
		{
			for (int j = 1; j < (M + 1); j++)
			{
				Uold[i][j] = Unew[i][j];
			}
		}
	} while (diff > eps* eps);

	for (int i = 0; i < (N + 2); i++) {
		double x = xmin + i * hx;
		for (int j = 0; j < (M + 2); j++) {
			double y = ymin + j * hy;
			Uexact[i][j] = y * (1. - y) * x * x * x;
			cout << "x = " << x << "  " << "y = " << y << "  " << "Uold = " << Uold[i][j] << "  " << "Unew = " << Unew[i][j] << "  " << "Uexact = " << Uexact[i][j] << "  " << "Uexact - Unew = " << Uexact[i][j] - Unew[i][j] << endl;
		}
	}


	for (int ii = 0; ii < (N + 2); ii++)
	{
		free(Unew[ii]);
		free(Uold[ii]);
		free(Uexact[ii]);
		free(f[ii]);
	}

	free(Unew);
	free(Uold);
	free(Uexact);
	free(f);

}
