#include <iostream>
#include <cmath>
#include "fftw3.h"

using namespace std;
double calc_f(double x)
{
	return -x * (x + 3) * exp(x);
}
int main()
{
	double xmin = 0.0;
	double xmax = 1.0;

	int N = 128;
	double* u, * Uf;
	double* f, * Ff;
	fftw_plan A;
	fftw_plan B;

	double pi = 3.141592356589793;
	double delta = (xmax - xmin) / (N - 1) ;

	cout << "Hello world!" << endl;

	u= (double*)fftw_malloc(sizeof(double) * N);
	f= (double*)fftw_malloc(sizeof(double) * N);
	Uf= (double*)fftw_malloc(sizeof(double) * N);
	Ff= (double*)fftw_malloc(sizeof(double) * N);

	for (int i = 0; i < N; i++)
	{
		double xn = xmin + i * delta;
		f[i] = calc_f(xn);
	}


	A = fftw_plan_r2r_1d(N, f, Ff, FFTW_RODFT00, FFTW_ESTIMATE);
	B = fftw_plan_r2r_1d(N, Uf, u, FFTW_RODFT00, FFTW_ESTIMATE);
	fftw_execute(A);

	for (int j = 0; j < N; j++)
	{
		Uf[j] = Ff[j] * delta * delta / (2 * (1. - cos(pi * (j + 1.) / (N - 1))));
	}

	for (int k = 0; k < N; k++)
	{
		u[k] = u[k] / N;
	}

	fftw_execute(B);

	for (int i = 0; i < N; i++)
	{
		double xn = xmin + i * delta;
		cout << "x= " << xn << "u[i]= " << u[i] << endl;
	}

	 fftw_free(u);
	  fftw_free(Uf);
	 fftw_free(f);
	 fftw_free(Ff);

	fftw_destroy_plan(A);
	fftw_destroy_plan(B);


	return 0;
}