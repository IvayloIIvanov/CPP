#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <complex>
#include <Lapacke/lapacke.h>

using namespace std;
int main()
{
	int indx, k, n, steps, cnt, info;

	complex<double>* Psi[n + 1];
	complex<double>* PsiPrime[n];
	complex<double>* D[n];
	complex<double>* DL[n - 1];
	complex<double>* DU[n - 1];
	complex<double>* Pot[n + 1];
	double* TrapCoef[n];

	char filename[32];

	double xmin, xmax, I0, T0, x, dx, t, dt, tout, sig, norm, xpmin, xpmax, Factor;
	double Ekin, Epot, V0, k0, sig;
	double PI = 3.14159265359;
	xmin = -25.;
	xmax = 25.;
	xpmin = 4.;
	xpmax = 6.;
	I0 = -5.;
	sig = 0.5;
	const complex<double> One(1.0, 0.0), Zero(0.0, 0.0), Im(0.0, 1.0);
	cout << "write grid resolution - dx=" << endl;
	cin >> dx;
	norm = (dx * dx) / (1 + V0 * 2 * dx * dx);
	cout << "write time integration step - dt=" << endl;
	cin >> dt;
	cout << "write Kinetic Energy" << endl;
	cin << Ekin;
	k0 = sqrt(2 * Ekin - 1.);
	cout << "write potential well - V0=" << endl;
	cin >> V0;
	cout << "write integration period - T0=" << endl;
	cin >> T0;
	cout << "write out period - tout=" << endl;
	cin >> tout;
	n = (int)(xmax - xmin / dx) - 1;
	steps = tout / dt;
	Factor = (0.5 * dt) / (dx * dx);
	for (int i = 2; i < n; i++) {
		TrapCoef[i] = 2.;
		TrapCoef[1] = 1.;
		TrapCoef[n] = 1.;
	};
	for (int i = 0; i < n + 1; i++)
	{
		double PI = 3.14159265359;
		x = xmin + i * dx;
		Psi[i] = (exp(-(x - I0) * (x - I0) / (4 * sig * sig)) / (sqrt(sqrt(2 * PI) * sig))) * complex(cos(k0 * x - I0), sin(k0 * (x - I0)));
		if (x >= xpmin && x <= xpmax)
		{
			Pot[i] = complex(V0, 0.);
		}
		else
		{
			Pot[i] = complex(0., 0.);
		}

	}
	Psi[0] = 0.;
	Psi[n + 1] = 0.;
	for (int i = 0; i < n + 1; i++)
	{
		complex<double> mycomplex;
		mycomplex = Psi[i] * std::conj(Psi[i]);
		norm = sqrt(dx * std::real(mycomplex);
		Psi[i] = Psi[i] / norm;
	}
	t = 0.;
	cnt = 0.;
	do {
		for (int i = 1; i < n + 1; i++) {
			Psi[i] = (One - Im * Factor - Im * 0.5 * dt * Pot[i]) * Psi[i] + 0.5 * Im * Factor * (Psi[i - 1] + Psi[i + 1]);
			D[i] = One + Im * Factor * (Psi[i - 1] + Psi[i + 1])
		};
		for (int i = 1; i < n; i++) {
			DL[i] = -0.5 * Im * Factor;
			DU[i] = -0.5 * Im * Factor;
		}
		;
		//call function?
		for (int i = 1; i < n + 1; i++)
		{
			PsiPrime[i] = 0.5 * ((Psi[i + 1] - Psi[i - 1]) / dx);
			Ekin = 0.25 * dx * (TrapCoef[i] * PsiPrime[i] * std::conj(PsiPrime[i]));
			Epot = 0.5 * dx * (TrapCoef[i] * std::conj(Psi[i]) * Pot[i] * Psi[i]);
			norm = 0.5 * dx * (TrapCoef[i] * std::conj(Psi[i]) * Psi[i]);
		}
		;
	} while (t < T0);
	if (cnt % steps == 0) {
		for (i = 1; i < n + 1; i++) {
			cout << "Ekin = " << Ekin << "   " << "Epot = " << Epot << "   " << "Ekin+Epot = " << Ekin + Epot << "   " << "norm = " << norm << endl;
		};
		for (i = 0; i < n + 1)
		{
			double PI = 3.14159265359;
			x = xmin + i * dx;
			Psi[i] = (exp(-(x - I0) * (x - I0) / (4 * sig * sig)) / (sqrt(sqrt(2 * PI) * sig))) * complex(cos(k0 * x - I0), sin(k0 * (x - I0)));
			if (x >= xpmin && x <= xpmax)
			{
				Pot[i] = complex(V0, 0.);
			}
			else
			{
				Pot[i] = complex(0., 0.);
			}
			cout << "real Psi = " << std::real(Psi[i]) << "   " << "imaginary Psi = " << std::imag(Psi[i]) << "   " << "abs Psi = " << abs(Psi[i]) << "   " << "real Pot = " << std::real(Pot[i]) << endl;
		}
	}
	else
	{
		cout << "Stop" << endl;
	}

}
return 0;
}

