#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <fstream>
// the code is slow - pbc was implemented without halo - cache misses

// ������: �������� ����������� ��������� ���, ���� �� ������������ �� 
// ������� ��������� � ��������� ������ �� ����� �� ����� 
// � ������� ������������ ��������

// �������� ����� [0,1] - ��������� �������������
inline double uniform()
{
	return double(rand()) / RAND_MAX;
}

// �������� ���� ����� [low,up] - ��������� �������������
inline int uniform_int(int low, int up)
{
	return rand() % (up - low + 1) + low;
}

class grid
{
public:
	const int up = 1;
	const int down = -1;
	int* m_data;
	int m_Nx, m_Ny;

	grid(int nx, int ny) : m_Nx(nx), m_Ny(ny)
	{
		int n = nx * ny;
		m_data = new int[n];
	};

	virtual ~grid()
	{
		delete[] m_data;
	};

	int& get(int row, int col)
	{
		auto element_index = row * m_Ny + col;
		return m_data[element_index];
	};

	int& get_pbc(int row, int col)
	{
		return get(row % m_Nx, col % m_Ny);
	};

	void init(int val)
	{
		// slower but more intuitive
		//     for (size_t i = 0; i < m_Nx; i++)
		//     {
		//         for (size_t j = 0; j < m_Ny; j++)
		//         {
		//             get(i,j) = val;
		//         }
		//     }

		// faster
		memset(m_data, val, m_Nx * m_Ny);
	}

	void init(int state, double probability)
	{
		for (size_t i = 0; i < m_Nx; i++)
		{
			for (size_t j = 0; j < m_Ny; j++)
			{
				if (uniform() < probability)
				{
					if (state == up)
					{
						get(i, j) = up;
					}
					else
					{
						get(i, j) = down;
					}
				}
				else
				{
					if (state == up)
					{
						get(i, j) = down;
					}
					else
					{
						get(i, j) = up;
					}
				}
			}
		}
	}
};

// eq. 13.68
double calculate_grid_E(grid& sigma, double J = 1.0)
{
	double E = 0.;

	for (size_t i = 0; i < sigma.m_Nx; i++)
	{
		for (size_t j = 0; j < sigma.m_Ny; j++)
		{
			float local_sum = sigma.get_pbc(i + 1, j) + sigma.get_pbc(i - 1, j);
			local_sum += sigma.get_pbc(i, j + 1) + sigma.get_pbc(i, j - 1);

			E += J * sigma.get(i, j) * local_sum;

		}
	}

	return E / (sigma.m_Nx * sigma.m_Ny);
}

// eq. 13.69
double calculate_grid_M(grid& sigma, double mu = 1.0)
{
	double M = 0.;

	for (size_t i = 0; i < sigma.m_Nx; i++)
	{
		for (size_t j = 0; j < sigma.m_Ny; j++)
		{

			M += mu * sigma.get(i, j);

		}
	}

	return M / (sigma.m_Nx * sigma.m_Ny);
}

////////////////////////////////////////////////
///                 MAIN                     ///
////////////////////////////////////////////////

int main()
{
	int N;         // ������ �� ��������� , grid dim.
	double T;      // �����������, temperature
	int NiterTerm; // ���� ������ �� ������������� �� ���������, number of the quilibration steps
	int Niter;     // ���� ����������� ������, number of the production sreps
	double Tstart, Tstop, dT;
	double Jconst = 1.;

	std::cout << "�������� ������ �� ���������: ";
	std::cin >> N;
	std::cout << "�������� ������� �����������: ";
	std::cin >> Tstart;
	std::cout << "�������� ������ �����������: ";
	std::cin >> Tstop;
	std::cout << "�������� ������ �� ��������� �� �����������: ";
	std::cin >> dT;
	std::cout << "�������� ���� ������ �� ������������� �� ���������: ";
	std::cin >> NiterTerm;
	std::cout << "�������� ���� ����������� ������: ";
	std::cin >> Niter;

	grid sigma(N, N);

	srand(time(NULL));

	std::ofstream outfile("output.dat");
	outfile << "T "
		<< " E"
		<< " devE"
		<< " M "
		<< "devM" << std::endl;
	for (T = Tstart; T <= Tstop; T += dT)
	{
		Jconst = 1.;

		sigma.init(sigma.up, 1.);

		double E = 0.0, devE = 0.0, M = 0.0, devM = 0.0;

		for (int k = 0; k < NiterTerm + Niter; ++k)
		{
			{ // MC ������, MC step

				int i, j;
				i = uniform_int(0, sigma.m_Nx - 1);
				j = uniform_int(0, sigma.m_Ny - 1);
				double dE = 2 * sigma.get_pbc(i, j) * (sigma.get_pbc(i - 1, j) + sigma.get_pbc(i - 1, j) + sigma.get_pbc(i, j + 1) + sigma.get_pbc(i, j + 1));
				double p = exp(-dE / T);
				if (p>1)
				{
					sigma.get(i, j) = -sigma.get(i, j);
				}
				else
				{
					double w = uniform();
					if (w < p)
					{
						sigma.get(i, j) = -sigma.get(i, j);
					}
				}

			} // ���� �� �� ��������, end MC step

			if (k >= NiterTerm)
			{ // �������� �� ����������, collect statistics
				/* (��������) double energy = calculate_grid_E(sigma, Jconst);
				E = E + energy;
				devE = devE + energy*energy; */

				double Einst = calculate_grid_E(sigma, Jconst);
				double Minst = calculate_grid_M(sigma, 1.0);
				E += Einst;
				M += Minst;
				devE += Einst * Einst;
				devM += Minst * Minst;

			} // ���� �� ���������� �� ����������, end statistics collection
		}

		{ // ���������� �� �������� ���������, averaging

			E /= Niter;
			M /= Niter;
			devE = sqrt(devE / Niter - E * E);
			devM = sqrt(devM / Niter - M * M);

			/* (��������) E = E / Niter;
			devE = sqrt(devE/Niter - E*E); */
		} // ���� �� ������������ �� �������� ���������

		std::cout << "T " << T << " E " << E << " +/- " << devE << " M " << M << " +/- " << devM << std::endl;
		outfile << T << " " << E << " " << devE << " " << M << " " << devM << std::endl;
	}
	outfile.close();
	return 0;
}
