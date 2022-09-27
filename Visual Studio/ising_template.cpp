#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <fstream>
// the code is slow - pbc was implemented without halo - cache misses

// Задача: Добавете необходимия програмен код, така че приложението да 
// сканира енергията и магнитния момент на модел на Изинг 
// в зададен температурен диапазон

// случайни числа [0,1] - равномено разпределение
inline double uniform()
{
    // добавете съответния програмен код
}

// случайни чели числа [low,up] - равномено разпределение
inline int uniform_int(int low, int up)
{
    // добавете съответния програмен код
}

class grid
{
public:
    const int up = 1;
    const int down = -1;
    int *m_data;
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

    int &get(int row, int col)
    {
        auto element_index = row * m_Ny + col;
        return m_data[element_index];
    };

    int &get_pbc(int row, int col)
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
double calculate_grid_E(grid &sigma, double J = 1.0)
{
    double E = 0.;

    for (size_t i = 0; i < sigma.m_Nx; i++)
    {
        for (size_t j = 0; j < sigma.m_Ny; j++)
        {
            float local_sum = sigma.get_pbc(i + 1, j) + sigma.get_pbc(i - 1, j);
            local_sum += sigma.get_pbc(i, j + 1) + sigma.get_pbc(i, j - 1);

            // добавете съответния програмен код

        }
    }

    return E / (sigma.m_Nx * sigma.m_Ny);
}

// eq. 13.69
double calculate_grid_M(grid &sigma, double mu = 1.0)
{
    double M = 0.;

    for (size_t i = 0; i < sigma.m_Nx; i++)
    {
        for (size_t j = 0; j < sigma.m_Ny; j++)
        {

            // добавете съответния програмен код

        }
    }

    return M / (sigma.m_Nx * sigma.m_Ny);
}

////////////////////////////////////////////////
///                 MAIN                     ///
////////////////////////////////////////////////

int main()
{
    int N;         // размер на решетката , grid dim.
    double T;      // температура, temperature
    int NiterTerm; // брой стъпки за термализиране на системата, number of the quilibration steps
    int Niter;     // брой продуктивни стъпки, number of the production sreps
    double Tstart, Tstop, dT;
    double Jconst = 1.;

    std::cout << "Въведете размет на решетката: ";
    std::cin >> N;
    std::cout << "Въведете начална температура: ";
    std::cin >> Tstart;
    std::cout << "Въведете крайна температура: ";
    std::cin >> Tstop;
    std::cout << "Въведете стъпка за сканиране по температура: ";
    std::cin >> dT;
    std::cout << "Въведете брой стъпки за термализиране на системата: ";
    std::cin >> NiterTerm;
    std::cout << "Въведете брой продуктивни стъпки: ";
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
            { // MC стъпка, MC step

                // добавете съответния програмен код

            } // край на МС стъпката, end MC step

            if (k >= NiterTerm)
            { // събиране на статистика, collect statistics
				/* (примерно) double energy = calculate_grid_E(sigma, Jconst);
				E = E + energy;
				devE = devE + energy*energy; */
                // добавете съответния програмен код

            } // край на събирането на статистика, end statistics collection
        }

        { // пресматане на средните стойности, averaging
            // добавете съответния програмен код
			/* (примерно) E = E / Niter;
			devE = sqrt(devE/Niter - E*E); */
        } // край на пресмятането на средните стойности

        std::cout << "T " << T << " E " << E << " +/- " << devE << " M " << M << " +/- " << devM << std::endl;
        outfile << T << " " << E << " " << devE << " " << M << " " << devM << std::endl;
    }
    outfile.close();
    return 0;
}
