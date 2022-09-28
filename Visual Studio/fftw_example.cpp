#include <iostream>
#include <cmath>
#include "fftw3.h"

using namespace std;

int main()
{
    cout << "Hello world!" << endl;   //I have made changes

    int N=12;

    double *in, *out;
    fftw_plan p;

    in = (double*) fftw_malloc(sizeof(double) * N);
    out = (double*) fftw_malloc(sizeof(double) * N);

    p = fftw_plan_r2r_1d(N,in,out,FFTW_RODFT00,FFTW_ESTIMATE);

    for (int i=0; i<N; i++)
    {
        in[i]=sin(i*8*180/N/3.14);
    }

    fftw_execute(p);

    for (int i=0; i<N; i++)
    {
        std::cout<<in[i]<<" "<<out[i]<<endl;
    }

    fftw_destroy_plan(p);

    fftw_free(in);
    fftw_free(out);



    return 0;
}
