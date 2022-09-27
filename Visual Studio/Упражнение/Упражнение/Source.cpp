#include <iostream>
#include <stdlib.h>
using namespace std;
float uniform(float umin, float umax)
{
	float rnd = 0;
	int ir = rand();
	float r = static_cast <float> (ir) / RAND_MAX;
	rnd = umin + (umax - umin) * r;
	return rnd;
}
int main() {

	cout << "Hello random!" << endl;
	for (int ii = 0; ii < 100; ++ii)
	{
		cout << static_cast <int>(uniform(1,7)) << "" << endl;
	}
		return 0;
}