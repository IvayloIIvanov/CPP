#include<iostream>
#include<complex>

using namespace std;

int main() {
	complex<double> a, b, c, D, x1, x2;
	cout << "a = ";
	cin >> a;
	cout << "b = ";
	cin >> b;
	cout << "c = ";
	cin >> c;
	D = b * b - 4. * a * c;
	x1 = (-b + sqrt(D)) /( 2. * a);
	x2 = (-b - sqrt(D))/( 2. * a);

	cout << "D = " << D << endl;
	cout << "x1 " << x1 << endl;
	cout << "x2 = " << x2 << endl;
}