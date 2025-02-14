#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <Windows.h>

using namespace std;

class Vector {
private:
	double* x;
	unsigned int Lenght;

public:

	Vector() {
		Lenght = 1;
		x = new double[1];
		x[0] = 0;
	}

	Vector(unsigned int N) {
		Lenght = N;
		x = new double[N];
		for (unsigned int i = 0; i < N; i++) {
			x[i] = 0;
		}
	}

	~Vector() {
		delete[] x;
	}

	Vector(const Vector& A) {
		Lenght = A.Lenght;
		x = new double[Lenght];
		for (unsigned int i = 0; i < Lenght; i++) {
			x[i] = A.x[i];
		}
	}

	Vector(unsigned int N, double arrgs, ...) {
		Lenght = N;
		x = new double[Lenght];
		double* p = &arrgs;
		for (unsigned int i = 0; i < N; i++) {
			x[i] = *p;
			p++;
		}
	}

	unsigned int GetLenght() {
		return Lenght;
	}

	double GetComponent(unsigned int i) {
		return x[i];
	}

	void SetLenght(unsigned int N) {
		double* tmp = x;
		x = new double[N];
		for (unsigned int i = 0; i < Lenght && i < N; i++) {
			x[i] = tmp[i];
		}
		if (N > Lenght) {
			for (unsigned int i = Lenght; i < N; i++) {
				x[i] = 0;
			}
		}
		Lenght = N;
		delete[] tmp;
	}

	void SetComponent(unsigned int i, double newX) {
		x[i] = newX;
	}

	friend ostream& operator<< (ostream& out, const Vector& vec);
	friend istream& operator>> (istream& in, Vector& vec);

	Vector& operator=(const Vector& newVec) {
		if (this == &newVec) {
			return *this;
		}
		this->SetLenght(newVec.Lenght);
		for (unsigned int i = 0; i < Lenght; i++) {
			this->SetComponent(i, newVec.x[i]);
		}
		return *this;
	}

	friend Vector operator+(const Vector& left, const Vector& right);
	friend Vector operator-(const Vector& left, const Vector& right);
	friend Vector operator*(const Vector& left, const double right);
	friend Vector operator*(const double left, const Vector& right);

	Vector& operator+=(const Vector& right) {
		for (unsigned int i = 0; i < Lenght; i++) {
			x[i] += right.x[i];
		}
		return *this;
	}

	Vector& operator-=(const Vector& right) {
		for (unsigned int i = 0; i < Lenght; i++) {
			x[i] -= right.x[i];
		}
		return *this;
	}

	Vector& operator*=(const double right) {
		for (unsigned int i = 0; i < Lenght; i++) {
			x[i] *= right;
		}
		return *this;
	}
};

Vector operator+(const Vector& left, const Vector& right) {
	Vector res(left.Lenght);
	for (unsigned int i = 0; i < res.Lenght; i++) {
		res.SetComponent(i, (left.x[i] + right.x[i]));
	}
	return res;
}

Vector operator-(const Vector& left, const Vector& right) {
	Vector res(left.Lenght);
	for (unsigned int i = 0; i < res.Lenght; i++) {
		res.SetComponent(i, (left.x[i] - right.x[i]));
	}
	return res;
}

Vector operator*(const Vector& left, const double right) {
	Vector res(left.Lenght);
	for (unsigned int i = 0; i < res.Lenght; i++) {
		res.SetComponent(i, (left.x[i] * right));
	}
	return res;
}

Vector operator*(const double left, const Vector& right) {
	Vector res(right.Lenght);
	for (unsigned int i = 0; i < res.Lenght; i++) {
		res.SetComponent(i, (right.x[i] * left));
	}
	return res;
}

ostream& operator<< (ostream& out, const Vector& vec) {
	out << "(";
	for (unsigned int i = 0; i < vec.Lenght - 1; i++) {
		out << vec.x[i] << ", ";
	}
	out << vec.x[vec.Lenght - 1] << ")";
	return out;
}

istream& operator>> (istream& in, Vector& vec) {
	unsigned int Lenght;
	in >> Lenght;
	vec.SetLenght(Lenght);
	double tmp;
	for (unsigned int i = 0; i < Lenght; i++) {
		in >> tmp;
		vec.SetComponent(i, tmp);
	}
	return in;
}

//optimal h (3M0Em/2M3)^0.33

Vector CalculateGradient(double(*f)(Vector&), Vector& x, double h) {
	Vector res(x.GetLenght());
	Vector H(x.GetLenght());
	Vector Hminus;
	Vector Hplus;
	for (unsigned int i = 0; i < res.GetLenght(); i++) {
		H.SetComponent(i, h);
		Hminus = x - H;
		Hplus = x + H;
		res.SetComponent(i, (f(Hplus) - f(Hminus)) * (1 / (2 * h)));
		H.SetComponent(i, 0);
	}
	return res;
}

struct Test {
public:
	double (*f)(Vector&);
	const unsigned int FunctionDim;
	Vector startX;
	string functionName;

	Test(double (*Func)(Vector&), const unsigned int FunctionDimensionCount, Vector startX, string funcName) : FunctionDim(FunctionDimensionCount) {
		f = Func;
		this->startX = startX;
		functionName = funcName;
	}
};

Test tests[] = {
	{ [](Vector& x) { return sin(0.5 * pow(x.GetComponent(0),2) - 0.25 * pow(x.GetComponent(1), 2) + 3) * cos(2 * x.GetComponent(0) + 1 + pow(M_E,x.GetComponent(1))); },
	  2, Vector(2, 0.1, 0.1), "sin(0.5x^2 - 0.25y^2+3)cos(2*x + 1 + e^y)"},
	{ [](Vector& x) {return pow(x.GetComponent(0),2) + x.GetComponent(0) + 1; },
	  1, Vector(1, 0), "x^2 + x + 1"},
	{ [](Vector& x) {return pow(x.GetComponent(0) - 100,2); },
	  1, Vector(1, 0), "(x-100)^2"}
};

int main()
{
	int counter;
	const double h = 0.01; //for calc derivatives
	const double eps = 0.00001; //accuracy
	const double lambda = 0.5; //speed param

	Vector Grad;
	Vector PrevX;
	Vector CurX;

	int chose;

	do {
		cout << " chose test function 0 - 2 	exit - 9" << endl;
		for (int i = 0; i < 3; i++) {
			cout << tests[i].functionName << endl;
		}
		cin >> chose;
		cout << " Start X " << endl;
		cout << tests[chose].startX << endl;
		if (chose != 9 && chose >= 0 && chose <= 2) {
			counter = 0;
			CurX = tests[chose].startX;
			do {
				cout << " Iter " << counter << endl;
				cout << " Cur X " << CurX << endl;
				PrevX = CurX;
				Grad = CalculateGradient(tests[chose].f, CurX, h);
				cout << " Grad " << Grad << endl;
				CurX -= lambda * Grad;
				counter++;
			} while (counter < 1000 && abs(tests[chose].f(PrevX) - tests[chose].f(CurX) > eps));

			if (abs(tests[chose].f(PrevX) - tests[chose].f(CurX) < eps)) {
				cout << " Res: " << CurX << endl;
				cout << " F = " << tests[chose].f(CurX) << endl;
			}
			else {
				cout << " Answer not found in 1000 iterations " << endl;
			}
		}
		cout << endl;
	} while (chose != 9);
	return 0;
}