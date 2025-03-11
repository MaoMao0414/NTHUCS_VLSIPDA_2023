#include "ExampleFunction.h"
#include "NumericalOptimizer.h"
#include <iostream>

using namespace std;

ExampleFunction::ExampleFunction()
{
}

void ExampleFunction::evaluateFG(const vector<double> &x, double &f, vector<double> &g)
{	
																 // minimize 3*x^2 + 2*x*y + 2*y^2 + 7
	f = 3 * x[0] * x[0] + 2 * x[0] * x[1] + 2 * x[1] * x[1] + 7; // objective function
	g[0] = 6 * x[0] + 2 * x[1];									 // gradient function of x[0]
	g[1] = 2 * x[0] + 4 * x[1];									 // gradient function of x[1]

	cout << x[0] << " " << x[1] << endl;

}

void ExampleFunction::evaluateF(const vector<double> &x, double &f)
{																 // minimize 3*x^2 + 2*x*y + 2*y^2 + 7
	f = 3 * x[0] * x[0] + 2 * x[0] * x[1] + 2 * x[1] * x[1] + 7; // objective function
}

unsigned ExampleFunction::dimension()
{
	return 2; // number of variables
}

int main()
{						// minimize 3*x^2 + 2*x*y + 2*y^2 + 7
	ExampleFunction ef; // require to define the object function and gradient function

	vector<double> x(2); // solution vector, number of variables
	x[0] = 100;			 // initialize the solution vector
	x[1] = 100;

	NumericalOptimizer no(ef);
	no.setX(x);				// set initial solution
	no.setNumIteration(35); // user-specified parameter
	no.setStepSizeBound(10); // user-specified parameter
	no.solve();				// Conjugate Gradient solver

	cout << "Current solution:\n";
	for (unsigned i = 0; i < no.dimension(); i++)
	{
		cout << "x[" << i << "] = " << no.x(i) << "\n";
	}
	cout << "Objective: " << no.objective() << "\n";
	////////////////////////////////////////////////////////////////
	return 0;
}
