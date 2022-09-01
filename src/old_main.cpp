#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <cmath>
#include <openacc.h>
#include "lambda_lanczos.hpp"

#include "grid.hpp"
#include "lattice_qcd.hpp"

using std::cout;
using std::endl;
using std::setprecision;
using std::setw;
using lambda_lanczos::LambdaLanczos;

using grid::Grid;
using grid::Grid2D;
using grid::Grid3D;

template<typename T>
using vector = std::vector<T>;

int main() {

	const int vlen = 27;
	const int num = 3;

	Grid3D<std::complex<double>> mv_mul(vlen);

	LambdaLanczos<std::complex<double>> engine(mv_mul, vlen, false); // false means to calculate the smallest eigenvalue.

	std::vector<lambda_lanczos::util::real_t<std::complex<double>>> eigvalues(num, 0.0);
	std::vector<std::vector<std::complex<double>>> eigvecs(num);

	auto res = engine.run(eigvalues, eigvecs);

	cout << "\n\nIteration count: " << res <<"\n\n";

	for(unsigned int c=0; c<std::min(res,static_cast<long unsigned int>(num)); c++)
	{
		cout << "Eigen value "<< setw(5) << c+1 << " :" << setprecision(16) << eigvalues[c] << endl;
		cout << endl;
	}

	for(unsigned int c=0; c<std::min(res,static_cast<long unsigned int>(num)); c++)
	{
		cout << "\nEigenvector "<< setw(5) << c+1 << " :" << endl;
		cout<<"[ ";
		for(unsigned int d=0; d<vlen; ++d)
			std::cout<<eigvecs[c][d]<<" ";
		cout<<"]";
		cout << endl;
	}

	return EXIT_SUCCESS;
}
