#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <cmath>
#include <openacc.h>
#include "lambda_lanczos.hpp"

#include "acc_alg/ll_acc_run.hpp"
#include "grid.hpp"
#include "lattice_qcd.hpp"
#include "ll_util_extended.hpp"

using std::cout;
using std::endl;
using std::setprecision;
using std::setw;
using lambda_lanczos::LambdaLanczos;

using grid::Grid;
using grid::Grid2D;
using grid::Grid3D;
using lattice_qcd::Complex3D;

template<typename T>
using vector = std::vector<T>;

//Uncomment code in TEST, SMALL_INPUT_SAMPLE & LARGE_INPUT_SAMPLE as required.

int main() {

//For checking operator and eigenvalues [TEST]
//	const int vlen = 27;
//	const int num = 27*3;
//SMALL_INPUT_SAMPLE
//	const int vlen = 343;
//	const int num = 10;
//LARGE_INPUT_SAMPLE
	const int vlen = 1000;
	const int num = 10;	//Uses around 16G RAM (for 10) with eigenvectors for
				//ACC version, specifically for su3 data

//	const int num = 25;	//Uses around 10G RAM for ACC version, specifically for U(-100,100) data

	Grid3D<Complex3D<double>> mv_mul(vlen);

//TEST
//	mv_mul.read_input(std::string("config_files/sample_su3_27.txt"));
//SMALL_INPUT_SAMPLE
//	mv_mul.read_input(std::string("config_files/sample_su3_343.txt"));
//LARGE_INPUT_SAMPLE
	mv_mul.read_input(std::string("config_files/sample_su3_1000.txt"));

//START: TESTING OPERATOR [TEST]
//Don't use in PARALLEL_2, will throw error as mv_mul expects data in GPU
/*
	std::random_device rnd_device;
	std::mt19937 mersenne_engine {rnd_device()};

	std::uniform_real_distribution<double> dist { -100.0, 100.0 };

	auto gen = [&dist, &mersenne_engine]()
	{
		return dist(mersenne_engine);
	};

	vector<double> testxd(vlen*3*2), testyd(vlen*3*2);
	vector<std::complex<double>> testx(vlen*3), testy(vlen*3);
	vector<Complex3D<double>> testx_3(vlen), testy_3(vlen);
	vector<Complex3D<double>> w2_3(vlen), w1_3(vlen);

	std::generate(std::begin(testxd), std::end(testxd), gen);
	std::generate(std::begin(testyd), std::end(testyd), gen);
	
	for(unsigned int i=0; i<vlen*3; i++)
	{
		testx[i].real(testxd[i*2]);
		testx[i].imag(testxd[i*2+1]);
		testy[i].real(testyd[i*2]);
		testy[i].imag(testyd[i*2+1]);
	}

	for(unsigned int i=0; i<vlen*3; i++)
	{
		testx_3[i/3].set_colour(i%3, testx[i]);
		testy_3[i/3].set_colour(i%3, testy[i]);
	}

	mv_mul(testy_3, w1_3);
	mv_mul(testx_3, w2_3);

	std::complex<double> inn_prod1(0.0, 0.0), inn_prod2(0.0, 0.0), diff(0.0, 0.0);

	for(unsigned int i=0; i<vlen; i++)
	{
		inn_prod1+=(lambda_lanczos::util::typed_conj(testx_3[i])*w1_3[i]);
		inn_prod2+=(lambda_lanczos::util::typed_conj(testy_3[i])*w2_3[i]);
	}

	diff = inn_prod1 - std::conj(inn_prod2);

	std::cout<<"\nDiff: "<<diff<<'\n';

	if(std::abs(diff.real()) < 1e-9 && std::abs(diff.imag()) < 1e-9)
		std::cout<<"Operator verified!\n";
*/
//END: TESTING OPERATOR [TEST]

	LambdaLanczos<Complex3D<double>> engine(std::ref(mv_mul), vlen, vlen*3, false);

	std::vector<lambda_lanczos::util::real_t<Complex3D<double>>> eigvalues(num, 0.0);
	std::vector<std::vector<Complex3D<double>>> eigvecs(num, {vlen});

	auto res = engine.run(eigvalues, eigvecs);
//For TEST vlen only!
//	double sum{0.0}, sum_sq{0};

	cout << "\n\nIteration count: " << res <<"\n\n";
	for(unsigned int c=0; c<std::min(res,static_cast<long unsigned int>(num)); c++)
	{
		cout << "Eigenvalue "<< setw(5) << c << " :" << setprecision(16) << eigvalues[c] << endl;
		cout << endl;

//		sum+=eigvalues[c];
//		sum_sq+=(eigvalues[c]*eigvalues[c]);
	}

//For TEST vlen only!
//	cout<<"\nSum: "<<sum<<", Sum of squares: "<<sum_sq<<'\n';

//Printing eigenvectors, for TEST vlen only!
/*
	for(unsigned int c=0; c<std::min(res,static_cast<long unsigned int>(num)); c++)
	{
		cout << "Eigenvector "<< setw(5) << c << " :" << endl;
		for(unsigned int d=0; d<vlen; ++d)
			std::cout<<eigvecs[c][d];
		cout << endl;
	}
*/

	return EXIT_SUCCESS;
}
