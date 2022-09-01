#ifndef _LL_ACC_RUN_HPP_
#define _LL_ACC_RUN_HPP_

#include <iostream>
#include <array>
#include <vector>
#include <tuple>
#include <functional>
#include <cassert>
#include <limits>
#include <cmath>
#include <numeric>
#include <random>
#include <algorithm>
#include <utility>
#include <unistd.h>
#include <string>
#include <chrono>

#include "grid.hpp"
#include "lambda_lanczos.hpp"
#include "lambda_lanczos_util.hpp"
#include "lambda_lanczos_tridiagonal.hpp"

#ifdef ACC_ALG
#include "ll_acc_util.hpp"
#include <openacc.h>
#endif

template <typename n_type>
using real_t = lambda_lanczos::util::real_t<n_type>;

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::minutes;
using std::chrono::seconds;
using std::chrono::milliseconds;

namespace lambda_lanczos {
#ifdef ACC_ALG
	template <typename T>
	size_t LambdaLanczos<T>:: acc_run(std::vector<real_t<T>>& eigvalues, std::vector<std::vector<T>>& eigvecs)
	{
		std::cout<<"\nIn accelerated code...\n";

		const size_t nroot = eigvecs.size();
		assert(eigvalues.size()==nroot);
		assert(0 < this->tridiag_eps_ratio && this->tridiag_eps_ratio < 1);

		std::vector<std::vector<T>> u; // Lanczos vectors
		std::vector<real_t<T>> alpha;  // Diagonal elements of an approximated tridiagonal matrix
		std::vector<real_t<T>> beta;   // Subdiagonal elements of an approximated tridiagonal matrix

		u.reserve(this->initial_vector_size);
		alpha.reserve(this->initial_vector_size);
		beta.reserve(this->initial_vector_size);

		const auto n = this->matrix_size;

		u.emplace_back(n);
		this->init_vector(u[0]);

		util::normalize(u[0]);

		std::vector<real_t<T>> evs(nroot); 						// Calculated eigenvalue
		std::vector<real_t<T>> pevs(nroot, std::numeric_limits<real_t<T>>::max()); 	// Previous eigenvalue

		size_t itern = this->max_iteration;

		real_t<T>	h_alpha, h_beta;
		std::vector<T>	ukm2(n, 0.0);
		std::vector<T>	ukm1(u[0]);
		std::vector<T>	uk(n, 0.0);							//This can au

		T *ukm2_ptr= 	&ukm2[0];
		T *ukm1_ptr= 	&ukm1[0];
		T *uk_ptr=	&uk[0];

		T *sliceA = nullptr;
		T *sliceB = nullptr;

		auto start_iter = high_resolution_clock::now();
#pragma acc data copyin(ukm1_ptr[0:n]) create(uk_ptr[0:n], ukm2_ptr[0:n])
		{	//DATA_REGION_1

		size_t arr_size=0;
		size_t arr_size_b=0;

		arr_size_b = acc_get_property(0, acc_device_current, acc_property_free_memory);

		arr_size_b/=2;

		if(arr_size_b>104857600)
			arr_size_b-=104857600;		// Spares 200MB in GPU global memory;
		else
			std::cout<<"\nNot enough memory on device\n";

		if(arr_size_b>0)
		{
			arr_size=arr_size_b/sizeof(T);
//			arr_size=2*u[0].size();
			sliceA=new T[arr_size];
			sliceB=new T[arr_size];
		}

		if(sliceA==nullptr || sliceB==nullptr)
			std::cout<<"\nNot enough memory on CPU\n";

#pragma acc data create(sliceA[0:arr_size], sliceB[0:arr_size])
		{	//DATA_REGION_2

		for(size_t k = 1; k <= this->max_iteration; ++k)
		{
//GPU_CODE_1
#pragma acc wait(1,2)
			this->mv_mul(ukm1, uk);

#pragma acc parallel loop present(ukm1_ptr[0:n], uk_ptr[0:n]) async(2)
			for(size_t i = 0; i < n; ++i)
				uk_ptr[i] += ukm1_ptr[i]*this->eigenvalue_offset;

			h_alpha = std::real(util::d_inner_prod(ukm1_ptr, uk_ptr, n, 2));
			//#pragma acc wait(2), inside d_inner_prod
			//As async ends there and requires a sync with
			//device for std::real

//CPU_CODE_1
			alpha.push_back(h_alpha);

			h_beta = (k==1)? 0 : beta[k-2];

//GPU_CODE_2
#pragma acc parallel loop present(uk_ptr[0:n], ukm1_ptr[0:n], ukm2_ptr[0:n]) copyin(h_alpha, h_beta) async(2)
			for(size_t i = 0; i < n; ++i)
				if(k == 1)
					uk_ptr[i] = uk_ptr[i] - h_alpha*ukm1_ptr[i];
				else
					uk_ptr[i] = uk_ptr[i] - h_beta*ukm2_ptr[i] - h_alpha*ukm1_ptr[i];

//DtoD: ukm2_ptr = ukm1_ptr
			acc_memcpy_device_async(acc_deviceptr(ukm2_ptr), acc_deviceptr(ukm1_ptr), sizeof(T)*n, 1);

			//Streams 2(data) & 3(compute) used inside
			if(k>1 && sliceA!=nullptr && sliceB!=nullptr)
				util::dualstream_schmidt_orth(	uk,
								sliceA, u.begin()		, u.begin()  + u.size()/2,
								sliceB, u.begin() + u.size()/2	, u.end()    - 1,
								arr_size/u[0].size()		, ((u.size()-1)%2==0)? (u.size()-1)/2: (u.size()-1)/2 + 1);

#pragma acc update host(uk_ptr[0:n]) wait(3)

//CPU_CODE_2

			u.push_back(uk);

      			beta.push_back(util::norm(u[k]));

			util::normalize(u[k]);

			T *t_uk = &u[k][0];

//HtoD: ukm1_ptr = uk_ptr
			acc_memcpy_to_device_async(acc_deviceptr(ukm1_ptr), t_uk, sizeof(T)*n, 1);

			for (size_t iroot = 0; iroot < nroot; ++iroot)
				evs[iroot] =	tridiagonal::find_mth_eigenvalue(alpha, beta,
							this->find_maximum ? alpha.size()-1-iroot : iroot,
							this->eps * this->tridiag_eps_ratio);

			const real_t<T> zero_threshold = util::minimum_effective_decimal<real_t<T>>()*1e-1;

			if(beta.back() < zero_threshold)
			{
				itern = k;
				break;
			}

			bool break_cond = true;
			for(size_t iroot = 0; iroot < nroot; ++iroot)
			{
				const auto& ev = evs[iroot];
				const auto& pev = pevs[iroot];
				if (std::abs(ev - pev) >= std::min(std::abs(ev), std::abs(pev)) * this->eps)
				{

					break_cond = false;
					break;
				}
			}

			if (break_cond)
			{
				itern = k;
				break;
			} 
			else 
				pevs = evs;
		}

		}	//DATA_REGION_2

		if(sliceA!=nullptr)
			delete []sliceA;
		if(sliceB!=nullptr)
			delete []sliceB;

		}	//DATA_REGION_1
		auto stop_iter = high_resolution_clock::now();

		auto mins = duration_cast<minutes>(stop_iter - start_iter);
		auto secs = duration_cast<seconds>(stop_iter - start_iter - mins);
		auto ms = duration_cast<milliseconds>(stop_iter - start_iter - mins - secs);
		std::cout<<"\nTime taken for Lanczos algorithm: "<<mins.count()<<"m"<<secs.count()<<"."<<ms.count()<<"s";

		eigvalues = evs;
		beta.back() = 0.0;

		for(size_t iroot = 0; iroot < nroot; ++iroot) {
			auto& eigvec = eigvecs[iroot];
			auto& ev = eigvalues[iroot];

			eigvec = eigenvector(alpha, beta, u, ev);
			ev -= this->eigenvalue_offset;
		}

		return itern;
	}
#endif	// ACC_ALG
}

#endif
