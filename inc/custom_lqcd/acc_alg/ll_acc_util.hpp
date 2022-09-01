#ifndef _LL_ACC_UTIL_HPP_
#define _LL_ACC_UTIL_HPP_

#include <complex>
#include "lambda_lanczos_util.hpp"

namespace lambda_lanczos {

	namespace util {

//Standard device versions of d_inner_prod and dualstream_schmidt_orth
//For Complex3D search -> COMPLEX3D

		template <typename T>
		inline T d_inner_prod(T *v1, T *v2, size_t n)
		{
			T ret{0.0};
#pragma acc parallel loop present(v1[0:n], v2[0:n]) reduction(+:ret) async(3)
			for(size_t i=0; i<n; ++i)
				ret+=(typed_conj(v1[i])*v2[i]);
#pragma acc wait(3)
			return ret;
		}

		template <typename ForwardIterator, typename T>
		inline void dualstream_schmidt_orth(std::vector<T>& uorth,
	                         T *sliceA, ForwardIterator first, ForwardIterator midp1,
        	                 T *sliceB, ForwardIterator midp1_2, ForwardIterator last,
				 size_t max_slice_size, size_t total_vec_half)
		{
  			const auto n = uorth.size();
			size_t cntA, cntB;
			size_t sizeA, sizeB;
			size_t item_size;
			bool break_check=false;

			ForwardIterator iterA_beg;
			ForwardIterator iterB_beg;
			ForwardIterator iterA_end;
			ForwardIterator iterB_end;

			item_size=sizeof(Complex3D<T>);

			T *uorth_ptr= &uorth[0];

			size_t outer_lim= static_cast<size_t>(std::ceil(total_vec_half/(double)max_slice_size));

			if(outer_lim>1)
			{
				iterA_beg = first;
				iterA_end = iterA_beg + max_slice_size;
				iterB_beg = iterA_end;
				iterB_end = iterB_beg + max_slice_size;
			}
			else
			{
				iterA_beg = first;
				iterA_end = midp1;
				iterB_beg = midp1_2;
				iterB_end = last;
			}


			for(size_t it=0; it<outer_lim; ++it)
			{
				if(it>0)
				{
					auto left= std::distance(iterB_end, last);
					if(left/2 >= max_slice_size)
					{
						iterA_beg = iterB_end;
						iterA_end = iterA_beg + max_slice_size;
						iterB_beg = iterA_end;
						iterB_end = iterB_beg + max_slice_size;
					}
					else if(left>1)
					{
						iterA_beg = iterB_end;
						iterA_end = iterA_beg + left/2;
						iterB_beg = iterA_end;
						iterB_end = last;
					}
					else
					{
						iterA_beg = iterB_end;
						iterA_end = last;

						break_check=true;
					}
				}

				cntA=0;
				sizeA=0;
				for(auto iter= iterA_beg; iter!=iterA_end; ++iter)
				{
					T *buf1 = &((*iter)[0]);
	
					acc_memcpy_to_device_async(acc_deviceptr(&sliceA[cntA]), buf1, n*item_size, 2);
	
					cntA+=n;
					sizeA++;
				}

#pragma acc wait(2,3)
			 	for(size_t iter = 0; iter < sizeA; ++iter)
				{
	    				T *uk = &sliceA[iter*n];
	    				T innprod = d_inner_prod(uk, uorth_ptr, n, 3);
	
#pragma acc parallel loop present(uk[0:n], uorth_ptr[0:n]) copyin(innprod) async(3)
	    				for(size_t i = 0; i < n; ++i)
						uorth_ptr[i] -= innprod * uk[i];
	  			}

				if(break_check)
					break;
	
				cntB=0;
				sizeB=0;
				for(auto iter= iterB_beg; iter!=iterB_end; ++iter)
				{
					T *buf2 = &((*iter)[0]);
	
					acc_memcpy_to_device_async(acc_deviceptr(&sliceB[cntB]), buf2, n*item_size, 2);
	
					cntB+=n;
					sizeB++;
				}

#pragma acc wait(2,3)
			 	for(size_t iter = 0; iter < sizeB; ++iter)
				{
	    				T *uk = &sliceB[iter*n];
	    				T innprod = d_inner_prod(uk, uorth_ptr, n, 3);
	
#pragma acc parallel loop present(uk[0:n], uorth_ptr[0:n]) copyin(innprod) async(3)
	    				for(size_t i = 0; i < n; ++i)
						uorth_ptr[i] -= innprod * uk[i];
	  			}
			}

		}


#ifdef COMPLEX3D
		template <typename T>
		inline std::complex<T> d_inner_prod(Complex3D<T> *v1, Complex3D<T> *v2, size_t n, int compute_stream)
		{
			std::complex<T> ret{0.0};
			T ret_r{0.0}, ret_i{0.0};	//Hack - As std::complex is not supported in OpenACC reduction

#pragma acc parallel loop present(v1[0:n], v2[0:n]) reduction(+:ret_r, ret_i) private(ret) async(compute_stream)
			for(size_t i=0; i<n; ++i)
			{
				ret= typed_conj(v1[i])*v2[i];

				ret_r += ret.real();
				ret_i += ret.imag();
			}
#pragma acc wait(compute_stream)
			ret.real(ret_r);
			ret.imag(ret_i);

			return ret;
		}

		template <typename ForwardIterator, typename T>
		inline void dualstream_schmidt_orth(std::vector<Complex3D<T>>& uorth,
	                         Complex3D<T> *sliceA, ForwardIterator first, ForwardIterator midp1,
        	                 Complex3D<T> *sliceB, ForwardIterator midp1_2, ForwardIterator last,
				 size_t max_slice_size, size_t total_vec_half)
		{
  			const auto n = uorth.size();
			size_t cntA, cntB;
			size_t sizeA, sizeB;
			size_t item_size;
			bool break_check=false;

			ForwardIterator iterA_beg;
			ForwardIterator iterB_beg;
			ForwardIterator iterA_end;
			ForwardIterator iterB_end;

			item_size=sizeof(Complex3D<T>);

			Complex3D<T> *uorth_ptr= &uorth[0];

			size_t outer_lim= static_cast<size_t>(std::ceil(total_vec_half/(double)max_slice_size));

			if(outer_lim>1)
			{
				iterA_beg = first;
				iterA_end = iterA_beg + max_slice_size;
				iterB_beg = iterA_end;
				iterB_end = iterB_beg + max_slice_size;
			}
			else
			{
				iterA_beg = first;
				iterA_end = midp1;
				iterB_beg = midp1_2;
				iterB_end = last;
			}


			for(size_t it=0; it<outer_lim; ++it)
			{
				if(it>0)
				{
					auto left= std::distance(iterB_end, last);
					if(left/2 >= max_slice_size)
					{
						iterA_beg = iterB_end;
						iterA_end = iterA_beg + max_slice_size;
						iterB_beg = iterA_end;
						iterB_end = iterB_beg + max_slice_size;
					}
					else if(left>1)
					{
						iterA_beg = iterB_end;
						iterA_end = iterA_beg + left/2;
						iterB_beg = iterA_end;
						iterB_end = last;
					}
					else
					{
						iterA_beg = iterB_end;
						iterA_end = last;

						break_check=true;
					}
				}

				cntA=0;
				sizeA=0;
				for(auto iter= iterA_beg; iter!=iterA_end; ++iter)
				{
					Complex3D<T> *buf1 = &((*iter)[0]);
	
					acc_memcpy_to_device_async(acc_deviceptr(&sliceA[cntA]), buf1, n*item_size, 2);
	
					cntA+=n;
					sizeA++;
				}

#pragma acc wait(2,3)
			 	for(size_t iter = 0; iter < sizeA; ++iter)
				{
	    				Complex3D<T> *uk = &sliceA[iter*n];
	    				std::complex<T> innprod = d_inner_prod(uk, uorth_ptr, n, 3);
	
#pragma acc parallel loop present(uk[0:n], uorth_ptr[0:n]) copyin(innprod) async(3)
	    				for(size_t i = 0; i < n; ++i)
						uorth_ptr[i] -= innprod * uk[i];
	  			}

				if(break_check)
					break;
	
				cntB=0;
				sizeB=0;
				for(auto iter= iterB_beg; iter!=iterB_end; ++iter)
				{
					Complex3D<T> *buf2 = &((*iter)[0]);
	
					acc_memcpy_to_device_async(acc_deviceptr(&sliceB[cntB]), buf2, n*item_size, 2);
	
					cntB+=n;
					sizeB++;
				}

#pragma acc wait(2,3)
			 	for(size_t iter = 0; iter < sizeB; ++iter)
				{
	    				Complex3D<T> *uk = &sliceB[iter*n];
	    				std::complex<T> innprod = d_inner_prod(uk, uorth_ptr, n, 3);
	
#pragma acc parallel loop present(uk[0:n], uorth_ptr[0:n]) copyin(innprod) async(3)
	    				for(size_t i = 0; i < n; ++i)
						uorth_ptr[i] -= innprod * uk[i];
	  			}
			}

		}
#endif //COMPLEX3D
	}
}

#endif
