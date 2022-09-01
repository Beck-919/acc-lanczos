#ifndef _LL_UTIL_EXTENDED_HPP_
#define _LL_UTIL_EXTENDED_HPP_

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

#include "lattice_qcd.hpp"
#include "lambda_lanczos_util.hpp"

using lattice_qcd::Complex3D;

namespace lambda_lanczos {

	namespace util {

		template <typename T>
		struct realTypeMap<Complex3D<T>> {
			typedef T type;
		};

		template<typename T>
		inline Complex3D<T> typed_conj(const Complex3D<T>& val)
		{
			Complex3D<T> ret1;
			for(auto i=0; i<3; i++)
				ret1.set_colour(i,TypedConjugate<std::complex<T>>::invoke(val.get_colour(i)));
			return ret1;
		}

	}
}

#endif
