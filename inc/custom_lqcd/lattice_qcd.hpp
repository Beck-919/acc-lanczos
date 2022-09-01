#ifndef _LATTICE_QCD_H_
#define _LATTICE_QCD_H_

#include <iostream>
#include <complex>
#include <vector>

namespace grid {
	template <typename T>
	class Grid3D;
}

namespace lattice_qcd{

	template <typename T> class Complex3D;

	//Addition
	template <typename T> Complex3D<T> operator+(const T&, const Complex3D<T>&);
	template <typename T> Complex3D<T> operator+(const Complex3D<T>&, const Complex3D<T>&);
	template <typename T> void operator+=(Complex3D<T>&, const Complex3D<T>&);

	//Subtraction
	template <typename T> Complex3D<T> operator-(const Complex3D<T>&, const Complex3D<T>&);
	template <typename T> void operator-=(Complex3D<T>&, const Complex3D<T>&);

	//Product
	template <typename T> Complex3D<T> operator*(const T&, const Complex3D<T>&);
	template <typename T> Complex3D<T> operator*(const Complex3D<T>&, const T&);
	template <typename T> Complex3D<T> operator*(const std::complex<T>&, const Complex3D<T>&);
	template <typename T> std::complex<T> operator*(const Complex3D<T>&, const Complex3D<T>&);
	template <typename T> void operator*=(Complex3D<T>&, const T&);
		
	//Out Stream
	template <typename T> std::ostream& operator<<(std::ostream& os, const Complex3D<T>&);
	
	template <typename T>
	class Complex3D
	{
		std::complex<T> colour[3];
		public:

		Complex3D()
		: colour{0.0}
		{}

		Complex3D(T val)
		: colour{val}
		{}

		Complex3D(std::complex<T> r, std::complex<T> g, std::complex<T> b)
		{
			colour[0]=r;
			colour[1]=g;
			colour[2]=b;
		}

		Complex3D(const Complex3D&)			=default;
		Complex3D(Complex3D &&)				=default;
		Complex3D &operator=(const Complex3D&)		=default;
		Complex3D &operator=(Complex3D &&)		=default;
		~Complex3D()					=default;

		std::complex<T> get_colour(int i) const
		{
			return colour[i];
		}

		void set_colour(int i, std::complex<T> val)
		{
			colour[i]= val;
		}

		//[SECTION1]
		//Complex3D= Complex3D [op] Complex3D
		friend Complex3D<T> operator+<T>(const Complex3D<T>&, const Complex3D<T>&);
		friend Complex3D<T> operator-<T>(const Complex3D<T>&, const Complex3D<T>&);
		friend std::complex<T> operator*<T>(const Complex3D<T>&, const Complex3D<T>&);

		//[SECTION2]
		//Complex3D= Complex3D [op] scalar *or* Complex3D= scalar [op] Complex3D
		friend Complex3D<T> operator+<T>(const T&, const Complex3D<T>&);
		friend Complex3D<T> operator*<T>(const T&, const Complex3D<T>&);
		friend Complex3D<T> operator*<T>(const Complex3D<T>&, const T&);
		friend Complex3D<T> operator*<T>(const std::complex<T>&, const Complex3D<T>&);

		//[SECTION3]
		//Complex3D [op]= scalar
		friend void operator*=<T>(Complex3D<T>&, const T&);

		//[SECTION4]
		//Complex3D [op]= Complex3D
		friend void operator+=<T>(Complex3D<T>&, const Complex3D<T>&);
		friend void operator-=<T>(Complex3D<T>&, const Complex3D<T>&);
		
		//[SECTION5]
		friend std::ostream& operator<<<T>(std::ostream& os, const Complex3D<T>&);
//		template<typename T>
//		friend Complex3D<T> typed_conj(const Complex3D<T>&);

		friend void grid::Grid3D<Complex3D<T>>:: operator() (const std::vector<Complex3D<T>>& , std::vector<Complex3D<T>>& );
	};

//[SECTION1]
	template<typename T>
	Complex3D<T> operator+(const Complex3D<T>& op1, const Complex3D<T>& op2)
	{
		Complex3D<T> ret1(op1);
		for(auto i=0; i<3; i++)
			ret1.colour[i]+=op2.colour[i];
		return ret1;
	}	

	template<typename T>
	Complex3D<T> operator-(const Complex3D<T>& op1, const Complex3D<T>& op2)
	{
		Complex3D<T> ret1(op1);
		for(auto i=0; i<3; i++)
			ret1.colour[i]-=op2.colour[i];
		return ret1;
	}	

	template<typename T>
	std::complex<T> operator*(const Complex3D<T>& op1, const Complex3D<T>& op2)
	{
		std::complex<T> ret1(0.0);
		for(auto i=0; i<3; i++)
			ret1=ret1 + op1.colour[i]*op2.colour[i];
		return ret1;
	}	

//[SECTION2]
	template<typename T>
	Complex3D<T> operator+(const T& scalar, const Complex3D<T>& op2)
	{
		Complex3D<T> ret1(op2);
		for(auto i=0; i<3; i++)
			ret1.colour[i]+=scalar;
		return ret1;
	}	

	template<typename T>
	Complex3D<T> operator*(const T& scalar, const Complex3D<T>& op2)
	{
		Complex3D<T> ret1(op2);
		for(auto i=0; i<3; i++)
			ret1.colour[i]*=scalar;
		return ret1;
	}

	template<typename T>
	Complex3D<T> operator*(const Complex3D<T>& op1, const T& scalar)
	{
		Complex3D<T> ret1(op1);
		for(auto i=0; i<3; i++)
			ret1.colour[i]*=scalar;
		return ret1;
	}	
	
	template<typename T>
	Complex3D<T> operator*(const std::complex<T>& scalar, const Complex3D<T>& op2)
	{
		Complex3D<T> ret1(op2);
		for(auto i=0; i<3; i++)
			ret1.colour[i]*=scalar;
		return ret1;
	}	

//[SECTION3]
	template<typename T>
	void operator*=(Complex3D<T>& op1, const T& scalar)
	{
		for(auto i=0; i<3; i++)
			op1.colour[i]*=scalar;
	}

//[SECTION4]
	template<typename T>
	void operator+=(Complex3D<T>& op1, const Complex3D<T>& op2)
	{
		for(auto i=0; i<3; i++)
			op1.colour[i]+=op2.colour[i];
	}

	template<typename T>
	void operator-=(Complex3D<T>& op1, const Complex3D<T>& op2)
	{
		for(auto i=0; i<3; i++)
			op1.colour[i]-=op2.colour[i];
	}

//[SECTION5]
	template<typename T>
	std::ostream& operator<<(std::ostream& os, const Complex3D<T>& op)
	{
		os<<"[ "<<op.colour[0]<<", "<<op.colour[1]<<", "<<op.colour[2]<<"]\n";
		return os;
	}

}

#endif
