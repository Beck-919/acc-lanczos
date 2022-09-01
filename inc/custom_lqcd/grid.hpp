#ifndef _GRID_HPP_
#define _GRID_HPP_

#include <cmath>
#include <vector>
#include <complex>
#include <string>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <fstream>
#include <iostream>
#include <iomanip>

#ifdef COMPLEX3D
#include "lattice_qcd.hpp"
#include "ll_util_extended.hpp"
using lattice_qcd::Complex3D;
#endif

namespace grid {

	template <typename T>
	class Grid
	{
		protected:

		unsigned int dim;
		unsigned int sq_root;
		unsigned int cu_root, cu_root_sq;

#pragma acc routine seq
		unsigned int neighbour2D(char direction, unsigned int pos)
		{
//			unsigned int n= (unsigned int)(std::sqrt(dim));

			switch(direction)
			{
				case 'l':
					return (sq_root+(pos%sq_root)-1)%sq_root + (pos/sq_root)*sq_root;
				case 'r':
					return (sq_root+(pos%sq_root)+1)%sq_root + (pos/sq_root)*sq_root;
				case 't':
					return (dim + pos + sq_root)%dim;
				case 'b':
					return (dim + pos - sq_root)%dim;
			}
			return 0;
		}

#pragma acc routine seq
		unsigned int neighbour3D(char direction, unsigned int pos)
		{
//			unsigned int n= (unsigned int)(std::cbrt(dim));
//			unsigned int n= std::lrint(std::cbrt(this->dim));
//			unsigned int nsq, quo, rem;
			unsigned int quo, rem;

//			nsq=		n*n;
			quo=		pos/cu_root_sq;
			rem=		pos%cu_root_sq;

			switch(direction)
			{
				case 'l':
					return (cu_root+(pos%cu_root)-1)%cu_root + (pos/cu_root)*cu_root;
				case 'r':
					return (cu_root+(pos%cu_root)+1)%cu_root + (pos/cu_root)*cu_root;
				case 't':
					return (cu_root_sq+rem+cu_root)%cu_root_sq + quo*cu_root_sq;
				case 'b':
					return (cu_root_sq+rem-cu_root)%cu_root_sq + quo*cu_root_sq;
				case 'f':
					return (dim+pos-cu_root_sq)%dim;
				case 'h':
					return (dim+pos+cu_root_sq)%dim;
			}
			return 0;
		}

#pragma acc routine seq
		unsigned int neighbour3D_num(int direction, unsigned int pos)
		{
//			unsigned int n= (unsigned int)(std::cbrt(dim));
//			unsigned int n= std::lrint(std::cbrt(this->dim));
//			unsigned int nsq, quo, rem;
			unsigned int quo, rem;

//			nsq=		n*n;
			quo=		pos/cu_root_sq;
			rem=		pos%cu_root_sq;

			switch(direction)
			{
				case 0:										//LEFT
					return (cu_root+(pos%cu_root)-1)%cu_root + (pos/cu_root)*cu_root;
				case 1:										//RIGHT
					return (cu_root+(pos%cu_root)+1)%cu_root + (pos/cu_root)*cu_root;
				case 2:										//TOP
					return (cu_root_sq+rem+cu_root)%cu_root_sq + quo*cu_root_sq;
				case 3:										//BOTTOM
					return (cu_root_sq+rem-cu_root)%cu_root_sq + quo*cu_root_sq;
				case 4:										//FRONT
					return (dim+pos-cu_root_sq)%dim;
				case 5:										//HIND
					return (dim+pos+cu_root_sq)%dim;
			}
			return 0;
		}

		Grid (unsigned int n)
		: dim {n}, sq_root{0}, cu_root{0}, cu_root_sq{0}
		{
			unsigned int c_sq_root = std::round(std::sqrt(this->dim));
			unsigned int c_cu_root = std::round(std::cbrt(this->dim));
			bool assigned=false;

			if((c_sq_root*c_sq_root)==this->dim)
			{
				sq_root=c_sq_root;
				assigned=true;
			}
			if((c_cu_root*c_cu_root*c_cu_root)==this->dim)
			{
				cu_root=c_cu_root;
				cu_root_sq=cu_root*cu_root;
				assigned=true;
			}

			if(!assigned)
			{
				std::cout<<"\nVector dimensions not compatible with a 2D or 3D grid.\nExiting...\n";
				exit(-1);
			}
		}

		Grid(const Grid&)			=default;
/*		Grid(const Grid& item)
		{
			std::cout<<"(Grid) Cpy ctor call for obj: "<<this<<'\n';

			this->dim = item.dim;
			this->sq_root = item.sq_root;
			this->cu_root = item.cu_root;
			this->cu_root_sq = item.cu_root_sq;
		}
*/
		Grid(Grid &&)				=default;
/*		Grid(Grid && item)
		{
			std::cout<<"(Grid) Mve ctor call for obj: "<<this<<'\n';

			this->dim = item.dim;
			this->sq_root = item.sq_root;
			this->cu_root = item.cu_root;
			this->cu_root_sq = item.cu_root_sq;

			item.dim=0;
			item.sq_root=0;
			item.cu_root=0;
			item.cu_root_sq=0;
		}
*/
		Grid &operator=(const Grid&)		=default;
		Grid &operator=(Grid &&)		=default;
		~Grid()					=default;

		public:

		virtual void operator() (const std::vector<T>& in, std::vector<T>& out) = 0;
	};

	template <typename T>
	class Grid2D: public Grid<T>
	{
		public:

		Grid2D(unsigned int n)
		: Grid<T>(n)
		{
			assert(this->sq_root!=0);
		}

		void operator() (const std::vector<T>& in, std::vector<T>& out) final
		{
			const T *in_ptr= &in[0];
			T *out_ptr= &out[0];

#ifdef PARALLEL_1
#pragma acc parallel loop copyin(in_ptr[0:dim]) copyout(out_ptr[0:dim])
#endif
			for(unsigned int i = 0; i < this->dim; ++i) 
			{
				out_ptr[i]= in_ptr[i] - 0.25*(
						in_ptr[this->neighbour2D('l', i)] +
						in_ptr[this->neighbour2D('r', i)] +
						in_ptr[this->neighbour2D('t', i)] +
						in_ptr[this->neighbour2D('b', i)]);
			}
		}

		Grid2D(const Grid2D&)			=default;
		Grid2D(Grid2D &&)			=default;
		Grid2D &operator=(const Grid2D&)	=default;
		Grid2D &operator=(Grid2D &&)		=default;
		~Grid2D()				=default;
	};

	template <typename T>
	class Grid3D: public Grid<T>
	{

		public:

		Grid3D(unsigned int n)
		: Grid<T>(n)
		{
			assert(this->cu_root!=0 && this->cu_root_sq!=0);
		}

		void operator() (const std::vector<T>& in, std::vector<T>& out) final
		{
			const T *in_ptr= &in[0];
			T *out_ptr= &out[0];

#ifdef PARALLEL_1
#pragma acc parallel loop copyin(in_ptr[0:dim]) copyout(out_ptr[0:dim])
#endif
			for(unsigned int i = 0; i < this->dim; ++i)
			{
				out_ptr[i]= in_ptr[i] - (1.0/6.0)*(
						in_ptr[this->neighbour3D('l', i)] +
						in_ptr[this->neighbour3D('r', i)] +
						in_ptr[this->neighbour3D('t', i)] +
						in_ptr[this->neighbour3D('b', i)] +
						in_ptr[this->neighbour3D('f', i)] +
						in_ptr[this->neighbour3D('h', i)]);
			}
		}

		Grid3D(const Grid3D&)			=default;
		Grid3D(Grid3D &&)			=default;
		Grid3D &operator=(const Grid3D&)	=default;
		Grid3D &operator=(Grid3D &&)		=default;
		~Grid3D()				=default;
	};

#ifdef COMPLEX3D
//MAIN-COMPUTE
	template <typename T>
	class Grid3D<Complex3D<T>>: public Grid<Complex3D<T>>
	{
		std::complex<T> *edges_ptr;

		public:

		Grid3D<Complex3D<T>>(unsigned int n)
		: Grid<Complex3D<T>>(n), edges_ptr(nullptr)
		{
			assert(this->cu_root!=0 && this->cu_root_sq!=0);

			edges_ptr = new std::complex<T>[3 * this->dim * 9];
			for(unsigned int i=0; i<(this->dim*3); ++i)
				for(unsigned int k=0; k<3; ++k)
					edges_ptr[i*9 + k*3+k]= std::complex<T>(1.0, 0.0);
#if defined(PARALLEL_1) || defined(PARALLEL_2)
#pragma acc enter data copyin(this)
#pragma	acc enter data copyin(this->edges_ptr[0:(3*this->dim*9)])
#endif
		}

		void operator()(const std::vector<Complex3D<T>>& in, std::vector<Complex3D<T>>& out) final
		{
			const Complex3D<T> *in_ptr= &in[0];
			Complex3D<T> *out_ptr= &out[0];

#ifdef PARALLEL_1
#pragma acc parallel loop collapse(2)								\
copyin(in_ptr[0:dim]) copyout(out_ptr[0:dim]) present(this->edges_ptr[0:(3*this->dim*9)])	\
private(j, l, direction, ref_node, row)
#endif
#ifdef PARALLEL_2
#pragma acc parallel loop collapse(2) 								\
present(in_ptr[0:dim], out_ptr[0:dim], this->edges_ptr[0:(3*this->dim*9)])			\
private(j, l, direction, ref_node, row) async(2)
#endif
			for(unsigned int i = 0; i < this->dim; ++i)			//VECTOR LENGTH
				for(unsigned int k=0; k < 3; ++k)			//MATRIX ROW
				{
					out_ptr[i].colour[k]=std::complex<T>(0.0, 0.0);
#if defined(PARALLEL_1) || defined(PARALLEL_2)
#pragma acc loop seq
#endif
					for(unsigned int j=0; j < 6; ++j)		//DIRECTIONS
					{
						unsigned int direction 	= (j/2)*this->dim*9;
						unsigned int ref_node 	= (((j+1)%2) 	* i 	+
									  (j%2) 	* this->neighbour3D_num(j, i))
									  * 9;
						unsigned int row	= k*3;

						for(unsigned int l=0; l<3; ++l)		//MATRIX COLUMN
							out_ptr[i].colour[k] += (std::complex<T>((j+1)%2)				  *
										edges_ptr[direction + ref_node + row + l]) 		  *
										(in_ptr[this->neighbour3D_num(j, i)].colour[l])  	  +
										(std::complex<T>((j%2))					  *
										std::conj(edges_ptr[direction + ref_node + l*3 + row/3])) *
										(in_ptr[this->neighbour3D_num(j, i)].colour[l]);
					}
					out_ptr[i].colour[k]= in_ptr[i].colour[k] - (1.0/6.0)*out_ptr[i].colour[k];
				}
		}

		Grid3D<Complex3D<T>>(const Grid3D<Complex3D<T>>&)			=delete; 	//To prevent copy contruction
/*		Grid3D<Complex3D<T>>(const Grid3D<Complex3D<T>>& item)
		: Grid<Complex3D<T>>(item)
		{
			std::cout<<"(Grid3D<Complex3D>) Cpy ctor call for obj: "<<this<<", Src: "<<&item<<'\n';

			auto complex_size=sizeof(std::complex<T>);
			this->edges_ptr = new std::complex<T>[3 * this->dim * 9];

			memcpy(this->edges_ptr, item.edges_ptr, complex_size * 3 * this->dim * 9);
#ifdef PARALLEL_1
#pragma	acc exit data delete(item.edges_ptr[0:3*dim*9])
#pragma acc exit data delete(item)
#pragma acc enter data copyin(this)
#pragma	acc enter data copyin(this->edges_ptr[0:(3*this->dim*9)])
#endif

		}
*/
		Grid3D<Complex3D<T>>(Grid3D<Complex3D<T>> &&)				=delete; 	//To prevent move contruction
/*		Grid3D<Complex3D<T>>(Grid3D<Complex3D<T>> && item)
		: Grid<Complex3D<T>>(std::move(item))
		{
			std::cout<<"(Grid3D<Complex3D>) Mve ctor call for obj: "<<this<<", Src: "<<&item<<'\n';

		 	edges_ptr=item.edges_ptr;
			item.edges_ptr=nullptr;
#ifdef PARALLEL_1
#pragma	acc exit data delete(item.edges_ptr[0:3*dim*9])
#pragma acc exit data delete(item)
#pragma acc enter data copyin(this)
#pragma	acc enter data copyin(this->edges_ptr[0:(3*this->dim*9)])
#endif
		}
*/
		Grid3D<Complex3D<T>> &operator=(const Grid3D<Complex3D<T>> &)		=delete;	//To prevent copy assignment
		Grid3D<Complex3D<T>> &operator=(Grid3D<Complex3D<T>> &&)		=delete;	//To prevent move assignment

		~Grid3D<Complex3D<T>>()
		{
//			std::cout<<"(Grid3D<Complex3D>) Dtor call for obj: "<<this<<'\n';
#if defined(PARALLEL_1) || defined(PARALLEL_2)
#pragma	acc exit data delete(this->edges_ptr[0:3*dim*9])
#pragma acc exit data delete(this)
#endif
			if(this->edges_ptr!=nullptr)
			{
				delete []this->edges_ptr;
				this->edges_ptr=nullptr;
			}
		}

		void read_input(const std::string &inp_name)
		{
			std::ifstream inp_file;
			inp_file.open(inp_name.c_str(), std::ios::in);

			unsigned int vlen;
			char sep;
			int node, direction;

			inp_file>>vlen;

			if(vlen == this->dim)
			{
				for(unsigned int i=0; i<vlen; i++)			//VECTOR LENGTH
					for(int j=0; j<3; j++)				//DIRECTIONS
					{
	
						inp_file>>node>>sep>>direction;
						for(int k=0; k<3; k++)			//MATRIX ROW
							for(int l=0; l<3; l++)		//MATRIX COLUMN
								inp_file>>std::setprecision(15)
									>>edges_ptr[	(direction*vlen*9) +
											(node*9) +
											(k*3)+
											l				];
					}
#if defined(PARALLEL_1) || defined(PARALLEL_2)
#pragma	acc update device(this->edges_ptr[0:(3*this->dim*9)])
#endif
			}
			else
			{
				std::cout<<"Dimensions in file don't match with the vlen. Not reading...\n";
				std::cout<<"vlen: "<<vlen<<", dim: "<<this->dim<<'\n';
			}

			inp_file.close();

		}
	};
#endif		//COMPLEX3D
}
#endif		//_GRID_HPP_
