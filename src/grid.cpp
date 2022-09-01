#include "grid.hpp"
#include "lattice_qcd.hpp"

using grid::Grid;
using grid::Grid2D;
using grid::Grid3D;

#ifdef COMPLEX3D
using lattice_qcd::Complex3D;
#endif

template class Grid <double>;
template class Grid <std::complex<double>>;
template class Grid2D <double>;
template class Grid2D <std::complex<double>>;
template class Grid3D <double>;
template class Grid3D <std::complex<double>>;

#ifdef COMPLEX3D
template class Grid3D <Complex3D<double>>;
#endif

