// Copyright (C) 2012 Ivan B. Ivanov. All rights reserved.

// Contact author: www.ivi.com or navi.adler@gmail.com.

// This file is part of ZeroSolver C++ library.

// ZeroSolver is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// ZeroSolver is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with ZeroSolver. If not, see <http://www.gnu.org/licenses/>.

/*! \file
    \brief The declaration and definition of Vertex class.
*/

#ifndef FZF_VERTEX_HPP
#define FZF_VERTEX_HPP

#include <cmath>
#include <complex>
#include <iostream>

/********************************************************************/

namespace type
{
/********************************************************************/

template <typename T = double> class Vertex
{
    typedef std::complex<T> Complex;

    public:

        Vertex (const Complex & Z, const Complex & F) : z(Z), f(F), a(arg(F))
        {
        }

        ~Vertex (void)
        {
        }

        friend std::ostream & operator << (std::ostream & output, const Vertex<T> & V)
        {
            output << std::scientific     << std::setprecision(16)
                   << real(V.z) << "    " << imag(V.z) << "\t"
                   << real(V.f) << "    " << imag(V.f) << "\t"
                   << V.a;
            return output;
        }

        friend std::ostream & operator << (std::ostream & output, const Vertex<T> * V)
        {
            output << std::scientific      << std::setprecision(16)
                   << real(V->z) << "    " << imag(V->z) << "\t"
                   << real(V->f) << "    " << imag(V->f) << "\t"
                   << V->a;
            return output;
        }

        const Complex & get_z (void) const { return z; }
        const Complex & get_f (void) const { return f; }
        const       T & get_a (void) const { return a; }

    private:

        const Complex z;
        const Complex f;
        const T       a;
};

/********************************************************************/
} //end of namespace

/********************************************************************/

#endif
