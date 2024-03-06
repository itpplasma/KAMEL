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
    \brief The declaration and definition of Interval class.
*/

#ifndef FZF_INTERVAL_HPP
#define FZF_INTERVAL_HPP

#include <iostream>

#include "vertex.hpp"
#include "list.hpp"

/********************************************************************/

namespace winding
{
/********************************************************************/

using namespace type;

/********************************************************************/

template <typename T = double> class Interval
{
    typedef Vertex<T> * pVertex;

    public:

        Interval (List_Node<pVertex> * Nl, List_Node<pVertex> * Nr, int Delta = -1024)
        :
        left(Nl), right(Nr), delta(Delta)
        {
        }

        ~Interval (void)
        {
        }

        List_Node<pVertex> * get_left  (void)       { return left;  }
        List_Node<pVertex> * get_right (void)       { return right; }
                       int   get_delta (void) const { return delta; }

        void set_delta (int Delta) { delta = Delta; }

        friend std::ostream & operator << (std::ostream & output, const Interval<T> & I)
        {
            //output << *(I.left) << "\t" << *(I.right) << "\t"<< I.delta;
            output << std::scientific << I.delta;
            return output;
        }

    private:

        List_Node<pVertex> * left;
        List_Node<pVertex> * right;

        int delta;
};

/********************************************************************/
} //end of namespace

/********************************************************************/

#endif
