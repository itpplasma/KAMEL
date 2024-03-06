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
    \brief The declaration and definition of Rectangle class.
*/

#ifndef FZF_RECTANGLE_HPP
#define FZF_RECTANGLE_HPP

#include <complex>
#include <iostream>
#include <fstream>

#include "list.hpp"
#include "vertex.hpp"
#include "edge.hpp"
#include "function.hpp"
#include "exceptions.hpp"

/********************************************************************/

namespace winding
{
/********************************************************************/

template <typename T = double> class Rectangle
{
    typedef Vertex<T>    * pVertex;
    typedef Edge<T>      * pEdge;
    typedef Rectangle<T> * pRectangle;

    public:

        Rectangle (pEdge E0_, pEdge E1_, pEdge E2_, pEdge E3_,
                   int split_,
                   List<pVertex> * zeros_,
                   const Function<T> & F_,
                   const Settings<T> & S_
                  )
        :
        E0(E0_), E1(E1_), E2(E2_), E3(E3_), MID(0),
        winding(-1024), split(split_), zeros(zeros_), F(F_), S(S_)
        {
        }

        ~Rectangle (void)
        {
            delete E0;
            delete E1;
            delete E2;
            delete E3;
            delete MID;
            delete zeros;
        }

        friend std::ostream & operator << (std::ostream & output, const Rectangle<T> & R)
        {
            output << "\n\nThe rectangle vertices:"
                   << "\n#\tRe{vertex}\t\t  Im{vertex})\t\t\tRe{f(vertex)}\t\t  Im{f(vertex)})\t\tArg{f(vertex)}"
                   << "\n" << 0 << "\t" << R.get_lb_vertex()
                   << "\n" << 1 << "\t" << R.get_rb_vertex()
                   << "\n" << 2 << "\t" << R.get_rt_vertex()
                   << "\n" << 3 << "\t" << R.get_lt_vertex();

            output << "\n\nThe rectangle winding number = " << R.winding << "\n";

            output << "\nThe list of zeros found up to the rectangle level:"
                   << "\n#\tRe{zero}\t\t  Im{zero})\t\t\tRe{f(zero)}\t\t  Im{f(zero)})\t\t\tArg{f(zero)}"
                   << *(R.zeros);

            output << "\nThe rectangle bottom edge:"; R.E0->write_vertices(output);
            output << "\nThe rectangle right edge:";  R.E1->write_vertices(output);
            output << "\nThe rectangle top edge:";    R.E2->write_vertices(output);
            output << "\nThe rectangle left edge:";   R.E3->write_vertices(output);

            return output;
        }

        friend std::ostream & operator << (std::ostream & output, const Rectangle<T> * R)
        {
            output << "\n\nThe rectangle vertices:"
                   << "\n#\tRe{vertex}\t\t  Im{vertex})\t\t\tRe{f(vertex)}\t\t  Im{f(vertex)})\t\tArg{f(vertex)}"
                   << "\n" << 0 << "\t" << R->get_lb_vertex()
                   << "\n" << 1 << "\t" << R->get_rb_vertex()
                   << "\n" << 2 << "\t" << R->get_rt_vertex()
                   << "\n" << 3 << "\t" << R->get_lt_vertex();

            output << "\n\nThe rectangle winding number = " << R->winding << "\n";

            output << "\nThe list of zeros found up to the rectangle level:"
                   << "\n#\tRe{zero}\t\t  Im{zero})\t\t\tRe{f(zero)}\t\t  Im{f(zero)})\t\t\tArg{f(zero)}"
                   << *(R->zeros);

            output << "\nThe rectangle bottom edge:"; R->E0->write_vertices(output);
            output << "\nThe rectangle right edge:";  R->E1->write_vertices(output);
            output << "\nThe rectangle top edge:";    R->E2->write_vertices(output);
            output << "\nThe rectangle left edge:";   R->E3->write_vertices(output);

            return output;
        }

        int calc_winding_number (void)
        {
            winding = E0->eval_delta_argument() +
                      E1->eval_delta_argument() -
                      E2->eval_delta_argument() -
                      E3->eval_delta_argument();

            return winding;
        }

        int get_winding_number (void) const
        {
            return winding;
        }

        List<pVertex> * get_zeros (void)
        {
            return zeros;
        }

        const Vertex<T> * get_lb_vertex (void) const
        {
            return E0->get_first_vertex();
        }

        const Vertex<T> * get_rb_vertex (void) const
        {
            return E0->get_last_vertex();
        }

        const Vertex<T> * get_rt_vertex (void) const
        {
            return E1->get_last_vertex();
        }

        const Vertex<T> * get_lt_vertex (void) const
        {
            return E3->get_last_vertex();
        }

        std::complex<T> get_center (void) const
        {
            return (E0->get_first_vertex()->get_z() + E1->get_last_vertex()->get_z()) / T(2.0);
        }

        T get_diagonal (void) const
        {
            return std::abs(E1->get_last_vertex()->get_z() - E0->get_first_vertex()->get_z());
        }

        void partition (pRectangle & R1, pRectangle & R2);

    private:

        pEdge E0;
        pEdge E1;
        pEdge E2;
        pEdge E3;
        pEdge MID;

        int winding;

        int split;

        List<pVertex> * zeros; //list of all zeros that have been found in the given rectangle up to this level of partition

        const Function<T> & F;
        const Settings<T> & S;
};

/********************************************************************/

template <typename T> inline void select_zeros (const Vertex<T> * V1, const Vertex<T> * V2,
                                                List<Vertex<T> *> * zeros, List<Vertex<T> *> * nulls);

/********************************************************************/

template <typename T> inline void Rectangle<T> :: partition (Rectangle<T> * & R1, Rectangle<T> * & R2)
{
pEdge           edge0 = 0;
pEdge           edge1 = 0;
pEdge           edge2 = 0;
pEdge           edge3 = 0;
List<pVertex> * nulls = 0;

if (split == 0)
{
    //rectangle 0:
    edge1 = E1->make_copy(1);
    edge3 = E3->make_copy(1);

    pVertex V1 = new Vertex<T>(*(E3->get_central_vertex()));
    pVertex V2 = new Vertex<T>(*(E1->get_central_vertex()));

    MID = new Edge<T>(V1, V2, F, S);

    edge0 =  E0->make_copy(0);
    edge2 = MID->make_copy(0);

    nulls = new List<pVertex>;

    select_zeros(edge0->get_first_vertex(), edge1->get_last_vertex(), zeros, nulls);

    R1 = new Rectangle<T>(edge0, edge1, edge2, edge3, 1, nulls, F, S);

    //rectangle 1:
    edge0 = MID->make_copy(0);
    edge1 =  E1->make_copy(2);
    edge2 =  E2->make_copy(0);
    edge3 =  E3->make_copy(2);

    nulls = new List<pVertex>;

    select_zeros(edge0->get_first_vertex(), edge1->get_last_vertex(), zeros, nulls);

    R2 = new Rectangle<T>(edge0, edge1, edge2, edge3, 1, nulls, F, S);
}
else
{
    //rectangle 0:
    edge0 = E0->make_copy(1);
    edge2 = E2->make_copy(1);

    pVertex V1 = new Vertex<T>(*(E0->get_central_vertex()));
    pVertex V2 = new Vertex<T>(*(E2->get_central_vertex()));

    MID = new Edge<T>(V1, V2, F, S);

    edge1 = MID->make_copy(0);
    edge3 =  E3->make_copy(0);

    nulls = new List<pVertex>;

    select_zeros(edge0->get_first_vertex(), edge1->get_last_vertex(), zeros, nulls);

    R1 = new Rectangle<T>(edge0, edge1, edge2, edge3, 0, nulls, F, S);

    //rectangle 1:
    edge0 =  E0->make_copy(2);
    edge1 =  E1->make_copy(0);
    edge2 =  E2->make_copy(2);
    edge3 = MID->make_copy(0);

    nulls = new List<pVertex>;

    select_zeros(edge0->get_first_vertex(), edge1->get_last_vertex(), zeros, nulls);

    R2 = new Rectangle<T>(edge0, edge1, edge2, edge3, 0, nulls, F, S);
}

//check that all known zeros are correctly distributed:
if (R1->get_zeros()->get_size() + R2->get_zeros()->get_size() != zeros->get_size())
{
    throw_Algorithm_Error_Exception("failed to distribute the zeros over the subrectangles");
}
}

/********************************************************************/

template <typename T> inline void select_zeros (const Vertex<T> * V1, const Vertex<T> * V2,
                                                List<Vertex<T> *> * zeros, List<Vertex<T> *> * nulls)
{
typedef std::complex<T> Complex;

const T xmin = real(V1->get_z());
const T ymin = imag(V1->get_z());
const T xmax = real(V2->get_z());
const T ymax = imag(V2->get_z());

List_Node<Vertex<T> *> * node = zeros->get_first();

while (node)
{
    const Complex & Z = node->get_data()->get_z();

    if (real(Z) > xmin && real(Z) < xmax && imag(Z) > ymin && imag(Z) < ymax) //check if Z lies inside the rectangle
    {
        nulls->push_back(node->get_data());
    }

    node = node->get_next();
}
}

/********************************************************************/
} //end of namespace

/********************************************************************/

#endif
