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
    \brief The declaration and definition of Box class.
*/

#ifndef TYPE_BOX_H
#define TYPE_BOX_H

/********************************************************************/

#include <iostream>
#include <fstream>

/********************************************************************/

namespace type
{
/********************************************************************/
/*!
\brief The class describes a 2D rectangular box.
*/
template <typename T = double> class Box
{
    public:

        /*!
         *  \brief A constructor of the Box object.
         *  \param xmin_ minimum x-coordinate of the box
         *  \param xmax_ maximum x-coordinate of the box
         *  \param ymin_ minimum y-coordinate of the box
         *  \param ymax_ maximum y-coordinate of the box
         */
        Box (const T & xmin_ = -1.0, const T & xmax_ = 1.0,
             const T & ymin_ = -1.0, const T & ymax_ = 1.0
            )
            :
            xmin(xmin_), xmax(xmax_),
            ymin(ymin_), ymax(ymax_)
        {
        }

        /*!
         *  \brief A destructor of the Box object.
         */
        ~Box (void)
        {
        }

        /*!
         *  \brief An insertion operator to print content of the object.
         *  \param output an output stream for the insertion
         *  \param B a Box object to be inserted (printed)
         *  \return output stream
         */
        friend std::ostream & operator << (std::ostream & output, const Box<T> & B)
        {
            output << "\nBox parameters:"
                   << "\nxmin=" << B.xmin
                   << "\nxmax=" << B.xmax
                   << "\nymin=" << B.ymin
                   << "\nymax=" << B.ymax
                   << std::endl;

            return output;
        }

        const T & get_xmin (void) const { return xmin; }
        const T & get_xmax (void) const { return xmax; }
        const T & get_ymin (void) const { return ymin; }
        const T & get_ymax (void) const { return ymax; }

    private:

        T xmin, xmax; //!< min and max x-coordinates
        T ymin, ymax; //!< min and max y-coordinates
};

/********************************************************************/
} //end of namespace

#endif
