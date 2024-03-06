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
    \brief The declaration and definition of the Newton's iterations on 1D complex plane.
*/

#ifndef ALGO_NEWTON_HPP
#define ALGO_NEWTON_HPP

#include <cmath>
#include <complex>
#include <iostream>
#include <limits>

#include "box.hpp"
#include "function.hpp"
#include "vertex.hpp"
#include "list.hpp"

/********************************************************************/

namespace alg
{
/*!
\brief The namespace contains classes and algorithms for the Newton's iterations.
*/
namespace newton
{
using namespace type;

/********************************************************************/
/*!
\brief The class contains settings for the Newton's iterations.
*/
template <typename T = double> class Settings
{
    public:

        typedef std::complex<T> Complex; //!<typedef for complex numbers

        /*!
         *  \brief A constructor of the Settings object.
         *  \param n_target_     desired number of the function zeros
         *  \param n_start_      number of a user supplied starting values
         *  \param start_        a user supplied starting values
         *  \param n_split_x_    number of automatic starting points along Re{z}
         *  \param n_split_y_    number of automatic starting points along Im{z}
         *  \param max_iter_num_ maximum number of the Newton's iterations for each starting point
         *  \param multiplicity_ expected multiplicity of a zero
         *  \param abs_eps_z_    absolute tolerance for the convergence condition on z
         *  \param rel_eps_z_    relative tolerance for the convergence condition on z
         *  \param abs_eps_f_    absolute tolerance for the convergence condition on f(z)
         *  \param rel_eps_f_    relative tolerance for the convergence condition on f(z)
         *  \param debug_level_  debug level (determines the amount of internal checks)
         *  \param print_level_  print level (determines the amount of details printed)
         */
        Settings (int           n_target_ = std::numeric_limits<int>::max(),
                  int            n_start_ = 0,
                  Complex *        start_ = 0,
                  int          n_split_x_ = 4,
                  int          n_split_y_ = 4,
                  int       max_iter_num_ = 16,
                  int       multiplicity_ = 1,
                  T            abs_eps_z_ = 2.0 * std::numeric_limits<T>::epsilon(),
                  T            rel_eps_z_ = 1.0e-6,
                  T            abs_eps_f_ = 2.0 * std::numeric_limits<T>::epsilon(),
                  T            rel_eps_f_ = 1.0e-12,
                  int        debug_level_ = 0,
                  int        print_level_ = 2
                 )
                 :
                         n_target(n_target_),
                           n_start(n_start_),
                               start(start_),
                       n_split_x(n_split_x_),
                       n_split_y(n_split_y_),
                 max_iter_num(max_iter_num_),
                 multiplicity(multiplicity_),
                       abs_eps_z(abs_eps_z_),
                       rel_eps_z(rel_eps_z_),
                       abs_eps_f(abs_eps_f_),
                       rel_eps_f(rel_eps_f_),
                   debug_level(debug_level_),
                   print_level(print_level_)
        {
        }

        /*!
         *  \brief A destructor of the Settings object.
         */
        ~Settings (void)
        {
        }

        friend std::ostream & operator << (std::ostream & output, const Settings<T> & S)
        {
            output << "\nThe Newton iterations settings:\n"
                   << "n_target="     << S.n_target     << "\n"
                   << "n_start="      << S.n_start      << "\n"
                   << "n_split_x="    << S.n_split_x    << "\n"
                   << "n_split_y="    << S.n_split_y    << "\n"
                   << "max_iter_num=" << S.max_iter_num << "\n"
                   << "multiplicity=" << S.multiplicity << "\n"
                   << "abs_eps_z="    << S.abs_eps_z    << "\n"
                   << "rel_eps_z="    << S.rel_eps_z    << "\n"
                   << "abs_eps_f="    << S.abs_eps_f    << "\n"
                   << "rel_eps_f="    << S.rel_eps_f    << "\n"
                   << "debug_level="  << S.debug_level  << "\n"
                   << "print_level="  << S.print_level  << std::endl;

            return output;
        }

        //set functions:
        /*!
         *  \brief Set a desired number of the function zeros.
         *  \param n_target_ target number of zeros
         */
        void set_n_target (const int & n_target_) const
        {
            n_target = n_target_;
        }

        /*!
         *  \brief Set the array of a user supplied starting values for the Newton's iterations.
         *  \param n_start_ number of starting values
         *  \param start_   array of starting values
         */
        void set_start_array (const int & n_start_, Complex * start_)
        {
            n_start = n_start_;
            start   = start_;
        }

        /*!
         *  \brief Set the number of automatic starting points along Re{z}.
         *  \param n_split_x_ number of starting points
         */
        void set_n_split_x (const int & n_split_x_) const
        {
            n_split_x = n_split_x_;
        }

        /*!
         *  \brief Set the number of automatic starting points along Im{z}.
         *  \param n_split_y_ number of starting points
         */
        void set_n_split_y (const int & n_split_y_) const
        {
            n_split_y = n_split_y_;
        }

        /*!
         *  \brief Set the maximum number of the Newton's iterations for each starting point.
         *  \param max_iter_num_ maximum number of iterations
         */
        void set_max_iter_num (const int & max_iter_num_)
        {
            max_iter_num = max_iter_num_;
        }

        /*!
         *  \brief Set the multiplicity of a zero.
         *  \param multiplicity_ expected multiplicity of a zero
         */
        void set_multiplicity (const int & multiplicity_) const
        {
            multiplicity = multiplicity_;
        }

        /*!
         *  \brief Set the absolute and relative tolerances for the convergence condition on z.
         *  \param abs_eps_z_ absolute tolerance
         *  \param rel_eps_z_ relative tolerance
         */
        void set_eps_for_arg (const T & abs_eps_z_, const T & rel_eps_z_)
        {
            abs_eps_z = abs_eps_z_;
            rel_eps_z = rel_eps_z_;
        }

        /*!
         *  \brief Set the absolute and relative tolerances for the convergence condition on f(z).
         *  \param abs_eps_f_ absolute tolerance
         *  \param rel_eps_f_ relative tolerance
         */
        void set_eps_for_func (const T & abs_eps_f_, const T & rel_eps_f_)
        {
            abs_eps_f = abs_eps_f_;
            rel_eps_f = rel_eps_f_;
        }

        void set_debug_level (const int & debug_level_)
        {
            debug_level = debug_level_;
        }

        void set_print_level (const int & print_level_)
        {
            print_level = print_level_;
        }

        //get functions:
        /*!
         *  \brief Get the target number of the function zeros.
         *  \return target number of zeros
         */
        const int & get_n_target (void) const
        {
            return n_target;
        }

        /*!
         *  \brief Get the array of a user supplied starting values for the Newton's iterations.
         *  \param n_start_ number of starting values
         *  \param start_   array of starting values
         */
        void get_start_array (int & n_start_, Complex * & start_) const
        {
            n_start_ = n_start;
            start_   = start;
        }

        /*!
         *  \brief Get the number of automatic starting points along Re{z}.
         *  \return number of starting points
         */
        const int & get_n_split_x (void) const
        {
            return n_split_x;
        }

        /*!
         *  \brief Get the number of automatic starting points along Im{z}.
         *  \return number of starting points
         */
        const int & get_n_split_y (void) const
        {
            return n_split_y;
        }

        /*!
         *  \brief Get the maximum number of the Newton's iterations for each starting point.
         *  \return maximum number of iterations
         */
        const int & get_max_iter_num (void) const
        {
            return max_iter_num;
        }

        /*!
         *  \brief Get the multiplicity of a zero.
         *  \return multiplicity of a zero
         */
        const int & get_multiplicity (void) const
        {
            return multiplicity;
        }

        /*!
         *  \brief Get the absolute and relative tolerances for the convergence condition on z.
         *  \param abs_eps_z_ absolute tolerance
         *  \param rel_eps_z_ relative tolerance
         */
        void get_eps_for_arg (T & abs_eps_z_, T & rel_eps_z_) const
        {
            abs_eps_z_ = abs_eps_z;
            rel_eps_z_ = rel_eps_z;
        }

        /*!
         *  \brief Get the absolute and relative tolerances for the convergence condition on f(z).
         *  \param abs_eps_f_ absolute tolerance
         *  \param rel_eps_f_ relative tolerance
         */
        void get_eps_for_func (T & abs_eps_f_, T & rel_eps_f_) const
        {
            abs_eps_f_ = abs_eps_f;
            rel_eps_f_ = rel_eps_f;
        }

        const int & get_debug_level (void) const
        {
            return debug_level;
        }

        const int & get_print_level (void) const
        {
            return print_level;
        }

    private:

        //The Newton iterations options:
mutable int       n_target;             //!< target number of the function zeros
        int       n_start;              //!< number of a user supplied starting values
        Complex * start;                //!< a user supplied starting values
mutable int       n_split_x;            //!< number of automatic starting points along Re{z}
mutable int       n_split_y;            //!< number of automatic starting points along Im{z}
        int       max_iter_num;         //!< maximum number of the Newton's iterations for each starting point
mutable int       multiplicity;         //!< multiplicity of a zero
        T         abs_eps_z, rel_eps_z; //!< absolute and relative tolerances for the convergence condition on z
        T         abs_eps_f, rel_eps_f; //!< absolute and relative tolerances for the convergence condition on f(z)
        int       debug_level;          //!< debug level (determines the amount of internal checks)
        int       print_level;          //!< print level (determines the amount of details printed)
};

/********************************************************************/

template <typename T> inline void
print_iteration_info (const int & iter,
                      const std::complex<T> & z_prev, const std::complex<T> & f_prev,
                      const std::complex<T> & z_next, const std::complex<T> & f_next)
{
std::clog << "\n" << "iter=" << iter << ": "
                  << "z_prev=" << z_prev << "\t" << "z_next=" << z_next << "\t"
                  << "f_prev=" << f_prev << "\t" << "f_next=" << f_next << std::endl;
}

/********************************************************************/

template <typename T> inline void
print_summary_info (const int & n_new_zeros, List<Vertex<T> *> * all_zeros)
{
std::clog << "\nThe Newton iterations summary:\n"
          << " n_new_zeros="  << n_new_zeros << std::endl;
}

/********************************************************************/

template <typename T> inline int
newton_iterations (const Function<T> & F,
                   const Box<T> & B,
                   const Settings<T> & S,
                   List< Vertex<T> *> * all_zeros
                  )
{
typedef std::complex<T> Complex;

//user supplied starting values:
int n_user;
Complex * user;
S.get_start_array(n_user, user);

List<Complex> start;
for (int k = 0; k < n_user; ++k) start.push_back(user[k]);

const T & xmin = B.get_xmin();
const T & xmax = B.get_xmax();
const T & ymin = B.get_ymin();
const T & ymax = B.get_ymax();

const int & n_split_x = S.get_n_split_x();
const int & n_split_y = S.get_n_split_y();

const T dx = (xmax - xmin) / n_split_x;
const T dy = (ymax - ymin) / n_split_y;
const T half = 0.5;

const Complex I(0.0e0, 1.0e0);
const Complex zlb(xmin, ymin);

for (int i = 0; i < n_split_x; ++i)
{
    for (int j = 0; j < n_split_y; ++j)
    {
        start.push_back(zlb + (i + half) * dx + (j + half) * dy * I);
    }
}

T abs_eps_z, rel_eps_z;
S.get_eps_for_arg(abs_eps_z, rel_eps_z);

T abs_eps_f, rel_eps_f;
S.get_eps_for_func(abs_eps_f, rel_eps_f);

const int & max_iter_num = S.get_max_iter_num();
const T   & multiplicity = S.get_multiplicity();
const int & n_target     = S.get_n_target();
const int & debug_level  = S.get_debug_level();

int n_new_zeros = 0;

List_Node<Complex> * guess = start.get_first();

while (guess)
{
    Complex z_prev = guess->get_data();
    Complex f_prev = F.eval_value(z_prev);
    Complex z_next;
    Complex f_next;

    const T mof = std::abs(f_prev); //modulo of function at the starting guess point

    int iter = -1;
    while (++iter < max_iter_num)
    {
        Complex df_prev = F.eval_deriv(z_prev);

        z_next = z_prev - multiplicity * f_prev / df_prev;

        //check if the point is still inside the rectangle:
        if (!(real(z_next) > xmin && real(z_next) < xmax && imag(z_next) > ymin && imag(z_next) < ymax)) break;

        f_next = F.eval_value(z_next);

        if (debug_level > 1) print_iteration_info(iter, z_prev, f_prev, z_next, f_next);

        //check convergence:
        if (std::abs(z_next - z_prev)  <  std::max(abs_eps_z, rel_eps_z * std::abs(z_next))  &&
                     std::abs(f_next)  <  std::max(abs_eps_f, rel_eps_f * mof)
           )
        {
            //check that z_next zero is new:
            bool is_new = true;

            List_Node< Vertex<T> *> * node = all_zeros->get_first();
            while (node)
            {
                if (std::abs(z_next - node->get_data()->get_z()) < std::max(abs_eps_z, rel_eps_z * std::abs(z_next)))
                {
                    is_new = false;
                    break;
                }
                node = node->get_next();
            }

            if (!is_new) break;

            n_new_zeros += multiplicity;

            for (int i = 0; i < multiplicity; ++i)
            {
                all_zeros->push_back(new Vertex<T>(z_next, f_next)); //add the new zero to the list of all zeros
            }

            break;
        }

        z_prev = z_next;
        f_prev = f_next;
    }

    if (n_new_zeros == n_target) break;

    guess = guess->get_next();
}

if (debug_level > 0) print_summary_info(n_new_zeros, all_zeros);

return n_new_zeros;
}

/********************************************************************/

} //end of namespace
} //end of namespace

/********************************************************************/

#endif
