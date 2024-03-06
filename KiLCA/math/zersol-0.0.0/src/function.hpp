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
    \brief The declaration and definition of Function class.
*/

#ifndef FZF_FUNCTION_HPP
#define FZF_FUNCTION_HPP

#include <complex>
#include <cstring>
#include <limits>

/*!
\brief The namespace contains the basic types needed for ZerSol algorithms.
*/
namespace type
{
/********************************************************************/
/*!
\brief The class implements properties of a function of one complex variable.
*/
template <typename T = double> class Function
{
        typedef std::complex<T> Complex; //!<typedef for complex numbers

    public:

        typedef Complex (*cmplx_func) (const Complex &, void *); //!<typedef for a function of one complex argument (if needed, additional data are passed via void * pointer)

        /*!
         *  \brief A constructor of the Function object.
         *  \param function name of a C++ function that evaluates the function value
         *  \param parameter pointer to a C++ structure that can be used to pass any data into the function and its derivative
         *  \param description string that indicates the function purpose
         */
        Function (cmplx_func function, void * parameter = 0, const std::string & description = "unknown")
                 :
                 f(function), df(0),
                 p(parameter),
                 name(description),
                 n_feval(0), n_deval(0),
                 method(0), h(pow(std::numeric_limits<T>::epsilon(), 1.0/3.0))
        {
        }

        /*!
         *  \brief Another constructor of the Function object.
         *  \param function    name of a C++ function that evaluates the function value
         *  \param dfunction   name of a C++ function that evaluates the function derivative
         *  \param parameter   pointer to a C++ structure that can be used to pass any data into the function and its derivative
         *  \param description string that indicates the function purpose
         */
        Function (cmplx_func function, cmplx_func dfunction, void * parameter = 0, const std::string & description = "unknown")
                 :
                 f(function), df(dfunction),
                 p(parameter),
                 name(description),
                 n_feval(0), n_deval(0),
                 method(0), h(pow(std::numeric_limits<T>::epsilon(), 1.0/3.0))
        {
        }

        /*!
         *  \brief A destructor of the Function object.
         */
        ~Function (void)
        {
        }

        /*!
         *  \brief An insertion operator to print content of the object.
         *  \param output an output stream for the insertion
         *  \param F a Function object to be inserted (printed)
         *  \return output stream
         */
        friend std::ostream & operator << (std::ostream & output, const Function<T> & F)
        {
            output << "\nComplex function info:"
                   << "\nname="     << F.name
                   << "\nfeval="    << F.n_feval
                   << "\ndeval="    << F.n_deval
                   << "\nndmethod=" << F.method
                   << "\nndstep="   << F.h
                   << std::endl;

            return output;
        }

        /*!
         *  \brief Set a finite difference method for numerical evaluation of the function derivative.
         *  \param method_ finite difference method: pass method_ = 0 for 3-point stencil approximation (second order) and method_ = 1 for 5-point stencil approximation (fourth order)
         */
        void set_nd_method (int method_) { method = method_; }

        /*!
         *  \brief Set a step size (dz) for numerical evaluation of the function derivative.
         *  \param h_ step size
         */
        void set_nd_step (const T & h_) { h = h_; }

        /*!
         *  \brief Get a finite difference method used for numerical evaluation of the function derivative.
         *  \return finite difference method: 0 for 3-point stencil approximation (second order) and 1 for 5-point stencil approximation (fourth order)
         */
        const int & get_nd_method (void) const { return method; }

        /*!
         *  \brief Get a step size (dz) used for numerical evaluation of the function derivative.
         *  \return step size
         */
        const T   & get_nd_step   (void) const { return h; }

        /*!
         *  \brief Get a pointer to the function data.
         *  \return data pointer
         */
        void * get_data_ptr (void) const { return p; }

        /*!
         *  \brief Evaluate the function value.
         *  \param z independent variable
         *  \return  function value
         */
        Complex eval_value (const Complex & z)           const { ++n_feval; return f(z, (void *)p); }

        /*!
         *  \brief Evaluate the function value for the q data.
         *  \param z independent variable
         *  \param q pointer to a C++ structure that can be used to pass any data into the function
         *  \return  function value
         */
        Complex eval_value (const Complex & z, void * q) const { ++n_feval; return f(z, (void *)q); }

        /*!
         *  \brief Evaluate the function derivative either by an analytic expression (if it was provided in constructor) or by a numerical approximation.
         *  \param z independent variable
         *  \return  function derivative
         */
        Complex eval_deriv (const Complex & z) const
        {
            ++n_deval;
            if (df) return        df(z, (void *)p);
            else    return num_deriv(z, (void *)p);
        }

        /*!
         *  \brief Evaluate the function derivative either by an analytic expression (if it was provided in constructor) or by a numerical approximation.
         *  \param z independent variable
         *  \param q pointer to a C++ structure that can be used to pass any data into the function
         *  \return  function value
         */
        Complex eval_deriv (const Complex & z, void * q) const
        {
            ++n_deval;
            if (df) return        df(z, q);
            else    return num_deriv(z, q);
        }

        /*!
         *  \brief Evaluate the function derivative by a numerical approximation.
         *  \param z independent variable
         *  \param s pointer to a C++ structure that can be used to pass any data into the function
         *  \return  function value
         */
        Complex num_deriv (const Complex & z, void * s) const
        {
            const Complex dz(h, h);

            Complex dFdz(0.0, 0.0);

            if (method == 0)
            {
                Complex fm  = eval_value(z - dz, s);
                Complex fp  = eval_value(z + dz, s);

                dFdz = (fp - fm) / (T(2.0)*dz);
            }
            else if (method > 0)
            {
                Complex fm  = eval_value(z - dz, s);
                Complex fp  = eval_value(z + dz, s);
                Complex fmm = eval_value(z - T(2.0)*dz, s);
                Complex fpp = eval_value(z + T(2.0)*dz, s);

                dFdz = ( fmm - fpp + T(8.0) * (fp - fm) ) / (T(12.0)*dz);
            }
            return dFdz;
        }

    private:

        cmplx_func f;         //!< pointer to a C++ function that evaluates the function value
        cmplx_func df;        //!< pointer to a C++ function that evaluates the function derivative

        void * p;             //!< pointer to a C++ structure that can be used to pass any data into the function and its derivative

        std::string name;     //!< string that indicates the function purpose

        mutable int n_feval;  //!< counter for number of the function evaluations
        mutable int n_deval;  //!< counter for number of the derivative evaluations

        int method;           //!< finite difference method for numerical evaluation of the function derivative
        T   h;                //!< step size for numerical evaluation of the function derivative
};

/********************************************************************/
} //end of namespace

#endif
