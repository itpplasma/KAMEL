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
    \brief The declaration and definition of main ZeroSolver class.
*/

#ifndef FZF_ZEROSOLVER_HPP
#define FZF_ZEROSOLVER_HPP

#include <complex>
#include <iostream>
#include <iomanip>

#include "box.hpp"
#include "tree.hpp"
#include "list.hpp"
#include "function.hpp"
#include "rectangle.hpp"
#include "newton.hpp"
#include "exceptions.hpp"

/********************************************************************/
/*!
\brief The namespace contains algorithms frequently used in scientific computations.
*/
namespace alg
{
/*!
\brief The namespace contains classes and algorithms for the zeros solver.
*/
namespace zersol
{
/********************************************************************/

using namespace type;
using namespace winding;
using namespace newton;

/********************************************************************/
/*!
\brief The class contains settings used by the zeros search algorithm.
*/
template <typename T = double> class Settings : public winding::Settings<T>, public newton::Settings<T>
{
    public:

        /*!
         *  \brief A constructor of the Settings object.
         *  \param use_winding_    define whether the solver will use the winding number evaluation or just simple Newton iterations
         *  \param max_part_level_ maximum allowed level of the rectangle partition
         *  \param debug_level_    debug level (determines the amount of internal checks)
         *  \param print_level_    print level (determines the amount of details printed)
         */
        Settings (bool use_winding_ = true, int max_part_level_ = 128, int debug_level_ = 0, int print_level_ = 2)
        :
        use_winding(use_winding_),
        max_part_level(max_part_level_),
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

        /*!
         *  \brief An insertion operator to print content of the object.
         *  \param output output stream for the insertion
         *  \param S Settings object to be inserted (printed)
         *  \return output stream
         */
        friend std::ostream & operator << (std::ostream & output, const Settings<T> & S)
        {
            output << "\nZero solver settings:\n"
                   << "use_winding="             << S.use_winding    << "\n"
                   << "max_partition_level="     << S.max_part_level << "\n"
                   << "debug_level="             << S.debug_level    << "\n"
                   << "print_level="             << S.print_level    << std::endl;

            output << *((const winding::Settings<T> *) &S);
            output << *((const  newton::Settings<T> *) &S);

            return output;
        }

        //set functions for recursive partitioning:
        /*!
         *  \brief Set a flag that defines whether the solver will use the winding number evaluation or just Newton iterations.
         *  \param use_winding_ flag: true - use the winding number evaluation, false - use the Newton's iterations only
         */
        void set_use_winding (const bool & use_winding_)
        {
            use_winding = use_winding_;
        }

        /*!
         *  \brief Set the maximum allowed level of the rectangle partition.
         *  \param max_part_level_ partition level
         */
        void set_max_part_level (const int & max_part_level_)
        {
            max_part_level = max_part_level_;
        }

        /*!
         *  \brief Set the debug level (determines the amount of internal checks).
         *  \param debug_level_ debug level
         */
        void set_debug_level (const int & debug_level_)
        {
            winding::Settings<T>::set_debug_level(debug_level_);
             newton::Settings<T>::set_debug_level(debug_level_);
              debug_level = debug_level_;
        }

        /*!
         *  \brief Set the print level (determines the amount of details printed).
         *  \param print_level_ print level
         */
        void set_print_level (const int & print_level_)
        {
            winding::Settings<T>::set_print_level(print_level_);
             newton::Settings<T>::set_print_level(print_level_);
              print_level = print_level_;
        }

        //get functions for recursive partitioning:
        /*!
         *  \brief Get a flag that defines whether the solver will use the winding number evaluation or just Newton iterations.
         *  \return flag: true - use the winding number evaluation, false - use the Newton's iterations only
         */
        const bool & get_use_winding (void) const
        {
            return use_winding;
        }

        /*!
         *  \brief Get the maximum allowed level of the rectangle partition.
         *  \return partition level
         */
        const int & get_max_part_level (void) const
        {
            return max_part_level;
        }

        /*!
         *  \brief Get the debug level (determines the amount of internal checks).
         *  \return debug level
         */
        const int & get_debug_level (void) const
        {
            return debug_level;
        }

        /*!
         *  \brief Get the print level (determines the amount of details printed).
         *  \return print level
         */
        const int & get_print_level (void) const
        {
            return print_level;
        }

    private:

        bool     use_winding; //!< define whether the solver will use the winding number evaluation or just simple Newton iterations
        int   max_part_level; //!< maximum allowed level of the rectangle partition
        int      debug_level; //!< debug level (determines the amount of internal checks)
        int      print_level; //!< print level (determines the amount of details printed)
};

/********************************************************************/

/*!
\brief The class implements the search algorithm for all zeros of a complex analytic function contained within a given rectangle.
*/
template <typename T = double> class ZeroSolver
{
        typedef std::complex<T> Complex; //!<typedef for complex numbers
        typedef Vertex<T>    *  pVertex;
        typedef Edge<T>      *  pEdge;
        typedef Rectangle<T> *  pRectangle;

    public:

        /*!
         *  \brief A constructor of the ZeroSolver object.
         *  \param F_ a type::Function object that defines the analytic function
         *  \param B_ a type::Box object that defines the rectangular region of a complex plane where the zeros are searched
         *  \param S_ a alg::zersol::Settings object that defines all settings used by the underlying algorithms
         */
        ZeroSolver (const Function<T> & F_, const Box<T> & B_, const Settings<T> & S_)
        :
        F(F_), B(B_), S(S_), is_found(false), status(0), zeros(0), partition(0), average(0)
        {
        }

        /*!
         *  \brief A destructor of the ZeroSolver object.
         */
        ~ZeroSolver (void)
        {
            if (zeros) zeros->destruct_data();
            if (partition) partition->destruct_data();

            delete status;
            delete zeros;
            delete partition;
        }

        /*!
         *  \brief An insertion operator to print content of the object.
         *  \param output an output stream for the insertion
         *  \param ZS a ZeroSolver object to be inserted (printed)
         *  \return output stream
         */
        friend std::ostream & operator << (std::ostream & output, const ZeroSolver<T> & ZS)
        {
            if (!ZS.S.get_print_level()) return output;

            output << "This is a dump of information used and generated during zeros search procedure:\n";

            output << ZS.F;
            output << ZS.B;
            output << ZS.S;

            if (ZS.status) ZS.print_status(output);

            if (ZS.zeros)
            {
                output << "\nThe list of zeros found in the given rectangular region:"
                       << "\n#\tRe{zero}\t\t  Im{zero})\t\t\tRe{f(zero)}\t\t  Im{f(zero)})\t\t\tArg{f(zero)}"
                       << *(ZS.zeros);
            }

            if (ZS.partition && ZS.S.get_print_level() > 1)
            {
                output << "\nThe region partition generated during the search of zeros:\n" << *(ZS.partition);
            }

            return output;
        }

        /*!
         *  \brief Print the status of the zeros search.
         *  \param output an output stream for the printing
         */
        void print_status (std::ostream & output = std::cout) const
        {
            output << "\nThe solver status: " << zeros->get_size() << " zeros have been found.\n";

            List_Node<std::string> * node = status->get_first();

            while (node)
            {
                output << node->get_data() << "\n";
                node = node->get_next();
            }
        }

        /*!
         *  \brief Run the zeros search procedure.
         *  \param N       length of the passed Z and V arrays
         *  \param Z       array for the function zeros
         *  \param V       array for the function values V = F(Z)
         *  \param n_zeros number of the zeros found
         *  \return        search status
         */
        int FindZeros (const int & N, Complex * Z, Complex * V, int & n_zeros);

    private:

        int find_zeros_in_rectangle (Tree_Node<pRectangle> * node);

        int copy_zeros (const int & N, Complex * Z, Complex * V, int & n_zeros);

        void add_to_status (const std::string & message)
        {
            status->push_back(message);
        }

    private:

        const Function<T>   F;
        const Box<T>        B;
        const Settings<T>   S;

        bool                is_found;

        List<std::string> * status;

        List<pVertex>     * zeros;

        Tree<pRectangle>  * partition;

        T                   average;
};

/********************************************************************/

template <typename T> int ZeroSolver<T> :: FindZeros (const int & N,
                                                      std::complex<T> * Z,
                                                      std::complex<T> * V,
                                                      int & n_zeros)
{
if (is_found) return copy_zeros(N, Z, V, n_zeros);

try
{
    status = new List<std::string>;

    zeros = new List<pVertex>;

    const winding::Settings<T> & WS = *((const winding::Settings<T> *) &S);
    const  newton::Settings<T> & NS = *((const  newton::Settings<T> *) &S);

    if (!S.get_use_winding()) //pure Newton iterations wanted
    {
        newton_iterations(F, B, NS, zeros);

        return copy_zeros(N, Z, V, n_zeros);
    }

    const Complex I(0.0e0, 1.0e0);

    const Complex Z0 = B.get_xmin() + B.get_ymin() * I;
    const Complex Z1 = B.get_xmax() + B.get_ymin() * I;
    const Complex Z2 = B.get_xmax() + B.get_ymax() * I;
    const Complex Z3 = B.get_xmin() + B.get_ymax() * I;

    const pVertex V0 = new Vertex<T>(Z0, F.eval_value(Z0));
    const pVertex V1 = new Vertex<T>(Z1, F.eval_value(Z1));
    const pVertex V2 = new Vertex<T>(Z2, F.eval_value(Z2));
    const pVertex V3 = new Vertex<T>(Z3, F.eval_value(Z3));

    average = ( std::abs(V0->get_f()) + std::abs(V1->get_f()) + std::abs(V2->get_f()) + std::abs(V3->get_f()) ) / 4.0;

    const pVertex V0c = new Vertex<T>(*V0);
    const pVertex V1c = new Vertex<T>(*V1);
    const pVertex V2c = new Vertex<T>(*V2);
    const pVertex V3c = new Vertex<T>(*V3);

    pEdge E0 = new Edge<T>(V0,  V1,  F, WS);
    pEdge E1 = new Edge<T>(V1c, V2,  F, WS);
    pEdge E2 = new Edge<T>(V3,  V2c, F, WS);
    pEdge E3 = new Edge<T>(V0c, V3c, F, WS);

    List<pVertex> * nulls = new List<pVertex>;

    pRectangle R = new Rectangle<T>(E0, E1, E2, E3, 0, nulls, F, WS);

    partition = new Tree<pRectangle>(R);

    Tree_Node<pRectangle> * root = partition->get_root();

    find_zeros_in_rectangle(root);
}
catch (std::bad_alloc & e)
{
    add_to_status("Error: memory allocation in heap failed. Check for the available memory.");
    //std::clog << *this;
}
catch (exc::Algorithm_Error_Exception & e)
{
    add_to_status(e.what());
    //std::clog << *this;
}

return copy_zeros(N, Z, V, n_zeros);
}

/********************************************************************/

template <typename T> int ZeroSolver<T> :: copy_zeros (const int & N,
                                                       std::complex<T> * Z,
                                                       std::complex<T> * V,
                                                       int & n_zeros)
{
n_zeros = zeros->get_size();

if (n_zeros > N)
{
    add_to_status("Warning: the size of output array is less than the number of zeros found: some zeros will not be copied.");
}

int i = 0;

List_Node<pVertex> * node = zeros->get_first();

while (node && i < N)
{
    Z[i] = node->get_data()->get_z();
    V[i] = node->get_data()->get_f();
    node = node->get_next();
    i++;
}

if (is_found)
{
    add_to_status("Warning: the solver has been called more than once.\n"
                  "Be aware that only the first call to FindZeros() does actual calculations.\n"
                  "Subsequent calls only return the zeros that already found by the solver.\n"
                  "Please, use a new ZeroSolver object to run the search with another settings.");

    return status->get_size() - 1;
}

if (status->get_size())
{
    add_to_status("The solver has detected some failures and therefore the list of zeros may be incomplete.\n"
                  "Please, change the solver settings, check the function behaviour"
                  " and have a look at the solver dump to find/fix the issue.");
}
else
{
    add_to_status("The solver has finished successfully without any detected failures.\n"
                  "To verify the correctness of the solution, improve the solver settings and run the search again.");
}

is_found = true;

return status->get_size() - 1;
}

/********************************************************************/

template <typename T> inline int ZeroSolver<T> :: find_zeros_in_rectangle (Tree_Node<pRectangle> * node)
{
if (node->get_level() > S.get_max_part_level())
{
    throw_Algorithm_Error_Exception("maximum allowed partition level is exceeded during recursive rectangle processing");
}

pRectangle rect = node->get_data();

int winding = rect->calc_winding_number();

int n_wanted = winding - rect->get_zeros()->get_size();

if (n_wanted == 0) return 0;

if (winding < 0)
{
    throw_Algorithm_Error_Exception("winding number of the rectangle is negative: the function has poles within it");
}

if (n_wanted < 0)
{
    throw_Algorithm_Error_Exception("the rectangle winding number is less than the number of already detected zeros");
}

//case when n_wanted > 0:

const Vertex<T> * lbv = rect->get_lb_vertex();
const Vertex<T> * rtv = rect->get_rt_vertex();

const Box<T> RB(real(lbv->get_z()), real(rtv->get_z()), imag(lbv->get_z()), imag(rtv->get_z()));

if (S.get_print_level() > 0)
{
    std::clog << "\nThe winding number " << winding << " is found for:" << RB;
}

const newton::Settings<T> & NS = *((const newton::Settings<T> *) &S);

NS.set_n_target(n_wanted);

int n_new_zeros = 0;

if (n_wanted > 1) //check for zeros with multiplicity > 1
{
    int nx = NS.get_n_split_x();
    int ny = NS.get_n_split_y();

    NS.set_multiplicity(n_wanted);
    NS.set_n_split_x(1);
    NS.set_n_split_y(1);

    n_new_zeros = newton_iterations(F, RB, NS, rect->get_zeros());

    NS.set_multiplicity(1);
    NS.set_n_split_x(nx);
    NS.set_n_split_y(ny);
}

if (n_new_zeros == 0) //check for zeros with multiplicity = 1
{
    n_new_zeros = newton_iterations(F, RB, NS, rect->get_zeros());
}

//add new zeros found to list of all zeros:
int n = 0;

List_Node<pVertex> * list_node = rect->get_zeros()->get_last();

while (list_node  &&  n++ < n_new_zeros)
{
    zeros->push_back(list_node->get_data());
    list_node = list_node->get_prev();
}

if (n_new_zeros == n_wanted) return n_wanted;

if (n_new_zeros > n_wanted)
{
    throw_Algorithm_Error_Exception("the rectangle winding number is less than the number of detected zeros");
}

//case when n_new_zeros < n_wanted:

int n_rest = n_wanted - n_new_zeros; //number of rest zeros after Newton iterations

T abs_eps_z, rel_eps_z;
S.get_eps_for_arg(abs_eps_z, rel_eps_z);

Complex Zc = rect->get_center();

if (rect->get_diagonal() < std::max(abs_eps_z, rel_eps_z * std::abs(Zc)))
{
    T abs_eps_f, rel_eps_f;
    S.get_eps_for_func(abs_eps_f, rel_eps_f);

    Complex Fc = F.eval_value(Zc);

    if (std::abs(Fc) < std::max(abs_eps_f, rel_eps_f * average))
    {
        //add the rectangle center as new zeros to list of all zeros:
        for (int l = 0; l < n_rest; ++l)
        {
            pVertex V = new Vertex<T>(Zc, Fc);
            rect->get_zeros()->push_back(V);
            zeros->push_back(V);
        }
        std::clog << "\nrectangle center is added due to small diagonal!\n";
        return n_wanted;
    }
}

//need to partition the rectangle further by creating subrectangles with existing zeros added:
pRectangle R1 = 0, R2 = 0;

rect->partition(R1, R2);

Tree_Node<pRectangle> * child1 = node->add_child(R1);
Tree_Node<pRectangle> * child2 = node->add_child(R2);

n_new_zeros = find_zeros_in_rectangle(child1) + find_zeros_in_rectangle(child2);

if (n_new_zeros < n_rest)
{
    throw_Algorithm_Error_Exception("the rectangle partition provided less zeros than it should be");
}
else if (n_new_zeros > n_rest)
{
    throw_Algorithm_Error_Exception("the rectangle partition provided more zeros than it should be");
}

return n_wanted;
}

/********************************************************************/
} //end of namespace
} //end of namespace

/********************************************************************/

#endif
