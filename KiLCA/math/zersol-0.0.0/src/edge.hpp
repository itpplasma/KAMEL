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
    \brief The declaration and definition of Edge class.
*/

#ifndef FZF_EDGE_HPP
#define FZF_EDGE_HPP

#include <cmath>
#include <complex>
#include <iostream>

#include "function.hpp"
#include "vertex.hpp"
#include "interval.hpp"
#include "list.hpp"
#include "tree.hpp"
#include "exceptions.hpp"

/********************************************************************/
/*!
\brief The namespace contains classes and algorithms needed for calculation of winding numbers.
*/
namespace winding
{
/********************************************************************/

using namespace type;

/********************************************************************/
/*!
\brief The class contains settings needed for evaluation of the winding number.
*/
template <typename T = double> class Settings
{
    public:

        /*!
         *  \brief A constructor of the Settings object.
         *  \param min_rec_lev_ minimum recursion level allowed in the winding number evaluation procedure
         *  \param max_rec_lev_ maximum recursion level allowed in the winding number evaluation procedure
         *  \param interp_err_  tolerance of interpolation for the winding number evaluation procedure
         *  \param jump_err_    tolerance of test for discontinuities for the winding number evaluation procedure
         *  \param debug_level_ debug level (determines the amount of internal checks)
         *  \param print_level_ print level (determines the amount of details printed)
         */
        Settings (int min_rec_lev_ = 4,
                  int max_rec_lev_ = 32,
                  T    interp_err_ = 1.0e-2,
                  T      jump_err_ = 1.0e-2,
                  int debug_level_ = 0,
                  int print_level_ = 2
                 )
                 :
                 min_rec_lev(min_rec_lev_),
                 max_rec_lev(max_rec_lev_),
                   interp_err(interp_err_),
                       jump_err(jump_err_),
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

        //set functions:
        /*!
         *  \brief Set the minimum recursion level for the winding number evaluation procedure.
         *  \param min_rec_lev_ minimum recursion level
         */
        void set_min_rec_lev (const int & min_rec_lev_) { min_rec_lev = min_rec_lev_; }

        /*!
         *  \brief Set the maximum recursion level for the winding number evaluation procedure.
         *  \param max_rec_lev_ maximum recursion level
         */
        void set_max_rec_lev (const int & max_rec_lev_) { max_rec_lev = max_rec_lev_; }

        /*!
         *  \brief Set a tolerance of interpolation for the winding number evaluation procedure.
         *  \param interp_err_ tolerance of interpolation
         */
        void set_interp_err (const T & interp_err_) { interp_err = interp_err_; }

        /*!
         *  \brief Set a tolerance of test for discontinuities for the winding number evaluation procedure.
         *  \param jump_err_ tolerance of test for discontinuities
         */
        void set_jump_err   (const T & jump_err_)   { jump_err = jump_err_; }

        void set_debug_level (const int & debug_level_) { debug_level = debug_level_; }
        void set_print_level (const int & print_level_) { print_level = print_level_; }

        //get functions:
        /*!
         *  \brief Get the minimum recursion level for the winding number evaluation procedure.
         *  \return minimum recursion level
         */
        const int & get_min_rec_lev (void) const { return min_rec_lev; }

        /*!
         *  \brief Get the maximum recursion level for the winding number evaluation procedure.
         *  \return maximum recursion level
         */
        const int & get_max_rec_lev (void) const { return max_rec_lev; }

        /*!
         *  \brief Get a tolerance of interpolation for the winding number evaluation procedure.
         *  \return tolerance of interpolation
         */
        const T & get_interp_err (void) const { return interp_err; }

        /*!
         *  \brief Get a tolerance of test for discontinuities for the winding number evaluation procedure.
         *  \return tolerance of test for discontinuities
         */
        const T & get_jump_err   (void) const { return jump_err; }

        const int & get_debug_level (void) const { return debug_level; }
        const int & get_print_level (void) const { return print_level; }

        friend std::ostream & operator << (std::ostream & output, const Settings<T> & S)
        {
            output << "\nWinding settings:\n"
                   << "min_rec_lev=" << S.min_rec_lev << "\n"
                   << "max_rec_lev=" << S.max_rec_lev << "\n"
                   << "interp_err="  << S.interp_err  << "\n"
                   << "jump_err="    << S.jump_err    << "\n"
                   << "debug_level=" << S.debug_level << "\n"
                   << "print_level=" << S.print_level << std::endl;

            return output;
        }

    private:

        //adaptive grid options for winding number calculation:
        int   min_rec_lev;  //!< minimum recursion level allowed in the winding number evaluation procedure
        int   max_rec_lev;  //!< maximum recursion level allowed in the winding number evaluation procedure
        T     interp_err;   //!< tolerance of interpolation for the winding number evaluation procedure
        T     jump_err;     //!< tolerance of test for discontinuities for the winding number evaluation procedure
        int   debug_level;  //!< debug level (determines the amount of internal checks)
        int   print_level;  //!< print level (determines the amount of details printed)
};

/********************************************************************/

template <typename T = double> class Edge
{
    typedef std::complex<T> Complex;
    typedef Vertex<T>   *   pVertex;
    typedef Interval<T> *   pInterval;
    typedef Edge<T>     *   pEdge;

    public:

        Edge (pVertex V1, pVertex V2, const Function<T> & F_, const Settings<T> & S_)
        :
        is_alloc(true), list(0), tree(0), F(F_), S(S_)
        {
            try
            {
                list = new List<pVertex>;

                List_Node<pVertex> * Nl = list->push_back(V1);
                List_Node<pVertex> * Nr = list->push_back(V2);

                pInterval interval = new Interval<T>(Nl, Nr);

                tree = new Tree<pInterval>(interval);

                calc_delta_argument_on_interval(tree->get_root());
            }
            catch (std::bad_alloc & e)
            {
                if (list) list->destruct_data(); delete list;
                if (tree) tree->destruct_data(); delete tree;

                throw;
            }
            catch (exc::Algorithm_Error_Exception & e)
            {
                throw;
            }
        }

        Edge (List<pVertex> * list_, Tree<pInterval> * tree_, const Function<T> & F_, const Settings<T> & S_)
        :
        is_alloc(false), list(list_), tree(tree_), F(F_), S(S_)
        {
        }

        ~Edge (void)
        {
            if (is_alloc)
            {
                //clear data stored in list and tree:
                list->destruct_data();
                tree->destruct_data();
            }

            delete list;
            delete tree;
        }

        friend std::ostream & operator << (std::ostream & output, const Edge<T> & E)
        {
            output << *(E.list);
            output << *(E.tree);

            return output;
        }

        void write_vertices (std::ostream & output) const
        {
            output << "\nThe list of edge vertices:"
                   << "\n#\tRe{vertex}\t\t  Im{vertex})\t\t\tRe{f(vertex)}\t\t  Im{f(vertex)})\t\tArg{f(vertex)}";

            output << *list;
        }

        int eval_delta_argument (void) const
        {
            return tree->get_root()->get_data()->get_delta();
        }

        const Vertex<T> * get_first_vertex (void) const
        {
            return list->get_first()->get_data();
        }

        const Vertex<T> * get_central_vertex (void) const
        {
            return tree->get_root()->get_first_child()->get_data()->get_right()->get_data();
        }

        const Vertex<T> * get_last_vertex (void) const
        {
            return list->get_last()->get_data();
        }

        const List_Node<pVertex> * get_first_node (void) const
        {
            return list->get_first();
        }

        const List_Node<pVertex> * get_central_node (void) const
        {
            return tree->get_root()->get_first_child()->get_data()->get_right();
        }

        const List_Node<pVertex> * get_last_node (void) const
        {
            return list->get_last();
        }

        const Tree_Node<pInterval> * get_top_interval (void) const
        {
            return tree->get_root();
        }

        const Tree_Node<pInterval> * get_first_half (void) const
        {
            return tree->get_root()->get_first_child();
        }

        const Tree_Node<pInterval> * get_second_half (void) const
        {
            return tree->get_root()->get_last_child();
        }

        int calc_delta_argument_on_interval (Tree_Node<pInterval> * parent);

        pEdge make_copy (int flag);

    private:

        bool is_alloc;

        List<pVertex>   * list;
        Tree<pInterval> * tree;

        const Function<T> & F;
        const Settings<T> & S;
};

/********************************************************************/

template <typename T> inline int Edge<T> :: calc_delta_argument_on_interval (Tree_Node<pInterval> * parent)
{
if (parent->get_level() > S.get_max_rec_lev())
{
    throw_Algorithm_Error_Exception("maximum allowed recursion level is exceeded during the argument evaluation on the edge: a non-analytic function or a zero crossing are suspected");
}

pInterval I = parent->get_data();

List_Node<pVertex> * left  = I->get_left();
List_Node<pVertex> * right = I->get_right();

const pVertex V1 = left->get_data();
const pVertex V2 = right->get_data();

const T & arg1 = V1->get_a();
const T & arg2 = V2->get_a();

const int level_flag = (parent->get_level() >= S.get_min_rec_lev());

const T pi = 3.141592653589793238462643383279502884197;

const T & jump_err   = S.get_jump_err();
const T & interp_err = S.get_interp_err();

//check if there is a jump in argument from +pi -> -pi or -pi -> +pi
if
(level_flag  &&  std::abs(arg1-pi) < jump_err  &&  std::abs(arg2+pi) < jump_err)
{
    I->set_delta(1);
    return 1;
}
else if
(level_flag  &&  std::abs(arg1+pi) < jump_err  &&  std::abs(arg2-pi) < jump_err)
{
    I->set_delta(-1);
    return -1;
}

//if there is no jump then we add vertex and check accuracy of interpolation:
const Complex & z1 = V1->get_z();
const Complex & z2 = V2->get_z();

const Complex z = (z1 + z2) / T(2.0);

const Complex f = F.eval_value(z);

const T argi = arg1 + (arg2 - arg1) * real( (z - z1)/(z2 - z1) );

const T dev = std::abs(argi - arg(f));

if (level_flag  &&  dev < interp_err)
{
    I->set_delta(0);
    return 0;
}

//at this point we split the interval in 2 parts by adding a new vertex:
pVertex Vc = new Vertex<T>(z, f);

List_Node<pVertex> * cent = list->insert_after(Vc, left);

pInterval Il = new Interval<T>(left, cent);
pInterval Ir = new Interval<T>(cent, right);

Tree_Node<pInterval> * Nl = parent->add_child(Il);
Tree_Node<pInterval> * Nr = parent->add_child(Ir);

int delta = calc_delta_argument_on_interval(Nl) + calc_delta_argument_on_interval(Nr);

I->set_delta(delta);

return delta;
}

/********************************************************************/

template <typename T> inline Edge<T> * Edge<T> :: make_copy (int flag)
{
//flag = 0 - the whole edge
//flag = 1 - first half of the edge
//flag = 2 - second half of the edge

Tree_Node<pInterval> * node = tree->get_root();

if (flag && !node->get_first_child()) //the edge has no subintervals
{
    std::clog << "\nadditional splitting required!\n";

    pInterval I = node->get_data();

    List_Node<pVertex> * left  = I->get_left();
    List_Node<pVertex> * right = I->get_right();

    assert(left == list->get_first() && right == list->get_last());

    const Complex & z1 = left->get_data()->get_z();
    const Complex & z2 = right->get_data()->get_z();

    const Complex z = (z1 + z2) / T(2.0);
    const Complex f = F.eval_value(z);

    pVertex Vc = new Vertex<T>(z, f);

    List_Node<pVertex> * cent = list->insert_after(Vc, left);

    pInterval Il = new Interval<T>(left, cent);
    pInterval Ir = new Interval<T>(cent, right);

    Tree_Node<pInterval> * Nl = node->add_child(Il);
    Tree_Node<pInterval> * Nr = node->add_child(Ir);

    int delta = calc_delta_argument_on_interval(Nl) + calc_delta_argument_on_interval(Nr);

    if (delta != I->get_delta())
    {
        throw_Algorithm_Error_Exception("the argument deltas on subintervals do not fit to the whole interval argument delta");
    }
}

List_Node<pVertex>   * first = 0;
List_Node<pVertex>   * last  = 0;
size_t                 size  = 0;
Tree_Node<pInterval> * top   = 0;

if (flag < 2) first = list->get_first();
else          first = node->get_first_child()->get_data()->get_right();

if (flag != 1) last = list->get_last();
else           last = node->get_first_child()->get_data()->get_right();

if (flag == 0) size = list->get_size();
else           size = list->get_size() / 2 + 1;

if      (flag == 0) top = node;
else if (flag == 1) top = node->get_first_child();
else if (flag == 2) top = node->get_last_child();

List<pVertex> * list_ = new List<pVertex>(first, last, size);

Tree<pInterval> * tree_ = new Tree<pInterval>(top);

return new Edge<T>(list_, tree_, F, S);
}

/********************************************************************/
} //end of namespace

/********************************************************************/

#endif
