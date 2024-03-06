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

/*
The implementation of the Fortran bindings for the ZerSol C++ library functions.
*/

#include "zerosolver.hpp"

/********************************************************************/

#define _Real double /* float, double (recommended) */

/********************************************************************/

typedef _Real     _real_;
typedef _real_ _complex_;
typedef void (*complex_function) (_complex_ *, _complex_ *, void *);

/********************************************************************/

class Fparam
{
    public:

        Fparam (complex_function f_, complex_function df_, void * data_)
               :
               f(f_),
               df(df_),
               data(data_)
        {
        }

        complex_function f;
        complex_function df;

        void * data;
};

/********************************************************************/

template <typename T> class Solver
{
    public:

        Solver (void) : D(0), F(0), B(0), S(0), Z(0), Cstart(0), Pstart(0)
        {
        }

        ~Solver (void)
        {
            delete D;
            delete F;
            delete B;
            delete S;
            delete Z;
            delete Cstart;
            delete Pstart;
        }

        Fparam                     * D;
        type::Function<T>          * F;
        type::Box<T>               * B;
        alg::zersol::Settings<T>   * S;
        alg::zersol::ZeroSolver<T> * Z;

        _complex_                  * Cstart;
        std::complex<T>            * Pstart;
};

/********************************************************************/

template <typename T> std::complex<T> func (const std::complex<T> & z, void * p)
{
Fparam * P = static_cast<Fparam *>(p);

_complex_ Z[2], F[2];

Z[0] = z.real();
Z[1] = z.imag();

P->f(Z, F, P->data);

return std::complex<T>(F[0], F[1]);
}

/********************************************************************/

template <typename T> std::complex<T> dfunc (const std::complex<T> & z, void * p)
{
Fparam * P = static_cast<Fparam *>(p);

_complex_ Z[2], dF[2];

Z[0] = z.real();
Z[1] = z.imag();

P->df(Z, dF, P->data);

return std::complex<T>(dF[0], dF[1]);
}

/********************************************************************/

template std::complex<_real_>  func<_real_> (const std::complex<_real_> & z, void * p);
template std::complex<_real_> dfunc<_real_> (const std::complex<_real_> & z, void * p);

/********************************************************************/

extern "C"
{
void create_solver_(void ** sol,
                    complex_function * f,
                    complex_function * df,
                    void * data,
                    _real_ * xmin,
                    _real_ * xmax,
                    _real_ * ymin,
                    _real_ * ymax
                   )
{
Solver<_real_> * solver = 0;

try
{
    solver = new Solver<_real_>();

    solver->D = new Fparam(*f, *df, data);

    type::Function<_real_>::cmplx_func F  = func<_real_>;
    type::Function<_real_>::cmplx_func dF = (*df ? dfunc<_real_> : 0);

    solver->F = new type::Function<_real_>(F, dF, static_cast<void *>(solver->D), std::string());

    solver->B = new type::Box<_real_>(*xmin, *xmax, *ymin, *ymax);

    solver->S = new alg::zersol::Settings<_real_>();
}
catch (std::bad_alloc & e)
{
    std::clog << "Error: memory allocation in heap failed. Check for the available memory.";
    delete solver;
    solver = 0;
}

*sol = static_cast<void *>(solver);
}

/********************************************************************/

void find_zeros_(void ** solver, int * N, _complex_ * Z, _complex_ * V, int * n_zeros, int * status)
{
Solver<_real_> * sol = static_cast<Solver<_real_> *>(*solver);

sol->Z = new alg::zersol::ZeroSolver<_real_>(*sol->F, *sol->B, *sol->S);

std::complex<_real_> * ZZ = new std::complex<_real_>[*N];
std::complex<_real_> * VV = new std::complex<_real_>[*N];

*status = sol->Z->FindZeros(*N, ZZ, VV, *n_zeros);

int count = std::min(*n_zeros, *N);

for (int i = 0; i < count; ++i)
{
    Z[2*i+0] = ZZ[i].real();
    Z[2*i+1] = ZZ[i].imag();
    V[2*i+0] = VV[i].real();
    V[2*i+1] = VV[i].imag();
}

delete [] ZZ;
delete [] VV;
}

/********************************************************************/

void print_status_(void ** solver, char * file_name, int * length)
{
file_name[*length] = '\0'; //end of string symbol

if (*length > 0)
{
    std::ofstream file(file_name, std::ofstream::out);
    static_cast<Solver<_real_> *>(*solver)->Z->print_status(file);
    file.close();
}
else
{
    static_cast<Solver<_real_> *>(*solver)->Z->print_status(std::cout);
}
}

/********************************************************************/

void print_dump_(void ** solver, char * file_name, int * length)
{
file_name[*length] = '\0'; //end of string symbol

if (*length > 0)
{
    std::ofstream file(file_name, std::ofstream::out);
    file << *(static_cast<Solver<_real_> *>(*solver)->Z);
    file.close();
}
else
{
    std::cout << *(static_cast<Solver<_real_> *>(*solver)->Z);
}
}

/********************************************************************/

void free_solver_(void ** solver)
{
delete static_cast<Solver<_real_> *>(*solver);
}

/********************************************************************/

void set_nd_method_(void ** solver, int * method)
{
static_cast<Solver<_real_> *>(*solver)->F->set_nd_method(*method);
}

/********************************************************************/

void set_nd_step_(void ** solver, _real_ * h)
{
static_cast<Solver<_real_> *>(*solver)->F->set_nd_step(*h);
}

/********************************************************************/

void set_min_rec_lev_(void ** solver, int * min_rec_lev)
{
static_cast<Solver<_real_> *>(*solver)->S->set_min_rec_lev(*min_rec_lev);
}

/********************************************************************/

void set_max_rec_lev_(void ** solver, int * max_rec_lev)
{
static_cast<Solver<_real_> *>(*solver)->S->set_max_rec_lev(*max_rec_lev);
}

/********************************************************************/

void set_interp_err_(void ** solver, _real_ * interp_err)
{
static_cast<Solver<_real_> *>(*solver)->S->set_interp_err(*interp_err);
}

/********************************************************************/

void set_jump_err_(void ** solver, _real_ * jump_err)
{
static_cast<Solver<_real_> *>(*solver)->S->set_jump_err(*jump_err);
}

/********************************************************************/

void set_n_target_(void ** solver, int * n_target)
{
static_cast<Solver<_real_> *>(*solver)->S->set_n_target(*n_target);
}

/********************************************************************/

void set_start_array_(void ** solver, int * n_start, _complex_ * start)
{
static_cast<Solver<_real_> *>(*solver)->Pstart = new std::complex<_real_>[*n_start];

for (int i = 0; i < *n_start; ++i)
{
    static_cast<Solver<_real_> *>(*solver)->Pstart[i] = std::complex<_real_>(start[2*i+0], start[2*i+1]);
}

static_cast<Solver<_real_> *>(*solver)->S->set_start_array(*n_start, static_cast<Solver<_real_> *>(*solver)->Pstart);
}

/********************************************************************/

void set_n_split_x_(void ** solver, int * n_split_x)
{
static_cast<Solver<_real_> *>(*solver)->S->set_n_split_x(*n_split_x);
}

/********************************************************************/

void set_n_split_y_(void ** solver, int * n_split_y)
{
static_cast<Solver<_real_> *>(*solver)->S->set_n_split_y(*n_split_y);
}

/********************************************************************/

void set_max_iter_num_(void ** solver, int * max_iter_num)
{
static_cast<Solver<_real_> *>(*solver)->S->set_max_iter_num(*max_iter_num);
}

/********************************************************************/

void set_multiplicity_(void ** solver, int * multiplicity)
{
static_cast<Solver<_real_> *>(*solver)->S->set_multiplicity(*multiplicity);
}

/********************************************************************/

void set_eps_for_arg_(void ** solver, _real_ * abs_eps_z, _real_ * rel_eps_z)
{
static_cast<Solver<_real_> *>(*solver)->S->set_eps_for_arg(*abs_eps_z, *rel_eps_z);
}

/********************************************************************/

void set_eps_for_func_(void ** solver, _real_ * abs_eps_f, _real_ * rel_eps_f)
{
static_cast<Solver<_real_> *>(*solver)->S->set_eps_for_func(*abs_eps_f, *rel_eps_f);
}

/********************************************************************/

void set_use_winding_(void ** solver, int * use_winding)
{
static_cast<Solver<_real_> *>(*solver)->S->set_use_winding(*use_winding);
}

/********************************************************************/

void set_max_part_level_(void ** solver, int * max_part_level)
{
static_cast<Solver<_real_> *>(*solver)->S->set_max_part_level(*max_part_level);
}

/********************************************************************/

void set_debug_level_(void ** solver, int * debug_level)
{
static_cast<Solver<_real_> *>(*solver)->S->set_debug_level(*debug_level);
}

/********************************************************************/

void set_print_level_(void ** solver, int * print_level)
{
static_cast<Solver<_real_> *>(*solver)->S->set_print_level(*print_level);
}

/********************************************************************/

void get_nd_method_(void ** solver, int * nd_method)
{
*nd_method = static_cast<Solver<_real_> *>(*solver)->F->get_nd_method();
}

/********************************************************************/

void get_nd_step_(void ** solver, _real_ * nd_step)
{
*nd_step = static_cast<Solver<_real_> *>(*solver)->F->get_nd_step();
}

/********************************************************************/

void get_min_rec_lev_(void ** solver, int * min_rec_lev)
{
*min_rec_lev = static_cast<Solver<_real_> *>(*solver)->S->get_min_rec_lev();
}

/********************************************************************/

void get_max_rec_lev_(void ** solver, int * max_rec_lev)
{
*max_rec_lev = static_cast<Solver<_real_> *>(*solver)->S->get_max_rec_lev();
}

/********************************************************************/

void get_interp_err_(void ** solver, _real_ * interp_err)
{
*interp_err = static_cast<Solver<_real_> *>(*solver)->S->get_interp_err();
}

/********************************************************************/

void get_jump_err_(void ** solver, _real_ * jump_err)
{
*jump_err = static_cast<Solver<_real_> *>(*solver)->S->get_jump_err();
}

/********************************************************************/

void get_n_target_(void ** solver, int * n_target)
{
*n_target = static_cast<Solver<_real_> *>(*solver)->S->get_n_target();
}

/********************************************************************/

void get_start_array_(void ** solver, int * n_start, _complex_ * start)
{
int N_start = 0;

std::complex<_real_> * ptr = 0;

static_cast<Solver<_real_> *>(*solver)->S->get_start_array(N_start, ptr);

if (N_start == 0) { *n_start = 0; return; }

int count = std::min(*n_start, N_start);

for (int i = 0; i < count; ++i)
{
    start[2*i+0] = ptr[i].real();
    start[2*i+1] = ptr[i].imag();
}

*n_start = N_start;
}

/********************************************************************/

void get_n_split_x_(void ** solver, int * n_split_x)
{
*n_split_x = static_cast<Solver<_real_> *>(*solver)->S->get_n_split_x();
}

/********************************************************************/

void get_n_split_y_(void ** solver, int * n_split_y)
{
*n_split_y = static_cast<Solver<_real_> *>(*solver)->S->get_n_split_y();
}

/********************************************************************/

void get_max_iter_num_(void ** solver, int * max_iter_num)
{
*max_iter_num = static_cast<Solver<_real_> *>(*solver)->S->get_max_iter_num();
}

/********************************************************************/

void get_multiplicity_(void ** solver, int * multiplicity)
{
*multiplicity = static_cast<Solver<_real_> *>(*solver)->S->get_multiplicity();
}

/********************************************************************/

void get_eps_for_arg_(void ** solver, _real_ * abs_eps_z, _real_ * rel_eps_z)
{
static_cast<Solver<_real_> *>(*solver)->S->get_eps_for_arg(*abs_eps_z, *rel_eps_z);
}

/********************************************************************/

void get_eps_for_func_(void ** solver, _real_ * abs_eps_f, _real_ * rel_eps_f)
{
static_cast<Solver<_real_> *>(*solver)->S->get_eps_for_func(*abs_eps_f, *rel_eps_f);
}

/********************************************************************/

void get_use_winding_(void ** solver, int * use_winding)
{
*use_winding = static_cast<Solver<_real_> *>(*solver)->S->get_use_winding();
}

/********************************************************************/

void get_max_part_level_(void ** solver, int * max_part_level)
{
*max_part_level = static_cast<Solver<_real_> *>(*solver)->S->get_max_part_level();
}

/********************************************************************/

void get_debug_level_(void ** solver, int * debug_level)
{
*debug_level = static_cast<Solver<_real_> *>(*solver)->S->get_debug_level();
}

/********************************************************************/

void get_print_level_(void ** solver, int * print_level)
{
*print_level = static_cast<Solver<_real_> *>(*solver)->S->get_print_level();
}

/********************************************************************/
}
