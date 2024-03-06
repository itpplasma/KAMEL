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

#include "zerosolver.hpp"

#include "zersolc.h"

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

_complex_ Z = z.real() + z.imag() * IU;

_complex_ F = P->f(Z, P->data);

return std::complex<T>(creal(F), cimag(F));
}

/********************************************************************/

template <typename T> std::complex<T> dfunc (const std::complex<T> & z, void * p)
{
Fparam * P = static_cast<Fparam *>(p);

_complex_ Z = z.real() + z.imag() * IU;

_complex_ dF = P->df(Z, P->data);

return std::complex<T>(creal(dF), cimag(dF));
}

/********************************************************************/

template std::complex<_real_>  func<_real_> (const std::complex<_real_> & z, void * p);
template std::complex<_real_> dfunc<_real_> (const std::complex<_real_> & z, void * p);

/********************************************************************/

void * create_solver (complex_function f,
                      complex_function df,
                      void * data,
                      char const * description,
                      _real_ xmin,
                      _real_ xmax,
                      _real_ ymin,
                      _real_ ymax
                     )
{
Solver<_real_> * solver = 0;

try
{
    solver = new Solver<_real_>();

    solver->D = new Fparam(f, df, data);

    type::Function<_real_>::cmplx_func F  = func<_real_>;
    type::Function<_real_>::cmplx_func dF = (df ? dfunc<_real_> : 0);

    solver->F = new type::Function<_real_>(F, dF, static_cast<void *>(solver->D), std::string(description));

    solver->B = new type::Box<_real_>(xmin, xmax, ymin, ymax);

    solver->S = new alg::zersol::Settings<_real_>();
}
catch (std::bad_alloc & e)
{
    std::clog << "Error: memory allocation in heap failed. Check for the available memory.";
    delete solver;
    solver = 0;
}

return static_cast<void *>(solver);
}

/********************************************************************/

int find_zeros (void * solver, int N, _complex_ * Z, _complex_ * V, int * n_zeros)
{
Solver<_real_> * sol = static_cast<Solver<_real_> *>(solver);

sol->Z = new alg::zersol::ZeroSolver<_real_>(*sol->F, *sol->B, *sol->S);

std::complex<_real_> * ZZ = new std::complex<_real_>[N];
std::complex<_real_> * VV = new std::complex<_real_>[N];

int status = sol->Z->FindZeros(N, ZZ, VV, *n_zeros);

int count = std::min(*n_zeros, N);

for (int i = 0; i < count; ++i)
{
    Z[i] = ZZ[i].real() + ZZ[i].imag() * IU;
    V[i] = VV[i].real() + VV[i].imag() * IU;
}

delete [] ZZ;
delete [] VV;

return status;
}

/********************************************************************/

void print_status (void * solver, const char * file_name)
{
if (std::string(file_name).length() > 0)
{
    std::ofstream file(file_name, std::ofstream::out);
    static_cast<Solver<_real_> *>(solver)->Z->print_status(file);
    file.close();
}
else
{
    static_cast<Solver<_real_> *>(solver)->Z->print_status(std::cout);
}
}

/********************************************************************/

void print_dump (void * solver, const char * file_name)
{
if (std::string(file_name).length() > 0)
{
    std::ofstream file(file_name, std::ofstream::out);
    file << *(static_cast<Solver<_real_> *>(solver)->Z);
    file.close();
}
else
{
    std::cout << *(static_cast<Solver<_real_> *>(solver)->Z);
}
}

/********************************************************************/

void free_solver (void * solver)
{
delete static_cast<Solver<_real_> *>(solver);
}

/********************************************************************/

void set_nd_method(void * solver, int method)
{
static_cast<Solver<_real_> *>(solver)->F->set_nd_method(method);
}

/********************************************************************/

void set_nd_step(void * solver, _real_ h)
{
static_cast<Solver<_real_> *>(solver)->F->set_nd_step(h);
}

/********************************************************************/

void set_min_rec_lev(void * solver, int min_rec_lev)
{
static_cast<Solver<_real_> *>(solver)->S->set_min_rec_lev(min_rec_lev);
}

/********************************************************************/

void set_max_rec_lev(void * solver, int max_rec_lev)
{
static_cast<Solver<_real_> *>(solver)->S->set_max_rec_lev(max_rec_lev);
}

/********************************************************************/

void set_interp_err(void * solver, _real_ interp_err)
{
static_cast<Solver<_real_> *>(solver)->S->set_interp_err(interp_err);
}

/********************************************************************/

void set_jump_err(void * solver, _real_ jump_err)
{
static_cast<Solver<_real_> *>(solver)->S->set_jump_err(jump_err);
}

/********************************************************************/

void set_n_target(void * solver, int n_target)
{
static_cast<Solver<_real_> *>(solver)->S->set_n_target(n_target);
}

/********************************************************************/

void set_start_array(void * solver, int n_start, _complex_ * start)
{
static_cast<Solver<_real_> *>(solver)->Pstart = new std::complex<_real_>[n_start];

for (int i = 0; i < n_start; ++i)
{
    static_cast<Solver<_real_> *>(solver)->Pstart[i] = std::complex<_real_>(creal(start[i]), cimag(start[i]));
}

static_cast<Solver<_real_> *>(solver)->S->set_start_array(n_start, static_cast<Solver<_real_> *>(solver)->Pstart);
}

/********************************************************************/

void set_n_split_x(void * solver, int n_split_x)
{
static_cast<Solver<_real_> *>(solver)->S->set_n_split_x(n_split_x);
}

/********************************************************************/

void set_n_split_y(void * solver, int n_split_y)
{
static_cast<Solver<_real_> *>(solver)->S->set_n_split_y(n_split_y);
}

/********************************************************************/

void set_max_iter_num(void * solver, int max_iter_num)
{
static_cast<Solver<_real_> *>(solver)->S->set_max_iter_num(max_iter_num);
}

/********************************************************************/

void set_multiplicity(void * solver, int multiplicity)
{
static_cast<Solver<_real_> *>(solver)->S->set_multiplicity(multiplicity);
}

/********************************************************************/

void set_eps_for_arg(void * solver, _real_ abs_eps_z, _real_ rel_eps_z)
{
static_cast<Solver<_real_> *>(solver)->S->set_eps_for_arg(abs_eps_z, rel_eps_z);
}

/********************************************************************/

void set_eps_for_func(void * solver, _real_ abs_eps_f, _real_ rel_eps_f)
{
static_cast<Solver<_real_> *>(solver)->S->set_eps_for_func(abs_eps_f, rel_eps_f);
}

/********************************************************************/

void set_use_winding(void * solver, int use_winding)
{
static_cast<Solver<_real_> *>(solver)->S->set_use_winding(use_winding);
}

/********************************************************************/

void set_max_part_level(void * solver, int max_part_level)
{
static_cast<Solver<_real_> *>(solver)->S->set_max_part_level(max_part_level);
}

/********************************************************************/

void set_debug_level(void * solver, int debug_level)
{
static_cast<Solver<_real_> *>(solver)->S->set_debug_level(debug_level);
}

/********************************************************************/

void set_print_level(void * solver, int print_level)
{
static_cast<Solver<_real_> *>(solver)->S->set_print_level(print_level);
}

/********************************************************************/

int get_nd_method(void * solver)
{
return static_cast<Solver<_real_> *>(solver)->F->get_nd_method();
}

/********************************************************************/

_real_ get_nd_step(void * solver)
{
return static_cast<Solver<_real_> *>(solver)->F->get_nd_step();
}

/********************************************************************/

int get_min_rec_lev(void * solver)
{
return static_cast<Solver<_real_> *>(solver)->S->get_min_rec_lev();
}

/********************************************************************/

int get_max_rec_lev(void * solver)
{
return static_cast<Solver<_real_> *>(solver)->S->get_max_rec_lev();
}

/********************************************************************/

_real_ get_interp_err(void * solver)
{
return static_cast<Solver<_real_> *>(solver)->S->get_interp_err();
}

/********************************************************************/

_real_ get_jump_err(void * solver)
{
return static_cast<Solver<_real_> *>(solver)->S->get_jump_err();
}

/********************************************************************/

int get_n_target(void * solver)
{
return static_cast<Solver<_real_> *>(solver)->S->get_n_target();
}

/********************************************************************/

void get_start_array(void * solver, int * n_start, _complex_ ** start)
{
*n_start = 0;  *start = 0;

std::complex<_real_> * ptr = 0;

static_cast<Solver<_real_> *>(solver)->S->get_start_array(*n_start, ptr);

if (*n_start == 0) return;

static_cast<Solver<_real_> *>(solver)->Cstart = new _complex_[*n_start];

for (int i = 0; i < *n_start; ++i)
{
    static_cast<Solver<_real_> *>(solver)->Cstart[i] = ptr[i].real() + ptr[i].imag() * IU;
}

*start = static_cast<Solver<_real_> *>(solver)->Cstart;
}

/********************************************************************/

int get_n_split_x(void * solver)
{
return static_cast<Solver<_real_> *>(solver)->S->get_n_split_x();
}

/********************************************************************/

int get_n_split_y(void * solver)
{
return static_cast<Solver<_real_> *>(solver)->S->get_n_split_y();
}

/********************************************************************/

int get_max_iter_num(void * solver)
{
return static_cast<Solver<_real_> *>(solver)->S->get_max_iter_num();
}

/********************************************************************/

int get_multiplicity(void * solver)
{
return static_cast<Solver<_real_> *>(solver)->S->get_multiplicity();
}

/********************************************************************/

void get_eps_for_arg(void * solver, _real_ * abs_eps_z, _real_ * rel_eps_z)
{
static_cast<Solver<_real_> *>(solver)->S->get_eps_for_arg(*abs_eps_z, *rel_eps_z);
}

/********************************************************************/

void get_eps_for_func(void * solver, _real_ * abs_eps_f, _real_ * rel_eps_f)
{
static_cast<Solver<_real_> *>(solver)->S->get_eps_for_func(*abs_eps_f, *rel_eps_f);
}

/********************************************************************/

int get_use_winding(void * solver)
{
return static_cast<Solver<_real_> *>(solver)->S->get_use_winding();
}

/********************************************************************/

int get_max_part_level(void * solver)
{
return static_cast<Solver<_real_> *>(solver)->S->get_max_part_level();
}

/********************************************************************/

int get_debug_level(void * solver)
{
return static_cast<Solver<_real_> *>(solver)->S->get_debug_level();
}

/********************************************************************/

int get_print_level(void * solver)
{
return static_cast<Solver<_real_> *>(solver)->S->get_print_level();
}

/********************************************************************/
