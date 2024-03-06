/*! \file core.h
    \brief The declaration of core_data class.
*/

#ifndef CORE_INCLUDE

#define CORE_INCLUDE

#include <cstring>

#include "code_settings.h"

#include "settings.h"
#include "background.h"
#include "mode.h"

using namespace :: std;

/*! \class core_data
    \brief The top level class of KiLCA library.
*/
class core_data
{
public:
    char *path2project; //!<project path

    settings *sd; //!<pointer to settings object

    background *bp; //!<pointer to background object

    int dim; //!<dimension of the mode_data array (number of modes)

    mode_data **mda; //!<array of mode_data objects

public:
    core_data (char *path); //!<constructor

    ~core_data (void); //!<destructor

    void delete_modes_array (void);

    void calc_and_set_mode_independent_core_data (void);

    void calc_and_set_mode_dependent_core_data_antenna (void);

    void calc_and_set_mode_dependent_core_data_eigmode (void);

    void calc_and_set_mode_dependent_core_data_antenna_interface (void);

    void calc_and_set_mode_dependent_core_data_antenna_interface (int m, int n, int flag = 0);
};

//!The functions below are needed for C - Fortran interface
extern "C"
{
void get_pointer_precision_ (int *);

void set_core_data_in_core_module_ (core_data **);

void set_background_in_core_module_ (background **);

void set_settings_in_core_module_ (settings **);

void clear_all_data_in_mode_data_module_ (void);
}

/*******************************************************************/

#endif
