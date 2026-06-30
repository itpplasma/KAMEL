/*! \file
    \brief The basic driver program to run KiLCA code.
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <climits>
#include <iostream>
#include <fstream>
#include <unistd.h>

#include "core.h"
#include "main_linear.h"

using namespace std;

/*******************************************************************/
//!The main program
/*!
Calls KiLCA library functions
*/
int main (int argc, char **argv)
{

//!gets the path to the project from the current dir or the command line
char *path = new char[1024];

if (argc == 1)
{
   getcwd (path, 1024); //current directory
}
else if (argc >= 2)
{
   strcpy (path, argv[1]); //command line argument
}

if (path[strlen(path)-1] != '/') strcat(path, "/");

intptr_t cd = core_data_create_ (path); //!<core data structure contains pointers to all important code data
set_core_data_in_core_module_ (&cd);  //!<stores the core data pointer in fortran module for later use

core_data_calc_and_set_mode_independent_ (cd); //!<evident from the function name

if (get_antenna_flag_eigmode_ () == 0)
{
    core_data_calc_and_set_mode_dependent_antenna_ (cd); //!<normal operation with antenna
}
else
{
    core_data_calc_and_set_mode_dependent_eigmode_ (cd); //!<special operation for instabilities hunting
}

//!frees the used memory
delete [] path;

core_data_destroy_ (cd);

return 0;
}

/*******************************************************************/
