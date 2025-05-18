/*! \file
    \brief The basic driver program to run KiLCA code.
*/

#include "core.h"

#include <climits>
#include <cmath>
#include <cstring>
#include <stdlib.h>
#include <unistd.h>

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

core_data *cd = new core_data (path); //!<core data structure contains pointers to all important code data
set_core_data_in_core_module_ (&cd);  //!<stores the core data pointer in fortran module for later use

cd->calc_and_set_mode_independent_core_data (); //!<evident from the function name

if (cd->sd->as->flag_eigmode == 0)
{
    cd->calc_and_set_mode_dependent_core_data_antenna (); //!<normal operation with antenna
}
else
{
    cd->calc_and_set_mode_dependent_core_data_eigmode (); //!<special operation for instabilities hunting
}

//!frees the used memory
delete [] path;

delete cd;

return 0;
}

/*******************************************************************/
