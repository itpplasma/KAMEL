#ifndef SETTINGS_INCLUDE
#define SETTINGS_INCLUDE

#include "antenna.h"
#include "back_sett.h"
#include "eigmode_sett.h"
#include "output_sett.h"

#include <cstring>

struct settings_t {
  char *path2project; // project path
  antenna *as; // antenna settings
  back_sett *bs; // background settings
  output_sett *os; // output settings
  eigmode_sett *es; // settings for eigmode search

  void read_settings(); // read settings from input config file

  settings_t(char *path) {
    path2project = new char[1024];
    strcpy(path2project, path);
  }

  ~settings_t(void) {
    delete[] path2project;

    delete as;
    delete bs;
    delete os;
    delete es;
  }
};

/*******************************************************************/

extern "C"
{
void copy_antenna_data_to_antenna_module_ (antenna **);

void copy_background_data_to_background_module_ (back_sett **);
}

/*******************************************************************/

#endif
