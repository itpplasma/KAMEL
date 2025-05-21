#ifndef SETTINGS_INCLUDE
#define SETTINGS_INCLUDE

#include "antenna_sett.h"
#include "back_sett.h"
#include "eigmode_sett.h"
#include "output_sett.h"

struct settings_t {
  char *path2project; // project path
  antenna_sett *as; // antenna settings
  back_sett *bs; // background settings
  output_sett *os; // output settings
  eigmode_sett *es; // settings for eigmode search
  bool flag_debug; // flag for debugging mode

  explicit settings_t(char *path);
  ~settings_t();
};

/*******************************************************************/

extern "C"
{
void copy_antenna_data_to_antenna_module_ (antenna_sett **);

void copy_background_data_to_background_module_ (back_sett **);
}

/*******************************************************************/

#endif
