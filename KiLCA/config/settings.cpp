#include "settings.h"

#include "constants.h"
#include "namelist.h"

#include <cstring>

settings_t::settings_t(char *path) {
  path2project = new char[256];
  strncpy(path2project, path, 256);

  bs = new back_sett;

  es = new eigmode_sett;
  es->fstart = new std::complex<double>[100];

  as = new antenna;
  os = new output_sett;

  read_namelist(
      &bs->rtor, &bs->rp, &bs->B0, bs->path2profiles, &bs->calc_back, bs->flag_back,
      &bs->N, &bs->V_gal_sys, &bs->V_scale, &bs->m_i, &bs->zele, &bs->zion,
      es->fname, &es->search_flag, &es->rdim, &es->rfmin, &es->rfmax,
      &es->idim, &es->ifmin, &es->ifmax, &es->stop_flag, &es->eps_res,
      &es->eps_abs, &es->eps_rel, &es->delta, &es->test_roots, &es->Nguess,
      &es->kmin, &es->kmax, es->fstart, &this->flag_debug);

  bs->mass[0] = bs->m_i * m_p; // ions mass

  as->read_settings(path2project);
  copy_antenna_data_to_antenna_module_(&as);
  copy_background_data_to_background_module_(&bs);
  os->read_settings(path2project);
}

settings_t::~settings_t() {
  delete[] es->fstart;

  delete as;
  delete bs;
  delete es;
  delete os;

  delete[] path2project;
}
