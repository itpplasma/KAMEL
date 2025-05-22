#include "settings.h"

#include "constants.h"
#include "namelist.h"

#include <cstdlib>
#include <cstring>

settings_t::settings_t(char *path) {
  using complx = std::complex<double>;

  path2project = new char[1024];
  strncpy(path2project, path, 1024);

  bs = new back_sett;

  es = new eigmode_sett;
  es->fstart = static_cast<complx*>(malloc(100 * sizeof(complx)));

  as = new antenna_sett;
  os = new output_sett;

  read_namelist(
      &bs->rtor, &bs->rp, &bs->B0, bs->path2profiles, &bs->calc_back, bs->flag_back,
      &bs->N, &bs->V_gal_sys, &bs->V_scale, &bs->m_i, &bs->zele, &bs->zion,
      es->fname, &es->search_flag, &es->rdim, &es->rfmin, &es->rfmax,
      &es->idim, &es->ifmin, &es->ifmax, &es->stop_flag, &es->eps_res,
      &es->eps_abs, &es->eps_rel, &es->delta, &es->test_roots, &es->Nguess,
      &es->kmin, &es->kmax, es->fstart, &this->flag_debug);

  // complete background settings
  bs->mass[0] = bs->m_i * m_p; // ions mass
  bs->flag_debug = this->flag_debug;

  as->read_settings(path2project);
  // complete eigmode settings
  es->fstart = static_cast<complx*>(realloc(es->fstart, es->Nguess));
  es->flag_debug = this->flag_debug;

  copy_antenna_data_to_antenna_module_(&as);
  copy_background_data_to_background_module_(&bs);
  os->read_settings(path2project);
}

settings_t::~settings_t() {
  free(es->fstart);

  delete as;
  delete bs;
  delete es;
  delete os;

  delete[] path2project;
}
