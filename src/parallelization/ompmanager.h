// -------------------------------------------------------------------------
// Copyright (C) Max Planck Institute of Biophysics - All Rights Reserved
// Unauthorized copying of this file, via any medium is strictly prohibited
// Proprietary and confidential
// The code comes without warranty of any kind
// Please refer to Kim and Hummer J.Mol.Biol. 2008
// -------------------------------------------------------------------------
#ifndef OMPMANAGER_H
#define OMPMANAGER_H

#if defined(_OPENMP)
#include <omp.h>
#else
inline int omp_get_max_threads() {
  return 1;
}
inline int omp_get_thread_num() {
  return 0;
}
inline int omp_get_num_threads() {
  return 1;
}

struct omp_lock_t {};

inline void omp_init_lock(omp_lock_t* /*lock*/) {
}

inline void omp_destroy_lock(omp_lock_t* /*lock*/) {
}

inline void omp_set_lock(omp_lock_t* /*lock*/) {
}

inline void omp_unset_lock(omp_lock_t* /*lock*/) {
}

#endif

#endif
