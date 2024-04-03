#ifndef _VOLE_OFFLINE_H_
#define _VOLE_OFFLINE_H_
#include "base_svole.h"
#include "lpn.h"
#include "mpfss_reg.h"
#include "vole_triple.h"

template <typename IO> class Offline {
public:
  IO *io;
  IO **ios;
  int party;
  int threads;
  int mat_sz;
  /// @brief 
  __uint128_t *Delta;
  __uint128_t *pre_a = nullptr;
  __uint128_t *pre_b = nullptr;
  __uint128_t *pre_c = nullptr;
  __uint128_t *pre_a_t = nullptr;
  __uint128_t *pre_r = nullptr;
  __uint128_t *pre_r_t = nullptr;
  __uint128_t *mac_a = nullptr;
  __uint128_t *mac_b = nullptr;
  __uint128_t *mac_c = nullptr;
  __uint128_t *mac_a_t = nullptr;
  __uint128_t *mac_r = nullptr;
  __uint128_t *mac_r_t = nullptr;

  ThreadPool *pool=nullptr;

  Offline(int party, int threads, IO **ios, int mat_sz) {
    this->io = ios[0];
    this->threads = threads;
    this->party = party;
    this->ios = ios;
    this->mat_sz=mat_sz;
    Delta = new __uint128_t[mat_sz];
    pool = new ThreadPool(threads);
  }

  ~Offline() {
    if (pre_a != nullptr)
      delete[] pre_a;
    if (pre_b != nullptr)
      delete[] pre_b;
    if (pre_c != nullptr)
      delete pre_c;
    if (pre_a_t != nullptr)
      delete pre_a_t;
    if (pre_r != nullptr)
      delete[] pre_r;
    if (pre_r_t != nullptr)
      delete[] pre_r_t;
    if (mac_a != nullptr)
      delete mac_a;
    if (mac_b != nullptr)
      delete mac_b;
    if (mac_c != nullptr)
      delete[] mac_c;
    if (mac_a_t != nullptr)
      delete[] mac_a_t;
    if (mac_r != nullptr)
      delete mac_r;
    if (mac_r_t != nullptr)
      delete mac_r_t;
    if (pool != nullptr)
      delete pool;
  }

// generate global key Delta
void initialize(){
  PRG prg;
  __int128_t delta = (__uint128_t)0;
  prg.random_data(&Delta, sizeof(__uint128_t));
  Delta = Delta & ((__uint128_t)0xFFFFFFFFFFFFFFFFLL);
  Delta = mod(Delta, pr);
  for(int i=0;i < mat_sz; i++){
    prg.random_data(&delta, sizeof(__uint128_t));
    Delta[i] = delta& ((__uint128_t)0xFFFFFFFFFFFFFFFFLL);
    Delta[i] = mod(Delta[i], pr);
  }
}





}