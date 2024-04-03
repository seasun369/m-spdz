#include "emp-tool/emp-tool.h"
#include "emp-zk/emp-zk.h"
#include "../emp-vole/base_vole.h"
#include "../emp-vole/vole_triple_new.h"
#include <iostream>
#include <fstream> // for file I/O，这个头文件包含了iostream头文件，所以此处可以不显示包含iostream头文件
#include <cstdlib> // support for exit()
#if defined(__linux__)
#include <sys/time.h>
#include <sys/resource.h>
#elif defined(__APPLE__)
#include <unistd.h>
#include <sys/resource.h>
#include <mach/mach.h>
#endif

using namespace emp;
using namespace std;

int party, port;
const int threads = 4;
int nparty = 2; //number of party
int tau = 4; 
int mat_sz = 128;

/*
void test_vole_triple(NetIO *ios[threads + 1], int party) {
  VoleTriple<NetIO> vtriple(party, threads, ios);

  __uint128_t Delta = (__uint128_t)0;
  if (party == ALICE) {
    PRG prg;
    prg.random_data(&Delta, sizeof(__uint128_t));
    Delta = Delta & ((__uint128_t)0xFFFFFFFFFFFFFFFFLL);
    Delta = mod(Delta, pr);
    auto start = clock_start();
    vtriple.setup(Delta);
    std::cout << "setup " << time_from(start) / 1000 << " ms" << std::endl;
    vtriple.check_triple(Delta, vtriple.pre_yz, vtriple.param.n_pre);
  } else {
    auto start = clock_start();
    vtriple.setup();
    std::cout << "setup " << time_from(start) / 1000 << " ms" << std::endl;
    vtriple.check_triple(0, vtriple.pre_yz, vtriple.param.n_pre);
  }

  int triple_need = vtriple.ot_limit;
  auto start = clock_start();
  __uint128_t *buf = new __uint128_t[triple_need];
  if (party == ALICE) {
    vtriple.extend(buf, triple_need);
    std::cout << triple_need << "\t" << time_from(start) / 1000 << " ms"
              << std::endl;
    vtriple.check_triple(Delta, buf, triple_need);
  } else {
    vtriple.extend(buf, triple_need);
    std::cout << triple_need << "\t" << time_from(start) / 1000 << " ms"
              << std::endl;
    vtriple.check_triple(0, buf, triple_need);
  }
  delete[] buf;

  // triple generation inplace
  uint64_t triple_need_inplace = vtriple.ot_limit;
  uint64_t memory_need = vtriple.byte_memory_need_inplace(triple_need_inplace);
  buf = new __uint128_t[memory_need];
  start = clock_start();
  vtriple.extend_inplace(buf, memory_need);
  std::cout << triple_need_inplace << "\tinplace\t" << memory_need << "\t"
            << time_from(start) / 1000 << " ms" << std::endl;
  if (party == ALICE)
    vtriple.check_triple(Delta, buf, memory_need);
  else
    vtriple.check_triple(0, buf, memory_need);

#if defined(__linux__)
  struct rusage rusage;
  if (!getrusage(RUSAGE_SELF, &rusage))
    std::cout << "[Linux]Peak resident set size: " << (size_t)rusage.ru_maxrss
              << std::endl;
  else
    std::cout << "[Linux]Query RSS failed" << std::endl;
#elif defined(__APPLE__)
  struct mach_task_basic_info info;
  mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
  if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info,
                &count) == KERN_SUCCESS)
    std::cout << "[Mac]Peak resident set size: "
              << (size_t)info.resident_size_max << std::endl;
  else
    std::cout << "[Mac]Query RSS failed" << std::endl;
#endif
}
*/

void gen_key(__uint128_t *Delta){
  PRG prg;
  __int128_t delta = (__uint128_t)0;
  for(int i=0;i < mat_sz; i++){
    prg.random_data(&delta, sizeof(__uint128_t));
    Delta[i] = delta& ((__uint128_t)0xFFFFFFFFFFFFFFFFLL);
    Delta[i] = mod(Delta[i], pr);
  }
}

void gen_tuple_b(__uint128_t *pre_b){
  for(int i=0;i < tau*mat_sz*mat_sz; i++){
    PRG prg;
    __int128_t delta = (__uint128_t)0;
    prg.random_data(&delta, sizeof(__uint128_t));
    pre_b[i] = delta& ((__uint128_t)0xFFFFFFFFFFFFFFFFLL);
    pre_b[i] = mod(Delta[i], pr);
  }
}

// w = vx + u
void rMVOPE(NetIO *ios[threads + 1], int party, __uint128_t *v, __uint128_t *u, __uint128_t *w, __uint128_t *x){
  for(int i = 0; i < mat_sz; i++){
  VoleTriple_N<NetIO> vtriple(party, threads, ios, x);

  if (party == ALICE) {
    // auto start = clock_start();
    vtriple.setup(v[i]);
    // std::cout << "setup " << time_from(start) / 1000 << " ms" << std::endl;
    vtriple.check_triple(v[i], vtriple.pre_yz, vtriple.param.n_pre);
  } else {
    // auto start = clock_start();
    vtriple.setup();
    // std::cout << "setup " << time_from(start) / 1000 << " ms" << std::endl;
    vtriple.check_triple(0, vtriple.pre_yz, vtriple.param.n_pre);
  }

  int triple_need = vtriple.ot_limit;
  // auto start = clock_start();
  __uint128_t *buf = new __uint128_t[triple_need];
  if (party == ALICE) {
    vtriple.extend(buf, triple_need);
    for (int j = 0; j < mat_sz*tau; ++j) {
        w[j] = mod(w[j] + buf[j], pr);
    }
    //std::cout << triple_need << "\t" << time_from(start) / 1000 << " ms"
    //          << std::endl;
    vtriple.check_triple(v[i], buf, triple_need);
  } else {
    vtriple.extend(buf, triple_need);
    for (int j = 0; j < mat_sz*tau; ++j) {
        x[i*mat_sz*tau+j] = mod((buf[i] >> 64), pr);
        u[j] = mod(u[j] + (buf[j] & 0xFFFFFFFFFFFFFFFFLL), pr); //haven't add negative
    }
    //std::cout << triple_need << "\t" << time_from(start) / 1000 << " ms"
    //          << std::endl;
    vtriple.check_triple(0, buf, triple_need);
  }
  delete[] buf;
  delete[] vtriple;
  }
}

// w = vx + u
void Auth(NetIO *io, int party, __uint128_t v, __uint128_t *u, __uint128_t *w, uint64_t *x){
    Base_vole<NetIO> *vole;
    __uint128_t *mac = new __uint128_t[mat_sz];

    if (party == ALICE) {
      vole = new Base_vole<NetIO>(party, io, v);
      vole->triple_gen_send(mac, mat_sz);
      for (int j = 0; j < mat_sz; ++j) {
        w[j] = mod(w[j] + mac[j], pr);
      }
    } else {
      vole = new Base_vole<NetIO>(party, io, x);
      vole->triple_gen_recv(mac, mat_sz);
      for (int j = 0; j < mat_sz; ++j) {
        u[j] = mod(u[j] + (mac[j] & 0xFFFFFFFFFFFFFFFFLL), pr);
      }
    }
}

int main(int argc, char **argv) {

__uint128_t *Delta = new __int128_t[mat_sz];
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

  parse_party_and_port(argv, &party, &port);

  // sample key
  gen_key(Delta);

  // sample [B]
  gen_tuple_b(pre_b);

  for (int i = 0; i < nparty; i++)
  {
    NetIO *ios[threads];
    for (int i = 0; i < threads; ++i)
    ios[i] = new NetIO(party == ALICE ? nullptr : "127.0.0.1", port + i);

    rMVOPE(ios, party, );

    for (int i = 0; i < threads; ++i) delete ios[i];

    party++;
    if (party > nparty)
    {
        party = 1;
    }
  }
  return 0;
}