#include "../emp-vole/base_vole.h"
#include "emp-tool/emp-tool.h"

using namespace emp;
using namespace std;

int party, port;

void check_triple(NetIO *io, __uint128_t *x, __uint128_t *y, int size) {
  if (party == ALICE) {
    io->send_data(x, sizeof(__uint128_t));
    io->send_data(y, size * sizeof(__uint128_t));
  } else {
    __uint128_t delta;
    __uint128_t *k = new __uint128_t[size];
    io->recv_data(&delta, sizeof(__uint128_t));
    io->recv_data(k, size * sizeof(__uint128_t));
    for (int i = 0; i < size; ++i) {
      __uint128_t tmp = mod(delta * (y[i] >> 64), pr);
      tmp = mod(tmp + k[i], pr);
      if (tmp != (y[i] & 0xFFFFFFFFFFFFFFFFLL)) {
        std::cout << "base_svole error" << std::endl;
        abort();
      }
    }
  }
}

void test_base_vole(NetIO *io, int party) {
  int test_n = 1024;
  __uint128_t *mac = new __uint128_t[test_n];
  __uint128_t *buf = new __uint128_t[test_n];
    PRG prg;
    uint64_t *x = new uint64_t[test_n + 1];
    prg.random_data(x, (test_n + 1) * sizeof(uint64_t));
    for (int i = 0; i < test_n + 1; ++i) {
      x[i] = mod(x[i]);
    }

  Base_vole<NetIO> *vole;

  __uint128_t Delta;
  if (party == ALICE) {
    PRG prg;
    prg.random_data(&Delta, sizeof(__uint128_t));
    Delta = Delta & ((__uint128_t)0xFFFFFFFFFFFFFFFFLL);
    Delta = mod(Delta, pr);

    vole = new Base_vole<NetIO>(party, io, Delta);
  } else {
    vole = new Base_vole<NetIO>(party, io, x);
  }

  // test single
  auto start = clock_start();
  if (party == ALICE) {
    vole->triple_gen_send(mac, test_n);
    vole->triple_gen_send(buf, test_n);
    check_triple(io, &Delta, mac, test_n);
  } else {
    //std::cout << x[0] << std::endl;
    vole->triple_gen_recv(mac, test_n);
    vole->triple_gen_recv(buf, test_n);
    //uint64_t tmp = (mac[0] >> 64);
    //std::cout << tmp << std::endl;
    std::cout << "base svole: " << time_from(start) * 1000 / test_n
              << " ns per entry" << std::endl;
    check_triple(io, nullptr, mac, test_n);
    std::cout << (uint64_t)(mac[0] >> 64) << std::endl;
    std::cout << (uint64_t)(buf[0] >> 64) << std::endl;
    if((__uint128_t)x[0] != (mac[0] >> 64)) abort();
  }
  std::cout << "pass check" << std::endl;

  delete[] mac;
}

int main(int argc, char **argv) {
  parse_party_and_port(argv, &party, &port);
  NetIO *io = new NetIO(party == ALICE ? nullptr : "127.0.0.1", port);

  std::cout << std::endl
            << "------------ BASE SVOLE ------------" << std::endl
            << std::endl;
  ;

  test_base_vole(io, party);

  delete io;
  return 0;
}