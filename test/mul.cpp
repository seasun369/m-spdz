#include "emp-tool/emp-tool.h"
#include "emp-ot/emp-ot.h"
#include "SPDZ/network/network.h"

using namespace emp;
using namespace std;

int party, port;
int matrix_sz=128;

__uint128_t* transform(__uint128_t* A){
  __uint128_t *A_T;
  for(int i=0; i < matrix_sz; ++i){
    for(int j =0; j < matrix_sz; ++j){
      __uint128_t tmp = A[i*matrix_sz+j];
      A_T[i*matrix_sz+j]=A[j*matrix_sz+i];
      A_T[j*matrix_sz+i]=tmp;
    }
  }
  return A_T;
}

void test_online_mul(NetIO *io, int party) {
    
  Online<NetIO> *online;

  online = new Online<NetIO>(party, io);
  
  long long test_n = matrix_sz * matrix_sz;

  __uint128_t *A, *B, *C, *R, *A_T, *R_T;
  A = new __uint128_t[test_n];
  B = new __uint128_t[test_n];
  C = new __uint128_t[test_n];
  R = new __uint128_t[test_n];

  get_sixtuple(NetIO *io, int party, uint128_t &A, uint128_t &B, uint128_t &C, uint128_t &R);

  A_T=transform(A);
  R_T=transform(R);

  __uint128_t *D, *E, *F;

  D = X - A;
  E = Y - B;
  F = transform(E)A_T-R_T;

  // test open
  auto start = clock_start();
  if (party == ALICE) {
    __uint128_t *k = new __uint128_t[test_n];
    io->recv_data(k, test_n * sizeof(__uint128_t));
  } else {
    io->send_data(D, test_n * sizeof(__uint128_t));
  }

  for (int i = 0; i < 10; ++i) {
    start = clock_start();
    if (party == ALICE) {
      online->mac_gen_send(mac, test_n);
      check_mac(io, &Delta, mac, test_n);
    } else {
      online->mac_gen_recv(mac, test_n);
      std::cout << "base mul: " << time_from(start) * 1000 / test_n
                << " ns per entry" << std::endl;
      check_mac(io, nullptr, mac, test_n);
    }
  }
  std::cout << "pass check" << std::endl;
}

int main(int argc, char **argv) {
  parse_party_and_port(argv, &party, &port);
  NetIO *io = new NetIO(party == ALICE ? nullptr : "127.0.0.1", port);

  std::cout << std::endl
            << "------------ BASE MUL ------------" << std::endl
            << std::endl;
  ;

  test_online_mul(io, party);

  delete io;
  return 0;
}