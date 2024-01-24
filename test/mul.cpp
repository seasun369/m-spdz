#include "emp-tool/emp-tool.h"
#include "emp-ot/emp-ot.h"
#include "spdz/network/network.h"
#include "spdz/online/spdz.h"

using namespace emp;
using namespace std;

int main(int argc, char **argv) {
  int port, party;
	parse_party_and_port(argv, &party, &port);

	const static int nP = 3;
	NetIOMP<nP> io(party, port);
	NetIOMP<nP> io2(party, port+2*(nP+1)*(nP+1)+1);
	NetIOMP<nP> *ios[2] = {&io, &io2};
	ThreadPool pool(4);	

  std::cout << std::endl
            << "------------ BASE MUL ------------" << std::endl
            << std::endl;
  ;

  test_online_mul(io, party);

  delete io;
  return 0;
}