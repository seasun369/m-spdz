//one element mul test
#include "emp-tool/emp-tool.h"
#include "emp-ot/emp-ot.h"
#include "../network/network.h"
#include "../online/spdz.h"

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

  uint64_t *x, *y, *mac_x, *mac_y, *output, *output_mac;
  int num_mul=1;
  x = new uint64_t[num_mul];
  mac_x = new uint64_t[num_mul];
  y = new uint64_t[num_mul];
  mac_y = new uint64_t[num_mul];
  output = new uint64_t[num_mul];
  output_mac = new uint64_t[num_mul];

  std::ifstream fin1,fin2;
	fin1.open("/home/seasun/spdz/pre_data/input_x.txt",std::ios::in);
	fin2.open("/home/seasun/spdz/pre_data/input_y.txt",std::ios::in);

  if(!(fin1.is_open() && fin2.is_open()))
		{
		    std::cerr<<"cannot open the file";
		}
	char line1[1024]={0};
	char line2[1024]={0};

  std::vector<triple> Inputx,Inputy;
		while(fin1.getline(line1,sizeof(line1)))
		{
		    triple t;
		    std::stringstream word(line1);
		    word>>t.party;
		    uint64_t num;
		    while(word>>num)
		        t.value.push_back(num);
		    Inputx.push_back(t);
		}

		while(fin2.getline(line2,sizeof(line2)))
		{
		    triple t;
		    std::stringstream word(line2);
		    word>>t.party;
		    uint64_t num;
		    while(word>>num)
		        t.value.push_back(num);
		    Inputy.push_back(t);
		}

    	for(int i=0; i< nP; i++){
			if(party == stoi(Inputx[i].party)){
				x[0] = Inputx[i].value[0];
				mac_x[0] = Inputx[i].value[1];
			}
      		if(party == stoi(Inputy[i].party)){
        		y[0] = Inputy[i].value[0];
				mac_y[0] = Inputy[i].value[1];
      		}
		}

  std::cout << std::endl
            << "------------ BASE MUL ------------" << std::endl
            << std::endl;
  ;

  SPDZ<nP>* mpc = new SPDZ<nP>(ios, &pool, party);
  cout <<"Setup:\t"<<party<<"\n";

  mpc->get_triple();

  auto start = clock_start();
  mpc->Online_mul(x,y,mac_x,mac_y,output,output_mac);
  mpc->mac_check(output,output_mac);
  std::cout << "BASE MUL: " << time_from(start) * 1000
                << " ns" << std::endl;

  delete mpc;
  return 0;
}