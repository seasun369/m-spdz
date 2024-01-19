#ifndef MSPDZ_H__
#define MSPDZ_H__

#include "spdz/network/network.h"
#include <emp-tool/emp-tool.h>
using namespace emp;

// Matrix SPDZ
template<int nP>
class MSPDZ { public:
	const static int SSP = 5;// 
	const block MASK = makeBlock(0x0ULL, 0xFFFFFULL);
	//Offline<nP>* fpre = nullptr; //not write
	block* mac[nP+1];
	block key;
	block* value[nP+1];

    //triple
    block* a,b,c;

	//block * labels;
	//BristolFormat * cf;
	NetIOMP<nP> * io;
	//int num_ands = 0, num_in;
	int party, total_pre, ssp;
	ThreadPool * pool;
	//block Delta;
		
	block* (*GTM)[4][nP+1];
	block* (*GTK)[4];
	block* (*GTv)[4][nP+1];
	//block* (*GT)[nP+1][4][nP+1];
	//block * eval_labels[nP+1];
	PRP prp;
	Online(NetIOMP<nP> * io[2], ThreadPool * pool, int party, bool * _delta = nullptr, int ssp = 40) {
		this->party = party;
		this->io = io[0];
		this->ssp = ssp;
		this->pool = pool;

		//num_in = cf->n1+cf->n2;
		//total_pre = num_in + num_ands + 3*ssp;
		//fpre = new FpreMP<nP>(io, pool, party, _delta, ssp);

		if(party == 1) {
			GTM = new block[num_ands][4][nP+1]; //...not down
			GTK = new block[num_ands][4][nP+1];
			GTv = new block[num_ands][4];
			GT = new block[num_ands][nP+1][4][nP+1];
		}

		//labels = new block[cf->num_wire];
		for(int i  = 1; i <= nP; ++i) {
			key[i] = new block[cf->num_wire];
			mac[i] = new block[cf->num_wire];
		}
		value = new block[cf->num_wire];
        a = new block[];
        b = new block[];
        c = new block[];
	}
	~Online() {
		delete fpre;
		if(party == 1) {
			delete[] GTM;
			delete[] GTK;
			delete[] GTv;
			delete[] GT;
		}
		for(int i = 1; i <= nP; ++i) {
			delete[] key[i];
			delete[] mac[i];
		}
		delete[] value;
	}
	PRG prg;

    // it should be implemented by offline.
    void gen_triple()

	void Online_mul (block * x, block * y, block * mac_x, block * mac_y, block * output, block *output_mac, int num_mul) {
		block *d = new block[num_mul];
        block *e = new block[num_mul];
        block *mac_d = new block[num_mul];
        block *mac_e = new block[num_mul];

        for(int i=0; i<num_mul; ++i) {
            d[i]=x[i]-a[i];
            e[i]=y[i]-b[i];
            mac_d[i]=mac_x[i]-mac_a[i];
            mac_e[i]=mac_y[i]-mac_b[i];
        }
        // partially open
		if(party != 1) {
			io->send_data(1, d, num_mul);
			io->flush(1);
			io->recv_data(1, d, num_mul);
		} else {
			block * tmp[nP+1];
			for(int i = 1; i <= nP; ++i) tmp[i] = new block[num_mul];
			vector<future<void>> res;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], num_mul);
				}));
			}
			joinNclean(res);
			for(int i = 0; i < num_mul; ++i)
				for(int j = 2; j <= nP; ++j)
					tmp[1][i] = add_mod(tmp[1][i] , tmp[j][i]);
            d=tmp[1];
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, d, party2]() {
					io->send_data(party2, d, num_mul);
					io->flush(party2);
				}));
			}
			joinNclean(res);
			for(int i = 1; i <= nP; ++i) delete[] tmp[i];
		}

        if(party != 1) {
			io->send_data(1, e, num_mul);
			io->flush(1);
			io->recv_data(1, e, num_mul);
		} else {
			block * tmp[nP+1];
			for(int i = 1; i <= nP; ++i) tmp[i] = new block[num_mul];
			vector<future<void>> res;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], num_mul);
				}));
			}
			joinNclean(res);
			for(int i = 0; i < num_mul; ++i)
				for(int j = 2; j <= nP; ++j)
					tmp[1][i] = add_mod(tmp[1][i] , tmp[j][i]);
            e=tmp[1];
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, e, party2]() {
					io->send_data(party2, e, num_mul);
					io->flush(party2);
				}));
			}
			joinNclean(res);
			for(int i = 1; i <= nP; ++i) delete[] tmp[i];
		}
        for(int i = 0; i < num_mul; ++i){
            output[i] = add_mod(add_mod(c[i] , mult_mod(d[i], b[i])) , add_mod(mult_mod(a[i], e[i]) , mult_mod(d[i],e[i])));
            output_mac[i] = add_mod(add_mod(mac_c[i] , mult_mod(d[i], mac_b[i])) , add_mod(mult_mod(mac_a[i], e[i]) , mult_mod(mult_mod(d[i],e[i]),key)));
        }
	}
	
    void mac_check(block * x, block * mac_x, int size){
        if(party != 1) {
			io->send_data(1, x, size);
			io->flush(1);
			io->recv_data(1, x, size);
		} else {
			block * tmp[nP+1];
			for(int i = 1; i <= nP; ++i) tmp[i] = new block[size];
			vector<future<void>> res;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], size);
				}));
			}
			joinNclean(res);
			for(int i = 0; i < size; ++i)
				for(int j = 2; j <= nP; ++j)
					tmp[1][i] = add_mod(tmp[1][i] , tmp[j][i]);
            x=tmp[1];
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, x, party2]() {
					io->send_data(party2, x, size);
					io->flush(party2);
				}));
			}
			joinNclean(res);
			for(int i = 1; i <= nP; ++i) delete[] tmp[i];
		}
        block * sigma_x = new block *[size];
        for(int i = 0; i < size; ++i) {
            sigma_x[i] = mac_x[i] - mult_mod(x[i], key);
        }
        if(party != 1) {
			io->send_data(1, sigma_x, size);
			io->flush(1);
		} else {
			block * tmp[nP+1];
			for(int i = 1; i <= nP; ++i) tmp[i] = new block[size];
			vector<future<void>> res;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], size);
				}));
			}
			joinNclean(res);
			for(int i = 0; i < size; ++i){
				for(int j = 2; j <= nP; ++j)
					tmp[1][i] = add_mod(tmp[1][i] , tmp[j][i]);
                if(tmp[1][i] != 0) cout <<"check error"<<endl; // problem?
            }
			for(int i = 1; i <= nP; ++i) delete[] tmp[i];
		}
    }
};