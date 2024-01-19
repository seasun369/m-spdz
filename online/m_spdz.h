#ifndef MSPDZ_H__
#define MSPDZ_H__

#include "spdz/network/network.h"
#include <emp-tool/emp-tool.h>
#include "utility.h"
using namespace emp;

// Matrix SPDZ
template<int nP>
class MSPDZ { public:
	const static int SSP = 5;// 
	const uint64_t MASK = makeuint64_t(0x0ULL, 0xFFFFFULL);
	//Offline<nP>* fpre = nullptr; //not write
	uint64_t* mac[nP+1];
	uint64_t key;
	uint64_t* value[nP+1];

    //triple
    uint64_t* a,b,c;

	//uint64_t * labels;
	//BristolFormat * cf;
	NetIOMP<nP> * io;
	//int num_ands = 0, num_in;
	int party, total_pre, ssp, mat_sz;
	ThreadPool * pool;
	//uint64_t Delta;
		
	uint64_t* (*GTM)[4][nP+1];
	uint64_t* (*GTK)[4];
	uint64_t* (*GTv)[4][nP+1];
	//uint64_t* (*GT)[nP+1][4][nP+1];
	//uint64_t * eval_labels[nP+1];
	PRP prp;
	Online(NetIOMP<nP> * io[2], ThreadPool * pool, int party, bool * _delta = nullptr, int ssp = 40, int mat_sz = 64) {
		this->party = party;
		this->io = io[0];
		this->ssp = ssp;
		this->pool = pool;

		//num_in = cf->n1+cf->n2;
		//total_pre = num_in + num_ands + 3*ssp;
		//fpre = new FpreMP<nP>(io, pool, party, _delta, ssp);

		if(party == 1) {
			GTM = new uint64_t[num_ands][4][nP+1]; //...not down
			GTK = new uint64_t[num_ands][4][nP+1];
			GTv = new uint64_t[num_ands][4];
			GT = new uint64_t[num_ands][nP+1][4][nP+1];
		}

		//labels = new uint64_t[cf->num_wire];
		for(int i  = 1; i <= nP; ++i) {
			key[i] = new uint64_t[cf->num_wire];
			mac[i] = new uint64_t[cf->num_wire];
		}
		value = new uint64_t[cf->num_wire];
        a = new uint64_t[];
        b = new uint64_t[];
        c = new uint64_t[];
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

	void Online_mul (uint64_t * x, uint64_t * y, uint64_t * mac_x, uint64_t * mac_y, uint64_t * output, uint64_t *output_mac) {
        int sz = mat_sz * mat_sz;
		uint64_t *d = new uint64_t[sz];
        uint64_t *e = new uint64_t[sz];
        uint64_t *f = new uint64_t[sz];
        uint64_t *mac_d = new uint64_t[sz];
        uint64_t *mac_e = new uint64_t[sz];
        uint64_t *mac_f = new uint64_t[sz];

        for(int i=0; i<sz; ++i) {
            d[i]=x[i]-a[i];
            e[i]=y[i]-b[i];
            mac_d[i]=mac_x[i]-mac_a[i];
            mac_e[i]=mac_y[i]-mac_b[i];
        }
        // partially open
		if(party != 1) {
			io->send_data(1, d, sz);
			io->flush(1);
			io->recv_data(1, d, sz);
		} else {
			uint64_t * tmp[nP+1];
			for(int i = 1; i <= nP; ++i) tmp[i] = new uint64_t[sz];
			vector<future<void>> res;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], sz);
				}));
			}
			joinNclean(res);
			for(int i = 0; i < sz; ++i)
				for(int j = 2; j <= nP; ++j)
					tmp[1][i] = add_mod(tmp[1][i] , tmp[j][i]);
            d=tmp[1];
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, d, party2]() {
					io->send_data(party2, d, sz);
					io->flush(party2);
				}));
			}
			joinNclean(res);
			for(int i = 1; i <= nP; ++i) delete[] tmp[i];
		}

        if(party != 1) {
			io->send_data(1, e, sz);
			io->flush(1);
			io->recv_data(1, e, sz);
		} else {
			uint64_t * tmp[nP+1];
			for(int i = 1; i <= nP; ++i) tmp[i] = new uint64_t[sz];
			vector<future<void>> res;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], sz);
				}));
			}
			joinNclean(res);
			for(int i = 0; i < sz; ++i)
				for(int j = 2; j <= nP; ++j)
					tmp[1][i] = add_mod(tmp[1][i] , tmp[j][i]);
            e=tmp[1];
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, e, party2]() {
					io->send_data(party2, e, sz);
					io->flush(party2);
				}));
			}
			joinNclean(res);
			for(int i = 1; i <= nP; ++i) delete[] tmp[i];
		}
        uint64_t *e_t = new uint64_t[sz];
        e_t = transform(e, mat_sz);
        f = mat_mult_mod(e_t,a_t,mat_sz);
        for(int i=0; i < sz ; ++i){
            f[i]=(f[i]-r_t[i] ) % pr;
        }

        if(party != 1) {
			io->send_data(1, f, sz);
			io->flush(1);
			io->recv_data(1, f, sz);
		} else {
			uint64_t * tmp[nP+1];
			for(int i = 1; i <= nP; ++i) tmp[i] = new uint64_t[sz];
			vector<future<void>> res;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], sz);
				}));
			}
			joinNclean(res);
			for(int i = 0; i < sz; ++i)
				for(int j = 2; j <= nP; ++j)
					tmp[1][i] = add_mod(tmp[1][i] , tmp[j][i]);
            f=tmp[1];
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, f, party2]() {
					io->send_data(party2, f, sz);
					io->flush(party2);
				}));
			}
			joinNclean(res);
			for(int i = 1; i <= nP; ++i) delete[] tmp[i];
		}
        uint64_t *f_t = new uint64_t[sz];
        uint64_t *db = new uint64_t[sz];
        uint64_t *de = new uint64_t[sz];
        uint64_t *mac_db = new uint64_t[mat_sz];
        uint64_t *mac_de = new uint64_t[mat_sz];
        f_t = transform(f, mat_sz);
        db = mat_mult_mod(d, b, mat_sz);
        de = mat_mult_mod(d, e, mat_sz);
        mac_db = mat_vec_mult(d, mac_b, mat_sz);
        mac_de = mat_vec_mult(de, key, mat_sz);

        for(int i = 0; i < sz; ++i){
            output[i] = mod(c[i]+db[i]+r[i]+de[i]+f_t[i]);
            output_mac[i] = mod(mac_c[i]+mac_db[i]+mac_r[i]+mac_de[i]+f_t[i]);
        }
	}
	
    void mac_check(uint64_t * x, uint64_t * mac_x){
        int size = mat_sz * mat_sz;
        if(party != 1) {
			io->send_data(1, x, size);
			io->flush(1);
			io->recv_data(1, x, size);
		} else {
			uint64_t * tmp[nP+1];
			for(int i = 1; i <= nP; ++i) tmp[i] = new uint64_t[size];
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
        uint64_t * sigma_x = new uint64_t *[mat_sz];
        sigma_x = mat_vec_mult(x,key,mat_sz);
        for(int i = 0; i < mat_sz; ++i) {
            sigma_x[i] = mod(mac_x[i] - sigma_x[i]);
        }
        if(party != 1) {
			io->send_data(1, sigma_x, mat_sz);
			io->flush(1);
		} else {
			uint64_t * tmp[nP+1];
			for(int i = 1; i <= nP; ++i) tmp[i] = new uint64_t[mat_sz];
			vector<future<void>> res;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], mat_sz);
				}));
			}
			joinNclean(res);
			for(int i = 0; i < mat_sz; ++i){
				for(int j = 2; j <= nP; ++j)
					tmp[1][i] = add_mod(tmp[1][i] , tmp[j][i]);
                if(tmp[1][i] != 0) cout <<"check error"<<endl; // problem?
            }
			for(int i = 1; i <= nP; ++i) delete[] tmp[i];
		}
    }
};