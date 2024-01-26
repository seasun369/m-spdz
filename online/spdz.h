#ifndef SPDZ_H__
#define SPDZ_H__

#include "../network/network.h"
#include "utility.h"
#include <emp-tool/emp-tool.h>
using namespace emp;

// base SPDZ protocol with scalar elements
template<int nP>
class SPDZ { public:
	//const static int SSP = 5; 
	const block MASK = makeBlock(0x0ULL, 0xFFFFFULL);
	//Offline<nP>* fpre = nullptr; //not write
	//uint64_t* mac[nP+1];
	uint64_t key;
	//uint64_t* value[nP+1];

    //triple
    uint64_t *a,*b,*c;
	uint64_t *mac_a,*mac_b,*mac_c;

	//block * labels;
	//BristolFormat * cf;
	NetIOMP<nP> * io;
	//int num_ands = 0, num_in;
	int num_mul=1;
	int party, total_pre, ssp;
	ThreadPool * pool;
	//block Delta;
		
	//uint64_t* (*GTM)[4][nP+1];
	//uint64_t* (*GTK)[4];
	//uint64_t* (*GTv)[4][nP+1];
	//uint64_t* (*GT)[nP+1][4][nP+1];
	//block * eval_labels[nP+1];
	PRP prp;
	SPDZ(NetIOMP<nP> * io[2], ThreadPool * pool, int party, bool * _delta = nullptr, int ssp = 40) {
		this->party = party;
		this->io = io[0];
		this->ssp = ssp;
		this->pool = pool;

		if(party == 1) key = 7;
		if(party == 2) key = 13;
		if(party == 3) key = 17;//key values are artificially set for test

		//num_in = cf->n1+cf->n2;
		//total_pre = num_in + num_ands + 3*ssp;
		//fpre = new FpreMP<nP>(io, pool, party, _delta, ssp);

        a = new uint64_t[num_mul];
        b = new uint64_t[num_mul];
        c = new uint64_t[num_mul];
		mac_a = new uint64_t[num_mul];
        mac_b = new uint64_t[num_mul];
        mac_c = new uint64_t[num_mul];
	}
	~SPDZ () {
		delete[] a;
		delete[] b;
		delete[] c;
		delete[] mac_a;
		delete[] mac_b;
		delete[] mac_c;
	}
	PRG prg;

    // it should be implemented by offline.
    void get_triple(){
		std::ifstream fin1,fin2,fin3;
		fin1.open("/home/seasun/spdz/pre_data/triple_a.txt",std::ios::in);
		fin2.open("/home/seasun/spdz/pre_data/triple_b.txt",std::ios::in);
		fin3.open("/home/seasun/spdz/pre_data/triple_c.txt",std::ios::in);

		if(!(fin1.is_open() && fin2.is_open() && fin3.is_open()))
		{
		    std::cerr<<"cannot open the file";
		}
		char line1[1024]={0};
		char line2[1024]={0};
		char line3[1024]={0};
		std::vector<triple> Triple1,Triple2,Triple3;
		//从文件中提取“行”
		while(fin1.getline(line1,sizeof(line1)))
		{
		    //定义局部变量
		    triple t;
		    //从“行”中提取“单词”
		    std::stringstream word(line1);
		    word>>t.party;
		    uint64_t num;
		    while(word>>num)
		        t.value.push_back(num);
		    Triple1.push_back(t);
		}

		while(fin2.getline(line2,sizeof(line2)))
		{
		    //定义局部变量
		    triple t;
		    //从“行”中提取“单词”
		    std::stringstream word(line2);
		    word>>t.party;
		    uint64_t num;
		    while(word>>num)
		        t.value.push_back(num);
		    Triple2.push_back(t);
		}

		while(fin3.getline(line3,sizeof(line3)))
		{
		    //定义局部变量
		    triple t;
		    //从“行”中提取“单词”
		    std::stringstream word(line3);
		    word>>t.party;
		    uint64_t num;
		    while(word>>num)
		        t.value.push_back(num);
		    Triple3.push_back(t);
		}

		for(int i=0; i< nP; i++){
			if(party == stoi(Triple1[i].party)){
				a[0] = Triple1[i].value[0];
				mac_a[0] = Triple1[i].value[1];
			}
			if(party == stoi(Triple2[i].party)){
				b[0] = Triple2[i].value[0];
				mac_b[0] = Triple2[i].value[1];
			}
			if(party == stoi(Triple3[i].party)){
				c[0] = Triple3[i].value[0];
				mac_c[0] = Triple3[i].value[1];
			}
		}
	}

	void Online_mul (uint64_t * x, uint64_t * y, uint64_t * mac_x, uint64_t * mac_y, uint64_t * output, uint64_t *output_mac) {
		uint64_t *d = new uint64_t[num_mul];
        uint64_t *e = new uint64_t[num_mul];
        uint64_t *mac_d = new uint64_t[num_mul];
        uint64_t *mac_e = new uint64_t[num_mul];

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
			uint64_t * tmp[nP+1];
			for(int i = 1; i <= nP; ++i) tmp[i] = new uint64_t[num_mul];
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
			uint64_t * tmp[nP+1];
			for(int i = 1; i <= nP; ++i) tmp[i] = new uint64_t[num_mul];
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
	
    void mac_check(uint64_t * x, uint64_t * mac_x){
        if(party != 1) {
			io->send_data(1, x, num_mul);
			io->flush(1);
			io->recv_data(1, x, num_mul);
		} else {
			uint64_t * tmp[nP+1];
			for(int i = 1; i <= nP; ++i) tmp[i] = new uint64_t[num_mul];
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
            x=tmp[1];
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, x, party2]() {
					io->send_data(party2, x, num_mul);
					io->flush(party2);
				}));
			}
			joinNclean(res);
			for(int i = 1; i <= nP; ++i) delete[] tmp[i];
		}
        uint64_t * sigma_x = new uint64_t[num_mul];
        for(int i = 0; i < num_mul; ++i) {
            sigma_x[i] = mac_x[i] - mult_mod(x[i], key);
        }
        if(party != 1) {
			io->send_data(1, sigma_x, num_mul);
			io->flush(1);
		} else {
			uint64_t * tmp[nP+1];
			for(int i = 1; i <= nP; ++i) tmp[i] = new uint64_t[num_mul];
			vector<future<void>> res;
			for(int i = 2; i <= nP; ++i) {
				int party2 = i;
				res.push_back(pool->enqueue([this, tmp, party2]() {
					io->recv_data(party2, tmp[party2], num_mul);
				}));
			}
			joinNclean(res);
			for(int i = 0; i < num_mul; ++i){
				for(int j = 2; j <= nP; ++j)
					tmp[1][i] = add_mod(tmp[1][i] , tmp[j][i]);
                if(tmp[1][i] != 0) std::cout <<"check error"<< std::endl; 
            }
			for(int i = 1; i <= nP; ++i) delete[] tmp[i];
		}
    }
};
#endif// SPDZ_H__