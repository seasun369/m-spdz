// matrix mul test
#include "emp-tool/emp-tool.h"
#include "emp-ot/emp-ot.h"
#include "../network/network.h"
#include "../online/m_spdz_mpir_256.h"
#include "../mpir/mpir.h"

using namespace emp;
using namespace std;

int main(int argc, char **argv)
{
	int port, party;
	parse_party_and_port(argv, &party, &port);

	const static int nP = 2;
	NetIOMP<nP> io(party, port);
	NetIOMP<nP> io2(party, port + 2 * (nP + 1) * (nP + 1) + 1);
	NetIOMP<nP> *ios[2] = {&io, &io2};
	ThreadPool pool(4);

	mpz_t *x, *y, *mac_x, *mac_y, *output, *output_mac;
	int mat_sz = 256;
	int sz = mat_sz * mat_sz;
	x = new mpz_t[sz];
	y = new mpz_t[sz];
	output = new mpz_t[sz];
	mac_x = new mpz_t[mat_sz];
	mac_y = new mpz_t[mat_sz];
	output_mac = new mpz_t[mat_sz];

	std::ifstream fin1, fin2, fin3, fin4;
	fin1.open("/home/jackie/spdz/pre_data/predata_mspdz/256/input_x.txt", std::ios::in);
	fin2.open("/home/jackie/spdz/pre_data/predata_mspdz/256/input_y.txt", std::ios::in);
	fin3.open("/home/jackie/spdz/pre_data/predata_mspdz/256/mac_x.txt", std::ios::in);
	fin4.open("/home/jackie/spdz/pre_data/predata_mspdz/256/mac_y.txt", std::ios::in);

	if (!(fin1.is_open() && fin2.is_open() && fin3.is_open() && fin4.is_open()))
	{
		std::cerr << "cannot open the file";
	}

	char line1[4096000] = {0};
	std::vector<triple> Triple;
	while (fin1.getline(line1, sizeof(line1)))
	{
		triple t;
		std::stringstream word(line1);
		word >> t.party;
		uint64_t num;
		while (word >> num)
			t.value.push_back(num);
		Triple.push_back(t);
	}

	for (int i = 0; i < nP; i++)
	{
		if (party == stoi(Triple[i].party))
		{
			for (int j = 0; j < sz; j++)
			{
				mpz_init_set_ui(x[j], Triple[i].value[j]);
			}
		}
	}
	
	Triple.clear();
	memset(line1, 0, sizeof(line1));

	while (fin2.getline(line1, sizeof(line1)))
	{
		triple t;
		std::stringstream word(line1);
		word >> t.party;
		uint64_t num;
		while (word >> num)
			t.value.push_back(num);
		Triple.push_back(t);
	}

	for (int i = 0; i < nP; i++)
	{
		if (party == stoi(Triple[i].party))
		{
			for (int j = 0; j < sz; j++)
			{
				mpz_init_set_ui(y[j], Triple[i].value[j]);
			}
		}
	}

	Triple.clear();
	memset(line1, 0, sizeof(line1));

	while (fin3.getline(line1, sizeof(line1)))
	{
		triple t;
		std::stringstream word(line1);
		word >> t.party;
		uint64_t num;
		while (word >> num)
			t.value.push_back(num);
		Triple.push_back(t);
	}

	for (int i = 0; i < nP; i++)
	{
		if (party == stoi(Triple[i].party))
		{
			for (int j = 0; j < mat_sz; j++)
			{
				mpz_init_set_ui(mac_x[j], Triple[i].value[j]);
			}
		}
	}

	Triple.clear();
	memset(line1, 0, sizeof(line1));

	while (fin4.getline(line1, sizeof(line1)))
	{
		triple t;
		std::stringstream word(line1);
		word >> t.party;
		uint64_t num;
		while (word >> num)
			t.value.push_back(num);
		Triple.push_back(t);
	}

	for (int i = 0; i < nP; i++)
	{
		if (party == stoi(Triple[i].party))
		{
			for (int j = 0; j < mat_sz; j++)
			{
				mpz_init_set_ui(mac_y[j], Triple[i].value[j]);
			}
		}
	}

	Triple.clear();
	memset(line1, 0, sizeof(line1));

	std::cout << std::endl
			  << "------------ MATRIX MUL ------------" << std::endl
			  << std::endl;
	;
	
	double sum_time = 0;
	for (int i = 0; i < 30; i++)
	{
		MSPDZ<nP> *mpc = new MSPDZ<nP>(ios, &pool, party);
		mpc->get_triple();

		auto start = clock_start();
		mpc->Online_mul(x, y, mac_x, mac_y, output, output_mac);
		mpc->mac_check(output, output_mac);

		sum_time += time_from(start) * 1000;
		delete mpc;
	}

	for (int i = 0; i < sz; i++)
	{
		mpz_clear(x[i]);
		mpz_clear(y[i]);
		mpz_clear(output[i]);
	}
	for (int i = 0; i < mat_sz; i++)
	{
		mpz_clear(mac_x[i]);
		mpz_clear(mac_y[i]);
		mpz_clear(output_mac[i]);
	}
	delete[] x;
	delete[] y;
	delete[] output;
	delete[] mac_x;
	delete[] mac_y;
	delete[] output_mac;

	std::cout << "MATRIX MUL: " << sum_time / 30
			  << " ns" << std::endl;

	return 0;
}
