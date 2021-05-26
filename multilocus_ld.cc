#include <ctime>
#include <random>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>

using namespace std;

map<int, vector<vector<int> > > recombine(map<int, vector<vector<int> > >&, double&);
void print_data(int g, ofstream&, map<int, vector<vector<int> > >&, double&);

double ab_r, bc_r, ef_r;
uniform_int_distribution<int> randomind(0,9999);
uniform_int_distribution<int> random_chromatid(0,3);
uniform_real_distribution<double> randomnum(0,1);
mt19937 engine(time(0));

map<int, vector<vector<int> > > recombine(map<int, vector<vector<int> > > &p, double &n) {
	map<int, vector<vector<int> > > recombined;
	for (int i=0; i < n; ++i){
		vector<vector<int> > parent1 = p[randomind(engine)]; //e.g., {0,0,0,0,0,1} & {0,1,1,0,0,1}
		vector<vector<int> > parent2 = p[randomind(engine)];

		// ab crossing-over
		// first parent
		if (randomnum(engine) <=  ab_r)
			swap(parent1[0][1], parent1[1][1]);
		// second parent
		if (randomnum(engine) <= ab_r)
			swap(parent2[0][1], parent2[1][1]);

		// bc crossing-over
		if (randomnum(engine) <= bc_r)
			swap(parent1[0][2], parent1[1][2]);
		if (randomnum(engine) <= bc_r)
			swap(parent2[0][2], parent2[1][2]);

		// ef crossing-over
		if (randomnum(engine) <= ef_r)
			swap(parent1[0][5], parent1[1][5]);
		if (randomnum(engine) <= ef_r)
			swap(parent2[0][5], parent2[1][5]);

		vector<int> chr1{randomnum(engine) < 0.5 ? 0 : 1, randomnum(engine) < 0.5 ? 0 : 1 };
		vector<int> chr3{randomnum(engine) < 0.5 ? 0 : 1, randomnum(engine) < 0.5 ? 0 : 1 };

		vector<int> variants1;
		vector<int> variants2;

		for (int j=0; j<3; ++j) {
			variants1.push_back(parent1[chr1[0]][j]);
			variants2.push_back(parent2[chr1[1]][j]);
		}

		variants1.push_back(parent1[randomnum(engine) < 0.5 ? 0 : 1][3]);
		variants2.push_back(parent2[randomnum(engine) < 0.5 ? 0 : 1][3]);

		for (int j=4; j<6; ++j) {
			variants1.push_back(parent1[chr3[0]][j]);
			variants2.push_back(parent2[chr3[1]][j]);
		}

		recombined[i].push_back(variants1);
		recombined[i].push_back(variants2);
	}
	return recombined;
}

void print_data(int g, ofstream &o, map<int, vector<vector<int> > > &p, double &n) {
	vector<vector<double> > coefs;
	vector<double> zeroes;
	zeroes.assign(6,0.);

	double nn = 2.*n; // number of alleles

	for (int i=0; i<6; ++i)
		coefs.push_back(zeroes);

	for (int i=0; i<5; ++i) {
		for (int j=i+1; j<6; ++j) {
			double freq1 = 0.;
			double freq2 = 0.;
			for (int k = 0; k<n; ++k) {
				for (int l=0; l<2; ++l) {
					freq1 += p[k][l][i];
					freq2 += p[k][l][j];
					if (p[k][l][i] ==1 && p[k][l][j] == 1)
						coefs[i][j]++;
				}
			}
			coefs[i][j] /= nn;
			coefs[i][j] -= (freq1/nn * freq2/nn);
		}
	}

	o << g << "\t";
	for (int i=0; i<5; ++i) {
		for (int j=i+1; j<6; ++j) {
			o << coefs[i][j] << "\t";
		}
	}
	o << endl;
}

int main (int argc, char *argv[] ) {
	int gen = atoi(argv[1]);
	double numind = atof(argv[2]);
	ab_r = atof(argv[3]);
	bc_r = atof(argv[4]);
	ef_r = atof(argv[5]);

	randomind.param(uniform_int_distribution<int>::param_type(0,numind-1));

	string fname("ld_coefs");
	ofstream ofile;
	ofile.open(fname.c_str());

	ofile << "gen\tDab\tDac\tDad\tDae\tDaf\tDbc\tDbd\tDbe\tDbf\tDcd\tDce\tDcf\tDde\tDdf\tDef" << endl;

	map<int, vector<vector<int> > > population;
   vector<int> p1{1,1,1,1,1,1};
	vector<int> p2{0,0,0,0,0,0};

	// parent population 1
	for (int i=0; i < numind/2; ++i) {
		population[i].push_back(p1);
		population[i].push_back(p1);
	}
	// parent population 2
	for (int i=numind/2; i<numind; ++i) {
		population[i].push_back(p2);
		population[i].push_back(p2);
	}

	for (int g=0; g<gen; ++g) {
		print_data(g, ofile, population, numind);
		population = recombine(population, numind);
	}

	ofile.close();

	return 0;
}
