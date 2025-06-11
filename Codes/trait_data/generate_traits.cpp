#include <fstream>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <map>
#include <algorithm>
#include <iterator>
#include <random>
#include <cstring>
#include <string>
#include <set>
#include <omp.h>
#include <cmath>
#include "functions.h"

using namespace std;
typedef vector<double>   vec;

int main(int argc, char *argv[]) {
	vec Vs_dia;
        vec Vs_euk;
        vec Vs_syne;
        vec Vs_pico;
        vec mumax_dia;
        vec mumax_euk;
        vec mumax_syne;
        vec mumax_pico;
        vec Qns_dia;
        vec Qns_euk;
        vec Qns_syne;
        vec Qns_pico;

	int n_species=100; // 100 species per phytoplankton groups

	//dutkiewicz et al 2020 allometric coefficents
        const double a_diatom =3.9;
        const double a_euk =1.4;
        const double a_cyano=0.9;
        const double b =-0.08;
        const double b_cyano =0.08;

	// random seed (5)
	default_random_engine generator;
        string se;
        se += argv[1];
        cout << "seed:" << endl;
        cout << se << endl;

	int seed_r = stoi(se);
        generator.seed(seed_r);
        // volume distribution of the different phytoplankton groups
        uniform_real_distribution<double> distribution_d(2.0,6.0); // diatom
        uniform_real_distribution<double> distribution_e(0.6,5.0); // eukaryote
        uniform_real_distribution<double> distribution_s(-0.2,2.0); // Synechococcus
        uniform_real_distribution<double> distribution_p(-1.0,1.0); // Prochlorococcus
        for (int i=0; i<n_species; ++i) {
                double d = distribution_d(generator);
                double e = distribution_e(generator);
                double sy = distribution_s(generator);
                double p = distribution_p(generator);
                Vs_dia.push_back(exp(d*log(10.0)));
                Vs_euk.push_back(exp(e*log(10.0)));
                Vs_syne.push_back(exp(sy*log(10.0)));
                Vs_pico.push_back(exp(p*log(10.0)));
            
                mumax_dia.push_back(a_diatom*pow(Vs_dia[i], b));
                mumax_euk.push_back(a_euk*pow(Vs_euk[i], b));
                mumax_syne.push_back(a_cyano*pow(Vs_syne[i], b_cyano));
                mumax_pico.push_back(a_cyano*pow(Vs_pico[i], b_cyano));

                auto Q_dia=Qn_diatom(Vs_dia[i]);
                auto Q_euk=Qn_euk(Vs_euk[i]);
                auto Q_syne=Qn_cyano(Vs_syne[i]);
                auto Q_pico=Qn_cyano(Vs_pico[i]);
                Qns_dia.push_back(Q_dia);
                Qns_euk.push_back(Q_euk);
                Qns_syne.push_back(Q_syne);
                Qns_pico.push_back(Q_pico);
        }

        vec Vs_all;
        vec mumax_all;
        vec Qns_all;
	
	Vs_all.insert(Vs_all.end(), Vs_dia.begin(), Vs_dia.end());
        Vs_all.insert(Vs_all.end(), Vs_euk.begin(), Vs_euk.end());
        Vs_all.insert(Vs_all.end(), Vs_syne.begin(), Vs_syne.end());
        Vs_all.insert(Vs_all.end(), Vs_pico.begin(), Vs_pico.end());
        mumax_all.insert(mumax_all.end(), mumax_dia.begin(), mumax_dia.end());
        mumax_all.insert(mumax_all.end(), mumax_euk.begin(), mumax_euk.end());
        mumax_all.insert(mumax_all.end(), mumax_syne.begin(), mumax_syne.end());
        mumax_all.insert(mumax_all.end(), mumax_pico.begin(), mumax_pico.end());
        Qns_all.insert(Qns_all.end(), Qns_dia.begin(), Qns_dia.end());
        Qns_all.insert(Qns_all.end(), Qns_euk.begin(), Qns_euk.end());
        Qns_all.insert(Qns_all.end(), Qns_syne.begin(), Qns_syne.end());
        Qns_all.insert(Qns_all.end(), Qns_pico.begin(), Qns_pico.end());

	const double aK = 0.17;
        const double bK = 0.27;
        const double aQ = 0.07;
        const double bQ = -0.17;
        const double aVmax = 0.51;
        const double bVmax = -0.27;
        vec K;
        vec Qmin;
        vec Vmax;
        vec Nc_all;
        const int n_species_tot=Vs_all.size();
        double n_species_tot1=(double) n_species_tot;
        for (int i=0; i<n_species_tot; ++i) {
                K.push_back(aK*pow(Vs_all[i], bK));
                Qmin.push_back(aQ*pow(Vs_all[i], bQ));
                Vmax.push_back(aVmax*pow(Vs_all[i], bVmax));
                Nc_all.push_back(K[i]*mumax_all[i]*Qmin[i]/Vmax[i]);
        }
        ofstream file_V;
        file_V.open("Vs_"+se+".txt");
        for (int i=0; i<Vs_all.size(); ++i){
                file_V << Vs_all[i];
                file_V << " ";
        }
        file_V.close();

        ofstream file_mu;
        file_mu.open("mumax_dutkiewicz_"+se+".txt");
        for (int i=0; i<mumax_all.size(); ++i){
                file_mu << mumax_all[i];
                file_mu << " ";
        }
        file_mu.close();

        ofstream file_Nc;
        file_Nc.open("Nc_dutkiewicz_"+se+".txt");
        for (int i=0; i<Nc_all.size(); ++i){
                file_Nc << Nc_all[i];
                file_Nc << " ";
        }
        file_Nc.close();
}

