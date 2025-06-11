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

using namespace std;
typedef vector<double>   vec;

auto Qn_diatom(double V){
  double Qn=pow(10, -0.541 + 0.811*log10(V));
  Qn=Qn*pow(10,-6);
  Qn=Qn/12;
  Qn=Qn*16/106;
  return Qn;
}

auto Qn_euk(double V){
  double Qn=pow(10, -0.665 + 0.939*log10(V));
  Qn=Qn*pow(10,-6);
  Qn=Qn/12;
  Qn=Qn*16/106;
  return Qn;
}

auto Qn_cyano(double V){
  double Qn=470*V;
  Qn=Qn*pow(10,-9);
  Qn=Qn/12;
  Qn=Qn*16/106;
  return Qn;
}

