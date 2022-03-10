#include <iostream>
#include <math.h>
#include <bit>
#include <cstdint>
#include <sdsl/bit_vectors.hpp>
using namespace std;
//using namespace sdsl;

int main()
{
  sdsl::bit_vector b = {0, 1, 0}, a = {1, 1, 1} , c = {0, 0} , d;
  ifstream fin("in.txt");
  ofstream fout("in.txt");
  //fout<<2<<" "<<3<<endl;
  a.serialize(fout);
  b.serialize(fout);
  fout.close();
  //cout<<a<<endl;
  //int t1,t2;
  //fin>>t1>>t2;
  c.load(fin);
  cout<<c<<endl;
  d.load(fin);
  cout<<d<<endl;
  //cout<<t1<<" "<<t2<<" "<<a<<endl;
}
