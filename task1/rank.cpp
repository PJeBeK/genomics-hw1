#include <iostream>
#include <math.h>
#include <sdsl/bit_vectors.hpp>
using namespace std;

struct rank_support
{
  sdsl::bit_vector b, rs, rb;
  uint64_t n, k, rs_width, rb_width;
  rank_support(sdsl::bit_vector *b)
  {
    this->b = *b;
    init();
  }

  uint64_t operator () (uint64_t i)
  {
    return rank1(i);
  }

  void compute_lengths()
  {
    n = b.size();
    rs_width = k = ceil(log2(b.size()));
    rb_width = ceil(log2(k * k));
  }
  
  void init()
  {
    compute_lengths();
    uint64_t rs_len = rs_width * ceil(n / (k * k));
    uint64_t rb_len = rb_width * ceil(n / k);
    rs = sdsl::bit_vector(rs_len);
    rb = sdsl::bit_vector(rb_len);
    for (int i = 0; i*k*k < n; i ++)
      {
	write_on_rs(i, (i==0 ? 0 : read_from_rs(i - 1) + read_from_rb(i*k - 1) + count_ones_on_b(i*k*k - k, k)));
	for (int j = 0; j < k && i*k*k + j*k < n; j ++)
	  write_on_rb(i*k + j, (j==0 ? 0 : read_from_rb(i*k + j - 1) + count_ones_on_b(i*k*k + (j-1)*k, k)));
      }
  }

  uint64_t count_ones_on_b(uint64_t st, uint64_t len)
  {
    return __popcount(b.get_int(st, len));
  }
  
  uint64_t read_from_rs(uint64_t idx)
  {
    return rs.get_int(idx * rs_width, rs_width);
  }

  uint64_t read_from_rb(uint64_t idx)
  {
    return rb.get_int(idx * rb_width, rb_width);
  }

  void write_on_rs(uint64_t idx, uint64_t value)
  {
    rs.set_int(idx * rs_width, value, rs_width);
  }

  void write_on_rb(uint64_t idx, uint64_t value)
  {
    rb.set_int(idx * rb_width, value, rb_width);
  }
  
  uint64_t rank1(uint64_t i)
  {
    i ++;
    uint64_t rs_idx = i / (k * k);
    uint64_t rb_idx = i / k;
    return read_from_rs(rs_idx) + read_from_rb(rb_idx) + count_ones_on_b(i - (i % k), i % k);
  }

  uint64_t overhead()
  {
    uint64_t rs_len = rs_width * ceil(n / (k * k));
    uint64_t rb_len = rb_width * ceil(n / k);
    return rs_len + rb_len;
  }

  void save(string& fname)
  {
    ofstream fout(fname);
    b.serialize(fout);
    rs.serialize(fout);
    rb.serialize(fout);
  }

  void load(string& fname)
  {
    ifstream fin(fname);
    b.load(fin);
    rs.load(fin);
    rb.load(fin);
    compute_lengths();
  }
};

int main()
{
  /*
    int N=100;
    sdsl::bit_vector b(N, 0);
    for(int i=0;i<N;i++)
    b[i] = rand()%2;
    rank_support r(&b);
    string fname = "in.txt";
    r.save(fname);
    r.load(fname);
    cout<<r.overhead()<<endl;
    int sum=0;
    for(int i=0;i<N;i++){
    sum+=b[i];
    cout<<r(i)<<" "<<sum<<endl;
    }
  */
  const int n = 7;
  int N[] = {0, 1000000, 10000000, 100000000, 200000000, 500000000, 1000000000};
  int overhead[7];
  cout<<"CONSUMED TIME\nFormat: size_of_bit_vector number_of_rank_query consumed_time"<<endl;
  for(int i=1;i<n;i++)
    {
      sdsl::bit_vector b(N[i], 0);
      for(int j=0;j<N[i];j++)
	b[i] = rand()%2;
      rank_support r(&b);

      const clock_t start_time = clock();
      for(int j=1;j<n;j++)
	{
	  for(int k=N[j-1]; k<N[j]; k++)
	    r(rand()%N[i]);
	  cout<<N[i]<<" "<<N[j]<<" "<<float(clock() - start_time) / CLOCKS_PER_SEC<<endl;
	}
      overhead[i] = r.overhead();
    }
  cout<<"\nOVERHEAD\nFormat: size_of_bit_vector overload"<<endl;
  for(int i=1;i<n;i++)
    cout<<N[i]<<" "<<overhead[i]<<endl;
}
