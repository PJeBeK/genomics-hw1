#include <iostream>
#include <math.h>
#include <sdsl/bit_vectors.hpp>
using namespace std;

struct rank_support
{
  sdsl::bit_vector b, rs, rb;
  uint64_t n, k, rs_width, rb_width;
  rank_support() {};

  rank_support(uint64_t n)
  {
    b = sdsl::bit_vector(n, 0);
    this->n = n;
  }
  
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

  void fsave(ofstream& fout)
  {
    b.serialize(fout);
    rs.serialize(fout);
    rb.serialize(fout);
  }

  void save(string& fname)
  {
    ofstream fout(fname);
    fsave(fout);
  }
  
  void fload(ifstream& fin)
  {
    b.load(fin);
    rs.load(fin);
    rb.load(fin);
    compute_lengths();
  }

  void load(string& fname)
  {
    ifstream fin(fname);
    fload(fin);
  }
};

struct sparse_array
{
  rank_support r;
  vector<string> elements;
  bool initiated = false;
  
  void create(uint64_t size)
  {
    r = rank_support(size);
    initiated = false;
  }

  void append(string elem, uint64_t pos)
  {
    if (pos >= r.n || r.b[pos] == 1)
      {
	cout<<"Warning: You can't append element at position: "<< pos <<endl;
	return;
      }
    initiated = false;
    elements.push_back(elem);
    r.b[pos] = 1;
  }

  void check_init()
  {
    if(!initiated)
      {
	r.init();
	initiated = true;
      }
  }
  
  bool get_at_rank(uint64_t i, string& elem)
  {
    check_init();
    if (i < 0 || elements.size() <= i)
      return false;
    elem = elements.at(i);
    return true;
  }

  bool get_at_index(uint64_t i, string& elem)
  {
    check_init();
    if (i < 0 || i >= r.n || r.b[i] == 0)
      return false;
    elem = elements.at(r(i) - 1);
    return true;
  }

  uint64_t num_elem_at(uint64_t i)
  {
    check_init();
    return r(i);
  }

  uint64_t size()
  {
    return r.n;
  }

  uint64_t num_elem()
  {
    return elements.size();
  }
  
  void save(string& fname)
  {
    ofstream fout(fname);
    r.fsave(fout);
    fout<<elements.size()<<" "<<initiated<<endl;
    for(int i=0;i<elements.size();i++)
      fout<<elements[i]<<" ";
    fout<<endl;
   }

  void load(string& fname)
  {
    ifstream fin(fname);
    r.fload(fin);
    int elements_len;
    fin>>elements_len>>initiated;
    elements.resize(elements_len);
    for(int i=0;i<elements_len;i++)
      fin>>elements[i];
   }
};

int main()
{
  /*
    sparse_array sa;
    sa.create(10);
    sa.append("foo", 1);
    sa.append("bar", 5);
    sa.append("baz", 9);
    string s1,s2,s3;
    cout<<sa.get_at_rank(1, s1)<<" "<<sa.get_at_index(3, s2)<<" "<<sa.get_at_index(5, s3)<<endl;
    cout<<s1<<endl<<s2<<endl<<s3<<endl;
    string fname = "in.txt";
    sa.save(fname);
    sa.load(fname);
    cout<<sa.size()<<" "<<sa.num_elem()<<endl;
    for(int i=0;i<=3;i++){
    string tmp;
    cout<<"rank "<<i<<" "<<sa.get_at_rank(i, tmp) <<" "<<tmp<<endl;
    }
    for(int i=0;i<=10;i++){
    string tmp;
    cout<<"index "<<i<<" "<<sa.num_elem_at(i)<<" "<<sa.get_at_index(i, tmp) <<" "<<tmp<<endl;
    }
  */
}
