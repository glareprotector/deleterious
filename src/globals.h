#ifndef GLOBALS_H
#define GLOBALS_H

#include <string>
#include <vector>
#include <ext/hash_map>
#include <cstring>
#include <utility>
#include <unordered_map>
#include <cmath>
#include <assert.h>

using namespace std;

using namespace __gnu_cxx;

namespace __gnu_cxx{                                                                
  template<> struct hash< std::string >{
    size_t operator()( const std::string& x ) const{ 
      return hash< const char* >()( x.c_str() );                                              
    }                                                                                         
  };                                                                                          
}     

template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std
{
  template<typename S, typename T> struct hash<pair<S, T>>
    {
      inline size_t operator()(const pair<S, T> & v) const
      {
	size_t seed = 0;
	::hash_combine(seed, v.first);
	::hash_combine(seed, v.second);
	return seed;
      }
    };
}

//typedef __gnu_cxx::hash_map<string, int> ih;
//typedef __gnu_cxx::hash_map<string, double> dh; 


typedef unordered_map<string, int> ih;
typedef unordered_map<string, double> dh; 
typedef unordered_map<char, double> cdmap;
typedef unordered_map<char, int> cimap;
typedef unordered_map<pair<char, char>, int> ccimap;
typedef unordered_map<pair<char, char>, double> ccdmap;

namespace globals{
  extern string base_folder;
  extern string protein_list_file;
  extern string results_folder;
  extern string all_aa;
  extern string all_aa_no_skip;
  extern cimap all_aa_to_i;
  void init_globals();
}

#endif
