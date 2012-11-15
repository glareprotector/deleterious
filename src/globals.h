#ifndef GLOBALS_H
#define GLOBALS_H

#include <string>
#include <vector>
#include <ext/hash_map>

namespace __gnu_cxx{                                                                                             
  template<> struct hash< std::string >{                                                                                           
    size_t operator()( const std::string& x ) const{                                                                                         
      return hash< const char* >()( x.c_str() );                                              
    }                                                                                         
  };                                                                                          
}     


using namespace std;

hash_map<string, int> asdf;

typedef hash_map<string, int> ih;
typedef hash_map<string, double> dh; 


namespace globals{
  extern string base_folder;
  extern string protein_list_file;
  extern string results_folder;
  extern string all_aa;
}

#endif
