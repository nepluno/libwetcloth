//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef __SPRING_UTILITIES_H__
#define __SPRING_UTILITIES_H__

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <fstream>

namespace Eigen { template <typename Derived> class DenseBase; }

namespace outputmod
{
  
  std::ostream& startred( std::ostream& stream );
  std::ostream& endred( std::ostream& stream );
  
  std::ostream& startgreen( std::ostream& stream );
  std::ostream& endgreen( std::ostream& stream );
  
  std::ostream& startpink( std::ostream& stream );
  std::ostream& endpink( std::ostream& stream );
  
  std::ostream& startblue( std::ostream& stream );
  std::ostream& endblue( std::ostream& stream );
  
}

#define FORCE_NO_PRINT
#define MAKE_REF(a) a, #a
//#define OVERRIDE

namespace stringutils
{
  
  template<class T>
  std::string convertToString( const T& tostring )
  {
    std::string out_string;
    std::stringstream ss;
    
    ss << tostring;
    ss >> out_string;
    
    return out_string;
  }
  
  template<class T>
  bool extractFromString( const std::string& in_string, T& output )
  {
    return (bool)(std::stringstream(in_string) >> output);
  }
  
  template<typename T>
  void print( const T& x, const char* name = NULL, bool override = false )
  {
#ifndef OVERRIDE
    if(override) {
#endif
      if(name) std::cout << std::string(name) << ":" << std::endl;
      std::cout << x << std::endl;
      std::cout << std::endl;
#ifndef OVERRIDE
    } else {
#ifndef NDEBUG
#ifndef FORCE_NO_PRINT
      if(name) std::cout << std::string(name) << ":" << std::endl;
      std::cout << x << std::endl;
      std::cout << std::endl;
#endif
#endif
    }
#endif
  }
  std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
  
  std::vector<std::string> split(const std::string &s, char delim);
  
  // Splits a string at given delimiter character
  void tokenize( const std::string& str, const char chr, std::vector<std::string>& tokens );
  std::vector<std::string> tokenize( const std::string& str, const char delimiter );
  void tokenize( const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters );

  // Reads a set number of values into an Eigen vector
  template<typename Derived>
  bool readList( const std::string& input_text, const char delimiter, Eigen::DenseBase<Derived>& list )
  {
    const std::vector<std::string> split_input = tokenize( input_text, delimiter );
    list.derived().resize( split_input.size() );

    for( unsigned entry_number = 0; entry_number < list.size(); ++entry_number )
    {
      if( !extractFromString( split_input[entry_number], list(entry_number) ) )
      {
        return false;
      }
    }
    return true;
  }  

  // Reads a set number of values into a std vector
  template<typename Derived>
  bool readList( const std::string& input_text, const char delimiter, std::vector<Derived>& list )
  {
    const std::vector<std::string> split_input = tokenize( input_text, delimiter );
    list.resize( split_input.size() );

    for( unsigned entry_number = 0; entry_number < list.size(); ++entry_number )
    {
      if( !extractFromString( split_input[entry_number], list[entry_number] ) )
      {
        return false;
      }
    }
    return true;
  }      
  
    
}

#endif
