//
// This file is part of the libWetCloth open source project
//
// The code is licensed under the same terms as a Clear BSD License but further
// restricted to academic and non-commercial use (commercial licenses may be
// obtained by contacting the faculty of the Columbia Computer Graphics Group
// or Columbia Technology Ventures).
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the disclaimer
// below) provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its contributors may be used
// to endorse or promote products derived from this software without specific
// prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY THIS
// LICENSE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
// GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.


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
