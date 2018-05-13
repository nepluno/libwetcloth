//
// This file is part of the libWetCloth open source project
//
// The code is licensed solely for academic and non-commercial use under the
// terms of the Clear BSD License. The terms of the Clear BSD License are
// provided below. Other licenses may be obtained by contacting the faculty
// of the Columbia Computer Graphics Group or a Columbia University licensing officer.
//
// We would like to hear from you if you appreciate this work.
//
// The Clear BSD License
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the disclaimer
// below) provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//  list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its contributors may be used
//  to endorse or promote products derived from this software without specific
//  prior written permission.
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


#include "StringUtilities.h"

namespace outputmod
{

std::ostream& startred( std::ostream& stream )
{
#ifndef __APPLE__
  stream << "\033[31;1m";
#endif
  return stream;
}

std::ostream& endred( std::ostream& stream )
{
#ifndef __APPLE__
  stream << "\033[m";
#endif
  return stream;
}

std::ostream& startgreen( std::ostream& stream )
{
#ifndef __APPLE__
  stream << "\033[32;1m";
#endif
  return stream;
}

std::ostream& endgreen( std::ostream& stream )
{
#ifndef __APPLE__
  stream << "\033[m";
#endif
  return stream;
}

std::ostream& startpink( std::ostream& stream )
{
#ifndef __APPLE__
  stream << "\033[35;1m";
#endif
  return stream;
}

std::ostream& endpink( std::ostream& stream )
{
#ifndef __APPLE__
  stream << "\033[m";
#endif
  return stream;
}

std::ostream& startblue( std::ostream& stream )
{
#ifndef __APPLE__
  stream << "\033[34;1m";
#endif
  return stream;
}

std::ostream& endblue( std::ostream& stream )
{
#ifndef __APPLE__
  stream << "\033[m";
#endif
  return stream;
}

}

namespace stringutils {
  
  std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
      elems.push_back(item);
    }
    return elems;
  }
  
  std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
  }

  void tokenize( const std::string& str, const char chr, std::vector<std::string>& tokens )
  {
    std::string::size_type substring_start = 0;
    std::string::size_type substring_end = str.find_first_of( chr, substring_start );
    while( substring_end != std::string::npos )
    {
      tokens.emplace_back( str.substr( substring_start, substring_end - substring_start ) );
      substring_start = substring_end + 1;
      substring_end = str.find_first_of( chr, substring_start );
    }
    // Grab the trailing substring, if present
    if( substring_start < str.size() )
    {
      tokens.emplace_back( str.substr( substring_start ) );
    }
    // Case of final character the delimiter
    if( str.back() == chr )
    {
      tokens.emplace_back( "" );
    }
  }

  std::vector<std::string> tokenize( const std::string& str, const char delimiter )
  {
    std::vector<std::string> tokens;
    tokenize( str, delimiter, tokens );
    return tokens;
  }  
  
}
