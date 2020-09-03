//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and
// Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "StringUtilities.h"

namespace outputmod {

std::ostream& startred(std::ostream& stream) {
#ifndef __APPLE__
  stream << "\033[31;1m";
#endif
  return stream;
}

std::ostream& endred(std::ostream& stream) {
#ifndef __APPLE__
  stream << "\033[m";
#endif
  return stream;
}

std::ostream& startgreen(std::ostream& stream) {
#ifndef __APPLE__
  stream << "\033[32;1m";
#endif
  return stream;
}

std::ostream& endgreen(std::ostream& stream) {
#ifndef __APPLE__
  stream << "\033[m";
#endif
  return stream;
}

std::ostream& startpink(std::ostream& stream) {
#ifndef __APPLE__
  stream << "\033[35;1m";
#endif
  return stream;
}

std::ostream& endpink(std::ostream& stream) {
#ifndef __APPLE__
  stream << "\033[m";
#endif
  return stream;
}

std::ostream& startblue(std::ostream& stream) {
#ifndef __APPLE__
  stream << "\033[34;1m";
#endif
  return stream;
}

std::ostream& endblue(std::ostream& stream) {
#ifndef __APPLE__
  stream << "\033[m";
#endif
  return stream;
}

}  // namespace outputmod

namespace stringutils {

std::vector<std::string>& split(const std::string& s, char delim,
                                std::vector<std::string>& elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

std::vector<std::string> split(const std::string& s, char delim) {
  std::vector<std::string> elems;
  split(s, delim, elems);
  return elems;
}

void tokenize(const std::string& str, const char chr,
              std::vector<std::string>& tokens) {
  std::string::size_type substring_start = 0;
  std::string::size_type substring_end =
      str.find_first_of(chr, substring_start);
  while (substring_end != std::string::npos) {
    tokens.emplace_back(
        str.substr(substring_start, substring_end - substring_start));
    substring_start = substring_end + 1;
    substring_end = str.find_first_of(chr, substring_start);
  }
  // Grab the trailing substring, if present
  if (substring_start < str.size()) {
    tokens.emplace_back(str.substr(substring_start));
  }
  // Case of final character the delimiter
  if (str.back() == chr) {
    tokens.emplace_back("");
  }
}

std::vector<std::string> tokenize(const std::string& str,
                                  const char delimiter) {
  std::vector<std::string> tokens;
  tokenize(str, delimiter, tokens);
  return tokens;
}

void tokenize(const std::string& str, std::vector<std::string>& tokens,
              const std::string& delimiters) {
  // Skip delimiters at beginning.
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  std::string::size_type pos = str.find_first_of(delimiters, lastPos);

  while (std::string::npos != pos || std::string::npos != lastPos) {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}
}  // namespace stringutils
