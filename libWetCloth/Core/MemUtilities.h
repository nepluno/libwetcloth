#ifndef MEM_UTILITIES_H
#define MEM_UTILITIES_H

//
// This file is part of the libWetCloth open source project
// Author:  David Robert Nadeau
// Site:    http://NadeauSoftware.com/
// License: Creative Commons Attribution 3.0 Unported License
//          http://creativecommons.org/licenses/by/3.0/deed.en_US
//

#if defined(_WIN32)
#include <Windows.h>
#endif

#if defined(_WIN32)
#include <psapi.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || \
    (defined(__APPLE__) && defined(__MACH__))
#include <sys/resource.h>
#include <unistd.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) || \
    (defined(__sun__) || defined(__sun) ||        \
     defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || \
    defined(__gnu_linux__)
#include <stdio.h>

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif

namespace memutils {
/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
size_t getPeakRSS();

/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
size_t getCurrentRSS();
}  // namespace memutils

#endif
