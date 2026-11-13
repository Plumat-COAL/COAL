// SPDX-License-Identifier: GPL-3.0-only
// A stand-in for <Rcpp.h>, only for gestimator_ref.cpp: orthologr's copies of the libsequence headers include
// Rcpp.h for Rcpp::Rcout (messages) and nothing else that matters to the numbers. With this file they compile
// without R, unmodified.
#pragma once
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
namespace Rcpp {
struct NullStream : std::ostream {
    NullStream() : std::ostream(nullptr) {}
};
static NullStream Rcout;  // the messages of the headers ("Wrong Nucleotide ...") are not wanted here
}  // namespace Rcpp
