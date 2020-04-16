// kentr.h
//
// Copyright (C) 2017 Kent Riemondy
//
// This file is part of kentr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#ifndef kentr__kentr_H
#define kentr__kentr_H

#include <algorithm>

#include <Rcpp.h>
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "htslib/faidx.h"
#include "htslib/tbx.h"
#include "ssw/ssw_cpp.h"
#include "ranks/ranker.h"

using namespace Rcpp ;

// class for handling bam file opening and closing
class BamReader {
public:
  samFile* in;
  hts_idx_t* idx;
  BGZF* bz;
  bam_hdr_t* header ;
  BamReader(const std::string& bampath,
            bool check_idx = true,
            int cache_size=10*BGZF_MAX_BLOCK_SIZE) ;

  ~BamReader(){
    hts_idx_destroy(idx);
    sam_close(in);
  }
};

class TabixReader {
public:
  htsFile* in;
  tbx_t* idx;
  BGZF* bz;
  TabixReader(const std::string& tbxpath,
              bool check_idx = true) ;

  ~TabixReader(){
    tbx_destroy(idx);
    hts_close(in);
  }
};

namespace kentr {
extern std::unordered_map<char, char> bpCompTable ;
} ;

CharacterVector revComp(CharacterVector vec) ;
std::string revComp(std::string seq) ;

struct MannWhitney {
  double pval ;
  int w ;
};

#endif
