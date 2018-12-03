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
#include "htslib/htslib/hts.h"
#include "htslib/htslib/sam.h"
#include "htslib/htslib/bgzf.h"
#include "htslib/htslib/faidx.h"
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

CharacterVector revComp(CharacterVector vec) ;
std::string revComp(std::string seq) ;

struct MannWhitney {
  double pval ;
  int w ;
};

#endif
