#include "kentr.h"

std::unordered_map<char, char> kentr::bpCompTable = {
  {'C', 'G'},
  {'G', 'C'},
  {'T', 'A'},
  {'A', 'T'},
  {'N', 'N'},
  {'n', 'n'},
  {'c', 'g'},
  {'g', 'c'},
  {'t', 'a'},
  {'a', 't'}
};
