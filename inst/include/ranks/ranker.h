// ranker.h
// source url: https://sites.google.com/site/jivsoft/Home/compute-ranks-of-elements-in-a-c---array-or-vector
// note: KR 2018, added the BSD text and a copyright notice
//
// Licensing: Anyone can use any of my code under the terms of the BSD license.
// Contact: I can be reached using my gmail account.  The account name is kjwilder
// Copyright (c) <2006>, <kjwilder>
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
// 3. All advertising materials mentioning features or use of this software
// must display the following acknowledgement:
// This product includes software developed by the <organization>.
// 4. Neither the name of the <organization> nor the
// names of its contributors may be used to endorse or promote products
// derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY <COPYRIGHT HOLDER> ''AS IS'' AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef RANKER_H
#define RANKER_H

#include <vector>
#include <string>
#include <algorithm>

using std::vector;
using std::string;

#ifndef uint
typedef unsigned int uint;
#endif

template <class T>
class lt { public: static int compare(T a, T b) { return(a < b); } };
template <class T>
class gt { public: static int compare(T a, T b) { return(a > b); } };

template <class T, class C>
class ranker
{
 private:
  const T* p;
  uint sz;

 public:
  ranker(const vector<T>& v) : p(&v[0]), sz(v.size()) { }
  ranker(const T* tp, uint s) : p(tp), sz(s) { }

  int operator()(uint i1, uint i2) const { return(C::compare(p[i1],p[i2])); }

  template <class S>
  void get_orders(vector<S>& w) const {
    w.resize(sz);
    w.front() = 0;
    for (typename vector<S>::iterator i = w.begin(); i != w.end() - 1; ++i)
      *(i + 1) = *i + 1;
    std::sort(w.begin(), w.end(), *this);
  }

  template <class S>
  void get_partial_orders(vector<S>& w, uint num) const {
    if (num > sz) num = sz;
    w.resize(sz);
    w.front() = 0;
    for (typename vector<S>::iterator i = w.begin(); i != w.end() - 1; ++i)
      *(i + 1) = *i + 1;
    std::partial_sort(w.begin(), w.begin() + num, w.end(), *this);
    w.resize(num);
  }

  template <class S>
  void get_ranks(vector<S>& w, const string& method) const {
    w.resize(sz);
    vector<uint> tmp(w.size());
    get_orders(tmp);
    if (method == "average") {
      for (uint c = 0, reps; c < w.size(); c += reps) { reps = 1;
        while (c + reps < w.size() && p[tmp[c]] == p[tmp[c + reps]]) ++reps;
	for (uint k = 0; k < reps; ++k)
          w[tmp[c + k]] = S(2 * c + reps - 1) / 2 + 1;
      }
    } else if (method == "min") {
      for (uint c = 0, reps; c < w.size(); c += reps) { reps = 1;
        while (c + reps < w.size() && p[tmp[c]] == p[tmp[c + reps]]) ++reps;
	for (uint k = 0; k < reps; ++k) w[tmp[c + k]] = c + 1;
      }
    } else if (method == "max") {
      for (uint c = 0, reps; c < w.size(); c += reps) { reps = 1;
        while (c + reps < w.size() && p[tmp[c]] == p[tmp[c + reps]]) ++reps;
	for (uint k = 0; k < reps; ++k) w[tmp[c + k]] = c + reps;
      }
    } else // default
      for (uint c = 0; c < w.size(); ++c) w[tmp[c]] = c + 1;
  }

  template <class S>
  void get_partial_ranks(vector<S>& w, const string& method, size_t num) const {
    if (num > sz) num = sz;
    vector<uint> tmp(sz);
    get_partial_orders(tmp, num);
    w.resize(sz);
    fill(w.begin(), w.end(), 0);
    if (method == "average") {
      for (uint c = 0, reps; c < num; c += reps) { reps = 1;
        while (c + reps < num && p[tmp[c]] == p[tmp[c + reps]]) ++reps;
	for (uint k = 0; k < reps; ++k)
          w[tmp[c + k]] = S(2 * c + reps - 1) / 2 + 1;
      }
    } else if (method == "min") {
      for (uint c = 0, reps; c < num; c += reps) { reps = 1;
        while (c + reps < num && p[tmp[c]] == p[tmp[c + reps]]) ++reps;
	for (uint k = 0; k < reps; ++k) w[tmp[c + k]] = c + 1;
      }
    } else if (method == "max") {
      for (uint c = 0, reps; c < num; c += reps) { reps = 1;
        while (c + reps < num && p[tmp[c]] == p[tmp[c + reps]]) ++reps;
	for (uint k = 0; k < reps; ++k) w[tmp[c + k]] = c + reps;
      }
    } else // default
      for (uint c = 0; c < num; ++c) w[tmp[c]] = c + 1;
  }

};

template <class T, class S>
inline void rank(const vector<T>& v, vector<S>& w,
          const string& method = "average")
  { ranker<T, lt<T> > r(v); r.get_ranks(w, method); }

template <class T, class S>
inline void rank(const T* d, uint size, vector<S>& w,
          const string& method = "average")
  { ranker<T, lt<T> > r(d, size); r.get_ranks(w, method); }

template <class T, class S>
inline void partial_rank(const vector<T>& v, vector<S>& w, uint num,
          const string& method = "average")
  { ranker<T, lt<T> > r(v); r.get_partial_ranks(w, method, num); }

template <class T, class S>
inline void partial_rank(const T* d, uint size, vector<S>& w, uint num,
          const string& method = "average")
  { ranker<T, lt<T> > r(d, size); r.get_partial_ranks(w, method, num); }

template <class T, class S>
inline void order(const vector<T>& v, vector<S>& w)
  { ranker<T, lt<T> > r(v); r.get_orders(w); }

template <class T, class S>
inline void order(const T* d, uint size, vector<S>& w)
  { ranker<T, lt<T> > r(d, size); r.get_orders(w); }

template <class T, class S>
inline void partial_order(const vector<T>& v, vector<S>& w, uint num)
  { ranker<T, lt<T> > r(v); r.get_partial_orders(w, num); }

template <class T, class S>
inline void partial_order(const T* d, uint size, vector<S>& w, uint num)
  { ranker<T, lt<T> > r(d, size); r.get_partial_orders(w, num); }

template <class T, class S>
inline void rankhigh(const vector<T>& v, vector<S>& w,
          const string& method = "average")
  { ranker<T, gt<T> > r(v); r.get_ranks(w, method); }

template <class T, class S>
inline void rankhigh(const T* d, uint size, vector<S>& w,
          const string& method = "average")
  { ranker<T, gt<T> > r(d, size); r.get_ranks(w, method); }

template <class T, class S>
inline void partial_rankhigh(const vector<T>& v, vector<S>& w, uint num,
          const string& method = "average")
  { ranker<T, gt<T> > r(v); r.get_partial_ranks(w, method, num); }

template <class T, class S>
inline void partial_rankhigh(const T* d, uint size, vector<S>& w, uint num,
          const string& method = "average")
  { ranker<T, gt<T> > r(d, size); r.get_partial_ranks(w, method, num); }

template <class T, class S>
inline void orderhigh(const vector<T>& v, vector<S>& w)
  { ranker<T, gt<T> > r(v); r.get_orders(w); }

template <class T, class S>
inline void orderhigh(const T* d, uint size, vector<S>& w)
  { ranker<T, gt<T> > r(d, size); r.get_orders(w); }

template <class T, class S>
inline void partial_orderhigh(const vector<T>& v, vector<S>& w, uint num)
  { ranker<T, gt<T> > r(v); r.get_partial_orders(w, num); }

template <class T, class S>
inline void partial_orderhigh(const T* d, uint size, vector<S>& w, uint num)
  { ranker<T, gt<T> > r(d, size); r.get_partial_orders(w, num); }

template <class T>
inline T quantile(const T* d, const uint size, const double q)
{
  if (size == 0) return T(0);
  if (size == 1) return d[0];
  if (q <= 0) return *std::min_element(d, d + size);
  if (q >= 1) return *std::max_element(d, d + size);

  double pos = (size - 1) * q;
  uint ind = uint(pos);
  double delta = pos - ind;
  vector<T> w(size); std::copy(d, d + size, w.begin());
  std::nth_element(w.begin(), w.begin() + ind, w.end());
  T i1 = *(w.begin() + ind);
  T i2 = *std::min_element(w.begin() + ind + 1, w.end());
  return i1 * (1.0 - delta) + i2 * delta;
}

template <class T>
inline T quantile(const vector<T>& v, const double q)
  { return quantile(&v[0], v.size(), q); }

#endif

