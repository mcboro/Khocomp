// Copyright 2025 Maciej Borodzik
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef PDCODE_H
#define PDCODE_H
/*
 * The plan of the computations is the following
 * 1. Read PD_CODE from file to the class pdcode;
 *
 * 2. Let n be the number of crossings. We compute
 *    all the 2^n resolutions. Each resolution is labelled
 *    by an integer in the range [0,2^n-1].
 *    The resolutions are stored in the vector "allresolutions"
 *    in the pdcode class
 *    Each of the resolutions contains informations about
 *    edges in each circle.
 *
 * 3. From within the class pdcode, we compute all the arrows.
 *    Note that in order to compte the arrows, we need to
 *    have access to the pdcode data.
 *
 * 4. Knowing the arrows, we store the data of fixed
 *    quantum grading and compute the differential.
 */

#include <cstring>
#include <map>
#include <set>
#include <string>
#include <vector>

using circle = std::set<unsigned long int>;
using local_arc = std::pair<unsigned long int, unsigned long int>;
using resolution_arcs = std::vector<local_arc>;

struct arrow
{
  // if ismerge, a,b are starting circles and c is the end circle
  // if not, a is the starting circle and b,c are the final circles.

  long int startpoint;
  long int endpoint;

  long int a;
  long int b;
  long int c;

  // A mapping for circles that don't change.
  std::vector<long int> maps;  // which circle corresponds to another

  bool ismerge;  // merge or split
};

// A collection of all circles.
class resolution
{
  // we represent a circle in a resolution via the set of edges
  // it consists of.
  // the circles are stored in a vector
  // given a single crossing data, we create
 public:
  std::vector<circle> w;
  resolution(std::vector<circle>& u) { w = u; };
  void show_all();
  unsigned long int size() { return w.size(); }
};


// A single crossing.
//
// The `singlecode` class represents a single crossing, i.e., a single element
// of the PD code (represented by the `pdcode` class).
struct singlecode
{
  unsigned long int a = 0;
  unsigned long int b = 0;
  unsigned long int c = 0;
  unsigned long int d = 0;
  singlecode() = default;
  singlecode(unsigned long int x, unsigned long int y, unsigned long int z,
             unsigned long int u)
  {
    a = x;
    b = y;
    c = z;
    d = u;
  }
  std::vector<unsigned long int> iter();
  resolution_arcs resolve(bool which_size);

  bool ispositive();

  bool compare(singlecode stwo);
};

class state
{
  /*
   A word of explanation. The state class does not contain data which is enough
   to compute all the arrows The reason is that the quantity ass assigns 1 or 0
   to circles, but this requires a consistent way of enumerating circles.

   This is hard to obtain when creating resolutions, and maybe even impossible
   sometimes So we artificially inject arrows.
   */
 private:
  // whether the state transitions started being calculated
  bool arrows_injected;

  // bool gone_with_arrows;
  std::vector<arrow>* state_arrows;

 public:
  unsigned long int resol;  // which resolution

  // a bitmap for the vector of circles
  unsigned long int ass;  // 1 assigns x to the circle, 0 assigns 1.

  unsigned long int numbercircles;

  void insert_arrows(std::vector<arrow>*);
  explicit state()
  {
    resol = 0;
    ass = 0;
    arrows_injected = false;
  }
  explicit state(unsigned long int a, unsigned long int b)
  {
    resol = a;
    ass = b;
    arrows_injected = false;
  }
  // we separate arrows in the Khovanov complex
  // from the arrows in the Bar-Natan complex
  // this is because we don't want to introduce polynomials here
  // and we have background to study Lee homology, not just BN.

  // this is the differential (1,1)->1, (1,0),(0,1)->0,
  // 0->(0,0), 1->(0,1)+(1,0)
  std::vector<state> go_with_arrow(arrow ar);

  // the split arrow 1->(1,1)
  std::vector<state> go_with_SplitBarNatan_arrow(arrow ar);

  // the merge arrow (0,0) -> 0
  std::vector<state> go_with_MergeBarNatan_arrow(arrow ar);

  // this is the new_differential.
  std::vector<state> go_with_LastBarNatan_arrow(arrow ar);

  // this gives all states with injected arrows.
  std::vector<state> go_with_injected_arrows();

  // this gives all states with injected arrows.
  std::vector<state> go_with_SplitBarNatan_injected_arrows();

  // this gives all states with injected arrows.
  std::vector<state> go_with_MergeBarNatan_injected_arrows();

  // this gives all states with injected arrows.
  std::vector<state> go_with_LastBarNatan_injected_arrows();

  // combines Merge+Split arrows into one set
  std::vector<state> go_with_BarNatan_injected_arrows();

  // the number or 1's in the binary representation of the number `resol`
  // returns homological degree
  unsigned long int homdeg() const;

  bool operator==(const state& rhs) const
  {
    return ((resol == rhs.resol) and (ass == rhs.ass));
  };

  // Returns quantum degree. (the number of q's or x's)
  long int quant();
};

using arrows_from_resolution = std::vector<arrow>;
using sorted_states = std::map<long int, std::vector<state> >;

// A representation of the PD code.
//
// The `pdcode` class contains allcodes: all the codes
// the constructor sets an empty vector
// the data is usually read from "myfile"
class pdcode
{
 private:
  std::vector<singlecode> allcodes;
  std::vector<resolution> allresolutions;
  std::vector<state> allstates;
  std::vector<arrows_from_resolution> allarrows;

  // Adds one resolution to allresolutions.
  resolution compute_single_resolution(unsigned long long int cube_element);

  arrow compute_single_arrow(unsigned long long int from,
                             unsigned long int byte);

  bool resolution_computed;
  bool arrows_computed;
  bool states_computed;

  // the number of q's (potentially the number of
  // q's - the number of 1's)
  std::set<long int> gradings;

 public:
  // contains signs of crossings.
  std::pair<unsigned long int, unsigned long int> signs();
  pdcode()
  {
    resolution_computed = false;
    arrows_computed = false;
    states_computed = false;
  };
  void insert(singlecode sig) { allcodes.push_back(sig); }
  void read_from_file(std::string myfile);
  unsigned long int length() { return allcodes.size(); }
  void compute_resolutions();
  void compute_arrows();
  void compute_all_states();
  void debug_one();
  std::vector<state> give_all_states();
  sorted_states compute_states();
  void trace_state(state sq);
  std::vector<resolution>* show_res() { return &allresolutions; };
  std::vector<singlecode>* show_code() { return &allcodes; };
};


class statecomplex
{
  // the class describes all states in a given quantum grading
 protected:
  long int quantdegree;

  // the min homological degree for all quantum grading states (for the given
  // quantum number)
  unsigned long int minimalhomdegree;

 public:
  std::vector<std::vector<state> > statesincomplex;
  explicit statecomplex()
  {
    quantdegree = 0;
    minimalhomdegree = 0;
  }
  explicit statecomplex(std::vector<state> all_with_given_quantum);
  std::vector<state> at(long int i) const
  {
    return statesincomplex[i];
  }  // needs a shift by mindegree
  long int dimension(long int i) const { return statesincomplex[i].size(); };
  long int size() const { return statesincomplex.size(); }
  long int mindeg() const { return minimalhomdegree; }
  std::vector<long int> trace_barnatan_map(const statecomplex& nextstatecomplex,
                                           long int statenumber,
                                           long int statehom);
};
long int where_is_my_edge(std::vector<circle> hiding, unsigned long int edge);

#endif /* PDCODE_H */
