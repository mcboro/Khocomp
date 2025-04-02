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

#include "pdcode.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

std::vector<unsigned long int> singlecode::iter()
{
  // this is syntactic sugar, turns the pdcode structure
  // into an iterable
  std::vector<unsigned long int> k = {a, b, c, d};
  return k;
}
bool singlecode::ispositive()
{
  // for the pdcode the crossing bd is the overcrossing.
  // if d is right after b, then the crossing is positive, if b is after d, the
  // crossing is negative
  if (b == d + 1) return true;
  if (d == b + 1) return false;
  // the remaining case is when we circle over, like (2,6,3,1). Then we treat
  // the smaller number as the next one hence (2,6,3,1) is a negative crossing,
  // because d=1 is actually regarded as 7.
  if (b < d) return true;
  return false;
}
bool singlecode::compare(singlecode st2)
{
  std::set<unsigned long int> setone{a, b, c, d};
  std::set<unsigned long int> settwo{st2.a, st2.b, st2.c, st2.d};
  if (setone != settwo)
    return false;
  else
    return true;
}
void ::resolution::show_all()
{
  int i = 0;
  for (auto& mycirc : w)
  {
    std::cout << "Circle number" << i << ":";
    for (auto& r : mycirc) std::cout << r << " ";
    std::cout << std::endl;
  }
}

void pdcode::read_from_file(std::string filename)
{
  // reads the file from the string.
  // it cuts out everything that is not a number,
  // so you can write PD[[3,1,4,5], [2,4,...
  // or 3,1,4,5,2,4
  // or even 3;1;4;5...
  std::ifstream myfile(filename);
  std::string handle_string;
  std::vector<std::string> my_split_string;
  std::vector<unsigned long int>
      number_stack;  // we will stack all numbers in this file.
  std::string aux_string;

  while (std::getline(myfile, handle_string))
  {
    for (const char& c : handle_string)
    {
      if (std::isdigit(c))
      {
        aux_string.push_back(c);
      }
      else
      {
        if (aux_string.size() > 0)
        {
          number_stack.push_back(std::stol(aux_string));
          aux_string.clear();
        }
      }
    }
  }
  for (int i = 0; i < number_stack.size(); i += 4)
  {
    allcodes.push_back(singlecode(number_stack[i], number_stack[i + 1],
                                  number_stack[i + 2], number_stack[i + 3]));
  }
}

std::pair<unsigned long int, unsigned long int> pdcode::signs()
{
  unsigned long int p(0), n(0);
  for (auto& cd : allcodes)
  {
    if (cd.ispositive())
      p++;
    else
      n++;
  }
  return std::pair<unsigned long int, unsigned long int>(p, n);
}

resolution_arcs singlecode::resolve(bool which_size)
{
  // gives two arcs in the resolution
  // the convention is that the crossing with SW-NE on the top
  // is resolved to cup-cap if 0 (which_size=false)
  // and to )( if 1 (which_size=true)
  // note that the a-entry in single code is the ingoing tunnel
  //
  local_arc la1, la2;
  la1.first = a;
  la2.first = c;
  resolution_arcs rs;
  if (which_size)
  {
    la1.second = d;
    la2.second = b;
  }
  else
  {
    la1.second = b;
    la2.second = d;
  }
  rs.push_back(la1);
  rs.push_back(la2);
  return rs;
}

// auxiliary function for "compute_single_resolution" procedure.
long int where_is_my_edge(std::vector<circle> hiding, unsigned long int edge)
{
  // returns the index of hiding for which edge is a subset
  // or -1, if there's no such index
  for (long int i = 0; i < hiding.size(); i++)
  {
    if (hiding[i].find(edge) != hiding[i].end()) return i;
  }
  return -1;
}

// warning! The procedure doesn't work if there is a loop in the diagram. That
// is, singlecode is forbidden to have two equal values.
resolution pdcode::compute_single_resolution(
    unsigned long long int cube_element)
{
  std::vector<circle> finaloutput;
  unsigned long long int powerindex = 1;
  bool crossresol = false;
  resolution_arcs myresol;
  local_arc p, q;
  singlecode s;
  long int a, b, c;
  std::vector<unsigned long int> t;
  circle aux_circle;
  for (long int index = 0; index < length(); index++)
  {
    s = allcodes[index];
    t = s.iter();
    // we check if we take the 0 or 1 resolution at step index
    if ((cube_element & powerindex) > 0)
      crossresol = true;
    else
      crossresol = false;
    myresol = s.resolve(crossresol);
    for (auto& arc : myresol)
    {
      // std::cout << arc.first << " " << arc.second << std::endl;
      a = where_is_my_edge(finaloutput, arc.first);
      b = where_is_my_edge(finaloutput, arc.second);
      // std::cout << a << " " << b << std::endl << std::endl;
      if ((a == -1) and (b == -1))
      {
        // no arcs are there, we need to start another circle.
        aux_circle.clear();
        aux_circle.insert(arc.first);
        aux_circle.insert(arc.second);
        finaloutput.push_back(aux_circle);
      }
      if ((a == -1) and (b >= 0))
      {
        // we add first edge to the b-th circle
        aux_circle = finaloutput[b];
        aux_circle.insert(arc.first);
        finaloutput[b] = aux_circle;
      }
      if ((a >= 0) and (b == -1))
      {
        aux_circle = finaloutput[a];
        aux_circle.insert(arc.second);
        finaloutput[a] = aux_circle;
      }
      if ((a >= 0) and (b >= 0))
      {
        // the two circles are there.
        if (a == b)
          continue;  // if the edges are there, we don't do anything. Otherwise
                     // we need to merge the two circles.
        if (b < a)
        {
          c = a;
          a = b;
          b = c;
        }  // ensure that b is larger than a
        aux_circle = finaloutput[a];
        for (auto& d : finaloutput[b])
          aux_circle.insert(d);  // merge b-th circle to a-th
        finaloutput[a] = aux_circle;
        finaloutput.erase(finaloutput.begin() +
                          b);  // hopefully, this removes the b-th circle
      }
    }
    powerindex = powerindex * 2;
  }
  return resolution(finaloutput);
}

void pdcode::compute_resolutions()
{
  // we compute the full cube of resolutions.
  for (unsigned long long int i = 0; i < std::pow(2, allcodes.size()); i++)
  {
    resolution res = compute_single_resolution(i);
    allresolutions.push_back(res);
  }
  resolution_computed = true;
}

arrow pdcode::compute_single_arrow(unsigned long long int from,
                                   unsigned long int byte)
{
  // from - index of the resolution (an element of the cube of resolutions)
  // byte indicates the crossing that is going to be modified
  long long int pw = std::pow(2, byte);
  arrow ar;
  ar.startpoint = from;     // resolution [from]
  ar.endpoint = from + pw;  // resolution [to]

  // std::cout << "Arrow from : "<<ar.startpoint << " To: " <<ar.endpoint;
  resolution r1 = allresolutions[from];
  resolution r2 = allresolutions[from + pw];
  // r1.show_all();
  // r2.show_all();
  singlecode sc = allcodes[byte];  // the crossing that is being modified
  // std::cout << sc.a <<" " << sc.c;
  //  we check, which circles meet at the crossing.
  unsigned long int k1 = where_is_my_edge(
      r1.w, sc.a);  // the index of the circle that 'a' belongs to
  unsigned long int k2 = where_is_my_edge(r1.w, sc.c);
  unsigned long int l1 = where_is_my_edge(r2.w, sc.a);
  unsigned long int l2 = where_is_my_edge(r2.w, sc.c);
  if ((k1 == k2) and (l1 == l2))
    std::cout << "Error! We merge circle to itself!" << std::endl;
  if ((k1 != k2) and (l1 != l2))
    std::cout << "Error! We split two circles into two" << std::endl;
  ar.a = k1;
  ar.c = l2;
  // split
  if (k1 == k2)
  {
    // there is one circle involved in the crossing
    // so we have a split.
    ar.b = l1;
    ar.ismerge = false;
    // merge
  }
  else
  {
    ar.b = k2;
    ar.ismerge = true;
  }
  // now we compute the map
  ar.maps.clear();
  for (auto& circ : r1.w)
  {
    ar.maps.push_back(where_is_my_edge(r2.w, *circ.begin()));
  }
  return ar;
}

void pdcode::compute_arrows()
{
  unsigned long long int pw = std::pow(2, allcodes.size());
  arrows_from_resolution ra;
  if (!resolution_computed) compute_resolutions();
  arrows_computed = true;
  unsigned long long int mybyte;
  for (unsigned long long int from = 0; from < pw; from++)
  {
    ra.clear();
    mybyte = 1;
    for (unsigned long int byte = 0; byte < allcodes.size(); byte++)
    {
      if ((from & mybyte) == 0) ra.push_back(compute_single_arrow(from, byte));
      mybyte = mybyte * 2;
    }
    allarrows.push_back(ra);
  }
}
void pdcode::debug_one()
{
  // shows the first circle
  resolution res = allresolutions[0];
  int i = 0;
  for (auto& circ : res.w)
  {
    i++;
    std::cout << "Circle number " << i << std::endl;
    for (auto& cross : circ) std::cout << cross << " ";
    std::cout << std::endl;
  }
  std::cout << std::endl << std::endl;
  std::cout << allarrows.size() << std::endl;
  arrows_from_resolution arr = allarrows[0];
  i = 0;
  for (auto& singlearrow : arr)
  {
    i++;
    std::cout << "Arrow number " << i << std::endl;
    std::cout << "From " << singlearrow.startpoint << " To "
              << singlearrow.endpoint << "  a=" << singlearrow.a
              << " b=" << singlearrow.b << " c=" << singlearrow.c;
    if (singlearrow.ismerge)
      std::cout << " Merges" << std::endl;
    else
      std::cout << " Splits" << std::endl;
    std::cout << std::endl;
  }
}

void pdcode::trace_state(state sq)
{
  // outputs info about the state.
  resolution res = allresolutions[sq.resol];
  std::cout << "Resolution number number " << sq.resol << " state number "
            << sq.ass << std::endl
            << "  Quantum grading " << sq.quant() << "Homological grading "
            << sq.homdeg() << " There are " << sq.numbercircles << "circles "
            << std::endl;
  unsigned long int i = 0;
  for (auto& mycirc : res.w)
  {
    std::cout << "Circle ";
    for (auto& r : mycirc) std::cout << r << " ";
    std::cout << "has label ";
    if ((sq.ass & (1 << i)) > 0)
      std::cout << "1";
    else
      std::cout << "0";
    std::cout << std::endl;
    i++;
  }
  std::cout << std::endl;
  std::cout << "There are arrows to states:" << std::endl;
  std::vector<state> vstate;
  for (auto& myarrow : allarrows[sq.resol])
  {
    vstate = sq.go_with_arrow(myarrow);
    for (auto& sq2 : vstate)
      std::cout << "(" << sq2.resol << "," << sq2.ass << ")" << std::endl;
  }
}

unsigned long int degree(unsigned long int i)
{
  // returns the degree of a number. If i is written in the binary form, we add
  // up all the ones.
  unsigned long int deg = 0;
  while (i > 0)
  {
    if ((i % 2) == 1) deg++;
    i = (i >> 1);
  }
  return deg;
}

std::vector<state> state::go_with_arrow(arrow ar)
{
  // we compute all the arrows in the Khovanov differential
  // we need to define the output as a vector, because we don't know a priori
  // how many arrows will be.
  // convention:
  // 1=v_+, 0=v_-
  // merge: (1,1) -> 1
  //        (1,0), (0,1) ->0
  //        (0,0) does not go anywhere
  // split: 0 -> (0,0)
  //        1 -> (0,1) + (1,0).
  std::vector<state> mystates;
  // long int myend = ar.endpoint;
  long int counter = 0;
  if (ar.ismerge)
  {
    // if the arrow is a merge, there is one vector in mystates.
    state mergestate;
    int gradingcount = 0;
    mergestate.resol = ar.endpoint;
    mergestate.ass = 0;
    mergestate.numbercircles =
        numbercircles - 1;  // merging decreases the number
    for (auto& j : ar.maps)
    {
      // copy states from other circles
      if ((counter != ar.a) and (counter != ar.b))
      {
        if ((ass & (1 << counter)) > 0) mergestate.ass += (1 << j);
      }
      counter++;
    }
    if ((ass & (1 << ar.a)) > 0) gradingcount++;
    if ((ass & (1 << ar.b)) > 0) gradingcount++;
    if (gradingcount == 0)
      return mystates;  // we have two times v_-, we don't add anything
    if ((mergestate.ass & (1 << ar.c)) > 0)
      std::cout << "Problem with merging" << std::endl;
    else if (gradingcount == 2)
      mergestate.ass += (1 << ar.c);
    mystates.push_back(mergestate);
    return mystates;
  }
  // now we deal with the case i is a split;
  //
  state firstsplit;
  state secondsplit;
  firstsplit.resol = ar.endpoint;
  secondsplit.resol = ar.endpoint;
  firstsplit.numbercircles = numbercircles + 1;
  secondsplit.numbercircles = numbercircles + 1;

  for (auto& j : ar.maps)
  {
    // copy states from other circles
    if (counter != ar.a)
    {
      if ((ass & (1 << counter)) > 0) firstsplit.ass += (1 << j);
    }
    counter++;
  }
  if ((firstsplit.ass & (1 << ar.c)) > 0)
    std::cout << "Problem with merging" << std::endl;
  if ((firstsplit.ass & (1 << ar.b)) > 0)
    std::cout << "Problem with merging" << std::endl;
  secondsplit.ass = firstsplit.ass;
  if ((ass & (1 << ar.a)) > 0)
  {
    firstsplit.ass += (1 << ar.b);
    secondsplit.ass += (1 << ar.c);
    mystates.push_back(firstsplit);
    mystates.push_back(secondsplit);
    return mystates;
  }
  else
  {
    mystates.push_back(firstsplit);  // it has 0,0 at positions b and c.
    return mystates;
  }
}

std::vector<state> state::go_with_LastBarNatan_arrow(arrow ar)
{
  // The arrow is the perturbation of Bar-Natan complex
  // with conventions
  // 1=v_+, 0=v_-
  // we set 0 -> (0,0).
  std::vector<state> mystates;
  mystates.clear();
  // long int myend = ar.endpoint;
  long int counter = 0;
  if (ar.ismerge) return mystates;
  state firstsplit;
  firstsplit.resol = ar.endpoint;
  firstsplit.numbercircles = numbercircles + 1;

  for (auto& j : ar.maps)
  {
    // copy states from other circles
    if (counter != ar.a)
    {
      if ((ass & (1 << counter)) > 0) firstsplit.ass += (1 << j);
    }
    counter++;
  }
  if ((firstsplit.ass & (1 << ar.c)) > 0)
    std::cout << "Problem with merging" << std::endl;
  if ((firstsplit.ass & (1 << ar.b)) > 0)
    std::cout << "Problem with merging" << std::endl;
  if ((ass & (1 << ar.a)) == 0) mystates.push_back(firstsplit);
  return mystates;
}

std::vector<state> state::go_with_SplitBarNatan_arrow(arrow ar)
{
  // The arrow is the perturbation of Bar-Natan complex
  // with conventions
  // 1=v_+, 0=v_-
  // as above, we set (0,0) -> 0
  // and 1->(1,1)?
  std::vector<state> mystates;
  mystates.clear();
  // long int myend = ar.endpoint;
  long int counter = 0;
  if (ar.ismerge) return mystates;
  state firstsplit;
  firstsplit.resol = ar.endpoint;
  firstsplit.numbercircles = numbercircles + 1;

  for (auto& j : ar.maps)
  {
    // copy states from other circles
    if (counter != ar.a)
    {
      if ((ass & (1 << counter)) > 0) firstsplit.ass += (1 << j);
    }
    counter++;
  }
  if ((firstsplit.ass & (1 << ar.c)) > 0)
    std::cout << "Problem with merging" << std::endl;
  if ((firstsplit.ass & (1 << ar.b)) > 0)
    std::cout << "Problem with merging" << std::endl;
  if ((ass & (1 << ar.a)) > 0)
  {
    firstsplit.ass += (1 << ar.b);
    firstsplit.ass += (1 << ar.c);
    mystates.push_back(firstsplit);
  }
  return mystates;
}
std::vector<state> state::go_with_MergeBarNatan_arrow(arrow ar)
{
  // if the arrow is a merge, there is one vector in mystates.
  std::vector<state> mystates;
  mystates.clear();
  // long int myend = ar.endpoint;
  long int counter = 0;
  if (!ar.ismerge) return mystates;
  state mergestate;
  int gradingcount = 0;
  mergestate.resol = ar.endpoint;
  mergestate.ass = 0;
  mergestate.numbercircles = numbercircles - 1;  // merging decreases the number
  for (auto& j : ar.maps)
  {
    // copy states from other circles
    if ((counter != ar.a) and (counter != ar.b))
    {
      if ((ass & (1 << counter)) > 0) mergestate.ass += (1 << j);
    }
    counter++;
  }
  if ((ass & (1 << ar.a)) > 0) gradingcount++;
  if ((ass & (1 << ar.b)) > 0) gradingcount++;
  if (gradingcount != 0) return mystates;
  mystates.push_back(mergestate);
  return mystates;
}

unsigned long int state::homdeg() const
{
  // naive homological degree
  return degree(resol);
}

long int state::quant()
{
  // naive quantum degree.
  // the quantinty 2*degree-numbercircles
  // is precisely what we get by assigning +1 to v_+ (state true) and -1 to v_-.
  return 2 * degree(ass) - numbercircles + degree(resol);
}
void pdcode::compute_all_states()
{
  if (!resolution_computed) compute_resolutions();
  state mystate;
  unsigned long int maxpow;
  for (unsigned long int n = 0; n < allresolutions.size(); n++)
  {
    maxpow = (1 << allresolutions[n].size());
    for (unsigned long int m = 0; m < maxpow; m++)
    {
      mystate = state(n, m);
      mystate.numbercircles =
          allresolutions[n].size();  // this is the size of circles in the state
      gradings.insert(mystate.quant());
      mystate.insert_arrows(&allarrows[n]);
      allstates.push_back(mystate);
    }
  }
  states_computed = true;
}
sorted_states pdcode::compute_states()
{
  if (!states_computed) compute_all_states();
  // returns a map that assigns to an integer all possible states with given
  // quantum grading
  sorted_states sort;
  std::vector<state> aux_states;
  // now allstates contains all states and gradings contains all the quantum
  // gradings. we start splitting allstates over all gradings.
  for (auto& grading : gradings)
  {
    // aux_states << all_states | std::views::filter([&grading](state r){return
    // r.quant()==grading;});
    aux_states.clear();
    for (auto& newstate : allstates)
    {
      if (grading == newstate.quant()) aux_states.push_back(newstate);
    }
    sort[grading] = aux_states;
  }
  return sort;
}
std::vector<state> pdcode::give_all_states()
{
  if (!states_computed) compute_all_states();
  return allstates;
}
void state::insert_arrows(std::vector<arrow>* a)
{
  state_arrows = a;
  arrows_injected = true;
}

std::vector<state> state::go_with_injected_arrows()
{
  std::vector<state> st;
  std::vector<state> uv;
  if (!arrows_injected) return st;
  for (auto& ar : *state_arrows)
  {
    uv = go_with_arrow(ar);
    for (auto& u : uv) st.push_back(u);
  }
  return st;
}

std::vector<state> state::go_with_SplitBarNatan_injected_arrows()
{
  std::vector<state> st;
  std::vector<state> uv;
  if (!arrows_injected) return st;
  for (auto& ar : *state_arrows)
  {
    uv = go_with_SplitBarNatan_arrow(ar);
    for (auto& u : uv) st.push_back(u);
  }
  return st;
}
std::vector<state> state::go_with_LastBarNatan_injected_arrows()
{
  std::vector<state> st;
  std::vector<state> uv;
  if (!arrows_injected) return st;
  for (auto& ar : *state_arrows)
  {
    uv = go_with_LastBarNatan_arrow(ar);
    for (auto& u : uv) st.push_back(u);
  }
  return st;
}

std::vector<state> state::go_with_MergeBarNatan_injected_arrows()
{
  std::vector<state> st;
  std::vector<state> uv;
  if (!arrows_injected) return st;
  for (auto& ar : *state_arrows)
  {
    uv = go_with_MergeBarNatan_arrow(ar);
    for (auto& u : uv) st.push_back(u);
  }
  return st;
}
std::vector<state> state::go_with_BarNatan_injected_arrows()
{
  std::vector<state> uv = go_with_MergeBarNatan_injected_arrows();
  std::vector<state> rs = go_with_SplitBarNatan_injected_arrows();
  uv.insert(uv.end(), rs.begin(), rs.end());
  return uv;
}
statecomplex::statecomplex(std::vector<state> all_with_given_quantum)
{
  if (all_with_given_quantum.size() == 0)
  {
    std::cout << "Error with statecomplex!" << std::endl;
    exit(-1);
  }
  quantdegree = all_with_given_quantum[0].quant();
  std::set<unsigned long int> ss;
  for (auto& st : all_with_given_quantum) ss.insert(st.homdeg());
  minimalhomdegree = *(ss.begin());  // the smallest hom grading
  std::vector<state> states_in_grading;
  states_in_grading.clear();
  for (unsigned long int i = minimalhomdegree; i <= *(ss.rbegin()); i++)
    statesincomplex.push_back(states_in_grading);
  // if minimal degree is, say, 5, we push all the states with given degree
  // to 5.
  for (auto& st : all_with_given_quantum)
    statesincomplex[st.homdeg() - minimalhomdegree].push_back(st);
}
std::vector<long int> statecomplex::trace_barnatan_map(
    const statecomplex& nextstatecomplex, long int statenumber,
    long int statehom)
{
  // takes the state number statenumber in statehom homological degree
  // goes with barnatan arrow and returns index of the element.
  // homological degree is second absolute. Namely it is subtracted
  // minimalhomdegree.
  int degreeshift = 1;
  std::vector<long int> allmaps;
  if (statehom < minimalhomdegree) return allmaps;
  if (statehom >= minimalhomdegree + statesincomplex.size()) return allmaps;
  if (statehom + degreeshift >=
      nextstatecomplex.mindeg() + nextstatecomplex.size())
    return allmaps;  // we have nowhere to go
  if (statenumber < 0) return allmaps;
  if (statenumber >= statesincomplex[statehom - minimalhomdegree].size())
    return allmaps;
  state initial_state =
      statesincomplex[statehom - minimalhomdegree].at(statenumber);
  std::vector<state> all_final_states =
      initial_state.go_with_BarNatan_injected_arrows();
  if (all_final_states.size() == 0) return allmaps;
  // if (all_final_states.size()>1) { std::cout << "Error! Too much states in
  // the differential" << std::endl; return -1;}
  long int indexstate = 0;
  for (auto& final_state : all_final_states)
  {
    indexstate = 0;
    // state final_state=all_final_states[0];
    for (auto& iterstate : nextstatecomplex.at(statehom + degreeshift -
                                               nextstatecomplex.mindeg()))
    {
      if (iterstate == final_state)
        allmaps.push_back(indexstate);  // return indexstate;
      indexstate++;
    }
  }
  // std::cout << "Found " << allmaps.size() << " states" << std::endl;
  return allmaps;
}
