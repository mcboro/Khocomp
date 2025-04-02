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

#include "symmetry.h"

#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

#include "algebra.h"

// we compute the vector translate_codes, which encodes resolutions.
action::action() {}

void action::codes(pdcode& pd)
{
  translate_codes.clear();
  std::vector<singlecode>* all_my_codes = pd.show_code();
  singlecode s2, s0;

  bool found;
  for (unsigned long int i = 0; i < all_my_codes->size(); i++)
  {
    s2 = all_my_codes->at(i);
    s2.a = edge_map[s2.a];
    s2.b = edge_map[s2.b];
    s2.c = edge_map[s2.c];
    s2.d = edge_map[s2.d];
    found = false;
    for (unsigned long int j = 0; j < all_my_codes->size(); j++)
    {
      s0 = all_my_codes->at(j);
      if (s2.compare(s0))
      {
        found = true;
        translate_codes.push_back(j);
      }
    }
    if (not found)
    {
      std::cout << "Error mapping crossings!" << std::endl;
      exit(-1);
    }
  }
}
void action::show_codes()
{
  for (auto& i : translate_codes) std::cout << i << " ";
  std::cout << std::endl;
}
void action::show_resolutions()
{
  for (auto& i : translate_resolutions) std::cout << i << " ";
  std::cout << std::endl;
}

void action::resolutions(pdcode& pd)
{
  translate_resolutions.clear();
  if (translate_codes.size() == 0) codes(pd);
  long int reslen = pd.show_res()->size();
  long int codelen = pd.show_code()->size();
  // std::cout << "codelen "<< codelen << " " << reslen << std::endl;

  unsigned long int b;
  for (long int a = 0; a < reslen; a++)
  {
    b = 0;
    for (long int j = 0; j < codelen; j++)
    {
      // this line assumes that the symmetry preserves the resolution. Note that
      // we *do not* verify this if you want, we need to change the "compare"
      // function in the pdcodes.
      if ((a & (1 << j)) > 0)  // xor (translate_codes[j]==j))
        b = b + (1 << translate_codes[j]);
    }
    translate_resolutions.push_back(b);
  }
}

// transforms a state (old) under the group action.
state action::newstate(state old, pdcode& pd)
{
  if (translate_resolutions.size() == 0) resolutions(pd);
  resolution oldres = pd.show_res()->at(old.resol);
  resolution newres = pd.show_res()->at(translate_resolutions[old.resol]);
  unsigned long int oldass = old.ass;
  unsigned long int newass = 0;
  unsigned long int oldcirc;
  unsigned long int newcirc;
  long int where;
  for (unsigned long int i = 0; i < oldres.size(); i++)
  {
    if ((oldass & (1 << i)) != 0)
    {
      // if the i-th Tcircle in the old resolution is marked with 1
      // we look at which Tcircle is the crossing.
      oldcirc = *oldres.w[i].begin();  // take an edge of the Tcircle
      newcirc =
          edge_map[oldcirc];  // it is mapped by the action to some other edge
      where = where_is_my_edge(newres.w, newcirc);
      if (where < 0)
      {
        std::cout << "Error mapping circles!" << std::endl;
        exit(-1);
      }
      newass = newass + (1 << where);
    }
  }
  state mystate = state(translate_resolutions[old.resol], newass);
  return mystate;
}

std::vector<unsigned long int> action::permute_states(std::vector<state> stvec,
                                                      pdcode& pd)
{
  std::vector<unsigned long int> ps;
  state st;
  for (auto& stold : stvec)
  {
    st = newstate(stold, pd);
    for (long int i = 0; i < stvec.size(); i++)
      if (st == stvec[i]) ps.push_back(i);
  }
  assert(stvec.size() == ps.size());
  return ps;
}

FMatrix action::permute_vector(const KhovanovChainComplex& ch,
                               const FMatrix& FMold, pdcode& pd,
                               long int degree)
{
  FMatrix FMnew(FMold.rows(), FMold.cols());
  std::vector<unsigned long int> translate =
      permute_states(ch.statesincomplex[degree], pd);
  for (long int i = 0; i < FMold.rows(); i++)
    for (long int j = 0; j < FMold.cols(); j++)
      if (FMold.get(i, j)) FMnew.set(translate[i], j, field(true));
  // std::cout << FMold << std::endl;
  // std::cout << FMnew << std::endl;

  return FMnew;
}

std::vector<std::vector<long int> > action::permute_homology(
    const KhovanovChainComplex& ch, pdcode& pd, long int degree)
{
  // works for homology of rank 1 or rank 2.
  std::vector<std::vector<long int> > results;
  FMatrix H = ch.homatdeg(degree);
  if (H.cols() == 0) return results;  // no columns, no results. Sorry.
  /*
  FMatrix Im;
  if (ch.size_of_maps()>0)
    Im=ch.imageatdeg(degree);
  else
    Im=FMatrix(H.rows(),0);
    */
  FMatrix Hnew = permute_vector(ch, H, pd, degree);
  FMatrix Hi;
  // unsigned long int myres;
  for (long int i = 0; i < H.cols(); i++)
  {
    Hi = Hnew.take_column(i);
    results.push_back(ch.trace_class_in_homology(Hi, degree));
  }
  // debug 1
  return results;
}

void action::read_from_file(std::string filename)
{
  std::ifstream myfile(filename);
  std::string handle_string;
  std::vector<std::string> my_split_string;
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
          edge_map.push_back(std::stol(aux_string));
          aux_string.clear();
        }
      }
    }
  }
}
