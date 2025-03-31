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

#ifndef SYMMETRY_H
#define SYMMETRY_H

#include <cstring>

#include "algebra.h"
#include "chaincomplex.h"
#include "pdcode.h"

class action
{
 private:
  std::vector<unsigned long int> edge_map;
  std::vector<unsigned long int> translate_codes;
  std::vector<unsigned long int> translate_resolutions;

  std::vector<unsigned long int> permute_states(std::vector<state> stvec,
                                                pdcode& pd);

 public:
  action();
  void resolutions(pdcode& pd);
  void codes(pdcode& pd);
  void show_codes();
  void show_resolutions();
  state newstate(state old, pdcode& pd);
  void read_from_file(std::string filename);
  FMatrix permute_vector(const KhovanovChainComplex& ch, const FMatrix& FMold,
                         pdcode& pd, long int degree);
  std::vector<std::vector<long int> > permute_homology(
      const KhovanovChainComplex& ch, pdcode& pd, long int degree);
  // new feature: permute kernel.
};

#endif
