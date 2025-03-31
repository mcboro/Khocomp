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

#ifndef KHOVANOV_H
#define KHOVANOV_H

#include <cstring>
#include <vector>

#include "algebra.h"
#include "chaincomplex.h"
#include "pdcode.h"
#include "symmetry.h"

class khovanov
{
 private:
  pdcode PD;
  std::map<long int, KhovanovChainComplex> chains;

  gradingshift mygradingshift;
  std::pair<long int, long int> crossingsigns;
  bool khovanovonly;
  long int khovanovthe;
  std::map<std::pair<long, long>, long int> homology_sizes;
  FMatrix move_bar_natan(FMatrix statevector, long int statehom,
                         long int statequant);
  std::vector<long int> barnatan_differential(long int homologyclass,
                                              long int statehom,
                                              long int statequant);
  std::map<std::pair<long, long>, std::vector<std::vector<long int> > >
      barnatan_maps;
  std::map<std::pair<long, long>, FMatrix> barnatan_chainmaps;
  std::map<std::pair<long, long>, SNF> barnatan_SNFs;
  std::map<std::pair<long, long>, FMatrix> barnatan_Kernels;
  std::map<std::pair<long, long>, FMatrix> barnatan_Images;
  std::map<std::pair<long, long>, FMatrix> barnatan_hom;
  void compute_all_vectors_of_barnatan_differentials();
  void compute_vector_of_barnatan_differentials_in_one_quantum(
      long int onequantum);
  void compute_bar_natan_homologies();
  FMatrix bar_natan_classes(
      std::pair<long,
                long>);  // returns vectors of bar_natan classes in homology
  FMatrix bar_natan_images_at_grading(std::pair<long, long>);
  // this is the action part
  action myaction;
  bool is_action_loaded;
  std::map<std::pair<long, long>, std::vector<std::vector<long int> > >
      action_maps;
  std::map<std::pair<long, long>, FMatrix> one_plus_tau_maps;
  std::map<std::pair<long, long>, SNF> one_plus_tau_SNFs;
  std::map<std::pair<long, long>, FMatrix> one_plus_tau_Kernels;
  std::map<std::pair<long, long>, bool> one_plus_tau_inter_d;
  void handle_action_at_bigrading(
      std::pair<std::pair<long int, long int>, long int>);
  void intersect_kernels_one_plus_tau_and_d();

 public:
  FMatrix show_kernel(std::pair<long int, long int>);
  FMatrix show_image_at_next(std::pair<long int, long int>);
  khovanov(pdcode P, bool khonly = false, long int khthe = 0);
  long int hom() const { return mygradingshift.homological; };
  long int quant() const { return mygradingshift.quantum; };
  void show_homology_at_grading(long int quantum);
  KhovanovChainComplex chain(long int deg) const { return chains.at(deg); };
  std::map<long int, KhovanovChainComplex> allchains() const { return chains; };
  void detailed_homology(long int quantum, long int homological);
  void show_all_homologies();
  void show_all_homologies_but_be_brief();
  long int size_of_hom(std::pair<long, long>);
  FMatrix show_map_at_bigrading(std::pair<long, long>);

  // bar natan part
  long int size_of_bn(std::pair<long, long>);
  std::vector<long int> show_bn_differential(std::pair<long int, long int>,
                                             long int);
  std::vector<long int> show_bn_differential(long int hom, long int quant,
                                             long int which)
  {
    return show_bn_differential(std::make_pair(hom, quant), which);
  };
  // action part
  void load_action_from_file(std::string filename);
  void print_action_at_bigrading(std::pair<long int, long int>);
};

#endif
