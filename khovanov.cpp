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

#include "khovanov.h"

#include <cassert>
#include <iostream>
#include <string>
#include <thread>

#include "algebra.h"
#include "chaincomplex.h"

#define ANSI_COLOR_RED "\x1b[31m"
#define ANSI_COLOR_GREEN "\x1b[32m"
#define ANSI_COLOR_YELLOW "\x1b[33m"
#define ANSI_COLOR_BLUE "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN "\x1b[36m"
#define ANSI_COLOR_RESET "\x1b[0m"

khovanov::khovanov(pdcode P, bool khonly, long int khthe)
{
  PD = P;
  // green = false; // don't compute invariance until we know it's safe

  sorted_states Q = PD.compute_states();
  std::pair<unsigned long int, unsigned long int> cross = PD.signs();
  // std::cout << cross.first << "  AA " << cross.second << std::endl;
  mygradingshift = gradingshift(cross.first, cross.second);
  khovanovonly = khonly;
  khovanovthe =
      khthe -
      mygradingshift.homological;  // we replace absolute shift by homological
  std::vector<state> w;
  for (auto& ss : Q)
  {
    std::cout << "Handling quantum grading "
              << ss.first + mygradingshift.quantum << std::endl;
    w = ss.second;
    chains[ss.first] = KhovanovChainComplex(w, khovanovonly, khovanovthe);
  }
  std::vector<std::thread> mythreads;
  for (auto& ss : chains)
  {
    KhovanovChainComplex* chain = &ss.second;
    int quantum_grading = ss.first + mygradingshift.quantum;
    mythreads.emplace_back(
        [chain, quantum_grading]()
        {
          chain->initialize();
          std::cout << "Quantum grading " + std::to_string(quantum_grading) +
                           " initialized\n";
        });
  }
  for (auto& th : mythreads) th.join();
  long int homdeg;
  long int quantdeg;
  for (auto& ss : chains)
  {
    quantdeg = ss.first;
    for (long int i = 0; i < ss.second.homlen(); i++)
    {
      homdeg = i + ss.second.mindeg();
      homology_sizes[std::make_pair(homdeg, quantdeg)] =
          ss.second.homatdeg(i).cols();
      // std::cout << "Homology sizes " << homdeg << "," << quantdeg << " = " <<
      // homology_sizes[std::make_pair(homdeg,quantdeg)] << std::endl;
    }
  }
  std::cout << "Computing Bar Natan differentials " << std::endl;
  if (!khovanovonly)
  {
    compute_all_vectors_of_barnatan_differentials();
    compute_bar_natan_homologies();
  }
  std::cout << "Done computing differentials " << std::endl;
}
long int khovanov::size_of_hom(std::pair<long, long> a)
{
  if (!homology_sizes.contains(a)) return 0;
  return homology_sizes[a];
}
long int khovanov::size_of_bn(std::pair<long, long> a)
{
  // std::cout << "Size of bn upoly with " << a.first << " " << a.second <<
  // std::endl;
  if (!barnatan_hom.contains(a)) return 0;
  return barnatan_hom[a].cols();
}

void khovanov::show_all_homologies()
{
  /* old version
     for (auto & ch:chains)
     {
     if (ch.second.homlen() > 0) {
     ch.second.export_homology(mygradingshift);
     }
     }
     */
  // new version
  std::pair<long, long> bigrading;
  std::vector<long int> kappa;
  long int homdegree;
  // long int genvector;
  // long int outvector;
  for (auto& chmap : chains)
  {
    homdegree = 0;
    for (auto& MyHom : chmap.second.homologies())
    {
      if (MyHom.cols() > 0)
      {
        bigrading =
            std::make_pair(homdegree + chmap.second.mindeg(), chmap.first);
        std::cout << "Homology at deg = "
                  << homdegree + mygradingshift.homological +
                         chmap.second.mindeg()
                  << " quantum " << chmap.first + mygradingshift.quantum
                  << " has rank " << MyHom.cols() << std::endl;
        if (!khovanovonly)
        {
          if (is_action_loaded and one_plus_tau_Kernels.contains(bigrading))
            std::cout << "Kernel of 1+tau has dimension "
                      << one_plus_tau_Kernels[bigrading].cols() << std::endl;
          std::cout << "Size of Bar Natan homology is " << size_of_bn(bigrading)
                    << std::endl;
          /* if (green and size_of_bn(bigrading)>0)
          {
              outvector = wye_complex(bigrading);
              if (outvector>0)
              {
                  std::cout << ANSI_COLOR_GREEN << " there are " << outvector <<
          " BarNatan generators that are not represented by invariant cycles."
          << ANSI_COLOR_RESET << std::endl;
              }
              else std::cout << ANSI_COLOR_MAGENTA << "all Bar-Natan generators
          are represented by invariant cycles " << ANSI_COLOR_RESET <<
          std::endl;
          } */
          for (long int i = 0; i < MyHom.cols(); i++)
          {
            std::cout << "Generator number " << i << " ";
            kappa = show_bn_differential(homdegree + chmap.second.mindeg(),
                                         chmap.first, i);
            if (kappa.size() == 0)
              std::cout << ANSI_COLOR_MAGENTA
                        << "Kernel of the Bar-Natan differential "
                        << ANSI_COLOR_RESET << std::endl;
            else
            {
              if (kappa[0] < 0)
                std::cout << ANSI_COLOR_MAGENTA
                          << "Kernel of the Bar-Natan differential ("
                          << kappa[0] << ")" << ANSI_COLOR_RESET << std::endl;
              else
              {
                std::cout << ANSI_COLOR_YELLOW
                          << "The Bar-Natan differential maps the generator to "
                             "generators ";
                for (auto& genvector : kappa) std::cout << genvector << " ";
                std::cout << ANSI_COLOR_RESET << std::endl;
              }
            }
          }
        }
        print_action_at_bigrading(bigrading);
        if (!khovanovonly and is_action_loaded)
        {
          if (one_plus_tau_inter_d.contains(bigrading))
          {
            std::cout << " The intersection of ker d_1 and 1+tau ";
            if (one_plus_tau_inter_d[bigrading])
              std::cout << ANSI_COLOR_BLUE << "contains";
            else
              std::cout << ANSI_COLOR_RED << "does not contain ";
            std::cout << ANSI_COLOR_RESET
                      << " a non-trivial element in bar-natan homology"
                      << std::endl;
          }
          else
          {
            std::cout << "The intersection of ker d_1 and 1+tau was not "
                         "computed or `obvious' reasons"
                      << std::endl;
          }
        }
        std::cout << std::endl;
      }
      homdegree++;
      // std::cout << homdegree << std::endl;
    }
  }
}

void khovanov::detailed_homology(long int quantum, long int homological)
{
  long int realquant = quantum - mygradingshift.quantum;
  long int realhomo = homological - mygradingshift.homological;
  // std::cout << "Homology input " << homological << " quant: "<< quantum <<
  // std::endl; std::cout << "Normalized " << realhomo << " quant: "<< realquant
  // << std::endl;
  KhovanovChainComplex ch = chains[realquant];
  // std::cout << "At this quantum grading homology is from " << ch.mindeg() <<
  // " to " << ch.mindeg()+ch.size() << std::endl;

  if (realhomo < ch.mindeg())
  {
    std::cout << "Homological degree too small" << std::endl;
    return;
  }
  if (realhomo >= ch.mindeg() + ch.size())
  {
    std::cout << "Homological degree too large" << std::endl;
    return;
  }
  FMatrix hom = ch.homatdeg(realhomo - ch.mindeg());
  if (hom.cols() == 0)
  {
    std::cout << " No homology at quantum " << quantum << " and homological "
              << homological << std::endl;
    return;
  }
  std::vector<state> mystates = ch.statesincomplex.at(realhomo - ch.mindeg());
  state r;
  std::cout << " Homology at quantum " << quantum << " and homological "
            << homological << " has rank " << hom.cols() << std::endl;
  assert(hom.rows() ==
         mystates.size());  // we make sure that we show correct states.
  for (long int i = 0; i < hom.cols(); i++)
  {
    std::cout << std::endl
              << "The generator number " << i
              << " is a sum of the following states" << std::endl;
    for (long int j = 0; j < hom.rows(); j++)
    {
      if (hom.get(j, i))
      {
        r = mystates[j];
        std::cout << "State number " << j << std::endl;
        // std::cout << "State at resolution " << r.resol << " with labelling "
        // << r.ass << std::endl;
        PD.trace_state(r);
      }
    }
  }
}
std::vector<long int> khovanov::show_bn_differential(
    std::pair<long int, long int> a, long int which)
{
  if (!barnatan_maps.contains(a))
  {
    std::vector<long int> w;
    return w;
  }
  if (barnatan_maps[a].size() <= which)
  {
    std::vector<long int> w;
    return w;
  }
  return barnatan_maps[a][which];
}
FMatrix khovanov::move_bar_natan(FMatrix statevector, long int statehom,
                                 long int statequant)
{
  // takes the statevector at quantum and hom degree and creates a new
  // statevector after barnatan map
  // int homdegreeshift = 1;
  int quantdegreeshift = 2;
  FMatrix empty(0, 0);
  if (!chains.contains(statequant))
  {
    // std::cout << "Empty grading " << statequant << std::endl;
    return empty;
  }  // this should never happen
  if (!chains.contains(statequant + quantdegreeshift))
  {
    // std::cout << "Empty second grading " << statequant+quantdegreeshift <<
    // std::endl;
    return empty;
  }
  // std::cout << statehom << " " << statequant << " not empty" << std::endl;
  FMatrix GG = chains[statequant].vector_of_barnatan_map(
      chains[statequant + quantdegreeshift], statevector, statehom);
  return GG;
}

std::vector<long int> khovanov::barnatan_differential(long int homologyclass,
                                                      long int statehom,
                                                      long int statequant)
{
  // std::cout << "barnatan_differential(" << homologyclass << "," << statehom
  // <<"," << statequant <<")" << std::endl; take homologyclass in statehom,
  // statequant degree
  FMatrix homologymatrix =
      chains[statequant].homatdeg(statehom - chains[statequant].mindeg());
  FMatrix source_of_barnatan = homologymatrix.take_column(homologyclass);
  //  std::cout <<"C" << std::endl;
  FMatrix image_of_barnatan =
      move_bar_natan(source_of_barnatan, statehom, statequant);
  //  std::cout <<"D" << std::endl;
  std::vector<long int> differ;
  if (image_of_barnatan.rows() == 0)
    return differ;  // no class for whatever reason
  int quantshift = 2;
  int homshift = 1;
  // we need to make sure that the matrix is not zero.
  bool has_nonzero = false;
  for (long int i = 0; i < image_of_barnatan.rows(); i++)
    if (image_of_barnatan.get(i, 0)) has_nonzero = true;
  if (!has_nonzero)
  {
    differ.push_back(-2);
    return differ;
  }
  differ = chains[statequant + quantshift].trace_class_in_homology(
      image_of_barnatan,
      statehom + homshift - chains[statequant + quantshift].mindeg());
  // std::cout <<"E" << std::endl;
  return differ;
}

void khovanov::compute_all_vectors_of_barnatan_differentials()
{
  std::vector<std::vector<long int> > bn_maps_fix_hom;
  std::vector<std::vector<std::vector<long int> > > bn_maps_all_hom;
  std::vector<std::thread> mythreads;
  long int qgrad;
  for (auto& quant : chains)
  {
    qgrad = quant.first;
    // std::cout << "Passing " << qgrad << std::endl;
    mythreads.push_back(std::thread(
        &khovanov::compute_vector_of_barnatan_differentials_in_one_quantum,
        this, qgrad));
  }
  for (auto& th : mythreads) th.join();
}
void khovanov::compute_vector_of_barnatan_differentials_in_one_quantum(
    long int onequantum)
{
  // std::cout << "Invoked with " << onequantum << std::endl;
  std::vector<std::vector<long int> > bn_maps_fix_hom;
  std::vector<long int> bn_diff;
  FMatrix homclass;
  FMatrix barnatan_diff;
  long int source_size;
  long int target_size;
  long int newhom;
  long int absolute_hom;
  std::vector<FMatrix> homologies_at_quantum = chains[onequantum].homologies();
  std::pair<long int, long int> a;
  bool empty = false;
  if (!chains.contains(onequantum + 2))
  {
    empty = true;
    // std::cout << "Cowardly refused computing differential at " << onequantum
    // << " quantum grading " << std::endl; return; // differential has nowhere
    // to go std::cout << "Bravely attempting to compute the differential at "
    // << onequantum << " quantum grading " << std::endl;
  }

  for (long int hc = 0; hc < homologies_at_quantum.size();
       hc++)  //(auto & homclass : homologies_at_quantum)
  {
    absolute_hom = hc + chains[onequantum].mindeg();
    a = std::make_pair(absolute_hom, onequantum);
    bn_maps_fix_hom.clear();
    if (chains[onequantum].homatdeg(hc).cols() == 0)
    {
      // std::cout << "Homology at " << absolute_hom << " quantum = " <<
      // onequantum << " appears to be trivial. " << std::endl;
      continue;
    }
    if (empty)
    {
      barnatan_maps[a] = bn_maps_fix_hom;
      continue;
    }
    newhom = absolute_hom + 1 -
             chains[onequantum + 2].mindeg();  // relative homological grading.
    if (newhom < 0 or newhom >= chains[onequantum + 2].homologies().size())
    {
      barnatan_maps[a] = bn_maps_fix_hom;
      // std::cout << "Homology at " << absolute_hom << " quantum = " <<
      // onequantum << ". Value of newhom = " << newhom +
      // chains[onequantum+2].mindeg() << " out of range " << std::endl;
      continue;
    }
    /*if (chains[onequantum + 2].homatdeg(newhom).cols() == 0) {
        //std::cout << "Homology at " << absolute_hom << " quantum = " <<
    onequantum << ". No homology at newhom = " << newhom +
    chains[onequantum+2].mindeg() << " out of range " << std::endl; continue;
    }*/
    homclass = homologies_at_quantum[hc];
    source_size = homclass.cols();
    target_size = chains[onequantum + 2].homatdeg(newhom).cols();
    barnatan_diff = FMatrix(target_size, source_size);
    for (long int i = 0; i < source_size; i++)
    {
      // std::cout << "Mindeg = " << chains[onequantum].mindeg() << std::endl;
      bn_diff = barnatan_differential(i, absolute_hom, onequantum);
      for (auto& j : bn_diff)
      {
        if (j >= 0)
        {
          barnatan_diff.set(j, i, field(true));
        }
      }
      bn_maps_fix_hom.push_back(bn_diff);
    }
    barnatan_maps[a] = bn_maps_fix_hom;
    barnatan_chainmaps[a] = barnatan_diff;
    barnatan_SNFs[a] = SNF(barnatan_diff);
  }
}

void khovanov::compute_bar_natan_homologies()
{
  long int homgrad;
  long int quantgrad;
  long int size;
  long int previous_size;
  long int next_size;
  std::pair<long int, long int> previous;
  std::pair<long int, long int> next;
  FMatrix F, I, K, E, Im;
  SNF S;
  for (auto& map : barnatan_maps)
  {
    size = size_of_hom(map.first);
    if (size == 0)
    {
      std::cout << "It happened something that should never happen!"
                << std::endl;
      continue;  // this should never happen
    }
    I = FMatrix(size, size);
    E = FMatrix(size, 0);
    for (long int i = 0; i < size; i++) I.set(i, i, field(true));
    homgrad = map.first.first;
    quantgrad = map.first.second;
    previous = std::make_pair(homgrad - 1, quantgrad - 2);
    next = std::make_pair(homgrad + 1, quantgrad + 2);
    previous_size = size_of_hom(previous);
    next_size = size_of_hom(next);
    // std::cout << " Case (" << homgrad << "," << quantgrad << ") previous_size
    // = " << previous_size << " next size = " << next_size << std::endl;
    if (previous_size == 0 and next_size == 0)
    {
      // std::cout << "Empty chain case " << std::endl;
      F = I;
      K = I;
      Im = E;
    }
    if (previous_size == 0 and next_size > 0)
    {
      // std::cout << "Empty previous case" << std::endl;
      S = barnatan_SNFs[map.first];
      F = ComputeStartingHomology(S);
      K = S.Kernel();
      Im = E;
    }
    if (previous_size > 0 and next_size == 0)
    {
      // std::cout << "Empty target case " << homgrad << " " << quantgrad <<
      // std::endl; here is an error
      // F = ComputeHomology(barnatan_SNFs[previous],SNF(FMatrix(size,size)));
      F = ComputeFinalHomology(barnatan_SNFs[previous]);
      K = I;
      Im = barnatan_SNFs[previous].Image();
    }
    if (previous_size > 0 and next_size > 0)
    {
      // std::cout << "All non-empty" << std::endl;
      S = barnatan_SNFs[map.first];
      F = ComputeHomology(barnatan_SNFs[previous], S);
      K = S.Kernel();
      Im = barnatan_SNFs[previous].Image();
    }
    barnatan_hom[map.first] = F;
    barnatan_Kernels[map.first] = K;
    barnatan_Images[map.first] = Im;
  }
}

void khovanov::handle_action_at_bigrading(
    std::pair<std::pair<long int, long int>, long int> myhom)
{
  std::vector<long int> avector;
  std::vector<std::vector<long int> > actionvector;
  long int dim = myhom.second;
  long int homdeg = myhom.first.first;
  long int quantdeg = myhom.first.second;
  long int relhomdeg = homdeg - chains[quantdeg].mindeg();
  if (dim > 0)
  {
    // std::cout << "We are here " << homdeg << "," << quantdeg << std::endl;
    actionvector = myaction.permute_homology(chains[quantdeg], PD, relhomdeg);
    action_maps[myhom.first] = actionvector;
    FMatrix F = FMatrix(dim, dim);
    for (long int j = 0; j < dim; j++)
    {
      for (auto& i : actionvector[j])
      {
        if (0 <= i && i < dim)
        {
          F.set(i, j, field(true));
        }
      }
      F.set(j, j, F.get(j, j) + field(true));  // add diagonal
    }
    SNF S = SNF(F);
    one_plus_tau_maps[myhom.first] = F;
    one_plus_tau_SNFs[myhom.first] = S;
    // std::cout << "We are here " << homdeg << "," << quantdeg << " " <<
    // S.Kernel() <<  std::endl;
    one_plus_tau_Kernels[myhom.first] = S.Kernel();
    // std::cout << "We are here " << homdeg << "," << quantdeg << " " <<
    // S.Kernel() <<  std::endl;
    //       std::cout << "Computed maps for " <<myhom.first.first << ","
    //       <<myhom.first.second << std::endl;
  }
}

void khovanov::load_action_from_file(std::string filename)
{
  std::cout << "Computing action " << std::endl;
  myaction.read_from_file(filename);
  myaction.codes(PD);
  myaction.resolutions(PD);
  std::cout << std::endl;
  // myaction.show_codes();
  // myaction.show_resolutions();
  std::vector<std::thread> mythreads;

  for (auto& myhom : homology_sizes)
  {
    mythreads.push_back(
        std::thread(&khovanov::handle_action_at_bigrading, this, myhom));
  }
  for (auto& th : mythreads) th.join();
  is_action_loaded = true;
  if (!khovanovonly) intersect_kernels_one_plus_tau_and_d();
}
void khovanov::print_action_at_bigrading(std::pair<long, long> bigrad)
{
  // std::cout << "Print action called with "<< bigrad.first << "," <<
  // bigrad.second << std::endl;
  if (!is_action_loaded) return;
  if (!one_plus_tau_maps.contains(bigrad))
  {
    std::cout << "But empty " << std::endl;
    return;
  }
  std::vector<std::vector<long int> > actionvector = action_maps[bigrad];
  std::vector<long int> myvector;
  for (long int i = 0; i < actionvector.size(); i++)
  {
    std::cout << " Generator number " << i;
    myvector = actionvector[i];
    if (myvector.empty())
      std::cout << " is mapped to nowhere" << std::endl;
    else
    {
      if (myvector.size() == 1 and myvector[0] == i)
        std::cout << " is mapped to itself" << std::endl;
      else
      {
        std::cout << ANSI_COLOR_RED << "is not fixed by the symmetry action "
                  << ANSI_COLOR_RESET << std::endl;
        std::cout << ANSI_COLOR_YELLOW
                  << "More precisely, it is mapped to generators ";
        for (auto& genvector : myvector) std::cout << genvector << " ";
        std::cout << ANSI_COLOR_RESET << std::endl;
      }
    }
  }
}
void khovanov::intersect_kernels_one_plus_tau_and_d()
{
  std::pair<long, long> bigrading;
  std::vector<long int> homology_class;
  // long int relhomdeg;
  // long int quantdeg;
  // long int homdeg;
  FMatrix Int;
  bool is_homology_non_trivial;
  for (auto& bn_ker : barnatan_Kernels)
  {
    bigrading = bn_ker.first;
    // homdeg = bigrading.first;
    // quantdeg = bigrading.second;
    if (barnatan_hom[bigrading].cols() > 0)
    {
      // relhomdeg = homdeg - chains[quantdeg].mindeg();

      if (one_plus_tau_Kernels.contains(bigrading))
      {
        //   std::cout << " passed" << std::endl;
        is_homology_non_trivial = false;
        Int = intersection_of_spaces(bn_ker.second,
                                     one_plus_tau_Kernels[bigrading]);
        for (long int j = 0; j < Int.cols(); j++)
        {
          homology_class = trace_homology_element(Int.take_column(j),
                                                  barnatan_hom[bigrading],
                                                  barnatan_Images[bigrading]);
          if (homology_class.size() > 0) is_homology_non_trivial = true;
        }
        one_plus_tau_inter_d[bn_ker.first] = is_homology_non_trivial;
      }
    }
    // else std::cout << " not passed" << std::endl;
  }
}