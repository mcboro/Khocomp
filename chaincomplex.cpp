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

#include "chaincomplex.h"

#include <iostream>

#include "algebra.h"

#define ANSI_COLOR_RED "\x1b[31m"
#define ANSI_COLOR_GREEN "\x1b[32m"
#define ANSI_COLOR_YELLOW "\x1b[33m"
#define ANSI_COLOR_BLUE "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN "\x1b[36m"
#define ANSI_COLOR_RESET "\x1b[0m"

void KhovanovChainComplex::initialize()
{
  isInitializedvalue = true;
  // std::cout << "Computing matrices " << std::endl;
  std::vector<state> Target;
  std::vector<state> Source;
  FMatrix M;
  std::vector<FMatrix> aux_maps;
  state s1;
  // state s2;
  std::vector<state> DSource;

  for (long int homgrad = 0; homgrad < statesincomplex.size() - 1; homgrad++)
  {
    Source = statesincomplex[homgrad];
    Target = statesincomplex[homgrad + 1];
    M = FMatrix(Target.size(), Source.size());
    for (long int i = 0; i < Source.size(); i++)
    {
      s1 = Source[i];
      DSource = s1.go_with_injected_arrows();
      for (auto& s2 : DSource)
      {
        for (long int j = 0; j < Target.size(); j++)
        {
          if (Target[j] == s2) M.set(j, i, field(true));
        }
      }
    }
    aux_maps.push_back(M);
  }
  if (aux_maps.size() > 0)
    Chain_At_Grading =
        ChainComplex(aux_maps, onlyone, theone - minimalhomdegree);
  else
  {
    ChainComplex C;
    C.from_empty(statesincomplex[0].size());
    Chain_At_Grading = C;
  }
}

void KhovanovChainComplex::export_homology(long int shift_homological,
                                           long int shift_quantum)
{
  // std::cout << "Grading shifts " << shift_homological << " quantum " <<
  // shift_quantum << std::endl; std::cout << "Homology at quantum grading " <<
  // quantdegree+shift_quantum << std::endl;
  long int homdeg = minimalhomdegree;
  for (auto& H : Chain_At_Grading.giveAllHomologies())
  {
    if (H.cols() > 0)
    {
      std::cout << "Homology groups at quantum grading "
                << quantdegree + shift_quantum << " homological grading "
                << homdeg + shift_homological << " has dimension " << H.cols()
                << std::endl;
      // std::cout << "Generators" << std::endl << H << std::endl << std::endl;
    }
    homdeg++;
  }
  //  std::cout << "Thanks for using our software" << std::endl;
}

bool KhovanovChainComplex::are_homology_equal(const FMatrix& A,
                                              const FMatrix& B, long int degree)
{
  // we take two vectors (not matrices) A and B.
  // we check if they are equal modulo image of d:degree-1 -> degree.
  // for this we take the matrix  (A | B | Image)
  // and see if after some moves, we get a zero column
  // we deduce that the matrix is not trivial.
  if (degree == 0)
  {
    // std::cout << "Degree zero case " << std::endl;
    return are_linearly_dependent(A, B);
    // if there's no previous homology, we don't need to compare.
  }
  FMatrix C = Chain_At_Grading.giveSNF(degree - 1).Image();
  if (C.cols() == 0)
  {
    std::cout << "Zero image case " << std::endl;
    return are_linearly_dependent(A, B);
  }
  /*
  std::cout << "A rank = " << A.rows() << "," << A.cols() << std::endl;
  std::cout << "B rank = " << B.rows() << "," << B.cols() << std::endl;
  std::cout << "Image rank = " << C.rows() << "," << C.cols() << std::endl <<
  std::endl;
  */

  return are_equal_mod_image(A, B,
                             Chain_At_Grading.giveSNF(degree - 1).Image());
}
std::vector<long int> KhovanovChainComplex::trace_class_in_homology(
    const FMatrix& A, long int degree) const
{
  // this function is essentially the same as the previous one
  // it takes a vector A
  // and returns a number b, like 101110, which means that A is equal to
  // first+second+third+fifth generator of homology. we take the relative degree
  // as the input
  FMatrix H = homatdeg(degree);
  if (H.cols() == 0)
  {
    std::vector<long int> a;
    return a;
  }  // no homology, no problem
  if (degree == 0)
  {
    // std::cout << "Degree zero case " << std::endl;
    return trace_homology_element(A, H);
    // if there's no previous homology, we don't need to compare.
  }
  FMatrix C = Chain_At_Grading.giveSNF(degree - 1).Image();
  if (C.cols() == 0)
  {
    std::cout << "Zero image case " << std::endl;
    return trace_homology_element(A, H);
  }
  return trace_homology_element(A, H, C);
}

FMatrix KhovanovChainComplex::imageatdeg(long int degree)
{
  if (degree >= 0) return (Chain_At_Grading.giveSNF(degree).Image());
  return FMatrix(Chain_At_Grading.giveChainMap(0).rows(), 0);
}

FMatrix KhovanovChainComplex::vector_of_barnatan_map(
    const KhovanovChainComplex& nextchaincomplex, FMatrix statevector,
    long int statehom)
{
  int degreeshift = 1;
  // takes a vector statevector of states (eventually: homology class) and
  // transforms it into a new vector
  FMatrix empty = FMatrix(0, 0);
  if (statehom < mindeg()) return empty;            // this should never happen
  if (statehom >= mindeg() + size()) return empty;  // this should never happen
  // the next two possibilities are actually important
  // the BarNatan differential can lead out of the part where it is defined.
  if (statehom + degreeshift < nextchaincomplex.mindeg())
  {
    // std::cout << "Size too small. " << statehom+degreeshift << " " <<
    // nextchaincomplex.mindeg() << std::endl;
    return empty;
  }
  if (statehom + degreeshift >=
      nextchaincomplex.mindeg() + nextchaincomplex.size())
  {
    // std::cout << "Size too large. " << statehom+degreeshift << " " <<
    // nextchaincomplex.mindeg()+nextchaincomplex.size() << std::endl;
    return empty;
  }
  // create the matrix with the number of rows equal to the number of elements
  // in the relevant homological grading
  FMatrix final_vector(
      nextchaincomplex.at(statehom + degreeshift - nextchaincomplex.mindeg())
          .size(),
      1);
  std::vector<long int> output;
  for (long int stateindex = 0; stateindex < statevector.rows(); stateindex++)
  {
    if (statevector.get(stateindex, 0))
    {
      output = trace_barnatan_map(nextchaincomplex, stateindex, statehom);
      for (auto& opk : output)
        final_vector.set(opk, 0, final_vector.get(opk, 0) + field(true));
    }
  }
  return final_vector;
}

gradingshift::gradingshift(long int a, long int b)
{
  homological = -b;
  quantum = a - 2 * b;
  std::cout << "Input grading data " << a << " and " << b << std::endl;
  std::cout << "Created grading shift. Homological = " << homological
            << " quantum=" << quantum << std::endl;
}
