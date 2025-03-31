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

#ifndef CHAINCOMPLEX_H
#define CHAINCOMPLEX_H

#include <vector>

#include "algebra.h"
#include "pdcode.h"

struct gradingshift
{
 public:
  long int quantum;
  long int homological;
  gradingshift(long int a, long int b);
  gradingshift(std::pair<long int, long int> crossings)
  {
    gradingshift(crossings.first, crossings.second);
  }
  gradingshift()
  {
    quantum = 0;
    homological = 0;
  }
};

class KhovanovChainComplex : public statecomplex
{
 private:
  ChainComplex Chain_At_Grading;
  bool onlyone;
  long int theone;
  bool isInitializedvalue;

 public:
  void initialize();
  bool isInitialized() { return isInitializedvalue; };
  // default parameters in the constructor are relative, that is, without
  // grading shift
  explicit KhovanovChainComplex(std::vector<state> all_with_given_quantum,
                                bool only = false, long int the = 0)
      : statecomplex(all_with_given_quantum)
  {
    onlyone = only;
    theone = the;
    isInitializedvalue = false;
  }

  KhovanovChainComplex() : statecomplex()
  {
    isInitializedvalue = false;
    onlyone = false;
    theone = 0;
  };
  std::vector<FMatrix> homologies() const
  {
    return Chain_At_Grading.giveAllHomologies();
  };
  FMatrix homatdeg(long int i) const
  {
    return Chain_At_Grading.giveHomology(i);
  };
  FMatrix imageatdeg(long int i);
  void export_homology(long int shift_homological, long int shift_quantum);
  void export_homology() { export_homology(0, 0); };
  void export_homology(gradingshift& mygradshift)
  {
    export_homology(mygradshift.homological, mygradshift.quantum);
  };
  FMatrix export_matrix(long int i) const
  {
    return Chain_At_Grading.giveChainMap(i);
  };
  FMatrix export_kernel(long int i) const
  {
    return Chain_At_Grading.giveKernel(i);
  };
  FMatrix export_image(long int i) const
  {
    return Chain_At_Grading.giveImage(i);
  };
  long int len() const { return Chain_At_Grading.mapsize(); }
  long int homlen() const { return Chain_At_Grading.homsize(); }
  bool are_homology_equal(const FMatrix& A, const FMatrix& B, long int degree);
  std::vector<long int> trace_class_in_homology(const FMatrix& A,
                                                long degree) const;
  FMatrix vector_of_barnatan_map(const KhovanovChainComplex& nextchaincomplex,
                                 FMatrix statevector, long int statehom);
};

#endif
