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

#ifndef ALGEBRA_H
#define ALGEBRA_H

#include <iostream>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "sparse_bool_matrix.h"

class field
{
 public:
  explicit field(bool value);

  field operator-() const;
  field operator+(const field&) const;
  field operator-(const field&) const;
  field operator*(const field&) const;

  operator bool() const;

  friend std::ostream& operator<<(std::ostream& os, const field& f);

 private:
  bool val;
};

class FMatrix
{
 public:
  FMatrix();

  FMatrix(const FMatrix& other);
  FMatrix& operator=(const FMatrix& other);
  FMatrix(FMatrix&&) = default;
  FMatrix& operator=(FMatrix&&) = default;

  // zero matrix
  FMatrix(long int rows, long int cols);

  FMatrix(std::unique_ptr<BoolMatrix> matrix) : matrix_(std::move(matrix)) {}

  long int rows() const;
  long int cols() const;
  bool isZero() const;

  field get(long int i, long int j) const;
  void set(long int i, long int j, field f);

  void sum_rows(long int which_row, long int to_which);
  void sum_cols(long int which_col, long int to_which);
  void switch_rows(long int which_row, long int with_which);
  void switch_cols(long int which_col, long int with_which);

  FMatrix take_column(long int which) const;

  // deletes column
  void kill_column(long int which);

  // adds columnwise Another after the matrix
  FMatrix Concatenate(const FMatrix& Another) const;

  FMatrix operator*(const FMatrix& f) const;

  // Convert to a matrix supporting a limited set of operations but using less
  // memory.
  void ReduceMemoryFootprint();

  // Convert to a matrix supporting all operations.
  void ConvertToCapableMatrix();

  // These accessors are necessary for some performance optimizations.
  BoolMatrix* matrix();
  const BoolMatrix* matrix() const;

  // Used for debugging and testing.
  std::vector<std::vector<bool>> GetContents() const;

  // This is a helper function for GoogleTest. It's not meant to provide
  // useful representation of the matrix, but to point out that below there
  // should be a readable comparison between the expected and actual matrices.
  // It's meant to be used with a custom matcher (FMatrixEq).
  // More: https://google.github.io/googletest/advanced.html.
  friend void PrintTo(const FMatrix&, std::ostream* os);

  bool operator==(const FMatrix&) const = default;

 private:
  std::unique_ptr<BoolMatrix> matrix_;
};

std::vector<long int> column_simplify(FMatrix& C);
void simplify_against_matrix(FMatrix& C, const FMatrix& D);

class SNF
{
 private:
  unsigned long int rank, cl, rw;
  FMatrix KernelMatrix;
  FMatrix ImageMatrix;
  bool isInitializedvalue;
  std::vector<long int> signature_columns;

 public:
  struct Composition
  {
    FMatrix left;
    FMatrix diagonal;
    FMatrix right_inverse;
  };

  bool isInitialized() { return isInitializedvalue; };
  std::optional<Composition> MatrixRecompute(FMatrix initial_matrix,
                                             bool return_composition = false);
  // Creates an SNF that has to be later initialized by a matrix.
  SNF() { isInitializedvalue = false; };
  SNF(FMatrix initial_matrix) { MatrixRecompute(std::move(initial_matrix)); };

  const FMatrix& Kernel() const { return KernelMatrix; }
  const FMatrix& Image() const { return ImageMatrix; }
  unsigned long int rows() const { return rw; }
};

FMatrix ComputeHomology(const SNF& D2, const SNF& D1);
FMatrix ComputeFinalHomology(const SNF& D2);
FMatrix ComputeStartingHomology(const SNF&);
bool are_equal_mod_image(const FMatrix& A, const FMatrix& B, const FMatrix& C);
bool are_linearly_dependent(const FMatrix& A, const FMatrix& B);
long int strip_zero_columns(FMatrix& A);
FMatrix intersection_of_spaces(const FMatrix& A, const FMatrix& B);

std::vector<long int> trace_homology_element(const FMatrix& A, const FMatrix& B,
                                             const FMatrix& C);
std::vector<long int> trace_homology_element(const FMatrix& A,
                                             const FMatrix& B);

class ChainComplex
{
 private:
  std::vector<FMatrix> ChainMaps;
  std::vector<SNF> AllMaps;
  std::vector<FMatrix> Homologies;
  void ComputeAllSNF();
  void ComputeHomologyGroup();
  bool onlyone;
  long int whichone;

 public:
  void from_empty(long int i);
  const FMatrix& giveChainMap(long i) const { return ChainMaps[i]; };
  const SNF& giveSNF(long i) const { return AllMaps[i]; };
  const FMatrix& giveKernel(long i) const { return AllMaps[i].Kernel(); };
  const FMatrix& giveImage(long i) const
  {
    return AllMaps[i].Image();
  };  // returns the image of the i-th differential as a subspace of C_{i+1}.
  const FMatrix& giveHomology(long i) const { return Homologies[i]; };
  const std::vector<FMatrix>& giveAllHomologies() const { return Homologies; };
  ChainComplex()
  {
    onlyone = false;
    whichone = 0;
  };
  ChainComplex(const std::vector<FMatrix>& Maps, bool reallyonlyone = false,
               long int whichistheone = 0)
  {
    onlyone = reallyonlyone;
    whichone = whichistheone;
    ChainMaps = Maps;
    for (long int i = 0; i < ChainMaps.size(); i++) AllMaps.push_back(SNF());
    ComputeAllSNF();
    ComputeHomologyGroup();
  }
  long int mapsize() const { return ChainMaps.size(); };
  long int homsize() const { return Homologies.size(); };
};

#endif