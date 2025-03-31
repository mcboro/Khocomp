// Copyright 2025 Jacek Olesiak
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

#ifndef SPARSE_BOOL_MATRIX_FAST_ROW_OPERATIONS_H
#define SPARSE_BOOL_MATRIX_FAST_ROW_OPERATIONS_H

#include <vector>

#include "sparse_bool_matrix.h"

class SparseBoolMatrixFastRowOperations
{
 public:
  static SparseBoolMatrixFastRowOperations CreateIdentityMatrix(
      int number_of_rows);

  void SumRows(Index which_row, Index to_which_row);

  void SwapRows(Index row_1, Index row_2);

  MatrixContentsAsRows ExtractContents() &&;

 private:
  SparseBoolMatrixFastRowOperations(int number_of_rows, int number_of_columns);

  std::vector<SparseBoolVector> rows_;
  int number_of_columns_;
};

#endif  // SPARSE_BOOL_MATRIX_FAST_ROW_OPERATIONS_H