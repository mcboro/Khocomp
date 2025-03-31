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

#include "sparse_bool_matrix_fast_row_operations.h"

#include "sparse_bool_matrix.h"
#include "sparse_bool_matrix_utils.h"

SparseBoolMatrixFastRowOperations
SparseBoolMatrixFastRowOperations::CreateIdentityMatrix(int number_of_rows)
{
  SparseBoolMatrixFastRowOperations matrix(number_of_rows, number_of_rows);
  for (Index row = 0; row < number_of_rows; ++row)
  {
    matrix.rows_[row].insert(row);
  }
  return matrix;
}

void SparseBoolMatrixFastRowOperations::SumRows(Index which_row,
                                                Index to_which_row)
{
  for (Index column : rows_[which_row])
  {
    auto [it, was_inserted] = rows_[to_which_row].insert(column);
    if (!was_inserted)
    {
      rows_[to_which_row].erase(it);
    }
  }
}

void SparseBoolMatrixFastRowOperations::SwapRows(Index row_1, Index row_2)
{
  using std::swap;
  swap(rows_[row_1], rows_[row_2]);
}

MatrixContentsAsRows
SparseBoolMatrixFastRowOperations::ExtractContents() &&
{
  return {.number_of_columns = number_of_columns_, .rows = std::move(rows_)};
}

SparseBoolMatrixFastRowOperations::SparseBoolMatrixFastRowOperations(
    int number_of_rows, int number_of_columns)
    : rows_(number_of_rows), number_of_columns_(number_of_columns)
{
  TerminateIfMatrixDimensionTooLarge(number_of_rows);
  TerminateIfMatrixDimensionTooLarge(number_of_columns);
}