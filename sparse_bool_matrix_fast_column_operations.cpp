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

#include "sparse_bool_matrix_fast_column_operations.h"

#include "sparse_bool_matrix.h"
#include "sparse_bool_matrix_utils.h"

SparseBoolMatrixFastColumnOperations
SparseBoolMatrixFastColumnOperations::CreateIdentityMatrix(
    int number_of_columns)
{
  SparseBoolMatrixFastColumnOperations matrix(number_of_columns,
                                              number_of_columns);
  for (Index column = 0; column < number_of_columns; ++column)
  {
    matrix.columns_[column].insert(column);
  }
  return matrix;
}

void SparseBoolMatrixFastColumnOperations::SumColumns(Index which_column,
                                                      Index to_which_column)
{
  for (Index row : columns_[which_column])
  {
    auto [it, was_inserted] = columns_[to_which_column].insert(row);
    if (!was_inserted)
    {
      columns_[to_which_column].erase(it);
    }
  }
}

void SparseBoolMatrixFastColumnOperations::SwapColumns(Index column_1,
                                                       Index column_2)
{
  using std::swap;
  swap(columns_[column_1], columns_[column_2]);
}

MatrixContentsAsColumns
SparseBoolMatrixFastColumnOperations::ExtractContents() &&
{
  return {.number_of_rows = number_of_rows_, .columns = std::move(columns_)};
}

SparseBoolMatrixFastColumnOperations::SparseBoolMatrixFastColumnOperations(
    int number_of_rows, int number_of_columns)
    : number_of_rows_(number_of_rows), columns_(number_of_columns)
{
  TerminateIfMatrixDimensionTooLarge(number_of_rows);
  TerminateIfMatrixDimensionTooLarge(number_of_columns);
}