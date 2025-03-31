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

#ifndef SPARSE_LOW_MEMORY_BOOL_MATRIX_H
#define SPARSE_LOW_MEMORY_BOOL_MATRIX_H

#include <vector>

#include "sparse_bool_matrix.h"

// An implementation using less memory that is not guaranteed to support all
// operations.
// The code should be refactored to use the proper interface, displaying only
// the supported methods, but we deemed it not worth it for now.
class SparseLowMemoryBoolMatrix : public BoolMatrix
{
 public:
  SparseLowMemoryBoolMatrix(MatrixContentsAsColumns contents);

  SparseLowMemoryBoolMatrix(const SparseLowMemoryBoolMatrix&) = default;
  SparseLowMemoryBoolMatrix& operator=(const SparseLowMemoryBoolMatrix&) =
      default;
  SparseLowMemoryBoolMatrix(SparseLowMemoryBoolMatrix&&) = default;
  SparseLowMemoryBoolMatrix& operator=(SparseLowMemoryBoolMatrix&&) = default;

  int GetNumberOfRows() const override;
  int GetNumberOfColumns() const override;

  bool Get(Index row, Index column) const override;
  void Set(Index row, Index column, bool value) override;

  std::unique_ptr<BoolMatrix> GetCopyWithFirstNColumns(
      int number_of_columns_to_copy) const override;

  bool IsZero() const override;

  const SparseBoolVector& GetRow(Index row) const override;
  const SparseBoolVector& GetColumn(Index column) const override;

  void SumRows(Index which_row, Index to_which_row) override;
  void SumColumns(Index which_column, Index to_which_column) override;

  void SwapRows(Index row_1, Index row_2) override;
  void SwapColumns(Index column_1, Index column_2) override;

  void DeleteColumn(Index column) override;

  void DeleteFirstNColumnsAndEmptyOnes(
      int number_of_first_columns_to_delete) override;

  void AppendColumn(SparseBoolVector column) override;

  void AddColumnToExistingOne(const SparseBoolVector& added_column,
                              Index to_which_column) override;

  void DeleteAllButOneElementInRow(
      Index row, std::optional<Index> column_to_keep) override;

  virtual std::unique_ptr<BoolMatrix> Clone() const override;
  virtual std::unique_ptr<BoolMatrix> MoveToCapableInstance() && override;
  virtual std::unique_ptr<BoolMatrix> MoveToLowMemoryInstance() && override;

 private:
  int number_of_rows_;
  std::vector<SparseBoolVector> columns_;
};

#endif  // SPARSE_LOW_MEMORY_BOOL_MATRIX_H