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

#ifndef SPARSE_BOOL_MATRIX_H
#define SPARSE_BOOL_MATRIX_H

#include <cassert>
#include <memory>
#include <optional>
#include <vector>

#include "absl/container/flat_hash_set.h"

// It is safe to use unsigned types for Index. Index should be a type
// that is safe to cast to int; otherwise, the code using SparseBoolMatrix
// may have to be adapted since it often casts to int.
// An example of a type that could limit memory usage is uint16_t.
// Most likely, Index will be one of (uint16_t, int32_t, int).
using Index = int32_t;

// A collection storing the indices of "true" elements in a bool vector.
using SparseBoolVector = absl::flat_hash_set<Index>;

struct MatrixContentsAsColumns
{
  int number_of_rows;
  std::vector<SparseBoolVector> columns;
};

struct MatrixContentsAsRows
{
  int number_of_columns;
  std::vector<SparseBoolVector> rows;
};

class BoolMatrix
{
 public:
  virtual int GetNumberOfRows() const = 0;
  virtual int GetNumberOfColumns() const = 0;

  virtual bool Get(Index row, Index column) const = 0;
  virtual void Set(Index row, Index column, bool value) = 0;

  virtual std::unique_ptr<BoolMatrix> GetCopyWithFirstNColumns(
      int number_of_columns_to_copy) const = 0;

  virtual bool IsZero() const = 0;

  // Indices in rows and columns are not sorted.
  virtual const SparseBoolVector& GetRow(Index row) const = 0;
  virtual const SparseBoolVector& GetColumn(Index column) const = 0;

  virtual void SumRows(Index which_row, Index to_which_row) = 0;
  virtual void SumColumns(Index which_column, Index to_which_column) = 0;

  virtual void SwapRows(Index row_1, Index row_2) = 0;
  virtual void SwapColumns(Index column_1, Index column_2) = 0;

  virtual void DeleteColumn(Index column) = 0;

  virtual void DeleteFirstNColumnsAndEmptyOnes(
      int number_of_first_columns_to_delete) = 0;

  virtual void AppendColumn(SparseBoolVector column) = 0;

  virtual void AddColumnToExistingOne(const SparseBoolVector& added_column,
                                      Index to_which_column) = 0;

  // Deletes all but one element in the given row. If `column_to_keep` is not
  // specified, it keeps the smallest element.
  virtual void DeleteAllButOneElementInRow(
      Index row, std::optional<Index> column_to_keep) = 0;

  virtual std::unique_ptr<BoolMatrix> Clone() const = 0;
  virtual std::unique_ptr<BoolMatrix> MoveToCapableInstance() && = 0;
  virtual std::unique_ptr<BoolMatrix> MoveToLowMemoryInstance() && = 0;

  // For debugging and testing.
  std::vector<std::vector<bool>> GetContents() const;

  // This is a helper function for GoogleTest. It's not meant to provide
  // useful representation of the matrix, but to point out that below there
  // should be a readable comparison between the expected and actual matrices.
  // It's meant to be used with a custom matcher (MatrixEq).
  // More: https://google.github.io/googletest/advanced.html.
  friend void PrintTo(const BoolMatrix&, std::ostream* os);

  virtual ~BoolMatrix() = default;
};

bool operator==(const BoolMatrix& lhs, const BoolMatrix& rhs);

class SparseBoolMatrix : public BoolMatrix
{
 public:
  // Constructs a zero matrix.
  SparseBoolMatrix(int number_of_rows, int number_of_columns);

  SparseBoolMatrix(std::vector<SparseBoolVector> rows,
                   std::vector<SparseBoolVector> columns);

  SparseBoolMatrix(MatrixContentsAsColumns contents);
  SparseBoolMatrix(MatrixContentsAsRows contents);

  SparseBoolMatrix(const SparseBoolMatrix&) = default;
  SparseBoolMatrix& operator=(const SparseBoolMatrix&) = default;
  SparseBoolMatrix(SparseBoolMatrix&&) = default;
  SparseBoolMatrix& operator=(SparseBoolMatrix&&) = default;

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
  std::vector<SparseBoolVector> rows_;
  std::vector<SparseBoolVector> columns_;
};

Index GetSmallestIndex(const SparseBoolVector& vector);

#endif  // SPARSE_BOOL_MATRIX_H