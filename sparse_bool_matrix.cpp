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

#include "sparse_bool_matrix.h"

#include <algorithm>
#include <memory>
#include <utility>

#include "absl/container/flat_hash_set.h"
#include "sparse_bool_matrix_utils.h"
#include "sparse_low_memory_bool_matrix.h"

namespace
{

void SwapSparseVectors(Index v1, Index v2,
                       std::vector<SparseBoolVector>& vectors,
                       std::vector<SparseBoolVector>& ortogonal_vectors)
{
  for (Index ortogonal_v : vectors[v1])
  {
    bool was_inserted = ortogonal_vectors[ortogonal_v].insert(v2).second;
    if (was_inserted)
    {
      ortogonal_vectors[ortogonal_v].erase(v1);
    }
  }
  for (Index ortogonal_v : vectors[v2])
  {
    bool was_inserted = ortogonal_vectors[ortogonal_v].insert(v1).second;
    if (was_inserted)
    {
      ortogonal_vectors[ortogonal_v].erase(v2);
    }
  }

  using std::swap;
  swap(vectors[v1], vectors[v2]);
}

void DeleteSparseVector(Index index, std::vector<SparseBoolVector>& vectors,
                        std::vector<SparseBoolVector>& ortogonal_vectors)
{
  for (SparseBoolVector& ov : ortogonal_vectors)
  {
    absl::erase_if(ov, [index](Index element) { return element >= index; });
  }
  vectors.erase(vectors.begin() + index);
  for (Index idx = index; idx < vectors.size(); ++idx)
  {
    for (Index ortogonal_idx : vectors[idx])
    {
      ortogonal_vectors[ortogonal_idx].insert(idx);
    }
  }
}

}  // namespace

std::vector<std::vector<bool>> BoolMatrix::GetContents() const
{
  std::vector<bool> empty_row(GetNumberOfColumns(), false);
  std::vector<std::vector<bool>> contents(GetNumberOfRows(), empty_row);
  for (int i = 0; i < GetNumberOfRows(); ++i)
  {
    for (int j : GetRow(i))
    {
      contents[i][j] = true;
    }
  }
  return contents;
}

void PrintTo(const BoolMatrix&, std::ostream* os)
{
  *os << "see the matrix below";
}

bool operator==(const BoolMatrix& lhs, const BoolMatrix& rhs)
{
  if (lhs.GetNumberOfRows() != rhs.GetNumberOfRows() ||
      lhs.GetNumberOfColumns() != rhs.GetNumberOfColumns())
  {
    return false;
  }
  for (Index i = 0; i < lhs.GetNumberOfColumns(); ++i)
  {
    if (lhs.GetColumn(i) != rhs.GetColumn(i))
    {
      return false;
    }
  }
  return true;
}

SparseBoolMatrix::SparseBoolMatrix(int number_of_rows, int number_of_columns)
    : rows_(number_of_rows), columns_(number_of_columns)
{
  TerminateIfMatrixDimensionTooLarge(rows_.size());
  TerminateIfMatrixDimensionTooLarge(columns_.size());
}

SparseBoolMatrix::SparseBoolMatrix(std::vector<SparseBoolVector> rows,
                                   std::vector<SparseBoolVector> columns)
    : rows_(std::move(rows)), columns_(std::move(columns))
{
  TerminateIfMatrixDimensionTooLarge(rows_.size());
  TerminateIfMatrixDimensionTooLarge(columns_.size());
}

SparseBoolMatrix::SparseBoolMatrix(MatrixContentsAsColumns contents)
    : rows_(contents.number_of_rows), columns_(std::move(contents.columns))
{
  TerminateIfMatrixDimensionTooLarge(rows_.size());
  TerminateIfMatrixDimensionTooLarge(columns_.size());
  for (Index column = 0; column < columns_.size(); ++column)
  {
    for (Index row : columns_[column])
    {
      rows_[row].insert(column);
    }
  }
}

SparseBoolMatrix::SparseBoolMatrix(MatrixContentsAsRows contents)
    : rows_(std::move(contents.rows)), columns_(contents.number_of_columns)
{
  TerminateIfMatrixDimensionTooLarge(rows_.size());
  TerminateIfMatrixDimensionTooLarge(columns_.size());
  for (Index row = 0; row < rows_.size(); ++row)
  {
    for (Index column : rows_[row])
    {
      columns_[column].insert(row);
    }
  }
}

int SparseBoolMatrix::GetNumberOfRows() const { return rows_.size(); }
int SparseBoolMatrix::GetNumberOfColumns() const { return columns_.size(); }

bool SparseBoolMatrix::Get(Index row, Index column) const
{
  return rows_[row].contains(column);
}

void SparseBoolMatrix::Set(Index row, Index column, bool value)
{
  if (value)
  {
    bool element_inserted = rows_[row].insert(column).second;
    if (element_inserted)
    {
      columns_[column].insert(row);
    }
  }
  else
  {
    bool element_removed = (rows_[row].erase(column) == 1);
    if (element_removed)
    {
      columns_[column].erase(row);
    }
  }
}

std::unique_ptr<BoolMatrix> SparseBoolMatrix::GetCopyWithFirstNColumns(
    int number_of_columns_to_copy) const
{
  // We don't expect many (or any) alterations to `rows` so we copy the whole
  // rows and remove the unnecessary elements.
  std::vector<SparseBoolVector> rows(rows_);
  std::vector<SparseBoolVector> columns(
      columns_.begin(), columns_.begin() + number_of_columns_to_copy);
  for (Index column = number_of_columns_to_copy; column < columns_.size();
       ++column)
  {
    for (Index row : columns_[column])
    {
      rows[row].erase(column);
    }
  }
  return std::make_unique<SparseBoolMatrix>(std::move(rows),
                                            std::move(columns));
}

bool SparseBoolMatrix::IsZero() const
{
  for (const SparseBoolVector& row : rows_)
  {
    if (!row.empty())
    {
      return false;
    }
  }
  return true;
}

const SparseBoolVector& SparseBoolMatrix::GetRow(Index row) const
{
  return rows_[row];
}

const SparseBoolVector& SparseBoolMatrix::GetColumn(Index column) const
{
  return columns_[column];
}

void SparseBoolMatrix::SumRows(Index which_row, Index to_which_row)
{
  for (Index column : rows_[which_row])
  {
    auto [it1, was_inserted1] = rows_[to_which_row].insert(column);
    auto [it2, was_inserted2] = columns_[column].insert(to_which_row);
    assert(was_inserted1 == was_inserted2);
    if (!was_inserted1)
    {
      rows_[to_which_row].erase(it1);
      columns_[column].erase(it2);
    }
  }
}

void SparseBoolMatrix::SumColumns(Index which_column, Index to_which_column)
{
  for (Index row : columns_[which_column])
  {
    auto [it1, was_inserted1] = columns_[to_which_column].insert(row);
    auto [it2, was_inserted2] = rows_[row].insert(to_which_column);
    assert(was_inserted1 == was_inserted2);
    if (!was_inserted1)
    {
      columns_[to_which_column].erase(it1);
      rows_[row].erase(it2);
    }
  }
}

void SparseBoolMatrix::SwapRows(Index row_1, Index row_2)
{
  SwapSparseVectors(row_1, row_2, rows_, columns_);
}

void SparseBoolMatrix::SwapColumns(Index column_1, Index column_2)
{
  SwapSparseVectors(column_1, column_2, columns_, rows_);
}

void SparseBoolMatrix::DeleteColumn(Index column)
{
  DeleteSparseVector(column, columns_, rows_);
}

void SparseBoolMatrix::DeleteFirstNColumnsAndEmptyOnes(
    int number_of_first_columns_to_delete)
{
  assert(number_of_first_columns_to_delete > 0 &&
         number_of_first_columns_to_delete <= GetNumberOfColumns());
  columns_.erase(columns_.begin(),
                 columns_.begin() + number_of_first_columns_to_delete);
  std::erase_if(columns_,
                [](const SparseBoolVector& column) { return column.empty(); });
  // Rebuild rows_ from scratch.
  rows_ = std::vector<SparseBoolVector>(rows_.size());
  for (Index column = 0; column < columns_.size(); ++column)
  {
    for (Index row : columns_[column])
    {
      rows_[row].insert(column);
    }
  }
}

void SparseBoolMatrix::AppendColumn(SparseBoolVector column)
{
  Index new_column_index = columns_.size();
  for (Index row : column)
  {
    rows_[row].insert(new_column_index);
  }
  columns_.push_back(std::move(column));
}

void SparseBoolMatrix::AddColumnToExistingOne(
    const SparseBoolVector& added_column, Index to_which_column)
{
  SparseBoolVector& recipient_column = columns_[to_which_column];
  for (Index row : added_column)
  {
    auto [it1, was_inserted1] = recipient_column.insert(row);
    auto [it2, was_inserted2] = rows_[row].insert(to_which_column);
    assert(was_inserted1 == was_inserted2);
    if (!was_inserted1)
    {
      recipient_column.erase(it1);
      rows_[row].erase(it2);
    }
  }
}

void SparseBoolMatrix::DeleteAllButOneElementInRow(
    Index row, std::optional<Index> column_to_keep)
{
  assert(!rows_[row].empty());
  Index first_column = columns_.size();
  for (Index column : rows_[row])
  {
    if (column < first_column)
    {
      first_column = column;
    }
    columns_[column].erase(row);
  }
  rows_[row].clear();
  Index column_to_restore =
      column_to_keep.has_value() ? *column_to_keep : first_column;
  rows_[row].insert(column_to_restore);
  columns_[column_to_restore].insert(row);
}

std::unique_ptr<BoolMatrix> SparseBoolMatrix::Clone() const
{
  return std::make_unique<SparseBoolMatrix>(*this);
}

std::unique_ptr<BoolMatrix> SparseBoolMatrix::MoveToCapableInstance() &&
{
  return std::make_unique<SparseBoolMatrix>(std::move(rows_),
                                            std::move(columns_));
}

std::unique_ptr<BoolMatrix> SparseBoolMatrix::MoveToLowMemoryInstance() &&
{
  return std::make_unique<SparseLowMemoryBoolMatrix>(
      MatrixContentsAsColumns{.number_of_rows = static_cast<int>(rows_.size()),
                              .columns = std::move(columns_)});
}

Index GetSmallestIndex(const SparseBoolVector& vector)
{
  // The smallest element could be cached in SparseBoolMatrix
  // if it becomes a bottleneck.
  // It would require updating the smallest element during matrix modifications.
  return *std::min_element(vector.begin(), vector.end());
}