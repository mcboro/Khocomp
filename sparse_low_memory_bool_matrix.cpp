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

#include "sparse_low_memory_bool_matrix.h"

#include <memory>

#include "absl/log/check.h"
#include "sparse_bool_matrix.h"

namespace
{

constexpr char kUnimplementedMessage[] =
    "This method is unimplemented. It could be supported in an unefficient "
    "way, but it's better to restructure the code to use a different matrix "
    "implementation for row operations.";

}  // namespace

SparseLowMemoryBoolMatrix::SparseLowMemoryBoolMatrix(
    MatrixContentsAsColumns contents)
    : number_of_rows_(contents.number_of_rows),
      columns_(std::move(contents.columns))
{
}

int SparseLowMemoryBoolMatrix::GetNumberOfRows() const
{
  return number_of_rows_;
}
int SparseLowMemoryBoolMatrix::GetNumberOfColumns() const
{
  return columns_.size();
}

bool SparseLowMemoryBoolMatrix::Get(Index row, Index column) const
{
  return columns_[column].contains(row);
}

void SparseLowMemoryBoolMatrix::Set(Index row, Index column, bool value)
{
  if (value)
  {
    columns_[column].insert(row);
  }
  else
  {
    columns_[column].erase(row);
  }
}

std::unique_ptr<BoolMatrix> SparseLowMemoryBoolMatrix::GetCopyWithFirstNColumns(
    int number_of_columns_to_copy) const
{
  std::vector<SparseBoolVector> columns(
      columns_.begin(), columns_.begin() + number_of_columns_to_copy);
  return std::make_unique<SparseLowMemoryBoolMatrix>(MatrixContentsAsColumns{
      .number_of_rows = number_of_rows_, .columns = std::move(columns)});
}

bool SparseLowMemoryBoolMatrix::IsZero() const
{
  for (const SparseBoolVector& column : columns_)
  {
    if (!column.empty())
    {
      return false;
    }
  }
  return true;
}

const SparseBoolVector& SparseLowMemoryBoolMatrix::GetRow(Index) const
{
  CHECK(false) << kUnimplementedMessage;
}

const SparseBoolVector& SparseLowMemoryBoolMatrix::GetColumn(Index column) const
{
  return columns_[column];
}

void SparseLowMemoryBoolMatrix::SumRows(Index, Index)
{
  CHECK(false) << kUnimplementedMessage;
}

void SparseLowMemoryBoolMatrix::SumColumns(Index which_column,
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

void SparseLowMemoryBoolMatrix::SwapRows(Index, Index)
{
  CHECK(false) << kUnimplementedMessage;
}

void SparseLowMemoryBoolMatrix::SwapColumns(Index column_1, Index column_2)
{
  using std::swap;
  swap(columns_[column_1], columns_[column_2]);
}

void SparseLowMemoryBoolMatrix::DeleteColumn(Index column)
{
  columns_.erase(columns_.begin() + column);
}

void SparseLowMemoryBoolMatrix::DeleteFirstNColumnsAndEmptyOnes(
    int number_of_first_columns_to_delete)
{
  assert(number_of_first_columns_to_delete > 0 &&
         number_of_first_columns_to_delete <= GetNumberOfColumns());
  columns_.erase(columns_.begin(),
                 columns_.begin() + number_of_first_columns_to_delete);
  std::erase_if(columns_,
                [](const SparseBoolVector& column) { return column.empty(); });
}

void SparseLowMemoryBoolMatrix::AppendColumn(SparseBoolVector column)
{
  columns_.push_back(std::move(column));
}

void SparseLowMemoryBoolMatrix::AddColumnToExistingOne(
    const SparseBoolVector& added_column, Index to_which_column)
{
  SparseBoolVector& recipient_column = columns_[to_which_column];
  for (Index row : added_column)
  {
    auto [it, was_inserted] = recipient_column.insert(row);
    if (!was_inserted)
    {
      recipient_column.erase(it);
    }
  }
}

void SparseLowMemoryBoolMatrix::DeleteAllButOneElementInRow(
    Index, std::optional<Index>)
{
  CHECK(false) << kUnimplementedMessage;
}

std::unique_ptr<BoolMatrix> SparseLowMemoryBoolMatrix::Clone() const
{
  return std::make_unique<SparseLowMemoryBoolMatrix>(*this);
}

std::unique_ptr<BoolMatrix>
SparseLowMemoryBoolMatrix::MoveToCapableInstance() &&
{
  return std::make_unique<SparseBoolMatrix>(MatrixContentsAsColumns{
      .number_of_rows = number_of_rows_, .columns = columns_});
}

std::unique_ptr<BoolMatrix>
SparseLowMemoryBoolMatrix::MoveToLowMemoryInstance() &&
{
  return std::make_unique<SparseLowMemoryBoolMatrix>(MatrixContentsAsColumns{
      .number_of_rows = number_of_rows_, .columns = std::move(columns_)});
}