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

#include <memory>
#include <optional>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "test_util.h"

namespace
{

using ::testing::UnorderedElementsAre;

// A custom matcher used to compare SparseBoolMatrix objects.
MATCHER_P(MatrixEq, expected, "")
{
  return testing::ExplainMatchResult(Vector2DimBoolEq(expected.GetContents()),
                                     arg.GetContents(), result_listener);
}

SparseBoolMatrix BuildMatrix(const std::vector<std::vector<bool>>& contents)
{
  int number_of_rows = contents.size();
  int number_of_columns = contents.empty() ? 0 : contents[0].size();
  SparseBoolMatrix matrix(number_of_rows, number_of_columns);
  for (int i = 0; i < contents.size(); ++i)
  {
    for (int j = 0; j < contents[i].size(); ++j)
    {
      if (contents[i][j])
      {
        matrix.Set(i, j, true);
      }
    }
  }
  return matrix;
}

TEST(SparseBoolMatrixTest, GetNumberOfRowsOK)
{
  SparseBoolMatrix matrix = BuildMatrix({
      {1, 1, 0, 1, 0},  //
      {1, 0, 1, 1, 0},  //
      {0, 0, 0, 0, 1},  //
      {1, 0, 0, 0, 0},  //
      {1, 0, 0, 1, 0},  //
      {0, 0, 0, 0, 0},  //
  });
  EXPECT_EQ(matrix.GetNumberOfRows(), 6);
}

TEST(SparseBoolMatrixTest, GetNumberOfColumnsOK)
{
  SparseBoolMatrix matrix = BuildMatrix({
      {1, 1, 0, 1, 0},  //
      {1, 0, 1, 1, 0},  //
      {0, 0, 0, 0, 1},  //
      {1, 0, 0, 0, 0},  //
      {1, 0, 0, 1, 0},  //
      {0, 0, 0, 0, 0},  //
  });
  EXPECT_EQ(matrix.GetNumberOfColumns(), 5);
}

TEST(SparseBoolMatrixTest, GetOK)
{
  SparseBoolMatrix matrix = BuildMatrix({
      {1, 1, 0},  //
      {1, 0, 1},  //
      {0, 0, 0},  //
      {1, 0, 0},  //
  });
  EXPECT_TRUE(matrix.Get(0, 0));
  EXPECT_TRUE(matrix.Get(0, 1));
  EXPECT_FALSE(matrix.Get(0, 2));

  EXPECT_TRUE(matrix.Get(1, 0));
  EXPECT_FALSE(matrix.Get(1, 1));
  EXPECT_TRUE(matrix.Get(1, 2));

  EXPECT_FALSE(matrix.Get(2, 0));
  EXPECT_FALSE(matrix.Get(2, 1));
  EXPECT_FALSE(matrix.Get(2, 2));

  EXPECT_TRUE(matrix.Get(3, 0));
  EXPECT_FALSE(matrix.Get(3, 1));
  EXPECT_FALSE(matrix.Get(3, 2));
}

TEST(SparseBoolMatrixTest, SetOK)
{
  SparseBoolMatrix matrix = BuildMatrix({
      {1, 1, 0},  //
      {1, 0, 1},  //
      {0, 0, 0},  //
      {1, 0, 0},  //
  });

  matrix.Set(0, 0, false);
  matrix.Set(1, 1, false);  // no change
  matrix.Set(3, 1, true);

  EXPECT_THAT(matrix, MatrixEq(BuildMatrix({
                          {0, 1, 0},  //
                          {1, 0, 1},  //
                          {0, 0, 0},  //
                          {1, 1, 0},  //
                      })));
}

TEST(SparseBoolMatrixTest, GetCopyWithFirstNColumnsGetsDeepCopy)
{
  SparseBoolMatrix matrix = BuildMatrix({
      {1, 1, 0},  //
      {1, 0, 1},  //
      {0, 0, 0},  //
  });
  std::unique_ptr<BoolMatrix> matrix_copy = matrix.GetCopyWithFirstNColumns(3);

  matrix.Set(0, 0, false);

  // Expect the copy to remain unchanged.
  EXPECT_THAT(*matrix_copy, MatrixEq(BuildMatrix({
                                {1, 1, 0},  //
                                {1, 0, 1},  //
                                {0, 0, 0},  //
                            })));
}

TEST(SparseBoolMatrixTest, GetCopyWithFirstNColumnsTrimsSuffixColumns)
{
  SparseBoolMatrix matrix = BuildMatrix({
      {1, 1, 0},  //
      {1, 0, 1},  //
      {0, 0, 0},  //
  });
  std::unique_ptr<BoolMatrix> matrix_copy = matrix.GetCopyWithFirstNColumns(2);

  EXPECT_THAT(*matrix_copy, MatrixEq(BuildMatrix({
                                {1, 1},  //
                                {1, 0},  //
                                {0, 0},  //
                            })));
}

TEST(SparseBoolMatrixTest, IsZeroForZeroMatrixReturnsTrue)
{
  SparseBoolMatrix matrix = BuildMatrix({
      {0, 0, 0},  //
      {0, 0, 0},  //
      {0, 0, 0},  //
  });
  EXPECT_TRUE(matrix.IsZero());
}

TEST(SparseBoolMatrixTest, IsZeroForNonzeroMatrixReturnsFalse)
{
  SparseBoolMatrix matrix = BuildMatrix({
      {1, 0, 0},  //
      {0, 0, 0},  //
      {0, 0, 0},  //
  });
  EXPECT_FALSE(matrix.IsZero());
}

TEST(SparseBoolMatrixTest, GetRowOK)
{
  SparseBoolMatrix matrix = BuildMatrix({
      {1, 0, 1},  //
      {0, 1, 0},  //
      {0, 1, 0},  //
  });
  EXPECT_THAT(matrix.GetRow(0), UnorderedElementsAre(0, 2));
}

TEST(SparseBoolMatrixTest, GetRowForEmptyRowsReturnsEmptyVector)
{
  SparseBoolMatrix matrix = BuildMatrix({
      {1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1},  //
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  //
      {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0},  //
  });
  EXPECT_THAT(matrix.GetRow(1), testing::IsEmpty());
}

TEST(SparseBoolMatrixTest, GetColumnOK)
{
  SparseBoolMatrix matrix = BuildMatrix({
      {1, 0, 1},  //
      {0, 1, 0},  //
      {0, 1, 0},  //
      {0, 1, 0},  //
      {0, 1, 0},  //
      {0, 1, 0},  //
  });
  EXPECT_THAT(matrix.GetColumn(1), UnorderedElementsAre(1, 2, 3, 4, 5));
}

TEST(SparseBoolMatrixTest, SwapRowsOK)
{
  SparseBoolMatrix matrix = BuildMatrix({
      {1, 0, 1, 1, 1},  //
      {0, 1, 0, 1, 0},  //
      {1, 1, 0, 0, 0},  //
      {0, 0, 0, 1, 0},  //
  });
  matrix.SwapRows(0, 1);
  EXPECT_THAT(matrix, MatrixEq(BuildMatrix({
                          {0, 1, 0, 1, 0},  //
                          {1, 0, 1, 1, 1},  //
                          {1, 1, 0, 0, 0},  //
                          {0, 0, 0, 1, 0},  //
                      })));
}

TEST(SparseBoolMatrixTest, SwapColumnsOK)
{
  SparseBoolMatrix matrix = BuildMatrix({
      {0, 1, 0, 1, 0},  //
      {1, 0, 1, 1, 1},  //
      {1, 1, 0, 0, 0},  //
      {0, 0, 0, 1, 0},  //
  });
  matrix.SwapColumns(1, 4);
  EXPECT_THAT(matrix, MatrixEq(BuildMatrix({
                          {0, 0, 0, 1, 1},  //
                          {1, 1, 1, 1, 0},  //
                          {1, 0, 0, 0, 1},  //
                          {0, 0, 0, 1, 0},  //
                      })));
}

TEST(SparseBoolMatrixTest, DeleteColumnOK)
{
  SparseBoolMatrix matrix = BuildMatrix({
      {0, 1, 0, 1, 0},  //
      {1, 0, 1, 1, 1},  //
      {1, 1, 0, 0, 0},  //
      {0, 0, 0, 1, 0},  //
  });
  matrix.DeleteColumn(3);
  EXPECT_THAT(matrix, MatrixEq(BuildMatrix({
                          {0, 1, 0, 0},  //
                          {1, 0, 1, 1},  //
                          {1, 1, 0, 0},  //
                          {0, 0, 0, 0},  //
                      })));
}

TEST(SparseBoolMatrixTest, DeleteLastColumnOK)
{
  SparseBoolMatrix matrix = BuildMatrix({
      {0},  //
      {1},  //
      {1},  //
      {0},  //
  });
  matrix.DeleteColumn(0);
  // It's a design choice: we've decided to allow matrices of sizes like
  // 0x4, 2x0, etc. to be different.
  EXPECT_THAT(matrix, MatrixEq(BuildMatrix({
                          {},  //
                          {},  //
                          {},  //
                          {},  //
                      })));
  EXPECT_EQ(matrix.GetNumberOfRows(), 4);
  EXPECT_EQ(matrix.GetNumberOfColumns(), 0);
}

TEST(SparseBoolMatrixTest, DeleteFirstNColumnsAndEmpyOnesOK)
{
  SparseBoolMatrix matrix = BuildMatrix({
      {0, 1, 0, 0, 0, 1, 0},  //
      {1, 0, 1, 0, 0, 1, 0},  //
      {1, 1, 0, 0, 0, 0, 0},  //
      {0, 0, 0, 0, 0, 1, 0},  //
  });
  matrix.DeleteFirstNColumnsAndEmptyOnes(2);
  EXPECT_THAT(matrix, MatrixEq(BuildMatrix({
                          {0, 1},  //
                          {1, 1},  //
                          {0, 0},  //
                          {0, 1},  //
                      })));
}

TEST(SparseBoolMatrixTest, AppendColumnOK)
{
  SparseBoolMatrix matrix = BuildMatrix({
      {0, 1, 0, 1, 0},  //
      {1, 0, 1, 1, 1},  //
      {1, 1, 0, 0, 0},  //
      {0, 0, 0, 1, 0},  //
  });
  matrix.AppendColumn({0, 1});
  EXPECT_THAT(matrix, MatrixEq(BuildMatrix({
                          {0, 1, 0, 1, 0, 1},  //
                          {1, 0, 1, 1, 1, 1},  //
                          {1, 1, 0, 0, 0, 0},  //
                          {0, 0, 0, 1, 0, 0},  //
                      })));
}

TEST(SparseBoolMatrixTest, AppendColumnToEmptyMatrixOK)
{
  SparseBoolMatrix matrix = BuildMatrix({
      {0},  //
      {1},  //
  });
  matrix.DeleteColumn(0);
  matrix.AppendColumn({0, 1});
  EXPECT_THAT(matrix, MatrixEq(BuildMatrix({
                          {1},  //
                          {1},  //
                      })));
}

TEST(SparseBoolMatrixTest, AddColumnToExisingOneOK)
{
  SparseBoolMatrix matrix = BuildMatrix({
      {0, 1, 0, 1, 0},  //
      {1, 0, 1, 1, 1},  //
      {1, 1, 0, 0, 0},  //
      {0, 0, 0, 1, 0},  //
  });
  matrix.AddColumnToExistingOne({0, 1, 2}, /*to_which_column=*/3);
  EXPECT_THAT(matrix, MatrixEq(BuildMatrix({
                          {0, 1, 0, 0, 0},  //
                          {1, 0, 1, 0, 1},  //
                          {1, 1, 0, 1, 0},  //
                          {0, 0, 0, 1, 0},  //
                      })));
}

TEST(SparseBoolMatrixTest, DeleteAllButOneElementInRowOK)
{
  SparseBoolMatrix matrix = BuildMatrix({
      {0, 1, 0, 1, 0},  //
      {1, 0, 1, 1, 1},  //
      {1, 1, 0, 0, 0},  //
      {0, 0, 0, 1, 0},  //
  });
  matrix.DeleteAllButOneElementInRow(1, /*column_to_keep=*/std::nullopt);
  EXPECT_THAT(matrix, MatrixEq(BuildMatrix({
                          {0, 1, 0, 1, 0},  //
                          {1, 0, 0, 0, 0},  //
                          {1, 1, 0, 0, 0},  //
                          {0, 0, 0, 1, 0},  //
                      })));
}

TEST(SparseBoolMatrixTest, DeleteAllButOneElementInRowKeepingNonfirstColumnOK)
{
  auto matrix = BuildMatrix({
      {0, 1, 0, 1, 0},  //
      {1, 0, 1, 1, 1},  //
      {1, 1, 0, 0, 0},  //
      {0, 0, 0, 1, 0},  //
  });
  matrix.DeleteAllButOneElementInRow(1, /*column_to_keep=*/2);
  EXPECT_THAT(matrix, MatrixEq(BuildMatrix({
                          {0, 1, 0, 1, 0},  //
                          {0, 0, 1, 0, 0},  //
                          {1, 1, 0, 0, 0},  //
                          {0, 0, 0, 1, 0},  //
                      })));
}

TEST(SparseBoolVectorTest, GetSmallestIndexOK)
{
  SparseBoolVector vector = {10, 5, 13, 1};
  EXPECT_EQ(GetSmallestIndex(vector), 1);

  vector = {2, 3, 4};
  EXPECT_EQ(GetSmallestIndex(vector), 2);

  vector = {1024, 555, 100};
  EXPECT_EQ(GetSmallestIndex(vector), 100);
}

TEST(SparseBoolVectorTest,
     GetSmallestIndexForVectorOfSizeOneReturnsCorrectIndex)
{
  SparseBoolVector vector = {42};
  EXPECT_EQ(GetSmallestIndex(vector), 42);
}

}  // namespace