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

#include "algebra.h"

#include <optional>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "test_util.h"

namespace
{

// A custom matcher used to compare FMatrix objects.
MATCHER_P(FMatrixEq, expected, "")
{
  static_assert(
      std::is_convertible_v<std::remove_cvref_t<arg_type>*, FMatrix*>);
  static_assert(
      std::is_convertible_v<std::remove_cvref_t<expected_type>*, FMatrix*>);
  return testing::ExplainMatchResult(Vector2DimBoolEq(expected.GetContents()),
                                     arg.GetContents(), result_listener);
}

FMatrix BuildMatrix(const std::vector<std::vector<bool>>& contents)
{
  assert(!contents.empty() && !contents[0].empty());
  FMatrix result(contents.size(), contents[0].size());
  for (int row_index = 0; row_index < contents.size(); ++row_index)
  {
    for (int column_index = 0; column_index < contents[0].size();
         ++column_index)
    {
      result.set(row_index, column_index,
                 field(contents[row_index][column_index]));
    }
  }
  return result;
}

TEST(SimplifyAgainstMatrixTest, SimplifyAgainstMatrixOK)
{
  FMatrix matrix = BuildMatrix({
      {1, 1, 0, 1},  //
      {1, 0, 1, 0},  //
      {0, 0, 0, 1},  //
  });
  FMatrix matrix_to_simplify_against = BuildMatrix({
      {1, 0, 0, 0},  //
      {0, 0, 1, 0},  //
      {0, 0, 0, 0},  //
  });

  simplify_against_matrix(matrix, matrix_to_simplify_against);

  EXPECT_THAT(matrix, FMatrixEq(BuildMatrix({
                          {0, 0, 0, 0},  //
                          {0, 0, 0, 0},  //
                          {0, 0, 0, 1},  //
                      })));
}

struct SNFTestParam
{
  FMatrix input_matrix;
  FMatrix expected_diagonal_matrix;
};

class SNFTest : public testing::TestWithParam<SNFTestParam>
{
};

TEST_P(SNFTest, MatrixRecomputeOK)
{
  const FMatrix& input_matrix = GetParam().input_matrix;
  const FMatrix& expected_diagonal_matrix = GetParam().expected_diagonal_matrix;

  SNF snf;
  std::optional<SNF::Composition> composition =
      snf.MatrixRecompute(input_matrix, /*return_composition=*/true);
  ASSERT_TRUE(composition.has_value());

  EXPECT_THAT(composition->diagonal, FMatrixEq(expected_diagonal_matrix));

  EXPECT_THAT(composition->left * composition->diagonal,
              FMatrixEq(input_matrix * composition->right_inverse));
}

INSTANTIATE_TEST_SUITE_P(
    TestMatrices, SNFTest,
    testing::Values(SNFTestParam{.input_matrix = BuildMatrix({
                                     {0, 0, 1},  //
                                     {0, 0, 0},  //
                                     {1, 0, 0},  //
                                 }),
                                 .expected_diagonal_matrix = BuildMatrix({
                                     {1, 0, 0},  //
                                     {0, 1, 0},  //
                                     {0, 0, 0},  //
                                 })},
                    SNFTestParam{.input_matrix = BuildMatrix({
                                     {1, 1, 0},  //
                                     {1, 0, 1},  //
                                     {0, 0, 0},  //
                                     {1, 1, 0},  //
                                 }),
                                 .expected_diagonal_matrix = BuildMatrix({
                                     {1, 0, 0},  //
                                     {0, 1, 0},  //
                                     {0, 0, 0},  //
                                     {0, 0, 0},  //
                                 })},
                    SNFTestParam{.input_matrix = BuildMatrix({
                                     {1, 1, 0, 1, 1, 0, 1, 1, 0, 1},  //
                                     {1, 0, 1, 1, 0, 1, 1, 0, 1, 0},  //
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},  //
                                     {1, 1, 0, 1, 1, 0, 1, 1, 0, 1},  //
                                     {1, 0, 1, 1, 0, 1, 1, 0, 1, 0},  //
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},  //
                                     {1, 1, 0, 1, 1, 0, 1, 1, 0, 1},  //
                                     {1, 0, 1, 1, 0, 1, 1, 0, 1, 0},  //
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},  //
                                 }),
                                 .expected_diagonal_matrix = BuildMatrix({
                                     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0},  //
                                     {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},  //
                                     {0, 0, 1, 0, 0, 0, 0, 0, 0, 0},  //
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  //
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  //
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  //
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  //
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  //
                                     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  //
                                 })}));

TEST(ComputeHomologyTest, ComputeHomologyOK)
{
  SNF snf1(BuildMatrix({
      {0},  //
      {0},  //
      {0},  //
      {0},  //
      {0},  //
      {0},  //
      {1},  //
      {1},  //
      {0},  //
      {0},  //
      {0},  //
      {0},  //
      {0},  //
      {0},  //
  }));
  SNF snf2(BuildMatrix({
      {0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  //
      {1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  //
      {0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},  //
      {0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0},  //
      {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0},  //
      {0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0},  //
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0},  //
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1},  //
      {1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1},  //
  }));
  EXPECT_THAT(ComputeHomology(snf1, snf2), FMatrixEq(BuildMatrix({
                                               {1, 0, 0, 0, 0},  //
                                               {0, 0, 1, 0, 0},  //
                                               {0, 0, 1, 0, 0},  //
                                               {1, 0, 1, 0, 0},  //
                                               {0, 0, 0, 1, 0},  //
                                               {0, 0, 0, 1, 0},  //
                                               {0, 0, 0, 0, 0},  //
                                               {1, 0, 1, 1, 0},  //
                                               {0, 0, 0, 0, 1},  //
                                               {0, 0, 0, 0, 1},  //
                                               {1, 0, 1, 1, 1},  //
                                               {0, 1, 0, 0, 0},  //
                                               {0, 1, 0, 0, 0},  //
                                               {1, 1, 1, 1, 1},  //
                                           })));
}

}  // namespace