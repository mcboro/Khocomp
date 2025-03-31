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

#include <set>
#include <type_traits>
#include <utility>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

// A matcher for std::vector<std::vector<bool>>.
// It prints in a readable way the differences in case of a mismatch.
MATCHER_P(Vector2DimBoolEq, expected, testing::PrintToString(expected))
{
  testing::StaticAssertTypeEq<std::vector<std::vector<bool>>,
                              std::remove_cvref_t<arg_type>>();
  testing::StaticAssertTypeEq<std::vector<std::vector<bool>>,
                              std::remove_cvref_t<expected_type>>();
  if (arg.size() != expected.size())
  {
    *result_listener << "\narg.size() != expected.size() (" << arg.size()
                     << " vs " << expected.size() << ")";
    return false;
  }

  std::set<std::pair<int, int>> errors;
  for (int i = 0; i < arg.size(); ++i)
  {
    if (arg[i].size() != expected[i].size())
    {
      *result_listener << "\narg[" << i << "].size() != expected[" << i
                       << "].size() (" << arg[i].size() << " vs "
                       << expected[i].size() << ")";
      return false;
    }
    for (int j = 0; j < arg[i].size(); ++j)
    {
      if (arg[i][j] != expected[i][j])
      {
        errors.insert({i, j});
      }
    }
  }
  if (errors.empty())
  {
    return true;
  }

  // Print the whole matrices. E.g.:
  //
  //              actual                            expected
  //  (  1  0! 0  0! 0! 0  0! 0! 0! )    (  1  1! 0  1! 1! 0  1! 1! 1! )
  //  (  1  1! 0! 0! 0  0! 0! 0  0  )    (  1  0! 1! 1! 0  1! 1! 0  0  )
  //  (  0  0  1! 0  0  0  0  0  0! )    (  0  0  0! 0  0  0  0  0  1! )
  //  (  1  0! 0  1  0! 0  0! 0! 0! )    (  1  1! 0  1  1! 0  1! 1! 1! )
  //  (  1  1! 0! 0! 1! 0! 0! 0  0  ) vs (  1  0! 1! 1! 0! 1! 1! 0  0  )
  //  (  0  0  1! 0  0  1! 0  0  0! )    (  0  0  0! 0  0  0! 0  0  1! )
  //  (  1  0! 0  0! 0! 0  1  0! 0! )    (  1  1! 0  1! 1! 0  1  1! 1! )
  //  (  1  1  0  0! 0! 0  0! 1  0! )    (  1  1  0  1! 1! 0  1! 1  1! )
  //  (  0  0  1! 0  0  0  0  0  1  )    (  0  0  0! 0  0  0  0  0  1  )
  //
  *result_listener << "\n";
  std::string first_matrix_label = "actual";
  std::string second_matrix_label = "expected";
  // Every element occupies 3 characters and each border occupies 2 characters.
  int matrix_width = 3 * arg[0].size() + 2 * 2;
  int separation_between_matrices_width = 4;
  int first_label_left_padding = (matrix_width - first_matrix_label.size()) / 2;
  int second_label_left_padding =
      (matrix_width - second_matrix_label.size()) / 2;
  for (int i = 0; i < first_label_left_padding; ++i)
  {
    *result_listener << " ";
  }
  *result_listener << first_matrix_label;
  int gap =
      (matrix_width - first_label_left_padding - first_matrix_label.size() +
       separation_between_matrices_width + second_label_left_padding);
  for (int i = 0; i < gap; ++i)
  {
    *result_listener << " ";
  }
  *result_listener << second_matrix_label << "\n";

  for (int i = 0; i < arg.size(); ++i)
  {
    *result_listener << "( ";
    for (int j = 0; j < arg[i].size(); ++j)
    {
      *result_listener << " " << arg[i][j];
      if (errors.contains({i, j}))
      {
        *result_listener << "!";
      }
      else
      {
        *result_listener << " ";
      }
    }
    *result_listener << " ) ";
    // Print "vs" somewhere in the middle.
    if (i == (arg.size() - 1) / 2)
    {
      *result_listener << "vs";
    }
    else
    {
      *result_listener << "  ";
    }
    *result_listener << " ( ";
    for (int j = 0; j < expected[i].size(); ++j)
    {
      *result_listener << " " << expected[i][j];
      if (errors.contains({i, j}))
      {
        *result_listener << "!";
      }
      else
      {
        *result_listener << " ";
      }
    }
    *result_listener << " )";
    if (i < arg.size() - 1)
    {
      *result_listener << "\n";
    }
  }
  return false;
}
