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

#include "sparse_bool_matrix_utils.h"

#include <iostream>
#include <limits>

#include "absl/log/check.h"
#include "sparse_bool_matrix.h"

void TerminateIfMatrixDimensionTooLarge(int dimension)
{
  CHECK(dimension < std::numeric_limits<Index>::max())
      << "The specified dimension (" << dimension
      << ") exceeds the index limit (" << std::numeric_limits<Index>::max() - 1
      << ")" << std::endl;
}