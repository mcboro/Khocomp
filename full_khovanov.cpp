// Copyright 2025 Maciej Borodzik
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

#include <cstring>

#include "khovanov.h"
#include "pdcode.h"

int main(int argc, char** argv)
{
  pdcode mycode;
  std::string st = "testfile";
  if (argc > 1) st = argv[1];
  long int quant, hom;
  if (argc > 2)
    hom = atol(argv[2]);
  else
    hom = 0;
  if (argc > 3)
    quant = atol(argv[3]);
  else
    quant = 0;
  mycode.read_from_file(st);
  mycode.compute_resolutions();
  mycode.compute_arrows();
  // khovanov kh(mycode,true,hom);
  khovanov kh(mycode);
  kh.detailed_homology(quant, hom);
}
