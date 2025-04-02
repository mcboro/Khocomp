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
#include <iostream>

#include "khovanov.h"
#include "pdcode.h"

#define ANSI_COLOR_RED "\x1b[31m"
#define ANSI_COLOR_GREEN "\x1b[32m"
#define ANSI_COLOR_YELLOW "\x1b[33m"
#define ANSI_COLOR_BLUE "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN "\x1b[36m"
#define ANSI_COLOR_RESET "\x1b[0m"

int main(int argc, char** argv)
{
  pdcode mycode;
  bool aconly = false;
  long int acwhich = 0;
  std::string st = "testfile";
  std::string at = "actionfile";
  if (argc > 1) st = argv[1];
  if (argc > 2) at = argv[2];
  if (argc > 3)
  {
    aconly = true;
    acwhich = atol(argv[3]);
  }
  mycode.read_from_file(st);
  mycode.compute_resolutions();
  mycode.compute_arrows();
  khovanov kh(mycode, aconly, acwhich);
  std::cout << "Loading " << std::endl;
  kh.load_action_from_file(at);
  kh.show_all_homologies();
}
