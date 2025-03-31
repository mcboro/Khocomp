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

#include "algebra.h"

#include <cassert>
#include <iostream>
#include <memory>
#include <optional>
#include <ostream>
#include <thread>
#include <utility>
#include <vector>

#include "sparse_bool_matrix.h"
#include "sparse_bool_matrix_fast_column_operations.h"
#include "sparse_low_memory_bool_matrix.h"

////////////////////////////////////////////////////////////////////////////////
// field                                                                      //
////////////////////////////////////////////////////////////////////////////////
field::field(bool value) : val(value) {}

field field::operator-() const { return field(not val); }

field field::operator+(const field& f) const { return field(val xor f.val); }

field field::operator-(const field& f) const { return field(val xor f.val); }

field field::operator*(const field& f) const { return field(val and f.val); }

field::operator bool() const { return val; }

std::ostream& operator<<(std::ostream& os, const field& f)
{
  os << static_cast<int>(f.val);
  return os;
}

////////////////////////////////////////////////////////////////////////////////
// FMatrix                                                                    //
////////////////////////////////////////////////////////////////////////////////
FMatrix::FMatrix() : FMatrix(/*rows=*/0, /*cols=*/0) {}

FMatrix::FMatrix(long int rows, long int cols)
    : matrix_(std::make_unique<SparseBoolMatrix>(rows, cols))
{
}

FMatrix::FMatrix(const FMatrix& other) : matrix_(other.matrix()->Clone()) {}

FMatrix& FMatrix::operator=(const FMatrix& other)
{
  matrix_ = other.matrix()->Clone();
  return *this;
}

long int FMatrix::rows() const { return matrix_->GetNumberOfRows(); }

long int FMatrix::cols() const { return matrix_->GetNumberOfColumns(); }

bool FMatrix::isZero() const { return matrix_->IsZero(); }

field FMatrix::get(long int i, long int j) const
{
  assert(0 <= i && i < rows() && 0 <= j && j < cols());
  return field(matrix_->Get(i, j));
}

void FMatrix::set(long int i, long int j, field f)
{
  assert(0 <= i && i < rows() && 0 <= j && j < cols());
  matrix_->Set(i, j, static_cast<bool>(f));
}

void FMatrix::sum_rows(long int which_row, long int to_which)
{
  matrix_->SumRows(which_row, to_which);
}

void FMatrix::sum_cols(long int which_col, long int to_which)
{
  matrix_->SumColumns(which_col, to_which);
}

void FMatrix::switch_rows(long int which_row, long int with_which)
{
  matrix_->SwapRows(which_row, with_which);
}

void FMatrix::switch_cols(long int which_col, long int with_which)
{
  matrix_->SwapColumns(which_col, with_which);
}

FMatrix FMatrix::take_column(long int which) const
{
  MatrixContentsAsColumns contents = {
      .number_of_rows = matrix_->GetNumberOfRows(),
      .columns = {matrix_->GetColumn(which)}};
  return FMatrix(std::make_unique<SparseBoolMatrix>(std::move(contents)));
}

void FMatrix::kill_column(long int which)
{
  assert(0 <= which && which < cols());
  matrix_->DeleteColumn(which);
}

FMatrix FMatrix::Concatenate(const FMatrix& Another) const
{
  assert(rows() == Another.rows());
  FMatrix M(matrix_->Clone());
  for (int i = 0; i < Another.cols(); ++i)
  {
    M.matrix_->AppendColumn(Another.matrix_->GetColumn(i));
  }
  return M;
}

FMatrix FMatrix::operator*(const FMatrix& f) const
{
  assert(cols() == f.rows());
  std::vector<SparseBoolVector> result_rows(rows());
  for (int i = 0; i < rows(); ++i)
  {
    SparseBoolVector& new_row = result_rows[i];
    for (int j = 0; j < f.cols(); ++j)
    {
      // Here, we calculate the (i, j) element of the resulting matrix.
      bool value = false;
      const SparseBoolVector& row = matrix_->GetRow(i);
      const SparseBoolVector& column = f.matrix_->GetColumn(j);
      for (int idx : row)
      {
        if (column.contains(idx))
        {
          value = !value;
        }
      }
      if (value)
      {
        new_row.insert(j);
      }
    }
  }
  MatrixContentsAsRows contents{.number_of_columns = static_cast<int>(f.cols()),
                                .rows = std::move(result_rows)};
  return FMatrix(std::make_unique<SparseBoolMatrix>(std::move(contents)));
}

void FMatrix::ReduceMemoryFootprint()
{
  matrix_ = std::move(*matrix_).MoveToLowMemoryInstance();
}

void FMatrix::ConvertToCapableMatrix()
{
  matrix_ = std::move(*matrix_).MoveToCapableInstance();
}

BoolMatrix* FMatrix::matrix() { return matrix_.get(); }
const BoolMatrix* FMatrix::matrix() const { return matrix_.get(); }

void PrintTo(const FMatrix&, std::ostream* os)
{
  *os << "see the matrix below";
}

std::vector<std::vector<bool>> FMatrix::GetContents() const
{
  std::vector<bool> empty_row(cols(), false);
  std::vector<std::vector<bool>> result(rows(), empty_row);
  for (int i = 0; i < rows(); ++i)
  {
    for (int j = 0; j < cols(); ++j)
    {
      if (get(i, j))
      {
        result[i][j] = true;
      }
    }
  }
  return result;
}

////////////////////////////////////////////////////////////////////////////////
// The rest                                                                   //
////////////////////////////////////////////////////////////////////////////////
std::vector<long int> column_simplify(FMatrix& C)
{
  std::vector<long int> signature_columns;
  signature_columns.reserve(C.cols());

  for (int column_index = 0; column_index < C.cols(); ++column_index)
  {
    const SparseBoolVector& column = C.matrix()->GetColumn(column_index);
    if (column.empty())
    {
      signature_columns.push_back(-1);
      continue;
    }
    int row_index = GetSmallestIndex(column);
    signature_columns.push_back(row_index);
    // This copy is intentional since we want to avoid modifying the
    // collection we're iterating over. The copy could be avoided, but
    // it's not worth it.
    SparseBoolVector row = C.matrix()->GetRow(row_index);
    for (int column_to_simplify : row)
    {
      if (column_to_simplify > column_index)
      {
        C.sum_cols(column_index, column_to_simplify);
      }
    }
  }
  return signature_columns;
}

std::optional<SNF::Composition> SNF::MatrixRecompute(FMatrix initial_matrix,
                                                     bool return_composition)
{
  initial_matrix.ConvertToCapableMatrix();
  isInitializedvalue = true;
  cl = initial_matrix.cols();
  rw = initial_matrix.rows();
  SparseBoolMatrixFastColumnOperations left_fast =
      SparseBoolMatrixFastColumnOperations::CreateIdentityMatrix(rw);
  SparseBoolMatrixFastColumnOperations right_inverse_fast =
      SparseBoolMatrixFastColumnOperations::CreateIdentityMatrix(cl);
  BoolMatrix* diagonal = initial_matrix.matrix();
  std::vector<int> mycolumns;
  mycolumns.reserve(cl);
  long int zerocolumns = 0;

  // first step: making matrix having 1 appearing 1 at most once in each column
  // and at most once in each row.
  for (long int pivotcolumn = 0; pivotcolumn < cl; pivotcolumn++)
  {
    // The copy is important here. We could avoid it, but we would have to be
    // careful, and it does not seem worth it.
    SparseBoolVector column = diagonal->GetColumn(pivotcolumn);
    if (column.empty())
    {
      mycolumns.push_back(-1);
      ++zerocolumns;
      continue;
    }

    int pivotrow = GetSmallestIndex(column);
    mycolumns.push_back(pivotrow);

    for (int row_index : column)
    {
      if (row_index == pivotrow)
      {
        continue;
      }

      diagonal->SumRows(pivotrow, row_index);
      left_fast.SumColumns(row_index, pivotrow);
    }

    const SparseBoolVector& row = diagonal->GetRow(pivotrow);
    for (int column_index : row)
    {
      if (column_index == pivotcolumn)
      {
        continue;
      }
      right_inverse_fast.SumColumns(pivotcolumn, column_index);
    }
    diagonal->DeleteAllButOneElementInRow(pivotrow, pivotcolumn);
  }

  if (!return_composition)
  {
    initial_matrix = FMatrix(0, 0);
  }

  rank = cl - zerocolumns;

  // int solution_1_swaps = 0;
  // now we start permuting. If there are TotalCols zero columns, then the
  // matrix is the zero matrix.
  if (zerocolumns < cl)
  {
    for (long int pivotcolumn = 0; pivotcolumn < cl; pivotcolumn++)
    {
      if (mycolumns[pivotcolumn] == -1)
      {
        // we replace it by the last column
        int currentcolumn = cl - 1;
        while (mycolumns[currentcolumn] == -1 and currentcolumn > pivotcolumn)
        {
          currentcolumn--;
        }
        if (currentcolumn > pivotcolumn)
        {
          mycolumns[pivotcolumn] = mycolumns[currentcolumn];
          mycolumns[currentcolumn] = -1;
          if (return_composition)
          {
            diagonal->SwapColumns(pivotcolumn, currentcolumn);
          }
          right_inverse_fast.SwapColumns(pivotcolumn, currentcolumn);
        }
      }

      if (mycolumns[pivotcolumn] >= 0)
      {
        // we want to get 1-s only on the diagonal
        // so if mycolumns[pivotcolumn] is not equal to pivotcolumn
        // we switch rows: pivotcolumn with currentcolumn
        // after that, the for loop updates the mycolumns vector
        if (mycolumns[pivotcolumn] != pivotcolumn)
        {
          int currentrow = mycolumns[pivotcolumn];
          if (return_composition)
          {
            diagonal->SwapRows(pivotcolumn, currentrow);
          }
          // Diagonal.switch_rows(pivotcolumn, currentrow);
          // if (isAllComputed)
          // {
          //   left_inverse_fast->SwapRows(pivotcolumn, currentrow);
          // }
          left_fast.SwapColumns(pivotcolumn, currentrow);
          mycolumns[pivotcolumn] = pivotcolumn;
          for (long int i = pivotcolumn + 1; i < cl; i++)
            if (mycolumns[i] == pivotcolumn) mycolumns[i] = currentrow;
        }
      }
    }
  }

  // compute kernel:
  FMatrix RightInverse = FMatrix(std::make_unique<SparseLowMemoryBoolMatrix>(
      std::move(right_inverse_fast).ExtractContents()));
  assert((RightInverse.rows() == cl) && (cl >= rank) &&
         (RightInverse.cols() == cl));
  KernelMatrix = FMatrix(
      std::make_unique<SparseLowMemoryBoolMatrix>(MatrixContentsAsColumns{
          .number_of_rows = static_cast<int>(cl),
          .columns = std::vector<SparseBoolVector>(cl - rank)}));
  for (int column_index = rank; column_index < cl; ++column_index)
  {
    const SparseBoolVector& column =
        RightInverse.matrix()->GetColumn(column_index);
    for (int row_index : column)
    {
      KernelMatrix.set(row_index, cl - column_index - 1, field(true));
    }
  }
  if (!return_composition)
  {
    RightInverse = FMatrix(0, 0);
  }

  // compute image:
  FMatrix Left = FMatrix(0, 0);
  if (return_composition)
  {
    SparseBoolMatrixFastColumnOperations left_copy = left_fast;
    Left = FMatrix(std::make_unique<SparseBoolMatrix>(
        std::move(left_copy).ExtractContents()));
  }
  MatrixContentsAsColumns left = std::move(left_fast).ExtractContents();
  assert(left.number_of_rows == rw && left.columns.size() >= rank);
  left.columns.erase(left.columns.begin() + rank, left.columns.end());
  ImageMatrix = FMatrix(std::make_unique<SparseBoolMatrix>(std::move(left)));
  signature_columns = column_simplify(ImageMatrix);
  ImageMatrix.ReduceMemoryFootprint();

  if (return_composition)
  {
    return Composition{.left = std::move(Left),
                       .diagonal = std::move(initial_matrix),
                       .right_inverse = std::move(RightInverse)};
  }
  return std::nullopt;
}

void simplify_against_matrix(FMatrix& C, const FMatrix& D)
{
  for (int column_index = 0; column_index < D.cols(); ++column_index)
  {
    const SparseBoolVector& column = D.matrix()->GetColumn(column_index);
    if (column.empty())
    {
      continue;
    }
    int row_index = GetSmallestIndex(column);
    // This copy is intentional since we want to avoid modifying the
    // collection we're iterating over. The copy could be avoided, but
    // it's not worth it.
    SparseBoolVector row = C.matrix()->GetRow(row_index);
    for (int column_to_simplify : row)
    {
      C.matrix()->AddColumnToExistingOne(column, column_to_simplify);
    }
  }
}

std::vector<long int> simplify_and_trace(FMatrix& C, const FMatrix& D)
{
  // the procedure is very similar to the one above
  // with the exception that we simplify only the first column of C
  // and we record all the classes that were used in the computations
  long int rownumber;
  std::vector<long int> columns;
  for (long int column_of_D = 0; column_of_D < D.cols(); column_of_D++)
  {
    rownumber = 0;
    while (rownumber < D.rows() and
           D.get(rownumber, column_of_D) == field(false))
      rownumber++;
    if (rownumber < D.rows())
    {
      if (C.get(rownumber, 0) == field(true))
      {
        columns.push_back(column_of_D);
        for (long int rowpivot = rownumber; rowpivot < C.rows(); rowpivot++)
        {
          C.set(rowpivot, 0, C.get(rowpivot, 0) + D.get(rowpivot, column_of_D));
        }
      }
    }
  }
  return columns;
}

void ChainComplex::ComputeAllSNF()
{
  std::vector<std::thread> mythreads;
  for (long int i = 0; i < ChainMaps.size(); i++)
  {
    if ((!onlyone) or (i = whichone) or (i == whichone - 1))
    {
      mythreads.emplace_back([this, i]()
                             { AllMaps[i].MatrixRecompute(ChainMaps[i]); });
    }
  }
  for (auto& th : mythreads) th.join();
}
void ChainComplex::ComputeHomologyGroup()
{
  // TODO : the chain complex is assumed to be "connected". no zero terms in
  // between.
  // TODO : need to implement the procedure.
  // TODO : maybe add dimensions as a vector to the constructor.
  FMatrix StupidMatrix(0, 0);
  if (AllMaps.size() == 0) return;
  if (AllMaps[0].isInitialized())
    Homologies.push_back(ComputeStartingHomology(AllMaps[0]));
  else
    Homologies.push_back(StupidMatrix);
  for (long int i = 1; i < AllMaps.size(); i++)
  {
    if (AllMaps[i - 1].isInitialized() and AllMaps[i].isInitialized())
    {
      /* Commented out a part.
      if ((ChainMaps[i].rows()==0 or ChainMaps[i].cols()==0) and
      (ChainMaps[i-1].rows()==0 or ChainMaps[i-1].cols()==0))
      Homologies.push_back(StupidMatrix);
      else
      {
          if (ChainMaps[i].rows()==0 or ChainMaps[i].cols()==0)
      Homologies.push_back(ComputeFinalHomology(AllMaps[i-1])); else
          {
              if (ChainMaps[i-1].rows()==0 or ChainMaps[i-1].cols()==0)
      Homologies.push_back(ComputeStartingHomology(AllMaps[i])); else {
                  //std::cout << "Source matrix " << ChainMaps[i-1].rows() <<
                  "
      rows and " << ChainMaps[i-1].cols() << " columns " << std::endl;
                  //std::cout << "Target matrix " << ChainMaps[i].rows() << "
      rows and " << ChainMaps[i].cols() << " columns " << std::endl;*/
      Homologies.push_back(ComputeHomology(AllMaps[i - 1], AllMaps[i]));
    }
    else
      Homologies.push_back(StupidMatrix);
  }
  if (AllMaps.back().isInitialized())
    Homologies.push_back(ComputeFinalHomology(AllMaps.back()));
  else
    Homologies.push_back(StupidMatrix);
}
void ChainComplex::from_empty(long i)
{
  // if there are no maps, we create homology at some level
  ChainMaps.clear();
  AllMaps.clear();
  Homologies.clear();
  FMatrix F(i, i);
  for (long int j = 0; j < i; j++) F.set(j, j, field(true));
  Homologies.push_back(F);
}

FMatrix ComputeHomology(const SNF& D2, const SNF& D1)
{
  const FMatrix& Im = D2.Image();
  const FMatrix& Ker = D1.Kernel();
  assert(Im.rows() == Ker.rows());

  long int dimker = Ker.cols();
  long int dimim = Im.cols();
  long int rw = Ker.rows();
  if (dimker == dimim)
    return FMatrix(rw, 0);  // no need to compute if it will be zero anyway.

  FMatrix All = Im.Concatenate(Ker);
  All.ConvertToCapableMatrix();
  column_simplify(All);
  All.ReduceMemoryFootprint();
  All.matrix()->DeleteFirstNColumnsAndEmptyOnes(dimim);
  return All;
}

FMatrix ComputeStartingHomology(const SNF& D)
{
  // returns the homology of the first entry in the chain complex
  FMatrix M = D.Kernel();
  M.ConvertToCapableMatrix();
  column_simplify(M);
  M.ReduceMemoryFootprint();
  return M;
}

FMatrix ComputeFinalHomology(const SNF& D2)
{
  // return the homology of the last entry in the chain complex
  FMatrix A1(D2.rows(), D2.rows());  // Rows not cols!
  SNF D1(A1);
  FMatrix H = ComputeHomology(D2, D1);
  return H;
}

bool has_zero_column(FMatrix T)
{
  // the name is misleading. We operate on matrix to make it `upper triangular'
  // and check for zero columns
  // TODO: there are some similar procedures in the code. Simplify to use just
  // one!
  long int rw = T.rows();
  long int cl = T.cols();
  long int startrw;
  for (long int col = 0; col < cl; col++)
  {
    startrw = 0;
    while ((startrw < rw) and (T.get(startrw, col) == false)) startrw++;
    if (startrw == rw)
      return true;  // here is a nice feature. If the column is zero, we just
                    // declare that A and B are linearly dependent mod the
                    // image.
    for (long int curcol = col + 1; curcol < cl; curcol++)
      if (T.get(startrw, curcol) == true) T.sum_cols(col, curcol);
  }
  // if the procedure works, we should have already found a zero column. So if
  // we haven't, the vectors are linearly independent.
  return false;
}

bool are_equal_mod_image(const FMatrix& A, const FMatrix& B, const FMatrix& C)
{
  //
  // we assume that A and B are vectors and C is a matrix
  // we check, whether A=B in the quotient space modulo C.
  // To this end, we create a matrix A|B|C
  // and after column operations we try to make kill as many columns as
  // possible.
  FMatrix T = A.Concatenate(B);
  T = T.Concatenate(C);
  return has_zero_column(T);
}

FMatrix intersection_of_spaces(const FMatrix& A, const FMatrix& B)
{
  // the procedure take subspaces generated by A and B.
  // and computes the intersection.
  // we assume that A and B have no kernel
  FMatrix C = A.Concatenate(B);
  // std::cout << "matrix C = " << C << std::endl;
  SNF S(C);
  FMatrix K = S.Kernel();
  long int newcols = K.cols();
  // std::cout << "newcols = " << newcols << std::endl;
  long int newrows = C.rows();
  FMatrix F(newrows, newcols);
  for (long int j = 0; j < newcols; j++)
  {
    for (int i = 0; i < A.cols(); i++)
    {
      if (K.get(i, j) == field(true))
      {
        for (long int k = 0; k < A.rows(); k++)
          F.set(k, j, F.get(k, j) + A.get(k, i));
      }
    }
  }
  strip_zero_columns(F);
  return F;
}

std::vector<long int> trace_homology_element(const FMatrix& A, const FMatrix& B,
                                             const FMatrix& C)
{
  // Takes a vector A. Checks if there exists a linear combination of elements
  // of B such that A is equal to this combination modulo C. if for none element
  // we have equality, we return -1
  if (C.cols() == 0) return trace_homology_element(A, B);
  FMatrix D = A;
  simplify_against_matrix(D, C);  // now we have reduced D mod C
  std::vector<long int> trace_vector = simplify_and_trace(D, B);
  bool is_non_zero = false;
  for (long int i = 0; i < D.rows(); i++)
    if (D.get(i, 0) == field(true)) is_non_zero = true;
  if (is_non_zero)
  {
    trace_vector.clear();
    trace_vector.push_back(-4);
  }
  return trace_vector;
}
std::vector<long int> trace_homology_element(const FMatrix& A, const FMatrix& B)
{
  FMatrix D = A;
  std::vector<long int> trace_vector = simplify_and_trace(D, B);
  bool is_non_zero = false;
  for (long int i = 0; i < D.rows(); i++)
    if (D.get(i, 0) == field(true)) is_non_zero = true;
  if (is_non_zero)
  {
    trace_vector.clear();
    trace_vector.push_back(-1);
  }
  return trace_vector;
  // this is an overloaded version of the previous procedure
  // in case C=0
}

bool are_linearly_dependent(const FMatrix& A, const FMatrix& B)
{
  // checks if the two vectors are dependent.
  //
  if (A.isZero()) return true;
  if (B.isZero()) return true;
  for (long int i = 0; i < A.rows(); i++)
    if (A.get(i, 0) != B.get(i, 0)) return false;
  return true;
}
long int strip_zero_columns(FMatrix& A)
{
  if (A.cols() == 0) return 0;
  long int number_of_zero_columns = 0;
  long int column_pivot = 0;
  bool isZero;
  while (column_pivot < A.cols())
  {
    isZero = true;
    for (long int row = 0; row < A.rows(); row++)
      if (A.get(row, column_pivot) == true) isZero = false;
    if (isZero)
    {
      A.kill_column(column_pivot);
      number_of_zero_columns++;
    }
    else
      column_pivot++;
  }
  return number_of_zero_columns;
}