#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() : rows_(1), cols_(1) { MemoryAllocate(); }

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows_ <= 0 || cols_ <= 0) {
    throw std::invalid_argument(
        "The number of columns or columns is less than 1");
  }
  MemoryAllocate();
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : rows_(other.rows_), cols_(other.cols_) {
  MemoryAllocate();
  CopyData(other);
}

S21Matrix::S21Matrix(S21Matrix&& other) noexcept
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.matrix_ = nullptr;
  other.cols_ = 0;
  other.rows_ = 0;
}

S21Matrix::~S21Matrix() {
  if (matrix_) {
    MemoryFree();
  }
}

void S21Matrix::MemoryAllocate() {
  matrix_ = new double*[rows_];
  matrix_[0] = new double[rows_ * cols_];
  std::memset(matrix_[0], 0, rows_ * cols_ * sizeof(double));
  for (int i = 1; i < rows_; ++i) {
    matrix_[i] = matrix_[i - 1] + cols_;
  }
}

void S21Matrix::MemoryFree() {
  delete[] matrix_[0];
  delete[] matrix_;
  rows_ = 0;
  cols_ = 0;
}

void S21Matrix::CopyData(S21Matrix const& other) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (cols_ != other.cols_ || rows_ != other.rows_) {
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (cols_ != other.cols_ || rows_ != other.rows_) {
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

bool S21Matrix::EqMatrix(const S21Matrix& other) const {
  double eps = 1e-7;
  bool result = (other.rows_ == rows_ && other.cols_ == cols_) ? true : false;
  for (int i = 0; i < rows_ && result; i++) {
    for (int j = 0; j < cols_; j++) {
      if (fabs(matrix_[i][j] - other.matrix_[i][j]) > eps) {
        result = false;
      }
    }
  }
  return result;
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.rows_) {
    throw std::out_of_range(
        "Incorrect input, the number of columns of the first matrix is not "
        "equal to the number of rows of the second matrix");
  }
  S21Matrix result_matrix(other.rows_, cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      long double tmp_sum = 0;
      for (int k = 0; k < cols_; k++) {
        tmp_sum += matrix_[i][k] * other.matrix_[k][j];
      }
      result_matrix.matrix_[i][j] = tmp_sum;
    }
  }
  *this = result_matrix;
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix result_matrix(cols_, rows_);
  for (int i = 0; i < cols_; i++) {
    for (int j = 0; j < rows_; j++) {
      result_matrix.matrix_[i][j] = matrix_[j][i];
    }
  }
  return result_matrix;
}

double S21Matrix::Determinant() {
  if (rows_ != cols_) {
    throw std::out_of_range("Incorrect input, the matrix must be square");
  }
  double result = HelperDeterminant(*this);
  return result;
}

double S21Matrix::HelperDeterminant(S21Matrix const& other) {
  int size = other.rows_;
  if (size == 1) {
    return other.matrix_[0][0];
  }
  double result = 0;
  S21Matrix* new_matrix = new S21Matrix(size - 1, size - 1);
  for (int k = 0; k < size; k++) {
    for (int i = 1, i1 = 0; i < size; i++, i1++) {
      for (int j = 0, j1 = 0; j < size; j++) {
        if (k != j) {
          new_matrix->matrix_[i1][j1] = other.matrix_[i][j];
          j1 += 1;
        }
      }
    }
    int sign = (k % 2) ? -1 : 1;
    result += sign * other.matrix_[0][k] * HelperDeterminant(*new_matrix);
  }
  delete new_matrix;
  return result;
}

void S21Matrix::CreateAddMatrix(S21Matrix& tmp_matrix, int l, int m) {
  for (int i = 0, i1 = 0, j1 = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (j != m && i != l) {
        tmp_matrix.matrix_[i1][j1] = matrix_[i][j];
        j1++;
      }
    }
    if (j1 != 0) i1++;
    j1 = 0;
  }
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_) {
    throw std::out_of_range("Incorrect input, the matrix must be square");
  }
  S21Matrix result(rows_, cols_);
  if (rows_ == 1) {
    result.matrix_[0][0] = matrix_[0][0];
  } else {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        int sign = pow(-1, i + j);
        S21Matrix* tmp_matrix = new S21Matrix(rows_ - 1, cols_ - 1);
        CreateAddMatrix(*tmp_matrix, i, j);
        double value = tmp_matrix->Determinant();
        result.matrix_[i][j] = sign * value;
        delete tmp_matrix;
      }
    }
  }
  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
  double determinant = Determinant();
  if (determinant == 0) {
    throw std::range_error("The determinat must not be zero");
  }
  S21Matrix result(rows_, cols_);
  S21Matrix tmp_matrix = CalcComplements();
  S21Matrix tmp_matrix_add = tmp_matrix.Transpose();
  long double tmp = 1.0 / determinant;
  result.CopyData(tmp_matrix_add);
  result.MulNumber(tmp);
  return result;
}

int S21Matrix::GetRows() const { return rows_; }

int S21Matrix::GetCols() const { return cols_; }

void S21Matrix::SetRows(int rows) {
  if (rows <= 0) {
    throw std::invalid_argument("Argument don't correct");
  }
  S21Matrix new_matrix(rows, cols_);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols_; j++) {
      if (i >= rows_) {
        new_matrix.matrix_[i][j] = 0;
      } else {
        new_matrix.matrix_[i][j] = matrix_[i][j];
      }
    }
  }
  *this = new_matrix;
}

void S21Matrix::SetCols(int cols) {
  if (cols <= 0) {
    throw std::invalid_argument("Argument don't correct");
  }
  S21Matrix new_matrix(rows_, cols);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols; j++) {
      if (j >= cols_) {
        new_matrix.matrix_[i][j] = 0;
      } else {
        new_matrix.matrix_[i][j] = matrix_[i][j];
      }
    }
  }
  *this = new_matrix;
}

S21Matrix& S21Matrix::operator=(S21Matrix const& other) {
  if (this != &other) {
    MemoryFree();
    cols_ = other.cols_;
    rows_ = other.rows_;
    MemoryAllocate();
    CopyData(other);
  }
  return *this;
}

S21Matrix& S21Matrix::operator=(S21Matrix&& other) noexcept {
  if (this != &other) {
    MemoryFree();
    cols_ = other.cols_;
    rows_ = other.rows_;
    matrix_ = other.matrix_;
    other.matrix_ = nullptr;
    other.cols_ = 0;
    other.rows_ = 0;
  }
  return *this;
}

S21Matrix S21Matrix::operator+(S21Matrix const& other) {
  S21Matrix result = *this;
  result.SumMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator-(S21Matrix const& other) {
  S21Matrix result = *this;
  result.SubMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(S21Matrix const& other) {
  S21Matrix result = *this;
  result.MulMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(double const number) {
  S21Matrix result = *this;
  result.MulNumber(number);
  return result;
}

S21Matrix operator*(double number, S21Matrix& other) { return other * number; }

bool S21Matrix::operator==(S21Matrix const& other) const {
  return EqMatrix(other);
}

S21Matrix S21Matrix::operator+=(S21Matrix const& other) {
  SumMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator-=(S21Matrix const& other) {
  SubMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator*=(S21Matrix const& other) {
  MulMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator*=(double const number) {
  MulNumber(number);
  return *this;
}

double& S21Matrix::operator()(int i, int j) const {
  if (i < 0 || i >= rows_ || j < 0 || j >= cols_) {
    throw std::out_of_range("The index isn't included in the range");
  }
  return matrix_[i][j];
}