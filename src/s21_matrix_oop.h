#ifndef CPP1_S21_MATRIXPLUS_SRC_S21_MATRIX_OOP_H_
#define CPP1_S21_MATRIXPLUS_SRC_S21_MATRIX_OOP_H_

#include <cmath>
#include <cstring>
#include <iostream>

class S21Matrix {
 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other) noexcept;
  ~S21Matrix();

  bool EqMatrix(const S21Matrix& other) const;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();

  S21Matrix& operator=(S21Matrix const& other);
  S21Matrix& operator=(S21Matrix&& other) noexcept;
  S21Matrix operator+(S21Matrix const& other);
  S21Matrix operator-(S21Matrix const& other);
  S21Matrix operator*(S21Matrix const& other);
  S21Matrix operator*(double const number);
  friend S21Matrix operator*(double number, S21Matrix& other);
  bool operator==(S21Matrix const& other) const;
  S21Matrix operator+=(S21Matrix const& other);
  S21Matrix operator-=(S21Matrix const& other);
  S21Matrix operator*=(S21Matrix const& other);
  S21Matrix operator*=(double const number);
  double& operator()(int i, int j) const;

  int GetRows() const;
  int GetCols() const;
  void SetRows(int rows);
  void SetCols(int cols);

 private:
  int rows_, cols_;
  double** matrix_;

  void MemoryAllocate();
  void MemoryFree();
  void CopyData(S21Matrix const& other);
  double HelperDeterminant(S21Matrix const& other);
  void CreateAddMatrix(S21Matrix& tmp_matrix, int l, int m);
};

#endif  // CPP1_S21_MATRIXPLUS_SRC_S21_MATRIX_OOP_H_