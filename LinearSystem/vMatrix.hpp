#pragma once

#define range(i, begin, end) for(size_t i = begin; i < (end); i++)

//#include <fstream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <string>

//#include "Frac.hpp"
using namespace std;

template <typename Field>
class Matrix {
    vector< vector<Field>* > _M;
    Matrix& self = *this;
    
    void clear() {
        for(auto& row: _M) delete row;
    }
    
public:
    
    size_t rows() const {
        return _M.size();
    }
    
    size_t cols() const {
        if(rows() == 0) return 0;
        return _M[0]->size();
    }
    
    
    Matrix(size_t rows, size_t cols) {
        _M.resize(rows);
        for(auto& row: _M) row = new vector<Field>(cols);
    }
    
    Matrix& operator = (const Matrix& A) {
        clear();
        _M.resize(A.rows());
        range(i, 0, A.rows()) _M[i] = new vector<Field>(*(A._M[i]));
        return self;
    }
    
    Matrix& operator = (Matrix&& A) {
        clear();
        _M = A._M;
        A._M.resize(0);
        return self;
    }
    
    Matrix(const Matrix &A) {
        self = A;
    }
    
    Matrix(Matrix&& A) {
        _M = A._M;
        A._M.resize(0);
    }
    
    ~Matrix() {
        clear();
    }
    
    
    
    Field& operator() (size_t i, size_t j) const {
        return (*_M[i])[j];
    }
    
    void add(size_t to, size_t from, Field k) {
        range(j, 0, cols()) self(to, j) += k * self(from, j);
    }
    
    void swap(size_t first, size_t second) {
        std::swap(_M[first], _M[second]);
    }
    
    void mul(size_t row, Field k) {
        range(j, 0, cols()) self(row,j) = self(row, j) * k;
    }
    
    void div(size_t row, Field k) {
        range(j, 0, cols()) self(row, j) = self(row,j) / k;
    }
    
    pair< bool, vector<bool> > row_reduce(vector<Field>& b) {
        if(b.size() != rows()) {
            throw "length of b don't match with quantity of rows";
        }
        vector<bool> pivots(cols(), false);
        
        size_t pivot = 0;
        range(k, 0, cols()) {
            range(i, pivot, rows()) {
                if(self(i,k) != Field()) {
                    pivots[k] = true;
                    
                    // swap
                    std::swap(b[pivot], b[i]);
                    swap(pivot, i);
                    
                    // div
                    b[pivot] /= self(pivot, k);
                    div(pivot, self(pivot, k));
                    
                    // nullify all elements in the row except the pivot
                    range(ii, 0, rows()) {
                        if(i != ii) {
                            b[ii] -= b[pivot]*self(ii, k);
                            add(ii, pivot, -self(ii,k));
                        }
                    }
                    pivot++;
                    break;
                }
            }
        }
        bool there_is_a_solution = true;
        range(i, pivot, rows()) {
            if(b[i] != Field()) {
                there_is_a_solution = false;
            }
        }
        
        return make_pair(there_is_a_solution, pivots);
    }
    
    friend istream& operator >> (istream& in, Matrix &M) {
        range(i, 0, M.rows())
            range(j, 0, M.cols())
                in >> M(i,j);
    }
    
    
    friend ostream& operator << (ostream& out, const Matrix &M) {
        out << "\\begin{bmatrix}\n";
        range(i, 0, M.rows()) {
            out << "    ";
            range(j, 0, M.cols()) {
                out << M(i, j) << " ";
            }
            out << (i < M.rows() ? "\\\\\n" : "\n");
        }
        out << "\\end{bmatrix}";
        return out;
    }
    
    friend Matrix operator * (const Matrix& A, const Matrix& B) {
        if(A.cols() != B.rows()) {
            throw "Sizes don't match";
        }
        Matrix C(A.rows(), B.cols());
        range(i, 0, A.rows()) {
            range(j, 0, B.cols()) {
                range(k, 0, A.cols()) {
                    C(i, j) += A(i, k) * B(k, j);
                }
            }
        }
        return C;
    }
    
    friend Matrix operator + (const Matrix& A, const Matrix& B) {
        if( (A.cols() != B.cols()) && (A.rows() != B.rows()) ) {
            throw "Sizes don't match";
        }
        Matrix C(A.rows(), B.cols());
        range(i, 0, A.rows()) {
            range(j, 0, A.cols()) {
                C(i, j) = A(i, j) + B(i, j);
            }
        }
        return C;
    }
    
    friend Matrix operator - (const Matrix& A) {
        Matrix B(A);
        range(i, 0, B.rows()) {
            range(j, 0, B.cols()) {
                B(i, j) = -B(i, j);
            }
        }
        return B;
    }
    
    friend Matrix operator - (const Matrix& A, const Matrix& B) {
        // return A + (-B);
        if( (A.cols() != B.cols()) || (A.rows() != B.rows()) ) {
            throw "Sizes don't match";
        }
        Matrix C(A.rows(), B.cols());
        range(i, 0, A.rows()) {
            range(j, 0, A.cols()) {
                C(i, j) = A(i, j) - B(i, j);
            }
        }
        return C;
        
    }
    
    friend bool operator == (const Matrix& A, const Matrix& B) {
        if( (A.cols() != B.cols()) || (A.rows() != B.rows()) ) {
            throw "Sizes don't match";
        }
        range(i, 0, A.rows()) {
            range(j, 0, A.cols()) {
                if(A(i, j) != B(i, j)) return false;
            }
        }
        return true;
    }
    
    friend bool operator != (const Matrix& A, const Matrix& B) {
        return !(A == B);
    }
    
    void set(size_t i, size_t j, const Field& val) {
        self(i, j) = val;
    }
    
    
    string tex_solution(vector<Field>& b) {
        Matrix A(self);
        
        stringstream out;
        
        bool is_there_a_solution;
        vector<bool> pivots;
        
        tie(is_there_a_solution, pivots) = A.row_reduce(b);
        
        if(!is_there_a_solution) {
            out << "There is no soulution\n";
            return out.str();
        }
        
        out << "X = \n";
        out << "\\begin{bmatrix}\n";
        
        range(j, 0, A.cols()) {
            out << "    " << (j < b.size() ? b[j] : Field())
            << ((j != A.cols()-1) ? "\\\\\n" : "\n");
        }
        
        out << "\\end{bmatrix}\n";
        
        size_t numeration = 0;
        range(k, 0, pivots.size()) {
            if(!pivots[k]) {
                out << " + \n" << "\\begin{bmatrix}\n";
                range(i, 0, A.cols()) {
                    out << "    " << ((i == k) ? -Field(1): ( i < A.rows() ? A(i, k) : Field()))
                    << ((i != A.cols()-1) ? "\\\\\n" : "\n");
                }
                numeration++;
                out << "\\end{bmatrix}C_{" << numeration << "}\n" ;
            }
        }
        return out.str();
    }
};
