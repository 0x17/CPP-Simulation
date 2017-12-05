//
// Created by André Schnabel on 03.12.17.
//

#include <vector>
#include <iostream>
#include "Matrix.h"

#pragma once

template<class T>
class Matrix3D {
    int l, m, n;
    std::vector<T> data;
public:
    template<class Func>
    Matrix3D(int _l, int _m , int _n, Func code) : l(_l), m(_m), n(_n), data(_l*_m*_n) {
        for(int h=0; h<l; h++)
            for(int i=0; i<m; i++)
                for(int j=0; j<n; j++)
                    data[indexOfElement(h,i,j)] = code(h,i,j);
    }

    Matrix3D(const Matrix3D& mx) : l(mx.l), m(mx.m), n(mx.n), data(mx.data) {}

    Matrix3D(const std::vector<std::vector<std::vector<T>>> &matrices)
    : l(static_cast<int>(matrices.size())),
    m(static_cast<int>(matrices[0].size())),
    n(static_cast<int>(matrices[0][0].size())),
    data(l*m*n) {
        for(int h=0; h<l; h++)
            for(int i=0; i<m; i++)
                for(int j=0; j<n; j++)
                    data[indexOfElement(h, i, j)] = matrices[h][i][j];
    }

    inline int indexOfElement(int h, int i, int j) const {
        return h * m * n + i * n + j;
    }

    Matrix3D() : l(0), m(0), n(0)  {}

    Matrix3D(int _l, int _m, int _n) : l(_l), m(_m), n(_n), data(_l*_m*_n) {}

    Matrix3D(int _l, int _m, int _n, int value) : Matrix3D(_l, _m, _n, [value](int h, int i, int j) { return value; }) {}

    ~Matrix3D() = default;

    int getL() const { return l; }
    int getM() const { return m; }
    int getN() const { return n; }

    inline T operator()(int h, int i, int j) const { return data[indexOfElement(h, i, j)]; }
    inline T &operator()(int h, int i, int j) { return data[indexOfElement(h, i, j)]; }

    T at(int h, int i, int j) const { return data[indexOfElement(h, i, j)]; }

    Matrix3D &operator=(const Matrix3D &mx) {
        data = mx.data;
        l = mx.l;
        m = mx.m;
        n = mx.n;
        return *this;
    }

    void resize(int _l, int _m, int _n) {
        data.resize(_l*_m*_n);
        l = _l;
        m = _m;
        n = _n;
    }

    Matrix<T> matrix(int h) {
        return Matrix<T>(m, n, [this,h](int i, int j) {
            return data[indexOfElement(h, i, j)];
        });
    }

    std::vector<T> row(int h, int i) const {
        std::vector<T> r(n);
        for(int j=0; j<n; j++)
            r[j] = data[indexOfElement(h, i, j)];
        return r;
    }

    void setRow(int h, int i, const std::vector<T> &row) {
        for(int j=0; j<n; j++)
            data[indexOfElement(h, i, j)] = row[j];
    }

    void setMatrix(int h, const Matrix<T> &mx) {
        for(int i=0; i<m; i++)
            for(int j=0; j<n; j++)
                data[indexOfElement(h, i, j)] = mx(i, j);
    }

    std::vector<T> column(int h, int j) const {
        std::vector<T> c(m);
        for (int i = 0; i<m; i++)
            c[i] = data[indexOfElement(h, i, j)];
        return c;
    }

    template<class Func>
    void foreach(Func f) const {
        for(int h=0; h<l; h++)
            for (int i = 0; i < m; i++)
                for (int j = 0; j < n; j++)
                    f(h, i, j, data[indexOfElement(h, i, j)]);
    }

    template<class Func>
    void foreach2(Func f) const {
        for (int h = 0; h < l; h++)
            for (int i = 0; i < m; i++)
                for (int j = 0; j < n; j++)
                    f(h, i, j);
    }

    template<class Func>
    void foreachAssign(Func f) {
        for (int h = 0; h < l; h++)
            for (int i = 0; i < m; i++)
                for (int j = 0; j < n; j++)
                    data[indexOfElement(h, i, j)] = f(h, i, j);
    }

    std::string toString() const {
        static auto commaExceptLast = [](int curIndex, int ub) {
            return curIndex + 1 == ub ? "" : ",";
        };

        std::stringstream out;
        out << "Matrix3D(l=" << l << ",m=" << m << ",n=" << n << "," << std::endl << "{";
        for(int h=0; h<l; h++) {
            out << "{";
            for (int i = 0; i < m; i++) {
                out << "{";
                for (int j = 0; j < n; j++) {
                    out << at(h, i, j) << commaExceptLast(j, n);
                }
                out << "}" << commaExceptLast(i, m) << std::endl;
            }
            out << "}" << commaExceptLast(h, l) << std::endl;
        }
        out << "}" << std::endl;
        return out.str();
    }

};