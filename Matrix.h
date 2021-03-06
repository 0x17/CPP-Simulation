//
// Created by André Schnabel on 05.11.15.
//

#pragma once

#include <vector>
#include <sstream>

template<class T>
class Matrix {
    int m, n;
    std::vector<T> data;
public:
    template<class Func>
    Matrix(int _m , int _n, Func code) : m(_m), n(_n), data(_m*_n) {
        for(int i=0; i<m; i++)
            for(int j=0; j<n; j++)
                data[i*n+j] = code(i,j);
    }

    Matrix(const Matrix& mx) : m(mx.m), n(mx.n), data(mx.data) {}

    Matrix(const std::vector<std::vector<T>> &rows)
    : m(static_cast<int>(rows.size())), n(static_cast<int>(rows[0].size())), data(m*n) {
        for(int i=0; i<m; i++)
            for(int j=0; j<n; j++)
                data[i*n+j] = rows[i][j];
    }

	enum class Mode {
	    COLUMN_VECTOR,
		ROW_VECTOR
    };

	Matrix(Mode m, const std::vector<T> &vec)
		: Matrix(
			m == Mode::COLUMN_VECTOR ? static_cast<int>(vec.size()) : 1,
			m == Mode::COLUMN_VECTOR ? 1 : static_cast<int>(vec.size()),
			[&](int i, int j) { return vec[m == Mode::COLUMN_VECTOR ? i : j]; }) {
    }

    Matrix() : m(0), n(0) {}

    Matrix(int _m, int _n) : m(_m), n(_n), data(_m*_n) {}

	Matrix(int _m, int _n, T value) : Matrix(_m, _n, [value](int i, int j) { return value; }) {}

    ~Matrix() = default;

	int getM() const { return m; }
    int getN() const { return n; }

    inline T operator()(int i, int j) const { return data[i*n+j]; }
    inline T &operator()(int i, int j) { return data[i*n+j]; }

	T at(int i, int j) const { return data[i*n + j]; }

    Matrix &operator=(const Matrix &mx) {
		data = mx.data;
		m = mx.m;
		n = mx.n;
		return *this;
    }

    void resize(int _m, int _n) {
        data.resize(_m*_n);
        m = _m;
        n = _n;
    }

    std::vector<T> row(int i) const {
        std::vector<T> r(n);
        for(int j=0; j<n; j++)
            r[j] = data[i*n+j];
        return r;
    }

	void setRow(int i, const std::vector<T> &row) {
		for(int j=0; j<n; j++)
			data[i*n+j] = row[j];
	}

	std::vector<T> column(int j) const {
		std::vector<T> c(m);
		for (int i = 0; i<m; i++)
			c[i] = data[i*n + j];
		return c;
    }

	template<class Func>
	void foreach(Func f) const {
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				f(i, j, data[i*n + j]);
    }

	template<class Func>
	void foreach2(Func f) const {
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				f(i, j);
	}

	template<class Func>
	void foreachAssign(Func f) {
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				data[i*n + j] = f(i, j);
    }

	std::string toString() const {
		std::stringstream out;
		out << "Matrix(m=" << m << ",n=" << n << "," << std::endl << "{";
		for (int i = 0; i < m; i++) {
			out << "{";
			for (int j = 0; j < n; j++) {
				out << at(i, j) << (j + 1 == n ? "" : ",");
			}
			out << "}" << (i + 1 == m ? "" : ",") << std::endl;
		}
		out << "}" << std::endl;
		return out.str();
	}

	void zeroOut() {
		for(int i=0; i<m; i++)
			for(int j=0; j<n; j++)
				data[i*n+j] = 0;
	}

	std::vector<std::vector<T>> nestedVectors() {
		std::vector<std::vector<T>> res(m);
		for(int i=0; i<m; i++) {
			std::vector<T> rrow(n);
			for (int j = 0; j < n; j++) {
				rrow[j] = data[i*n+j];
			}
			res[i] = rrow;
		}
		return res;
	}
};

template<class T>
using Matrix2D = Matrix<T>;
