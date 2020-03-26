//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef MATH_DEFS_H
#define MATH_DEFS_H

#ifdef WIN32
#define _USE_MATH_DEFINES
#endif
#include <cmath>

#include <Eigen/Core>
#include <Eigen/StdVector>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/SparseQR>

typedef double scalar;
typedef unsigned long long uint64;

struct int_scalar {
	int i;
	scalar v;

	inline bool operator() (const int_scalar& struct1, const int_scalar& struct2)
	{
		return (struct1.v < struct2.v);
	}
};

struct int5_scalar {
	int row_node;
	int row_dir;
	int col_bucket;
	int col_node;
	int col_dir;
	scalar v;

	int5_scalar(int rn, int rd, int cb, int cn, int cd, scalar value)
		: row_node(rn), row_dir(rd), col_bucket(cb), col_node(cn), col_dir(cd), v(value)
	{}

	int5_scalar()
		: row_node(0), row_dir(0), col_bucket(0), col_node(0), col_dir(0), v(0.0)
	{}
};

typedef Eigen::Matrix<scalar, 2, 1> Vector2s;
typedef Eigen::Matrix<float, 2, 1> Vector2f;
typedef Eigen::Matrix<int, 2, 1> Vector2i;
typedef Eigen::Matrix<scalar, 3, 1> Vector3s;
typedef Eigen::Matrix<float, 3, 1> Vector3f;
typedef Eigen::Matrix<int, 3, 1> Vector3i;
typedef Eigen::Matrix<scalar, 4, 1> Vector4s;
typedef Eigen::Matrix<float, 4, 1> Vector4f;
typedef Eigen::Matrix<int, 4, 1> Vector4i;
typedef Eigen::Matrix<scalar, 5, 1> Vector5s;
typedef Eigen::Matrix<float, 5, 1> Vector5f;
typedef Eigen::Matrix<int, 5, 1> Vector5i;
typedef Eigen::Matrix<scalar, 6, 1> Vector6s;
typedef Eigen::Matrix<float, 6, 1> Vector6f;
typedef Eigen::Matrix<int, 6, 1> Vector6i;
typedef Eigen::Matrix<scalar, 7, 1> Vector7s;
typedef Eigen::Matrix<float, 7, 1> Vector7f;
typedef Eigen::Matrix<int, 7, 1> Vector7i;
typedef Eigen::Matrix<scalar, 8, 1> Vector8s;
typedef Eigen::Matrix<float, 8, 1> Vector8f;
typedef Eigen::Matrix<int, 8, 1> Vector8i;
typedef Eigen::Matrix<scalar, 9, 1> Vector9s;
typedef Eigen::Matrix<float, 9, 1> Vector9f;
typedef Eigen::Matrix<int, 9, 1> Vector9i;
typedef Eigen::Matrix<scalar, 27, 1> Vector27s;
typedef Eigen::Matrix<float, 27, 1> Vector27f;
typedef Eigen::Matrix<int, 27, 1> Vector27i;
typedef Eigen::Matrix<scalar, 36, 1> Vector36s;
typedef Eigen::Matrix<float, 36, 1> Vector36f;
typedef Eigen::Matrix<int, 36, 1> Vector36i;
typedef Eigen::Matrix<scalar, 64, 1> Vector64s;
typedef Eigen::Matrix<float, 64, 1> Vector64f;
typedef Eigen::Matrix<int, 64, 1> Vector64i;
typedef Eigen::Matrix<scalar, 125, 1> Vector125s;
typedef Eigen::Matrix<float, 125, 1> Vector125f;
typedef Eigen::Matrix<int, 125, 1> Vector125i;
typedef Eigen::Matrix<scalar, Eigen::Dynamic, 1> VectorXs;
typedef Eigen::Matrix<float, Eigen::Dynamic, 1> VectorXf;
typedef Eigen::Matrix<int, Eigen::Dynamic, 1> VectorXi;
typedef Eigen::Matrix<char, Eigen::Dynamic, 1> VectorXc;
typedef Eigen::Matrix<unsigned char, Eigen::Dynamic, 1> VectorXuc;

typedef Eigen::Quaternion<scalar> Quat4s;
typedef Eigen::Quaternion<float> Quat4f;

typedef Eigen::Matrix<scalar, 1, 2> Vector2sT;
typedef Eigen::Matrix<float, 1, 2> Vector2fT;
typedef Eigen::Matrix<int, 1, 2> Vector2iT;
typedef Eigen::Matrix<scalar, 1, 3> Vector3sT;
typedef Eigen::Matrix<float, 1, 3> Vector3fT;
typedef Eigen::Matrix<int, 1, 3> Vector3iT;
typedef Eigen::Matrix<scalar, 1, 4> Vector4sT;
typedef Eigen::Matrix<float, 1, 4> Vector4fT;
typedef Eigen::Matrix<int, 1, 4> Vector4iT;
typedef Eigen::Matrix<scalar, 1, 5> Vector5sT;
typedef Eigen::Matrix<float, 1, 5> Vector5fT;
typedef Eigen::Matrix<int, 1, 5> Vector5iT;
typedef Eigen::Matrix<scalar, 1, 6> Vector6sT;
typedef Eigen::Matrix<float, 1, 6> Vector6fT;
typedef Eigen::Matrix<int, 1, 6> Vector6iT;
typedef Eigen::Matrix<scalar, 1, 7> Vector7sT;
typedef Eigen::Matrix<float, 1, 7> Vector7fT;
typedef Eigen::Matrix<int, 1, 7> Vector7iT;
typedef Eigen::Matrix<scalar, 1, 8> Vector8sT;
typedef Eigen::Matrix<float, 1, 8> Vector8fT;
typedef Eigen::Matrix<int, 1, 8> Vector8iT;
typedef Eigen::Matrix<scalar, 1, 9> Vector9sT;
typedef Eigen::Matrix<float, 1, 9> Vector9fT;
typedef Eigen::Matrix<int, 1, 9> Vector9iT;
typedef Eigen::Matrix<scalar, 1, 27> Vector27sT;
typedef Eigen::Matrix<float, 1, 27> Vector27fT;
typedef Eigen::Matrix<int, 1, 27> Vector27iT;
typedef Eigen::Matrix<scalar, 1, 36> Vector36sT;
typedef Eigen::Matrix<float, 1, 36> Vector36fT;
typedef Eigen::Matrix<int, 1, 36> Vector36iT;
typedef Eigen::Matrix<scalar, 1, 64> Vector64sT;
typedef Eigen::Matrix<float, 1, 64> Vector64fT;
typedef Eigen::Matrix<int, 1, 64> Vector64iT;
typedef Eigen::Matrix<scalar, 1, 125> Vector125sT;
typedef Eigen::Matrix<float, 1, 125> Vector125fT;
typedef Eigen::Matrix<int, 1, 125> Vector125iT;
typedef Eigen::Matrix<scalar, 1, Eigen::Dynamic> VectorXsT;
typedef Eigen::Matrix<int, 1, Eigen::Dynamic> VectorXiT;

typedef Eigen::Array<scalar, Eigen::Dynamic, Eigen::Dynamic> ArrayXs;
typedef Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic> ArrayXf;
typedef Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic> ArrayXi;

typedef Eigen::Block<float, 2, 2> Block2f;
typedef Eigen::Block<float, 2, 3> Block2x3f;
typedef Eigen::Block<float, 2, 4> Block2x4f;
typedef Eigen::Block<float, 2, 5> Block2x5f;
typedef Eigen::Block<float, 2, 6> Block2x6f;
typedef Eigen::Block<float, 2, 7> Block2x7f;
typedef Eigen::Block<float, 2, 8> Block2x8f;
typedef Eigen::Block<float, 2, 9> Block2x9f;
typedef Eigen::Block<float, 3, 2> Block3x2f;
typedef Eigen::Block<float, 3, 3> Block3f;
typedef Eigen::Block<float, 3, 4> Block3x4f;
typedef Eigen::Block<float, 3, 5> Block3x5f;
typedef Eigen::Block<float, 3, 6> Block3x6f;
typedef Eigen::Block<float, 3, 7> Block3x7f;
typedef Eigen::Block<float, 3, 8> Block3x8f;
typedef Eigen::Block<float, 3, 9> Block3x9f;
typedef Eigen::Block<float, 4, 2> Block4x2f;
typedef Eigen::Block<float, 4, 3> Block4x3f;
typedef Eigen::Block<float, 4, 4> Block4f;
typedef Eigen::Block<float, 4, 5> Block4x5f;
typedef Eigen::Block<float, 4, 6> Block4x6f;
typedef Eigen::Block<float, 4, 7> Block4x7f;
typedef Eigen::Block<float, 4, 8> Block4x8f;
typedef Eigen::Block<float, 4, 9> Block4x9f;
typedef Eigen::Block<float, 5, 2> Block5x2f;
typedef Eigen::Block<float, 5, 3> Block5x3f;
typedef Eigen::Block<float, 5, 4> Block5x4f;
typedef Eigen::Block<float, 5, 5> Block5f;
typedef Eigen::Block<float, 5, 6> Block5x6f;
typedef Eigen::Block<float, 5, 7> Block5x7f;
typedef Eigen::Block<float, 5, 8> Block5x8f;
typedef Eigen::Block<float, 5, 9> Block5x9f;
typedef Eigen::Block<float, 6, 2> Block6x2f;
typedef Eigen::Block<float, 6, 3> Block6x3f;
typedef Eigen::Block<float, 6, 4> Block6x4f;
typedef Eigen::Block<float, 6, 5> Block6x5f;
typedef Eigen::Block<float, 6, 6> Block6f;
typedef Eigen::Block<float, 6, 7> Block6x7f;
typedef Eigen::Block<float, 6, 8> Block6x8f;
typedef Eigen::Block<float, 6, 9> Block6x9f;
typedef Eigen::Block<float, 7, 2> Block7x2f;
typedef Eigen::Block<float, 7, 3> Block7x3f;
typedef Eigen::Block<float, 7, 4> Block7x4f;
typedef Eigen::Block<float, 7, 5> Block7x5f;
typedef Eigen::Block<float, 7, 6> Block7x6f;
typedef Eigen::Block<float, 7, 7> Block7f;
typedef Eigen::Block<float, 7, 8> Block7x8f;
typedef Eigen::Block<float, 7, 9> Block7x9f;
typedef Eigen::Block<float, 8, 2> Block8x2f;
typedef Eigen::Block<float, 8, 3> Block8x3f;
typedef Eigen::Block<float, 8, 4> Block8x4f;
typedef Eigen::Block<float, 8, 5> Block8x5f;
typedef Eigen::Block<float, 8, 6> Block8x6f;
typedef Eigen::Block<float, 8, 7> Block8x7f;
typedef Eigen::Block<float, 8, 8> Block8f;
typedef Eigen::Block<float, 8, 9> Block8x9f;
typedef Eigen::Block<float, 9, 2> Block9x2f;
typedef Eigen::Block<float, 9, 3> Block9x3f;
typedef Eigen::Block<float, 9, 4> Block9x4f;
typedef Eigen::Block<float, 9, 5> Block9x5f;
typedef Eigen::Block<float, 9, 6> Block9x6f;
typedef Eigen::Block<float, 9, 7> Block9x7f;
typedef Eigen::Block<float, 9, 8> Block9x8f;
typedef Eigen::Block<float, 9, 9> Block9f;

typedef Eigen::Matrix<float, 2, 2> Matrix2f;
typedef Eigen::Matrix<float, 2, 3> Matrix2x3f;
typedef Eigen::Matrix<float, 2, 4> Matrix2x4f;
typedef Eigen::Matrix<float, 2, 5> Matrix2x5f;
typedef Eigen::Matrix<float, 2, 6> Matrix2x6f;
typedef Eigen::Matrix<float, 2, 7> Matrix2x7f;
typedef Eigen::Matrix<float, 2, 8> Matrix2x8f;
typedef Eigen::Matrix<float, 2, 9> Matrix2x9f;
typedef Eigen::Matrix<float, 3, 2> Matrix3x2f;
typedef Eigen::Matrix<float, 3, 3> Matrix3f;
typedef Eigen::Matrix<float, 3, 4> Matrix3x4f;
typedef Eigen::Matrix<float, 3, 5> Matrix3x5f;
typedef Eigen::Matrix<float, 3, 6> Matrix3x6f;
typedef Eigen::Matrix<float, 3, 7> Matrix3x7f;
typedef Eigen::Matrix<float, 3, 8> Matrix3x8f;
typedef Eigen::Matrix<float, 3, 9> Matrix3x9f;
typedef Eigen::Matrix<float, 4, 2> Matrix4x2f;
typedef Eigen::Matrix<float, 4, 3> Matrix4x3f;
typedef Eigen::Matrix<float, 4, 4> Matrix4f;
typedef Eigen::Matrix<float, 4, 5> Matrix4x5f;
typedef Eigen::Matrix<float, 4, 6> Matrix4x6f;
typedef Eigen::Matrix<float, 4, 7> Matrix4x7f;
typedef Eigen::Matrix<float, 4, 8> Matrix4x8f;
typedef Eigen::Matrix<float, 4, 9> Matrix4x9f;
typedef Eigen::Matrix<float, 5, 2> Matrix5x2f;
typedef Eigen::Matrix<float, 5, 3> Matrix5x3f;
typedef Eigen::Matrix<float, 5, 4> Matrix5x4f;
typedef Eigen::Matrix<float, 5, 5> Matrix5f;
typedef Eigen::Matrix<float, 5, 6> Matrix5x6f;
typedef Eigen::Matrix<float, 5, 7> Matrix5x7f;
typedef Eigen::Matrix<float, 5, 8> Matrix5x8f;
typedef Eigen::Matrix<float, 5, 9> Matrix5x9f;
typedef Eigen::Matrix<float, 6, 2> Matrix6x2f;
typedef Eigen::Matrix<float, 6, 3> Matrix6x3f;
typedef Eigen::Matrix<float, 6, 4> Matrix6x4f;
typedef Eigen::Matrix<float, 6, 5> Matrix6x5f;
typedef Eigen::Matrix<float, 6, 6> Matrix6f;
typedef Eigen::Matrix<float, 6, 7> Matrix6x7f;
typedef Eigen::Matrix<float, 6, 8> Matrix6x8f;
typedef Eigen::Matrix<float, 6, 9> Matrix6x9f;
typedef Eigen::Matrix<float, 7, 2> Matrix7x2f;
typedef Eigen::Matrix<float, 7, 3> Matrix7x3f;
typedef Eigen::Matrix<float, 7, 4> Matrix7x4f;
typedef Eigen::Matrix<float, 7, 5> Matrix7x5f;
typedef Eigen::Matrix<float, 7, 6> Matrix7x6f;
typedef Eigen::Matrix<float, 7, 7> Matrix7f;
typedef Eigen::Matrix<float, 7, 8> Matrix7x8f;
typedef Eigen::Matrix<float, 7, 9> Matrix7x9f;
typedef Eigen::Matrix<float, 8, 2> Matrix8x2f;
typedef Eigen::Matrix<float, 8, 3> Matrix8x3f;
typedef Eigen::Matrix<float, 8, 4> Matrix8x4f;
typedef Eigen::Matrix<float, 8, 5> Matrix8x5f;
typedef Eigen::Matrix<float, 8, 6> Matrix8x6f;
typedef Eigen::Matrix<float, 8, 7> Matrix8x7f;
typedef Eigen::Matrix<float, 8, 8> Matrix8f;
typedef Eigen::Matrix<float, 8, 9> Matrix8x9f;
typedef Eigen::Matrix<float, 9, 2> Matrix9x2f;
typedef Eigen::Matrix<float, 9, 3> Matrix9x3f;
typedef Eigen::Matrix<float, 9, 4> Matrix9x4f;
typedef Eigen::Matrix<float, 9, 5> Matrix9x5f;
typedef Eigen::Matrix<float, 9, 6> Matrix9x6f;
typedef Eigen::Matrix<float, 9, 7> Matrix9x7f;
typedef Eigen::Matrix<float, 9, 8> Matrix9x8f;
typedef Eigen::Matrix<float, 9, 9> Matrix9f;

typedef Eigen::Block<scalar, 2, 2> Block2s;
typedef Eigen::Block<scalar, 2, 3> Block2x3s;
typedef Eigen::Block<scalar, 2, 4> Block2x4s;
typedef Eigen::Block<scalar, 2, 5> Block2x5s;
typedef Eigen::Block<scalar, 2, 6> Block2x6s;
typedef Eigen::Block<scalar, 2, 7> Block2x7s;
typedef Eigen::Block<scalar, 2, 8> Block2x8s;
typedef Eigen::Block<scalar, 2, 9> Block2x9s;
typedef Eigen::Block<scalar, 3, 2> Block3x2s;
typedef Eigen::Block<scalar, 3, 3> Block3s;
typedef Eigen::Block<scalar, 3, 4> Block3x4s;
typedef Eigen::Block<scalar, 3, 5> Block3x5s;
typedef Eigen::Block<scalar, 3, 6> Block3x6s;
typedef Eigen::Block<scalar, 3, 7> Block3x7s;
typedef Eigen::Block<scalar, 3, 8> Block3x8s;
typedef Eigen::Block<scalar, 3, 9> Block3x9s;
typedef Eigen::Block<scalar, 4, 2> Block4x2s;
typedef Eigen::Block<scalar, 4, 3> Block4x3s;
typedef Eigen::Block<scalar, 4, 4> Block4s;
typedef Eigen::Block<scalar, 4, 5> Block4x5s;
typedef Eigen::Block<scalar, 4, 6> Block4x6s;
typedef Eigen::Block<scalar, 4, 7> Block4x7s;
typedef Eigen::Block<scalar, 4, 8> Block4x8s;
typedef Eigen::Block<scalar, 4, 9> Block4x9s;
typedef Eigen::Block<scalar, 5, 2> Block5x2s;
typedef Eigen::Block<scalar, 5, 3> Block5x3s;
typedef Eigen::Block<scalar, 5, 4> Block5x4s;
typedef Eigen::Block<scalar, 5, 5> Block5s;
typedef Eigen::Block<scalar, 5, 6> Block5x6s;
typedef Eigen::Block<scalar, 5, 7> Block5x7s;
typedef Eigen::Block<scalar, 5, 8> Block5x8s;
typedef Eigen::Block<scalar, 5, 9> Block5x9s;
typedef Eigen::Block<scalar, 6, 2> Block6x2s;
typedef Eigen::Block<scalar, 6, 3> Block6x3s;
typedef Eigen::Block<scalar, 6, 4> Block6x4s;
typedef Eigen::Block<scalar, 6, 5> Block6x5s;
typedef Eigen::Block<scalar, 6, 6> Block6s;
typedef Eigen::Block<scalar, 6, 7> Block6x7s;
typedef Eigen::Block<scalar, 6, 8> Block6x8s;
typedef Eigen::Block<scalar, 6, 9> Block6x9s;
typedef Eigen::Block<scalar, 7, 2> Block7x2s;
typedef Eigen::Block<scalar, 7, 3> Block7x3s;
typedef Eigen::Block<scalar, 7, 4> Block7x4s;
typedef Eigen::Block<scalar, 7, 5> Block7x5s;
typedef Eigen::Block<scalar, 7, 6> Block7x6s;
typedef Eigen::Block<scalar, 7, 7> Block7s;
typedef Eigen::Block<scalar, 7, 8> Block7x8s;
typedef Eigen::Block<scalar, 7, 9> Block7x9s;
typedef Eigen::Block<scalar, 8, 2> Block8x2s;
typedef Eigen::Block<scalar, 8, 3> Block8x3s;
typedef Eigen::Block<scalar, 8, 4> Block8x4s;
typedef Eigen::Block<scalar, 8, 5> Block8x5s;
typedef Eigen::Block<scalar, 8, 6> Block8x6s;
typedef Eigen::Block<scalar, 8, 7> Block8x7s;
typedef Eigen::Block<scalar, 8, 8> Block8s;
typedef Eigen::Block<scalar, 8, 9> Block8x9s;
typedef Eigen::Block<scalar, 9, 2> Block9x2s;
typedef Eigen::Block<scalar, 9, 3> Block9x3s;
typedef Eigen::Block<scalar, 9, 4> Block9x4s;
typedef Eigen::Block<scalar, 9, 5> Block9x5s;
typedef Eigen::Block<scalar, 9, 6> Block9x6s;
typedef Eigen::Block<scalar, 9, 7> Block9x7s;
typedef Eigen::Block<scalar, 9, 8> Block9x8s;
typedef Eigen::Block<scalar, 9, 9> Block9s;

typedef Eigen::Matrix<scalar, 2, 2> Matrix2s;
typedef Eigen::Matrix<scalar, 2, 3> Matrix2x3s;
typedef Eigen::Matrix<scalar, 2, 4> Matrix2x4s;
typedef Eigen::Matrix<scalar, 2, 5> Matrix2x5s;
typedef Eigen::Matrix<scalar, 2, 6> Matrix2x6s;
typedef Eigen::Matrix<scalar, 2, 7> Matrix2x7s;
typedef Eigen::Matrix<scalar, 2, 8> Matrix2x8s;
typedef Eigen::Matrix<scalar, 2, 9> Matrix2x9s;
typedef Eigen::Matrix<scalar, 3, 2> Matrix3x2s;
typedef Eigen::Matrix<scalar, 3, 3> Matrix3s;
typedef Eigen::Matrix<scalar, 3, 4> Matrix3x4s;
typedef Eigen::Matrix<scalar, 3, 5> Matrix3x5s;
typedef Eigen::Matrix<scalar, 3, 6> Matrix3x6s;
typedef Eigen::Matrix<scalar, 3, 7> Matrix3x7s;
typedef Eigen::Matrix<scalar, 3, 8> Matrix3x8s;
typedef Eigen::Matrix<scalar, 3, 9> Matrix3x9s;
typedef Eigen::Matrix<scalar, 4, 2> Matrix4x2s;
typedef Eigen::Matrix<scalar, 4, 3> Matrix4x3s;
typedef Eigen::Matrix<scalar, 4, 4> Matrix4s;
typedef Eigen::Matrix<scalar, 4, 5> Matrix4x5s;
typedef Eigen::Matrix<scalar, 4, 6> Matrix4x6s;
typedef Eigen::Matrix<scalar, 4, 7> Matrix4x7s;
typedef Eigen::Matrix<scalar, 4, 8> Matrix4x8s;
typedef Eigen::Matrix<scalar, 4, 9> Matrix4x9s;
typedef Eigen::Matrix<scalar, 5, 2> Matrix5x2s;
typedef Eigen::Matrix<scalar, 5, 3> Matrix5x3s;
typedef Eigen::Matrix<scalar, 5, 4> Matrix5x4s;
typedef Eigen::Matrix<scalar, 5, 5> Matrix5s;
typedef Eigen::Matrix<scalar, 5, 6> Matrix5x6s;
typedef Eigen::Matrix<scalar, 5, 7> Matrix5x7s;
typedef Eigen::Matrix<scalar, 5, 8> Matrix5x8s;
typedef Eigen::Matrix<scalar, 5, 9> Matrix5x9s;
typedef Eigen::Matrix<scalar, 6, 2> Matrix6x2s;
typedef Eigen::Matrix<scalar, 6, 3> Matrix6x3s;
typedef Eigen::Matrix<scalar, 6, 4> Matrix6x4s;
typedef Eigen::Matrix<scalar, 6, 5> Matrix6x5s;
typedef Eigen::Matrix<scalar, 6, 6> Matrix6s;
typedef Eigen::Matrix<scalar, 6, 7> Matrix6x7s;
typedef Eigen::Matrix<scalar, 6, 8> Matrix6x8s;
typedef Eigen::Matrix<scalar, 6, 9> Matrix6x9s;
typedef Eigen::Matrix<scalar, 7, 2> Matrix7x2s;
typedef Eigen::Matrix<scalar, 7, 3> Matrix7x3s;
typedef Eigen::Matrix<scalar, 7, 4> Matrix7x4s;
typedef Eigen::Matrix<scalar, 7, 5> Matrix7x5s;
typedef Eigen::Matrix<scalar, 7, 6> Matrix7x6s;
typedef Eigen::Matrix<scalar, 7, 7> Matrix7s;
typedef Eigen::Matrix<scalar, 7, 8> Matrix7x8s;
typedef Eigen::Matrix<scalar, 7, 9> Matrix7x9s;
typedef Eigen::Matrix<scalar, 8, 2> Matrix8x2s;
typedef Eigen::Matrix<scalar, 8, 3> Matrix8x3s;
typedef Eigen::Matrix<scalar, 8, 4> Matrix8x4s;
typedef Eigen::Matrix<scalar, 8, 5> Matrix8x5s;
typedef Eigen::Matrix<scalar, 8, 6> Matrix8x6s;
typedef Eigen::Matrix<scalar, 8, 7> Matrix8x7s;
typedef Eigen::Matrix<scalar, 8, 8> Matrix8s;
typedef Eigen::Matrix<scalar, 8, 9> Matrix8x9s;
typedef Eigen::Matrix<scalar, 9, 2> Matrix9x2s;
typedef Eigen::Matrix<scalar, 9, 3> Matrix9x3s;
typedef Eigen::Matrix<scalar, 9, 4> Matrix9x4s;
typedef Eigen::Matrix<scalar, 9, 5> Matrix9x5s;
typedef Eigen::Matrix<scalar, 9, 6> Matrix9x6s;
typedef Eigen::Matrix<scalar, 9, 7> Matrix9x7s;
typedef Eigen::Matrix<scalar, 9, 8> Matrix9x8s;
typedef Eigen::Matrix<scalar, 9, 9> Matrix9s;
typedef Eigen::Matrix<scalar, 27, 3> Matrix27x3s;
typedef Eigen::Matrix<scalar, 27, 4> Matrix27x4s;
typedef Eigen::Matrix<scalar, 27, 5> Matrix27x5s;
typedef Eigen::Matrix<scalar, 36, 3> Matrix36x3s;
typedef Eigen::Matrix<scalar, 64, 3> Matrix64x3s;
typedef Eigen::Matrix<scalar, 64, 3> Matrix125x3s;

typedef Eigen::Block<int, 2, 2> Block2i;
typedef Eigen::Block<int, 2, 3> Block2x3i;
typedef Eigen::Block<int, 2, 4> Block2x4i;
typedef Eigen::Block<int, 2, 5> Block2x5i;
typedef Eigen::Block<int, 2, 6> Block2x6i;
typedef Eigen::Block<int, 2, 7> Block2x7i;
typedef Eigen::Block<int, 2, 8> Block2x8i;
typedef Eigen::Block<int, 2, 9> Block2x9i;
typedef Eigen::Block<int, 3, 2> Block3x2i;
typedef Eigen::Block<int, 3, 3> Block3i;
typedef Eigen::Block<int, 3, 4> Block3x4i;
typedef Eigen::Block<int, 3, 5> Block3x5i;
typedef Eigen::Block<int, 3, 6> Block3x6i;
typedef Eigen::Block<int, 3, 7> Block3x7i;
typedef Eigen::Block<int, 3, 8> Block3x8i;
typedef Eigen::Block<int, 3, 9> Block3x9i;
typedef Eigen::Block<int, 4, 2> Block4x2i;
typedef Eigen::Block<int, 4, 3> Block4x3i;
typedef Eigen::Block<int, 4, 4> Block4i;
typedef Eigen::Block<int, 4, 5> Block4x5i;
typedef Eigen::Block<int, 4, 6> Block4x6i;
typedef Eigen::Block<int, 4, 7> Block4x7i;
typedef Eigen::Block<int, 4, 8> Block4x8i;
typedef Eigen::Block<int, 4, 9> Block4x9i;
typedef Eigen::Block<int, 5, 2> Block5x2i;
typedef Eigen::Block<int, 5, 3> Block5x3i;
typedef Eigen::Block<int, 5, 4> Block5x4i;
typedef Eigen::Block<int, 5, 5> Block5i;
typedef Eigen::Block<int, 5, 6> Block5x6i;
typedef Eigen::Block<int, 5, 7> Block5x7i;
typedef Eigen::Block<int, 5, 8> Block5x8i;
typedef Eigen::Block<int, 5, 9> Block5x9i;
typedef Eigen::Block<int, 6, 2> Block6x2i;
typedef Eigen::Block<int, 6, 3> Block6x3i;
typedef Eigen::Block<int, 6, 4> Block6x4i;
typedef Eigen::Block<int, 6, 5> Block6x5i;
typedef Eigen::Block<int, 6, 6> Block6i;
typedef Eigen::Block<int, 6, 7> Block6x7i;
typedef Eigen::Block<int, 6, 8> Block6x8i;
typedef Eigen::Block<int, 6, 9> Block6x9i;
typedef Eigen::Block<int, 7, 2> Block7x2i;
typedef Eigen::Block<int, 7, 3> Block7x3i;
typedef Eigen::Block<int, 7, 4> Block7x4i;
typedef Eigen::Block<int, 7, 5> Block7x5i;
typedef Eigen::Block<int, 7, 6> Block7x6i;
typedef Eigen::Block<int, 7, 7> Block7i;
typedef Eigen::Block<int, 7, 8> Block7x8i;
typedef Eigen::Block<int, 7, 9> Block7x9i;
typedef Eigen::Block<int, 8, 2> Block8x2i;
typedef Eigen::Block<int, 8, 3> Block8x3i;
typedef Eigen::Block<int, 8, 4> Block8x4i;
typedef Eigen::Block<int, 8, 5> Block8x5i;
typedef Eigen::Block<int, 8, 6> Block8x6i;
typedef Eigen::Block<int, 8, 7> Block8x7i;
typedef Eigen::Block<int, 8, 8> Block8i;
typedef Eigen::Block<int, 8, 9> Block8x9i;
typedef Eigen::Block<int, 9, 2> Block9x2i;
typedef Eigen::Block<int, 9, 3> Block9x3i;
typedef Eigen::Block<int, 9, 4> Block9x4i;
typedef Eigen::Block<int, 9, 5> Block9x5i;
typedef Eigen::Block<int, 9, 6> Block9x6i;
typedef Eigen::Block<int, 9, 7> Block9x7i;
typedef Eigen::Block<int, 9, 8> Block9x8i;
typedef Eigen::Block<int, 9, 9> Block9i;

typedef Eigen::Matrix<int, 2, 2> Matrix2i;
typedef Eigen::Matrix<int, 2, 3> Matrix2x3i;
typedef Eigen::Matrix<int, 2, 4> Matrix2x4i;
typedef Eigen::Matrix<int, 2, 5> Matrix2x5i;
typedef Eigen::Matrix<int, 2, 6> Matrix2x6i;
typedef Eigen::Matrix<int, 2, 7> Matrix2x7i;
typedef Eigen::Matrix<int, 2, 8> Matrix2x8i;
typedef Eigen::Matrix<int, 2, 9> Matrix2x9i;
typedef Eigen::Matrix<int, 3, 2> Matrix3x2i;
typedef Eigen::Matrix<int, 3, 3> Matrix3i;
typedef Eigen::Matrix<int, 3, 4> Matrix3x4i;
typedef Eigen::Matrix<int, 3, 5> Matrix3x5i;
typedef Eigen::Matrix<int, 3, 6> Matrix3x6i;
typedef Eigen::Matrix<int, 3, 7> Matrix3x7i;
typedef Eigen::Matrix<int, 3, 8> Matrix3x8i;
typedef Eigen::Matrix<int, 3, 9> Matrix3x9i;
typedef Eigen::Matrix<int, 4, 2> Matrix4x2i;
typedef Eigen::Matrix<int, 4, 3> Matrix4x3i;
typedef Eigen::Matrix<int, 4, 4> Matrix4i;
typedef Eigen::Matrix<int, 4, 5> Matrix4x5i;
typedef Eigen::Matrix<int, 4, 6> Matrix4x6i;
typedef Eigen::Matrix<int, 4, 7> Matrix4x7i;
typedef Eigen::Matrix<int, 4, 8> Matrix4x8i;
typedef Eigen::Matrix<int, 4, 9> Matrix4x9i;
typedef Eigen::Matrix<int, 5, 2> Matrix5x2i;
typedef Eigen::Matrix<int, 5, 3> Matrix5x3i;
typedef Eigen::Matrix<int, 5, 4> Matrix5x4i;
typedef Eigen::Matrix<int, 5, 5> Matrix5i;
typedef Eigen::Matrix<int, 5, 6> Matrix5x6i;
typedef Eigen::Matrix<int, 5, 7> Matrix5x7i;
typedef Eigen::Matrix<int, 5, 8> Matrix5x8i;
typedef Eigen::Matrix<int, 5, 9> Matrix5x9i;
typedef Eigen::Matrix<int, 6, 2> Matrix6x2i;
typedef Eigen::Matrix<int, 6, 3> Matrix6x3i;
typedef Eigen::Matrix<int, 6, 4> Matrix6x4i;
typedef Eigen::Matrix<int, 6, 5> Matrix6x5i;
typedef Eigen::Matrix<int, 6, 6> Matrix6i;
typedef Eigen::Matrix<int, 6, 7> Matrix6x7i;
typedef Eigen::Matrix<int, 6, 8> Matrix6x8i;
typedef Eigen::Matrix<int, 6, 9> Matrix6x9i;
typedef Eigen::Matrix<int, 7, 2> Matrix7x2i;
typedef Eigen::Matrix<int, 7, 3> Matrix7x3i;
typedef Eigen::Matrix<int, 7, 4> Matrix7x4i;
typedef Eigen::Matrix<int, 7, 5> Matrix7x5i;
typedef Eigen::Matrix<int, 7, 6> Matrix7x6i;
typedef Eigen::Matrix<int, 7, 7> Matrix7i;
typedef Eigen::Matrix<int, 7, 8> Matrix7x8i;
typedef Eigen::Matrix<int, 7, 9> Matrix7x9i;
typedef Eigen::Matrix<int, 8, 2> Matrix8x2i;
typedef Eigen::Matrix<int, 8, 3> Matrix8x3i;
typedef Eigen::Matrix<int, 8, 4> Matrix8x4i;
typedef Eigen::Matrix<int, 8, 5> Matrix8x5i;
typedef Eigen::Matrix<int, 8, 6> Matrix8x6i;
typedef Eigen::Matrix<int, 8, 7> Matrix8x7i;
typedef Eigen::Matrix<int, 8, 8> Matrix8i;
typedef Eigen::Matrix<int, 8, 9> Matrix8x9i;
typedef Eigen::Matrix<int, 9, 2> Matrix9x2i;
typedef Eigen::Matrix<int, 9, 3> Matrix9x3i;
typedef Eigen::Matrix<int, 9, 4> Matrix9x4i;
typedef Eigen::Matrix<int, 9, 5> Matrix9x5i;
typedef Eigen::Matrix<int, 9, 6> Matrix9x6i;
typedef Eigen::Matrix<int, 9, 7> Matrix9x7i;
typedef Eigen::Matrix<int, 9, 8> Matrix9x8i;
typedef Eigen::Matrix<int, 9, 9> Matrix9i;
typedef Eigen::Matrix<int, 27, 2> Matrix27x2i;
typedef Eigen::Matrix<int, 36, 2> Matrix36x2i;
typedef Eigen::Matrix<int, 64, 2> Matrix64x2i;
typedef Eigen::Matrix<int, 125, 2> Matrix125x2i;
typedef Eigen::Matrix<int, 27, 3> Matrix27x3i;
typedef Eigen::Matrix<int, 36, 3> Matrix36x3i;
typedef Eigen::Matrix<int, 64, 3> Matrix64x3i;
typedef Eigen::Matrix<int, 125, 3> Matrix125x3i;

typedef Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;
typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> MatrixXf;
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> MatrixXi;

typedef Eigen::SparseMatrix<scalar, Eigen::ColMajor> SparseXs;
typedef Eigen::SparseMatrix<scalar, Eigen::RowMajor> SparseRXs;
typedef Eigen::Triplet<scalar> Triplets;
typedef std::vector< Triplets > TripletXs;
typedef Eigen::Triplet<int> Tripleti;
typedef std::vector< Tripleti > TripletXi;

template <unsigned N>
using Vectors = Eigen::Matrix<scalar, N, 1>;

template <unsigned N>
using Vectori = Eigen::Matrix<int, N, 1>;

template <unsigned N>
using Vectorf = Eigen::Matrix<float, N, 1>;

template <unsigned N>
using VectorsT = Eigen::Matrix<scalar, 1, N>;

template <unsigned N>
using VectoriT = Eigen::Matrix<int, 1, N>;

template <unsigned N>
using VectorfT = Eigen::Matrix<float, 1, N>;

template <unsigned N>
using Matrixs = Eigen::Matrix<scalar, N, N>;

template <unsigned N>
using Matrixi = Eigen::Matrix<int, N, N>;

template <unsigned N>
using Matrixf = Eigen::Matrix<float, N, N>;

template <unsigned N>
using Blocks = Eigen::Block<scalar, N, N>;

template <unsigned N>
using Blocki = Eigen::Block<int, N, N>;

template <unsigned N>
using Blockf = Eigen::Block<float, N, N>;
//typedef Matrix<int, 1, 2> RowVector2i;

template<int DIM>
struct int_Vectors {
	int i;
	Vectors<DIM> v;
	scalar d;

	int_Vectors(int i_, const Vectors<DIM>& v_) :
		i(i_), v(v_), d(v_.norm())
	{}

	inline bool operator() (const int_Vectors& struct1, const int_Vectors& struct2)
	{
		return (struct1.d < struct2.d);
	}
};

template<int DIM>
struct int_Vectors_scalar {
	int i;
	Vectors<DIM> v;
	scalar d;
	scalar eta;

	int_Vectors_scalar(int i_, const Vectors<DIM>& v_, const scalar& eta_) :
		i(i_), v(v_), d(v_.norm()), eta(eta_)
	{}

	bool operator<(const int_Vectors_scalar& b) const
	{
		return d < b.d;
	}
};

#endif
