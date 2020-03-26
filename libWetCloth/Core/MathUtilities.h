//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef MATH_UTILITIES_H
#define MATH_UTILITIES_H

#include <Eigen/Core>
#include "MathDefs.h"
#include <iostream>
#include <random>
#include <string>

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491298336733624406566430860213949463952247371907021798609437027705392171762931767523846748184676694051320005681271452635608277857713427577896091736371787214684409012249534301465495853710507922796892589235420199561121290219608640344181598136297747713099605187072113499999983729780499510597317328160963185950244594553469083026425223082533446850352619311881710100031378387528865875332083814206171776691473035982534904287554687311595628638823537875937519577818577805321712268066130019278766111959092164201989

#ifndef M_PI
const double M_PI = PI;
#endif

#ifdef WIN32
#undef min
#undef max
#endif

using std::min;
using std::max;
using std::swap;

namespace mathutils
{

template<typename T>
inline void check_isnan(const char* name, const std::vector< Eigen::Matrix<T, Eigen::Dynamic, 1> >& v)
{
	int i = 0;
	for (const auto& vec : v) {
		if (std::isnan((double) vec.sum()))
		{
			std::cerr << "NAN in " << name << " [" << i << "]:" << std::endl;
			std::cerr << vec << std::endl;
			exit(EXIT_FAILURE);
		}
		++i;
	}
}

template<typename T>
inline void check_isnan(const char* name, const std::vector< Eigen::Matrix<T, Eigen::Dynamic, 1> >& vx, const std::vector< Eigen::Matrix<T, Eigen::Dynamic, 1> >& vy, const std::vector< Eigen::Matrix<T, Eigen::Dynamic, 1> >& vz)
{
	std::string name0(name);
	name0 += ", x";
	check_isnan(name0.c_str(), vx);

	std::string name1(name);
	name1 += ", y";
	check_isnan(name1.c_str(), vy);

	std::string name2(name);
	name2 += ", z";
	check_isnan(name2.c_str(), vz);
}


template<typename T>
inline void check_isnan(const char* name, const Eigen::Matrix<T, Eigen::Dynamic, 1>& v)
{
	if (std::isnan((double) v.sum()))
	{
		std::cerr << "NAN in " << name << std::endl;
		std::cerr << v << std::endl;
		exit(EXIT_FAILURE);
	}
}

template<class S, class T>
inline S lerp(const S& value0, const S& value1, T f)
{ return (1 - f) * value0 + f * value1; }

template<class S, class T>
inline S bilerp(const S& v00, const S& v10,
                const S& v01, const S& v11,
                T fx, T fy)
{
	return lerp(lerp(v00, v10, fx),
	            lerp(v01, v11, fx),
	            fy);
}

template<class S, class T>
inline S trilerp(const S& v000, const S& v100,
                 const S& v010, const S& v110,
                 const S& v001, const S& v101,
                 const S& v011, const S& v111,
                 T fx, T fy, T fz)
{
	return lerp(bilerp(v000, v100, v010, v110, fx, fy),
	            bilerp(v001, v101, v011, v111, fx, fy),
	            fz);
}

template<class T>
void zero(std::vector<T>& v)
{ for (int i = (int)v.size() - 1; i >= 0; --i) v[i] = 0; }

template<class T>
T abs_max(const std::vector<T>& v)
{
	T m = 0;
	for (int i = (int)v.size() - 1; i >= 0; --i) {
		if (std::fabs(v[i]) > m)
			m = std::fabs(v[i]);
	}
	return m;
}

template<class T>
bool contains(const std::vector<T>& a, T e)
{
	for (unsigned int i = 0; i < a.size(); ++i)
		if (a[i] == e) return true;
	return false;
}

template<class T>
void add_unique(std::vector<T>& a, T e)
{
	for (unsigned int i = 0; i < a.size(); ++i)
		if (a[i] == e) return;
	a.push_back(e);
}

template<class T>
void insert(std::vector<T>& a, unsigned int index, T e)
{
	a.push_back(a.back());
	for (unsigned int i = (unsigned int)a.size() - 1; i > index; --i)
		a[i] = a[i - 1];
	a[index] = e;
}
template<class T>
void erase(std::vector<T>& a, unsigned int index)
{
	for (unsigned int i = index; i < a.size() - 1; ++i)
		a[i] = a[i + 1];
	a.pop_back();
}

template<class T>
void erase_swap(std::vector<T>& a, unsigned int index)
{
	for (unsigned int i = index; i < a.size() - 1; ++i)
		swap(a[i], a[i + 1]);
	a.pop_back();
}

template<class T>
void erase_unordered(std::vector<T>& a, unsigned int index)
{
	a[index] = a.back();
	a.pop_back();
}

template<class T>
void erase_unordered_swap(std::vector<T>& a, unsigned int index)
{
	swap(a[index], a.back());
	a.pop_back();
}

template<class T>
void find_and_erase_unordered(std::vector<T>& a, const T& doomed_element)
{
	for (unsigned int i = 0; i < a.size(); ++i)
		if (a[i] == doomed_element) {
			erase_unordered(a, i);
			return;
		}
}

template<class T>
void replace_once(std::vector<T>& a, const T& old_element, const T& new_element)
{
	for (unsigned int i = 0; i < a.size(); ++i)
		if (a[i] == old_element) {
			a[i] = new_element;
			return;
		}
}
template<class T>
inline Eigen::Matrix<T, 3, 1> grad_trilerp(const T& v000, const T& v100,
        const T& v010, const T& v110,
        const T& v001, const T& v101,
        const T& v011, const T& v111,
        T fx, T fy, T fz)
{
	return
	    Eigen::Matrix<T, 3, 1>(-(fy - 1.) * (fz - 1.), -(fx - 1.) * (fz - 1.), -(fx - 1.) * (fy - 1.)) * v000 +
	    Eigen::Matrix<T, 3, 1>((fy - 1.) * (fz - 1.), fx * (fz - 1.), fx * (fy - 1.)) * v100 +
	    Eigen::Matrix<T, 3, 1>(fy * (fz - 1.), (fx - 1.) * (fz - 1.), fy * (fx - 1.) ) * v010 +
	    Eigen::Matrix<T, 3, 1>(-fy * (fz - 1.), -fx * (fz - 1.), -fx * fy ) * v110 +
	    Eigen::Matrix<T, 3, 1>(fz * (fy - 1.), fz * (fx - 1.), (fx - 1.) * (fy - 1.)) * v001 +
	    Eigen::Matrix<T, 3, 1>(-fz * (fy - 1.), -fx * fz, -fx * (fy - 1.)) * v101 +
	    Eigen::Matrix<T, 3, 1>(-fy * fz, -fz * (fx - 1.), -fy * (fx - 1.) ) * v011 +
	    Eigen::Matrix<T, 3, 1>(fy * fz, fx * fz, fx * fy ) * v111;
}

template<class T>
inline void get_barycentric(T x, int& i, T& f, int i_low, int i_high)
{
	T s = std::floor(x);
	i = (int)s;
	if (i < i_low) {
		i = i_low;
		f = 0;
	} else if (i > i_high - 2) {
		i = i_high - 2;
		f = 1;
	} else
		f = (T)(x - s);
}

template<typename S, int K>
inline void swap(Eigen::Matrix<S, Eigen::Dynamic, 1>& v, int i, int j)
{
	Eigen::Matrix<S, K, 1> c = v.template segment<K>(i * K);
	v.template segment<K>(i * K) = v.template segment<K>(j * K);
	v.template segment<K>(j * K) = c;
}

template<typename S, int K>
inline void swap(Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic>& v, int i, int j)
{
	Eigen::Matrix<S, K, K> c = v.template block<K, K>(i * K, 0);
	v.template block<K, K>(i * K, 0) = v.template block<K, K>(j * K, 0);
	v.template block<K, K>(j * K, 0) = c;
}

inline void fisherYates( int n, std::vector<int>& indices )
{
	indices.resize(n);
	for (int i = 0; i < n; ++i) indices[i] = i;
	for (int i = n - 1; i >= 1; --i) {
		const int j = rand() % (i + 1);
		std::swap(indices[i], indices[j]);
	}
}

bool approxSymmetric( const MatrixXs& A, const scalar& eps );

scalar scalarRand(const scalar min, const scalar max);

inline scalar perimeter(const scalar& ra, const scalar& rb)
{
	// The second-Ramanujan approximation to ellipse perimeter
	const scalar h = (ra - rb) * (ra - rb) / ((ra + rb) * (ra + rb));
	return PI * (ra + rb) * (1. + 3. * h / (10. + sqrt(4. - 3. * h)));
}

inline int mod_floor(int a, int n) {
	return ((a % n) + n) % n;
}

inline scalar sgn(const scalar& x)
{
	return x == 0.0 ? 0.0 : (x > 0.0 ? 1.0 : -1.0);
}

inline scalar mean_curvature(const Vector9s& h, const scalar& dx)
{
	const scalar t2 = 1.0 / (dx * dx);
	const scalar t3 = h[1] - h[7];
	const scalar t7 = h[4] * 2.0;
	const scalar t4 = h[3] + h[5] - t7;
	const scalar t5 = 1.0 / (dx * dx * dx * dx);
	const scalar t6 = h[3] - h[5];
	const scalar t8 = dx * dx;
	const scalar t0 = pow(t8, 3.0 / 2.0) * (t2 * t4 + t2 * (h[1] - h[4] * 2.0 + h[7]) + (t3 * t3) * t4 * t5 * (1.0 / 4.0) + t5 * (t6 * t6) * (h[1] + h[7] - t7) * (1.0 / 4.0) - t3 * t5 * t6 * (h[0] - h[2] - h[6] + h[8]) * (1.0 / 8.0)) * 1.0 / pow(t8 * 4.0 - h[1] * h[7] * 2.0 - h[3] * h[5] * 2.0 + h[1] * h[1] + h[3] * h[3] + h[5] * h[5] + h[7] * h[7], 3.0 / 2.0) * 8.0;
	return t0;
}
//Given two signed distance values (line endpoints), determine what fraction of a connecting segment is "inside"
inline scalar fraction_inside(const scalar& phi_left, const scalar& phi_right) {
	if (phi_left < 0 && phi_right < 0)
		return 1;
	if (phi_left < 0 && phi_right >= 0)
		return phi_left / (phi_left - phi_right);
	if (phi_left >= 0 && phi_right < 0)
		return phi_right / (phi_right - phi_left);
	else
		return 0;
}

inline void cycle_array(scalar* arr, int size) {
	scalar t = arr[0];
	for (int i = 0; i < size - 1; ++i)
		arr[i] = arr[i + 1];
	arr[size - 1] = t;
}

inline scalar fraction_inside(const scalar& phi_bl, const scalar& phi_br,
                              const scalar& phi_tl, const scalar& phi_tr) {

	int inside_count = (phi_bl < 0 ? 1 : 0) + (phi_tl < 0 ? 1 : 0) + (phi_br < 0 ? 1 : 0) + (phi_tr < 0 ? 1 : 0);
	scalar list[] = { phi_bl, phi_br, phi_tr, phi_tl };

	if (inside_count == 4)
		return 1.;
	else if (inside_count == 3) {
		//rotate until the positive value is in the first position
		while (list[0] < 0) {
			cycle_array(list, 4);
		}

		//Work out the area of the exterior triangle
		scalar side0 = 1. - fraction_inside(list[0], list[3]);
		scalar side1 = 1. - fraction_inside(list[0], list[1]);
		return 1. - 0.5 * side0 * side1;
	}
	else if (inside_count == 2) {

		//rotate until a negative value is in the first position, and the next negative is in either slot 1 or 2.
		while (list[0] >= 0 || !(list[1] < 0 || list[2] < 0)) {
			cycle_array(list, 4);
		}

		if (list[1] < 0) { //the matching signs are adjacent
			scalar side_left = fraction_inside(list[0], list[3]);
			scalar side_right = fraction_inside(list[1], list[2]);
			return  0.5 * (side_left + side_right);
		}
		else { //matching signs are diagonally opposite
			//determine the centre point's sign to disambiguate this case
			scalar middle_point = 0.25 * (list[0] + list[1] + list[2] + list[3]);
			if (middle_point < 0) {
				scalar area = 0.;

				//first triangle (top left)
				scalar side1 = 1. - fraction_inside(list[0], list[3]);
				scalar side3 = 1. - fraction_inside(list[2], list[3]);

				area += 0.5 * side1 * side3;

				//second triangle (top right)
				scalar side2 = 1 - fraction_inside(list[2], list[1]);
				scalar side0 = 1 - fraction_inside(list[0], list[1]);
				area += 0.5 * side0 * side2;

				return 1. - area;
			}
			else {
				scalar area = 0.;

				//first triangle (bottom left)
				scalar side0 = fraction_inside(list[0], list[1]);
				scalar side1 = fraction_inside(list[0], list[3]);
				area += 0.5 * side0 * side1;

				//second triangle (top right)
				scalar side2 = fraction_inside(list[2], list[1]);
				scalar side3 = fraction_inside(list[2], list[3]);
				area += 0.5 * side2 * side3;
				return area;
			}

		}
	}
	else if (inside_count == 1) {
		//rotate until the negative value is in the first position
		while (list[0] >= 0) {
			cycle_array(list, 4);
		}

		//Work out the area of the interior triangle, and subtract from 1.
		scalar side0 = fraction_inside(list[0], list[3]);
		scalar side1 = fraction_inside(list[0], list[1]);
		return 0.5 * side0 * side1;
	}
	else
		return 0;

}

template<typename S, unsigned N, unsigned M>
inline Eigen::Matrix<S, N, M> frac(const Eigen::Matrix<S, N, M>& x)
{
	Eigen::Matrix<S, N, M> ret(x);
	for (unsigned i = 0; i < N; ++i) for (unsigned j = 0; j < M; ++j) ret(i, j) -= floor(ret(i, j));
	return ret;
}

inline scalar sqr(const scalar& s)
{
	return s * s;
}

inline Vector3s cross_x(const Vector3s& a, const scalar& b)
{
	return Vector3s(0.0, a(2) * b, -a(1) * b);
}

inline Vector3s cross_y(const Vector3s& a, const scalar& b)
{
	return Vector3s(-a(2) * b, 0.0, a(0) * b);
}

inline Vector3s cross_z(const Vector3s& a, const scalar& b)
{
	return Vector3s(a(1) * b, -a(0) * b, 0.0);
}

inline scalar cross_x_row(const Vector3s& a, const Vector3s& b)
{
	return -a(2) * b(1) + a(1) * b(2);
}

inline scalar cross_y_row(const Vector3s& a, const Vector3s& b)
{
	return a(2) * b(0) - a(0) * b(2);
}

inline scalar cross_z_row(const Vector3s& a, const Vector3s& b)
{
	return -a(1) * b(0) + a(0) * b(1);
}

template<typename T>
inline T smooth_kernel( const T& r2, const T& h ) {
	return std::max( (T) pow((T) 1.0 - r2 / (h * h), (T) 3.0), (T) 0.0 );
}

inline scalar quad_kernel(const scalar& x) {
	const scalar dx = fabs(x);
	if (dx < 0.5) return 0.75 - dx * dx;
	else if (dx < 1.5) return 0.5 * (1.5 - dx) * (1.5 - dx);
	else return 0.0;
}

inline scalar cubic_kernel(const scalar& x) {
	const scalar dx = fabs(x);
	if (dx < 1.0) return 0.5 * dx * dx * dx - dx * dx + 2.0 / 3.0;
	else if (dx < 2.0) return (2.0 - dx) * (2.0 - dx) * (2.0 - dx) / 6.0;
	else return 0.0;
}

inline scalar linear_kernel(const scalar& x) {
	const scalar dx = fabs(x);
	if (dx < 1.0) return 1.0 - dx;
	else return 0.0;
}

inline scalar grad_quad_kernel(const scalar& x) {
	const scalar dx = fabs(x);
	if (dx < 0.5) return -2.0 * dx * sgn(x);
	else if (dx < 1.5) return sgn(x) * (dx - 1.5);
	else return 0.0;
}

inline scalar grad_cubic_kernel(const scalar& x) {
	const scalar dx = fabs(x);
	if (dx < 1.0) return (x * sgn(x) * (3.0 * x - 4.0 * sgn(x))) / 2.0;
	else if (dx < 2.0) return -(sgn(x) * (dx - 2.0) * (dx - 2.0)) / 2.0;
	else return 0.0;
}

inline scalar grad_linear_kernel(const scalar& x) {
	const scalar dx = fabs(x);
	if (dx < 1.0) return -sgn(x);
	else return 0.0;
}

template<int order>
inline scalar N_kernel(const Vector3s& x) {
	if (order == 2)
		return quad_kernel(x(0)) * quad_kernel(x(1)) * quad_kernel(x(2));
	else if (order == 3)
		return cubic_kernel(x(0)) * cubic_kernel(x(1)) * cubic_kernel(x(2));
	else
		return linear_kernel(x(0)) * linear_kernel(x(1)) * linear_kernel(x(2));
}

template<int order>
inline scalar grad_N_kernel(const Vector3s& x, const scalar& h, Vector3s& g) {
	if (order == 1) {
		const scalar v0 = linear_kernel(x(0));
		const scalar v1 = linear_kernel(x(1));
		const scalar v2 = linear_kernel(x(2));
		const scalar g0 = grad_linear_kernel(x(0));
		const scalar g1 = grad_linear_kernel(x(1));
		const scalar g2 = grad_linear_kernel(x(2));

		g(0) = g0 * v1 * v2;
		g(1) = v0 * g1 * v2;
		g(2) = v0 * v1 * g2;
		g *= 1.0 / h;
		return v0 * v1 * v2;
	} else if (order == 2) {
		const scalar inv_h = 1. / h;
		const scalar dx0 = fabs(x(0));
		const scalar dx1 = fabs(x(1));
		const scalar dx2 = fabs(x(2));
		const scalar sx0 = sgn(x(0));
		const scalar sx1 = sgn(x(1));
		const scalar sx2 = sgn(x(2));

		const scalar v0 = dx0 < 0.5 ? (0.75 - dx0 * dx0) : (dx0 < 1.5 ? (0.5 * (1.5 - dx0) * (1.5 - dx0)) : 0.0);
		const scalar v1 = dx1 < 0.5 ? (0.75 - dx1 * dx1) : (dx1 < 1.5 ? (0.5 * (1.5 - dx1) * (1.5 - dx1)) : 0.0);
		const scalar v2 = dx2 < 0.5 ? (0.75 - dx2 * dx2) : (dx2 < 1.5 ? (0.5 * (1.5 - dx2) * (1.5 - dx2)) : 0.0);

		const scalar g0 = dx0 < 0.5 ? (-2.0 * dx0 * sx0) : (dx0 < 1.5 ? sx0 * (dx0 - 1.5) : 0.0);
		const scalar g1 = dx1 < 0.5 ? (-2.0 * dx1 * sx1) : (dx1 < 1.5 ? sx1 * (dx1 - 1.5) : 0.0);
		const scalar g2 = dx2 < 0.5 ? (-2.0 * dx2 * sx2) : (dx2 < 1.5 ? sx2 * (dx2 - 1.5) : 0.0);

		g(0) = g0 * v1 * v2 * inv_h;
		g(1) = v0 * g1 * v2 * inv_h;
		g(2) = v0 * v1 * g2 * inv_h;
		return v0 * v1 * v2;
	} else {
		const scalar v0 = cubic_kernel(x(0));
		const scalar v1 = cubic_kernel(x(1));
		const scalar v2 = cubic_kernel(x(2));
		const scalar g0 = grad_cubic_kernel(x(0));
		const scalar g1 = grad_cubic_kernel(x(1));
		const scalar g2 = grad_cubic_kernel(x(2));

		g(0) = g0 * v1 * v2;
		g(1) = v0 * g1 * v2;
		g(2) = v0 * v1 * g2;
		g *= 1.0 / h;
		return v0 * v1 * v2;
	}
}

inline scalar grad_N_kernel_x(const Vector3s& x, const scalar& h, const int order) {
	std::function<scalar(const scalar&)> val = (order == 1) ? linear_kernel : ((order == 3) ? cubic_kernel : quad_kernel);
	std::function<scalar(const scalar&)> grad = (order == 1) ? grad_linear_kernel : ((order == 3) ? grad_cubic_kernel : grad_quad_kernel);

	return grad(x(0)) * val(x(1)) * val(x(2)) / h;
}

inline scalar grad_N_kernel_y(const Vector3s& x, const scalar& h, const int order) {
	std::function<scalar(const scalar&)> val = (order == 1) ? linear_kernel : ((order == 3) ? cubic_kernel : quad_kernel);
	std::function<scalar(const scalar&)> grad = (order == 1) ? grad_linear_kernel : ((order == 3) ? grad_cubic_kernel : grad_quad_kernel);

	return val(x(0)) * grad(x(1)) * val(x(2)) / h;
}

inline scalar grad_N_kernel_z(const Vector3s& x, const scalar& h, const int order) {
	std::function<scalar(const scalar&)> val = (order == 1) ? linear_kernel : ((order == 3) ? cubic_kernel : quad_kernel);
	std::function<scalar(const scalar&)> grad = (order == 1) ? grad_linear_kernel : ((order == 3) ? grad_cubic_kernel : grad_quad_kernel);

	return val(x(0)) * val(x(1)) * grad(x(2)) / h;
}

template<typename S, unsigned N>
inline void QRDecompose(const Eigen::Matrix<S, N, N>& A, Eigen::Matrix<S, N, N>& Q, Eigen::Matrix<S, N, N>& R) {
	Eigen::HouseholderQR< Eigen::Matrix<S, N, N> > qr(A);
	R = qr.matrixQR().template triangularView<Eigen::Upper>();
	MatrixXs s = VectorXs(R.diagonal().array().sign()).asDiagonal();
//        std::cout<<"R: "<<R<<std::endl;
//        std::cout<<"s:"<<s<<std::endl;
	Q = qr.householderQ() * s;
	R = s * R;
}

template<class T>
inline T clamp(T a, T lower, T upper)
{
	if (a < lower) return lower;
	else if (a > upper) return upper;
	else return a;
}

inline scalar defaultRadiusMultiplier()
{
	return 1.0 / sqrt(3.0);
}

template<typename S, typename T>
inline S lerp_weno( const S value[], T f )
{
	S p1 = value[0] + (value[1] - value[0]) * (f + 2.0) + (value[2] - 2.0 * value[1] + value[0]) * (f + 2.0) * (f + 1.0) * 0.5
	       + (value[3] - 3.0 * value[2] + 3.0 * value[1] - value[0]) * (f + 2.0) * (f + 1.0) * f / 6.0;
	S p2 = value[1] + (value[2] - value[1]) * (f + 1.0) + (value[3] - 2.0 * value[2] + value[1]) * (f + 1.0) * f * 0.5
	       + (value[4] - 3.0 * value[3] + 3.0 * value[2] - value[1]) * (f + 1.0) * f * (f - 1.0) / 6.0;
	S p3 = value[2] + (value[3] - value[2]) * f + (value[4] - 2.0 * value[3] + value[2]) * f * (f - 1.0) * 0.5
	       + (value[5] - 3.0 * value[4] + 3.0 * value[3] - value[2]) * f * (f - 1.0) * (f - 2.0) / 6.0;

	T C1 = (2 - f) * (3 - f) / 20.0;
	T C2 = (3 - f) * (f + 2) / 10.0;
	T C3 = (f + 2) * (f + 1) / 20.0;

	T IS1 = (814.0 * value[3] * value[3] + 4326 * value[2] * value[2] + 2976 * value[1] * value[1] + 244 * value[0] * value[0] - 3579 * value[2] * value[3] - 6927 * value[2] * value[1]
	         + 1854 * value[2] * value[0] + 2634 * value[3] * value[1] - 683 * value[3] * value[0] - 1659 * value[1] * value[0])   / 180.0;
	T IS2 = (1986 * value[3] * value[3] + 1986 * value[2] * value[2] + 244 * value[1] * value[1] + 244 * value[4] * value[4] + 1074 * value[2] * value[4] - 3777 * value[2] * value[3]
	         - 1269 * value[2] * value[1] + 1074 * value[3] * value[1] - 1269 * value[4] * value[3] - 293 * value[4] * value[1])  / 180.0;
	T IS3 = (814 * value[2] * value[2] + 4326 * value[3] * value[3] + 2976 * value[4] * value[4] + 244 * value[5] * value[5] - 683 * value[2] * value[5] + 2634 * value[2] * value[4]
	         - 3579 * value[2] * value[3] - 6927 * value[3] * value[4] + 1854 * value[3] * value[5] - 1659 * value[4] * value[5]) / 180.0;

	const T epsilon = 1e-6;
	T alpha1 = C1 / ((IS1 + epsilon) * (IS1 + epsilon));
	T alpha2 = C2 / ((IS2 + epsilon) * (IS2 + epsilon));
	T alpha3 = C3 / ((IS3 + epsilon) * (IS3 + epsilon));

	T sumalpha = alpha1 + alpha2 + alpha3;
	T w1 = alpha1 / sumalpha;
	T w2 = alpha2 / sumalpha;
	T w3 = alpha3 / sumalpha;

	return p1 * w1 + p2 * w2 + p3 * w3;
}

inline scalar cosine_ease_function(const scalar& t, const scalar& t0, const scalar& t1, const scalar& ta, const scalar& tb, const scalar& amp, const scalar& freq)
{
	const scalar Ta = ta - t0;
	const scalar Tt = t - t0;
	const scalar Tb = t1 - tb;
	const scalar Te = t1 - t0;
	const scalar w = 2.0 * M_PI * freq;

	if (t < t0 || t > t1) return 0.0;
	else if (t < ta) {
		const scalar t2 = Ta * w;
		const scalar t3 = cos(t2);
		const scalar t4 = sin(t2);
		return 1.0 / (Ta * Ta * Ta) * (Tt * Tt) * (amp * t3 * 2.0 + amp * Ta * t4 * w) * -3.0 + amp * 1.0 / (Ta * Ta) * Tt * (t3 * 3.0 + Ta * t4 * w) * 2.0;
	} else if (t > tb) {
		const scalar t2 = Tb - Te;
		const scalar t3 = t2 * w;
		const scalar t4 = Te - Tt;
		const scalar t5 = 1.0 / (Tb * Tb * Tb);
		const scalar t6 = cos(t3);
		const scalar t7 = sin(t3);
		return -amp * t5 * (Te * 2.0 - Tt * 2.0) * (Tb * t6 * 3.0 - Te * t6 * 2.0 + t6 * Tt * 2.0 + (Tb * Tb) * t7 * w + Tb * t7 * w * Tt - Tb * Te * t7 * w) + amp * (t4 * t4) * t5 * (t6 * 2.0 + Tb * t7 * w);
	} else {
		return -amp * w * sin(w * Tt);
	}
}

inline scalar weno_ease_function(const scalar& t, const scalar& dt, const scalar& t0, const scalar& t1, const scalar& base_dt, const scalar& cur_pos, const std::vector<scalar>& bases)
{
	if (t < t0 || t > t1) return 0.0;

	scalar spos = (t - t0) / base_dt;
	int ipos = (int) spos;
	scalar fpos = spos - (scalar) ipos;

	int nb = (int) bases.size();

	scalar extracted_bases[] = {
		bases[clamp(ipos - 2, 0, nb - 1)],
		bases[clamp(ipos - 1, 0, nb - 1)],
		bases[clamp(ipos - 0, 0, nb - 1)],
		bases[clamp(ipos + 1, 0, nb - 1)],
		bases[clamp(ipos + 2, 0, nb - 1)],
		bases[clamp(ipos + 3, 0, nb - 1)]
	};

	scalar target_pos = lerp_weno(extracted_bases, fpos);
	return (target_pos - cur_pos) / dt;
}

inline scalar cubic_ease_function(const scalar& t, const scalar& t0, const scalar& t1, const scalar& ta, const scalar& tb, const scalar& L)
{
	scalar yh = (L * 2.0) / (t1 - t0 + tb - ta);
	if (t < t0 || t > t1) return 0.0;
	else {
		if (t < ta) return (yh * (t0 - t) * (t0 - t) * (t0 - 3.0 * ta + 2.0 * t)) / ((t0 - ta) * (t0 - ta) * (t0 - ta));
		else if (t > tb) return (yh * (t1 - t) * (t1 - t) * (t1 - 3.0 * tb + 2.0 * t)) / ((t1 - tb) * (t1 - tb) * (t1 - tb));
		else return yh;
	}
}

inline scalar inverse_D_coeff(const scalar& h, const int order) {
	switch (order) {
	case 1:
		return 1.0;
	case 3:
		return 3.0 / (h * h);
	default:
		return 4.0 / (h * h);
	}
}

inline scalar D_coeff(const scalar& h, const int order) {
	switch (order) {
	case 1:
		return 1.0;
	case 3:
		return (h * h) / 3.0;
	default:
		return (h * h) / 4.0;
	}
}

// [Botsch et al. 2010] (3.9)
inline void grad_triangle(const Vector3s& x0, const Vector3s& x1, const Vector3s& x2, Matrix3s& grad_coeff)
{
	Vector3s n = (x2 - x0).cross(x1 - x0);
	const scalar len = n.norm();
	if (len < 1e-63) {
		grad_coeff.setZero();
		return;
	}

	n /= len;

	Matrix3s mrot = Eigen::AngleAxis<scalar>(-M_PI / 2.0, n).toRotationMatrix();
	Vector3s v0 = mrot * (x0 - x2) / len;
	Vector3s v1 = mrot * (x1 - x0) / len;

	grad_coeff.block<3, 1>(0, 0) = -(v0 + v1);
	grad_coeff.block<3, 1>(0, 1) = v0;
	grad_coeff.block<3, 1>(0, 2) = v1;
}

inline void get_div_triangle(const scalar& va0, const scalar& va1, const scalar& va2, const scalar& thickness0, const scalar& thickness1, const scalar& thickness2, const Vector3s& x0, const Vector3s& x1, const Vector3s& x2, Vector3s& div0, Vector3s& div1, Vector3s& div2)
{
	Vector3s n = (x2 - x0).cross(x1 - x0);
	const scalar len = n.norm();
	if (len < 1e-63) {
		return;
	}

	n /= len;

	Matrix3s mrot = Eigen::AngleAxis<scalar>(-M_PI / 2.0, n).toRotationMatrix();

	div0 = -mrot * (x2 - x1) * 0.5 * thickness0 / va0;
	div1 = -mrot * (x0 - x2) * 0.5 * thickness1 / va1;
	div2 = -mrot * (x1 - x0) * 0.5 * thickness2 / va2;

}

inline scalar get_rotated_pore_coeff_x(const Vector3s& a, const scalar& sf)
{
	const scalar t2 = a(0) * a(0);
	const scalar t3 = a(1) * a(1);
	const scalar t4 = a(2) * a(2);
	return (t2 + t3 * 2.0 + t4 * 2.0 - a(0) * a(1) - a(0) * a(2) + sf * t2 - sf * t3 - sf * t4 + a(0) * a(1) * sf * 2.0 + a(0) * a(2) * sf * 2.0) / std::max(1e-63, t2 + t3 + t4);
}

inline scalar get_rotated_pore_coeff_y(const Vector3s& a, const scalar& sf)
{
	const scalar t2 = a(0) * a(0);
	const scalar t3 = a(1) * a(1);
	const scalar t4 = a(2) * a(2);
	return (t2 * 2.0 + t3 + t4 * 2.0 - a(0) * a(1) - a(1) * a(2) - sf * t2 + sf * t3 - sf * t4 + a(0) * a(1) * sf * 2.0 + a(1) * a(2) * sf * 2.0) / std::max(1e-63, t2 + t3 + t4);
}

inline scalar get_rotated_pore_coeff_z(const Vector3s& a, const scalar& sf)
{
	const scalar t2 = a(0) * a(0);
	const scalar t3 = a(1) * a(1);
	const scalar t4 = a(2) * a(2);
	return (t2 * 2.0 + t3 * 2.0 + t4 - a(0) * a(2) - a(1) * a(2) - sf * t2 - sf * t3 + sf * t4 + a(0) * a(2) * sf * 2.0 + a(1) * a(2) * sf * 2.0) / std::max(1e-63, t2 + t3 + t4);
}

inline scalar get_rotated_drag_x(const Vector3s& a, const Vector3s& d)
{
	const scalar t2 = a(0) * a(0);
	const scalar t3 = a(1) * a(1);
	const scalar t4 = a(2) * a(2);
	const scalar t5 = t2 + t3 + t4;
	const scalar t8 = sqrt(t5);
	const scalar t9 = t3 * t8;
	const scalar t10 = a(2) * t2;
	const scalar t6 = t9 + t10;
	const scalar t7 = 1.0 / std::max(1e-63, t5);
	const scalar t11 = t2 + t3;
	const scalar t12 = a(2) - t8;
	const scalar t13 = 1.0 / std::max(1e-63, t11 * t11);
	const scalar t14 = 1.0 / std::max(1e-63, t11);
	return d(2) * t2 * t7 + d(0) * (t6 * t6) * t7 * t13 + a(0) * a(1) * d(2) * t7 + a(0) * a(2) * d(2) * t7 - a(0) * d(0) * t6 * t7 * t14 + d(1) * t2 * t3 * t7 * (t12 * t12) * t13 - a(0) * d(1) * t3 * t7 * t12 * t14 + a(0) * a(1) * d(0) * t6 * t7 * t12 * t13 + a(0) * a(1) * d(1) * t7 * t12 * t13 * (a(2) * t3 + t2 * t8);
}

inline scalar get_rotated_drag_y(const Vector3s& a, const Vector3s& d)
{
	const scalar t2 = a(1) * a(1);
	const scalar t3 = a(0) * a(0);
	const scalar t4 = a(2) * a(2);
	const scalar t5 = t2 + t3 + t4;
	const scalar t8 = sqrt(t5);
	const scalar t9 = t3 * t8;
	const scalar t10 = a(2) * t2;
	const scalar t6 = t9 + t10;
	const scalar t7 = 1.0 / std::max(1e-63, t5);
	const scalar t11 = t2 + t3;
	const scalar t12 = a(2) - t8;
	const scalar t13 = 1.0 / std::max(1e-63, t11 * t11);
	const scalar t14 = 1.0 / std::max(1e-63, t11);
	return d(2) * t2 * t7 + d(1) * (t6 * t6) * t7 * t13 + a(0) * a(1) * d(2) * t7 + a(1) * a(2) * d(2) * t7 - a(1) * d(1) * t6 * t7 * t14 + d(0) * t2 * t3 * t7 * (t12 * t12) * t13 - a(1) * d(0) * t3 * t7 * t12 * t14 + a(0) * a(1) * d(1) * t6 * t7 * t12 * t13 + a(0) * a(1) * d(0) * t7 * t12 * t13 * (a(2) * t3 + t2 * t8);
}

inline scalar get_rotated_drag_z(const Vector3s& a, const Vector3s& d)
{
	const scalar t2 = a(0) * a(0);
	const scalar t3 = a(1) * a(1);
	const scalar t4 = a(2) * a(2);
	const scalar t5 = t2 + t3 + t4;
	const scalar t6 = sqrt(t5);
	return (d(0) * (t2 * t2) + d(1) * (t3 * t3) + d(0) * t2 * t3 + d(1) * t2 * t3 + d(2) * t2 * t4 + d(2) * t3 * t4 - a(0) * d(0) * t3 * t6 + a(1) * d(0) * t2 * t6 + a(0) * d(1) * t3 * t6 - a(1) * d(1) * t2 * t6 - a(0) * a(2) * d(0) * t2 - a(1) * a(2) * d(0) * t2 - a(0) * a(2) * d(1) * t3 - a(1) * a(2) * d(1) * t3 + a(0) * a(2) * d(2) * t2 + a(0) * a(2) * d(2) * t3 + a(1) * a(2) * d(2) * t2 + a(1) * a(2) * d(2) * t3) / std::max(1e-63, t5 * (t2 + t3));
}

template<typename Callable, typename Array3T>
inline scalar interpolate_2nd_order(const Vector3s& p, const Array3T& arr, Callable func)
{
	Vector3s frac_p = Vector3s(p(0) - floor(p(0)), p(1) - floor(p(1)), p(2) - floor(p(2)));
	Vector3i base_p = Vector3i((int) floor(p(0)) + (frac_p(0) < 0.5 ? -1 : 0),
	                           (int) floor(p(1)) + (frac_p(1) < 0.5 ? -1 : 0),
	                           (int) floor(p(2)) + (frac_p(2) < 0.5 ? -1 : 0));

	Vector3s local_p = frac_p + Vector3s( frac_p(0) < 0.5 ? 1 : 0,
	                                      frac_p(1) < 0.5 ? 1 : 0,
	                                      frac_p(2) < 0.5 ? 1 : 0);

	scalar sum_u = 0.0;

	for (int k = 0; k < 3; ++k) for (int j = 0; j < 3; ++j) for (int i = 0; i < 3; ++i)
			{
				Vector3i inode = base_p + Vector3i(i, j, k);
				if (!func(inode)) continue;

				Vector3s diff = local_p - Vector3s(i, j, k);

				sum_u += N_kernel<2>(diff) * arr(inode(0), inode(1), inode(2));
			}

	return sum_u;
}

template<typename Callable, typename Array3T>
inline scalar interpolate_2nd_order_with_grad(const Vector3s& p, Vector3sT& B, const Array3T& arr, const scalar& h, Callable func)
{
	Vector3s frac_p = Vector3s(p(0) - floor(p(0)), p(1) - floor(p(1)), p(2) - floor(p(2)));
	Vector3i base_p = Vector3i((int) floor(p(0)) + (frac_p(0) < 0.5 ? -1 : 0),
	                           (int) floor(p(1)) + (frac_p(1) < 0.5 ? -1 : 0),
	                           (int) floor(p(2)) + (frac_p(2) < 0.5 ? -1 : 0));

	Vector3s local_p = frac_p + Vector3s( frac_p(0) < 0.5 ? 1 : 0,
	                                      frac_p(1) < 0.5 ? 1 : 0,
	                                      frac_p(2) < 0.5 ? 1 : 0);

	scalar sum_u = 0.0;
	B.setZero();

	for (int k = 0; k < 3; ++k) for (int j = 0; j < 3; ++j) for (int i = 0; i < 3; ++i)
			{
				Vector3i inode = base_p + Vector3i(i, j, k);
				if (!func(inode)) continue;

				Vector3s diff = Vector3s(i, j, k) - local_p;
				const scalar w = N_kernel<2>(diff);
				const scalar arrval = arr(inode(0), inode(1), inode(2));

				sum_u += w * arrval;
				B += w * arrval * diff.transpose() * h;
			}

	return sum_u;
}

// swap so that a<b
template<class T>
inline void sort(T &a, T &b)
{
	if (a > b) std::swap(a, b);
}

// swap so that a<b<c
template<class T>
inline void sort(T &a, T &b, T &c)
{
	if (a > b) std::swap(a, b);
	if (a > c) std::swap(a, c);
	if (b > c) std::swap(b, c);
}

// swap so that a<b<c<d
template<class T>
inline void sort(T &a, T &b, T &c, T &d)
{
	if (a > b) std::swap(a, b);
	if (c > d) std::swap(c, d);
	if (a > c) std::swap(a, c);
	if (b > d) std::swap(b, d);
	if (b > c) std::swap(b, c);
}

inline Vector4s QuatVec4(const Quat4s& q)
{
	return Vector4s(q.w(), q.x(), q.y(), q.z());
}

inline void dhdr_yarn(double mu, double la, double r22, double r23, double r33, double *dhdr22, double *dhdr23, double *dhdr33) {
	Eigen::Matrix2d r(2, 2);
	r(0, 0) = r22;
	r(0, 1) = r23;
	r(1, 0) = 0.0;
	r(1, 1) = r33;

	Eigen::JacobiSVD<Eigen::Matrix2d> svd(r, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Vector2d sigm = svd.singularValues();
	Eigen::Vector2d sigm_inv = Eigen::Vector2d(1.0 / sigm(0), 1.0 / sigm(1));
	Eigen::Vector2d lnsigm = Eigen::Vector2d(log(sigm(0)), log(sigm(1)));

	if (lnsigm(0) + lnsigm(1) > 0.0) {
		*dhdr22 = 0.0;
		*dhdr23 = 0.0;
		*dhdr33 = 0.0;
	} else {
		Eigen::Matrix2d tmp = Eigen::Matrix2d(sigm_inv.asDiagonal()) * Eigen::Matrix2d(lnsigm.asDiagonal());
		Eigen::Matrix2d ret = svd.matrixU() * (2.0 * mu * tmp + la * lnsigm.sum() * Eigen::Matrix2d(sigm_inv.asDiagonal())) * svd.matrixV().transpose();

		*dhdr22 = ret(0, 0);
		*dhdr23 = ret(0, 1);
		*dhdr33 = ret(1, 1);
	}
}

inline void dgdr_cloth(double mu, double r13, double r23, double *dhdr13, double *dhdr23) {
	*dhdr13 = mu * r13;
	*dhdr23 = mu * r23;
}

inline void dhdr_cloth(double mu, double la, double r33, double *dhdr33) {
	if (r33 <= 1.0) *dhdr33 = -(2.0 * mu + la) * (1.0 - r33) * (1.0 - r33);
	else *dhdr33 = 0.0;
}

inline scalar twist_component( const Eigen::Quaternion<scalar>& rot, const Vector3s& dir )
{
	Vector3s ra = rot.vec();
	Vector3s p = dir * (ra.dot(dir));
	Eigen::Quaternion<scalar> twist(rot.w(), p(0), p(1), p(2));
	twist.normalize();
	return Eigen::AngleAxis<scalar>(twist).angle();
}

// build orthonormal basis from a unit vector
// Listing 2 in Frisvad 2012, Building an Orthonormal Basis from a 3D Unit Vector Without Normalization
inline void orthonormal_basis(const Vector3s& n, Vector3s& b1, Vector3s& b2)
{
	if (n(2) + 1.0 < 1e-12) {
		b1 = Vector3s(0.0, -1.0, 0.0);
		b2 = Vector3s(-1.0, 0.0, 0.0);
		return;
	}

	const scalar a = 1.0 / (1.0 + n(2));
	const scalar b = -n(0) * n(1) * a;
	b1 = Vector3s(1.0 - n(0) * n(0) * a, b, -n(0));
	b2 = Vector3s(b, 1.0 - n(1) * n(1) * a, -n(1));
}

void print_histogram_analysis(const std::vector< VectorXs >&, int, const std::string&, bool);
}

#endif
