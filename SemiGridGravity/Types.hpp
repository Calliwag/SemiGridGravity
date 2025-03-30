#pragma once

#include <cstring>
#include <algorithm>

constexpr double T_Epsilon = 1e-10;

template<typename T>
inline T lerp(T val1, T val2, double part)
{
	return (1 - part) * val1 + part * val2;
}

template<int N, typename T>
class Vec
{
	T vec[N];

	Vec() = default;
	Vec(const Vec<N, T>& other) { memcpy(vec, other.vec, N * sizeof(T)); };
	T& operator[](int n) { return vec[n]; };
	const T& operator[](int n) const { return vec[n]; };
};

typedef Vec<2, double> Vec2d;
typedef Vec<2, int> Vec2i;
template<typename T>
class Vec<2, T>
{
public:
	T x;
	T y;
	T& operator[](int n) { return n == 0 ? x : y; };
	const T& operator[](int n) const { return n == 0 ? x : y; };
};

template<typename T>
class Rect
{
public:
	Vec<2, T> c1;
	Vec<2, T> c2;

	inline bool CheckBound(Vec<2, T> point) const { return (point.x >= c1.x) && (point.y >= c1.y) && (point.x <= c2.x) && (point.y <= c2.y); }
};

typedef Rect<int> RectI;
typedef Rect<double> RectD;

template <typename T>
class Grid
{
public:
	T* arr;
	int X;
	int Y;

	Grid() { arr = new T[1](); X = 0; Y = 0; };
	Grid(int x, int y) : X(x), Y(y)
	{
		arr = new T[X * Y];
		std::fill(arr, arr + X * Y, T());
	};
	~Grid() { delete[] arr; };
	Grid(const Grid& other) : X(other.X), Y(other.Y)
	{
		arr = new T[X * Y];
		std::copy(other.arr, other.arr + X * Y, arr);
	};
	Grid& operator=(const Grid<T>& other)
	{
		if (this != &other)
		{
			delete[] arr;
			X = other.X;
			Y = other.Y;
			arr = new T[X * Y];
			std::copy(other.arr, other.arr + (X * Y), arr);
		}
		return *this;
	}
	T& operator[](int x, int y) { return arr[y * X + x]; };
	T& operator[](Vec2i vec) { return operator[](vec.x, vec.y); };
	T TryGet(Vec2i vec) const {
		if (CheckBound(vec)) return arr[vec.y * X + vec.x];
		else return T();
	};
	bool CheckBound(int x, int y) const { return (x >= 0) && (y >= 0) && (x < X) && (y < Y); };
	bool CheckBound(Vec2i vec) const { return CheckBound(vec.x, vec.y); };
};

template<typename T>
Grid<T> Convolution(const Grid<T>& image, Grid<T> kernel);

template<int N, typename T>
inline Vec<N, T> operator+(const Vec<N, T>& vec1, const Vec<N, T> vec2)
{
	Vec<N, T> ret;
	for (int n = 0; n < N; n++)
	{
		ret[n] = vec1[n] + vec2[n];
	}
	return ret;
}

template<typename T>
inline Vec<2, T> operator+(const Vec<2, T>& vec1, const Vec<2, T> vec2)
{
	return {vec1[0] + vec2[0],vec1[1] + vec2[1]};
}

template<int N, typename T>
inline Vec<N, T>& operator+=(Vec<N, T>& vec1, const Vec<N, T> vec2)
{
	vec1 = vec1 + vec2;
	return vec1;
}

template<typename T>
inline Vec<2, T>& operator+=(Vec<2, T>& vec1, const Vec<2, T> vec2)
{
	vec1 = { vec1[0] + vec2[0],vec1[1] + vec2[1] };
	return vec1;
}

template<int N, typename T>
inline Vec<N, T> operator-(const Vec<N, T>& vec1, const Vec<N, T> vec2)
{
	Vec<N, T> ret;
	for (int n = 0; n < N; n++)
	{
		ret[n] = vec1[n] - vec2[n];
	}
	return ret;
}

template<int N, typename T>
inline Vec<N, T>& operator-=(Vec<N, T>& vec1, const Vec<N, T> vec2)
{
	vec1 = vec1 - vec2;
	return vec1;
}

template<typename T>
inline Vec<2, T>& operator-=(Vec<2, T>& vec1, const Vec<2, T> vec2)
{
	vec1 = { vec1[0] - vec2[0],vec1[1] - vec2[1] };
	return vec1;
}

template<int N, typename T>
inline Vec<N, T> operator-(const Vec<N,T>& vec)
{
	Vec<N, T> ret;
	for (int n = 0; n < N; n++)
	{
		ret[n] = -vec[n];
	}
	return ret;
}

template<int N, typename T>
inline Vec<N, T> operator*(const Vec<N,T>& vec, double scale)
{
	Vec<N, T> ret;
	for (int n = 0; n < N; n++)
	{
		ret[n] = vec[n] * scale;
	}
	return ret;
}

template<typename T>
inline Vec<2, T> operator*(const Vec<2, T>& vec, double scale)
{
	return { (T)(vec[0] * scale),(T)(vec[1] * scale) };
}

template<int N, typename T>
inline Vec<N, T> operator*(double scale, const Vec<N, T>& vec)
{
	return vec * scale;
}

template<int N, typename T>
inline Vec<N, T>& operator*=(Vec<N, T>& vec, double scale)
{
	vec = vec * scale;
	return vec;
}

template<int N, typename T>
inline double Mag(Vec<N, T> vec)
{
	double mag = 0;
	for (int n = 0; n < N; n++)
	{
		mag += (double)(vec[n] * vec[n]);
	}
	mag = sqrt(mag);
	return mag;
}

template<typename T>
inline double Mag(Vec<2, T> vec)
{
	return sqrt(vec[0] * vec[0] + vec[1] * vec[1]);
}

template<int N, typename T>
inline Vec<N, T> Norm(Vec<N, T> vec)
{
	double mag = 0;
	for (int n = 0; n < N; n++)
	{
		mag += (double)(vec[n] * vec[n]);
	}
	mag = sqrt(mag);
	return vec * (1 / mag);
}

template<typename T>
inline Vec<2, T> Norm(Vec<2, T> vec)
{
	double mag = sqrt(vec[0] * vec[0] + vec[1] * vec[1]);
	if (mag < T_Epsilon) return { 0, 0 };
	return vec * (1 / mag);
}

template<int N, typename T>
inline double Dot(Vec<N, T> vec1, Vec<N, T> vec2)
{
	double dot = 0;
	for (int n = 0; n < N; n++)
	{
		dot += vec1[n] * vec2[n];
	}
	return dot;
}

template<typename T>
inline Grid<T> Convolution(const Grid<T>& image, Grid<T> kernel)
{
	Grid<T> out(image.X,image.Y);

	if (kernel.X % 2 != 1 || kernel.Y % 2 != 1)
		throw "Kernel dimension must be odd numbers";

	for (int X = 0; X < out.X; X++)
	{
		for (int Y = 0; Y < out.Y; Y++)
		{
			out[{X, Y}] = T();
			for (int x = 0; x < kernel.X; x++)
			{
				for (int y = 0; y < kernel.Y; y++)
				{
					out[{X, Y}] +=
						kernel[{x, y}]
						* image.TryGet({(int)(x + X - (kernel.X - 1) / 2.0), (int)(y + Y - (kernel.Y - 1) / 2.0)});
				}
			}
		}
	}

	return out;
}
