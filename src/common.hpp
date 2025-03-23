#ifndef COMMON_H
#define COMMON_H

#include "math.h"

static constexpr double GRAVITY_CONSTANT = 1;

struct Vec3 {
    double x, y, z;
};

inline void operator+=(Vec3 &a, const Vec3 &b) { a.x += b.x; a.y += b.y; a.z += b.z; }
inline void operator-=(Vec3 &a, const Vec3 &b) { a.x -= b.x; a.y -= b.y; a.z -= b.z; }
inline void operator*=(Vec3 &a, double b) { a.x *= b; a.y *= b; a.z *= b; }
inline void operator/=(Vec3 &a, double b) { a.x /= b; a.y /= b; a.z /= b; }

inline Vec3 operator+(const Vec3 &a, const Vec3 &b) { return Vec3 { a.x + b.x, a.y + b.y, a.z + b.z }; }
inline Vec3 operator-(const Vec3 &a, const Vec3 &b) { return Vec3 { a.x - b.x, a.y - b.y, a.z - b.z }; }
inline Vec3 operator*(const Vec3 &a, double b) { return Vec3 { a.x * b, a.y * b, a.z * b }; }
inline Vec3 operator/(const Vec3 &a, double b) { return Vec3 { a.x / b, a.y / b, a.z / b }; }

inline double length2(const Vec3 &a)
{
    return a.x*a.x + a.y*a.y + a.z*a.z;
}

inline double length(const Vec3 &a)
{
    return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

void M_InverseMatrix(double *mat, int n, double *inv);

#endif
