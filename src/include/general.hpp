/* Some general types used in multiple files.
 *
 * Currently only includes a (rough) 24-bit integer type.
 *
 * Author: Maxim van den Berg
 */

#pragma once

#include <cmath>

/* Wrapper of integer such that the bit depth is
 * only 24 bits. Used often for audio file storage.
 */
struct int24 {
    int data : 24; // TODO: for space optimisation, use char[3].

    // Constructors
    int24() { data = 0; }
    template<typename type>
    int24(type a) {
        data = (int) a & ((1 << 24) - 1);
    }

    // Conversions
    operator int() const { return ((int) (data << 8)) / (1 << 8); }
    operator double() const { return (double) (int) (*this); }
    operator float() const { return (float) (int) (*this); }
    operator unsigned() const { return (unsigned int) (int) (*this); }
    operator long() const { return (long) (int) (*this); }

    // Operators
    friend int24 operator*(const int24 a, const int24 b)  {
        return { (int) a * (int) b };
    }
    friend int24 operator+(const int24 a, const int24 b)  {
        return { (int) a.data + (int) b.data };
    }
    friend int24 operator-(const int24 a, const int24 b)  {
        return { (int) a - (int) b };
    }

    friend int operator<(const int24 a, const int b) {
        return ((int) a) < b;
    }
    friend int operator-(const int24 a, const int b)  {
        return (int) a - b;
    }
    friend int operator+(const int24 a, const int b)  {
        return (int) a + b;
    }

    friend double operator*(const double a, const int24 b)  {
        return a * (double) b;
    }
    friend double operator+(const double a, const int24 b)  {
        return a + (double) b;
    }
    friend double operator-(const int24 a, const double b)  {
        return (double) a - b;
    }

    friend void operator+=(long a, const int24 b)  {
        a = a + (long) b;
    }

    friend float operator*(const float a, const int24 b)  {
        return a * (float) b;
    }
    friend float operator+(const float a, const int24 b)  {
        return a + (float) b;
    }
    friend float operator-(const int24 a, const float b)  {
        return (float) a - b;
    }

    friend void operator+=(int24 a, const double b)  {
        a = { (double) a + b };
    }
    friend void operator-=(int24 a, const double b)  {
        a = { (double) a - b };
    }

    friend void operator+=(unsigned long a, const int24 b)  {
        a += (int) a;
    }


    friend int24 operator-(int24 a)  {
        return { -((int) a) };
    }
    friend int abs(int24 a) {
        return abs((int) a.data);
    }
};
