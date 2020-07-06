/* The Levinson-Durbin algorithm for solving a system of linear
 * equations based on the autocorrelation of a given input signal,
 * thus where the matrix in question is Toeplitz.
 *
 * Used in lpc.cpp for minimizing the residual error.
 *
 * Author: Maxim van den Berg
 */

#define _USE_MATH_DEFINES
#include <cmath>

#include "include/levinson.hpp"
#include "include/ringbuffer.hpp"
#include "include/general.hpp"

using namespace std;

template<typename type>
vector<double> levinsondurbin(RingBuffer<type>& input, int order, int N) {
    vector<long> R(order+1);

    // Apply the Hamming windowing.
    double x[N];
    for (int n = 0; n < N; n++)
        x[n] = (double) input[n] * (1 - cos(2 * M_PI * n / N)) / 2;

    // Setup the autocorrelation vector.
    for (int i = 0; i <= order; i++)
        for (int n = i; n < N; n++)
            R[i] += ((long) x[n-i] * (long) x[n]);

    vector<double> a(order+1);

    double backward[order];

    backward[0] = R[0] ? -1 / (double) R[0] : 1;
    a[0] = R[0] ? -R[1] / (double) R[0] : 1;

    for (int m = 1; m <= order; m++) {
        // Calculate of the forward and solution
        // vectors error at index n.
        double err_b = 0, err = 0;
        for (int i = 1; i <= m; i++) {
            err_b += backward[i-1] * R[i];
            err += a[i-1] * R[m+1-i];
        }

        // Update the forward vector.
        double backwardnew[order];
        backwardnew[0] = err_b * backward[m-1] / (err_b * err_b - 1);
        for (int i = 1; i <= m; i++)
            backwardnew[i] = (-backward[m-1-i] + err_b * backward[i-1]) /
                             (err_b * err_b - 1);

        // Update the backward and solution vectors.
        for (int i = 0; i <= m; i++) {
            backward[i] = backwardnew[i];
            a[i] -= (R[m] + err) * backward[i];
        }
    }

    input.reposition(input.position() + N);
    return a;
}


template<typename type>
vector<double> levinsondurbinmodified(RingBuffer<type>& input, int order, int N) {
    vector<double> R(order+1);

    // Apply the Hamming windowing.
    double x[N];
    for (int n = 0; n < N; n++)
        // x[n] = (double) input[n] * (1 - cos(2 * M_PI * n / N)) / 2;
        x[n] = (double) input[n];

    // Setup the autocorrelation vector.
    for (int i = 0; i <= order; i++)
        for (int n = i; n < N; n++)
            R[i] += x[n-i] * x[n];

    vector<vector<double>> a(order+1, vector<double>(order+1));

    // Zero block detection.
    if (R[0] == 0) {
        input.reposition(input.position() + N);

        return a[0];
    }

    // Levinson-Durbin setup.
    a[0][0] = 1;
    double V = R[0];

    // The Levinson-Durbin recursion.
    for (int m = 0; m < order; m++) {
        // Calculate the error term.
        double err = 0;
        for (int i = 0; i <= m; i++)
            err += a[m][i] * R[m+1-i];
        err /= V;

        // Calculate the next LPC coefficient set.
        a[m+1][0] = 1;
        a[m+1][m+1] = -err;
        for (int i = 1; i <= m; i++)
            a[m+1][i] = a[m][i] - err * a[m][m+1-i];

        // Update the correction value.
        V *= (1 - err*err);
    }

    input.reposition(input.position() + N);
    a[order].erase(a[order].begin());

    for (int i = 0; i < order; i++)
        a[order][i] *= -1;

    return a[order];
}

template vector<double> levinsondurbin(RingBuffer<char>& input,
                                       int N, int order);
template vector<double> levinsondurbin(RingBuffer<short>& input,
                                       int N, int order);
template vector<double> levinsondurbin(RingBuffer<int24>& input,
                                       int N, int order);
template vector<double> levinsondurbin(RingBuffer<int>& input,
                                       int N, int order);

template vector<double> levinsondurbinmodified(RingBuffer<char>& input,
                                       int N, int order);
template vector<double> levinsondurbinmodified(RingBuffer<short>& input,
                                       int N, int order);
template vector<double> levinsondurbinmodified(RingBuffer<int24>& input,
                                       int N, int order);
template vector<double> levinsondurbinmodified(RingBuffer<int>& input,
                                       int N, int order);
