/* Two linear prediction algorithms based on the Levinson-Durbin
 * method for solving linear systems of Toeplitz matrices.
 *
 * Author: Maxim van den Berg
 */

#pragma once

#include <vector>
#include <iostream>

#include "ringbuffer.hpp"

/* Standard Levinson-Durbin algorithm for linear predictive coding.
 * Computes the autocorrelation matrix and solves the linear Toeplitz system.
 * Returns the LPC coefficient vector.
 *
 * Also applies the Hamming window before analysing.
 */
template<typename T>
std::vector<double> levinsondurbin(RingBuffer<T>& input, int order, int N);

/* Modified version of the Levinson Durbin algorithm.
 */
template<typename T>
std::vector<double> levinsondurbinmodified(RingBuffer<T>& input, int order, int N);
