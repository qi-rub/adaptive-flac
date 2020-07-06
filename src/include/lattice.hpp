/* Lattice based methods for calculating the linear prediction
 * variables. Multiple methods are given.
 *
 * The methods use the reflection coefficients k_j, which will
 * need to be converted to the standard linear prediction
 * coefficients a_j at the end.
 *
 * Two alternative algorithms are given, which choose the final
 * blocksize in accordance with the analysed signal.
 *
 * Author: Maxim van den Berg
 */


#pragma once

#include <vector>
#include <iostream>

#include "general.hpp"
#include "plot.hpp"
#include "wav.hpp"

/* Calculates the standard LPC coefficients from the given
 * vector of reflection coefficients.
 */
std::vector<double> reflection2linear(std::vector<double>& k);


/* Burg algorithm for fixed blocksizes
 *
 * Uses the more efficient "covariance method", which means
 * the forward and backward predictions are never explicitly
 * computed, but estimated accurately from the covariance matrix.
 */
template<typename type>
std::vector<double> burgfixed(RingBuffer<type>& input, int order, int N);


/* Adaptive version of the covariance Burg algorithm above, which tries
 * to find a blocksizes which matches the given audio signal.
 *
 * Afterwards, updates the location in the given RingBuffer appropriately.
 *
 * See the report for more information.
 */
template<typename type>
std::vector<double> burgadaptive(RingBuffer<type>& input, int order,
                                 int Nmin, int Nmax);

/* Modification of the above adaptive version, which uses the autocorrelation
 * matrix instead of the covariance matrix. This makes the algorithm a lot
 * faster, at the cost of some blocksize prediction accuracy.
 */
template<typename type>
std::vector<double> burgadaptive_autocor(RingBuffer<type>& input, int order,
                                         int Nmin, int Nmax);


