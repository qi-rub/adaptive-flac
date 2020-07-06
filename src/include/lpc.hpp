/* Linear predictive coding using one of the following algorithms.
 *  - Fixed-blocksize Levinson-Durbin algorithm
 *  - Fixed-blocksize covariance Burg algorithm
 *  - Adaptive-blocksize covariance Burg algorithm
 *  - Adaptive-blocksize autocorrelation Burg algorithm
 *
 * See levinson.hpp and lattice.hpp for further information
 * about these algorithms.
 *
 * Also support quantization of the LPC coefficients if a
 * `precision` value is specified to the constructor.
 *
 * Author: Maxim van den Berg
 */

#pragma once

#include <vector>
#include <memory>

#include "levinson.hpp"
#include "lattice.hpp"
#include "ringbuffer.hpp"
#include "general.hpp"

enum LPCtype {
    LevinsonDurbinFixed,
    // LevinsonDurbinAdaptive,
    LatticeBurgFixed,
    LatticeBurgAdaptive,
    LatticeBurgAdaptiveAutocorrelation
};


/* A standard implementation for Linear Predictive Coding, where
 * - The parameters are computed using a linear least squares approach,
 *   using the Levinson-Durbin algorithm (see `levinson.hpp`).
 * - A windowing function is applied to the signal.
 *
 * The analysis computes all residuals, such that the original signal
 * can be reconstructed completely.
 */
template<typename T>
class LPC {

public:
    /* A struct for gathering of all data which the LPC analysis algorithm
     * outputs, which are required to reconstruct the signal losslessly.
     *
     * They include:
     * - parameters: The `order` parameters which specify the LPC model.
     * - warmup:     The `order` unencoded starting samples.
     * - residuals:  All calculated residuals which uniquely define the
     *               original signal.
     */
    struct Analysed {
        /* Initialise by making all vectors the correct size.  */
        Analysed(int ord, int N) {
            warmup.resize(ord, 0);
            residuals.resize(N - ord, 0);
            abs_sum = 0;

            order = ord;
            size = N;

            qlevel = 0;
            qprecision = -1;
        }

        // Size of analysed block and filter order.
        int size;
        int order;

        std::vector<T> warmup;
        std::vector<T> residuals;
        double abs_sum;

        std::vector<double> parameters;

        /* Quantizer variables. */
        std::vector<int> qparameters;
        int qlevel;
        int qprecision;
    };

    /* Performs a full analysis of the given input datastream. Puts the
     * analysed data into the `analysis` parameter. For more information
     * on the returned object and the computational method, see above.
     *
     * When precision is set, quantizes the parameters to `precision`
     * bits values.
     */
    std::shared_ptr<Analysed> analyse(RingBuffer<T>& input, int qlp_precision,
                                      LPCtype lpctype, int order_max,
                                      int Ndefault, int Nmin = -1,
                                      int Nmax = -1);

    /* Reassembles the original signal using the LPC model and the given
     * residuals.
     */
    std::shared_ptr<std::vector<T>> synthesize(Analysed& analysis);

private:
};
