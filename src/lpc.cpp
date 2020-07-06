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

#include "include/lpc.hpp"

#define Analysed typename LPC<type>::Analysed

using namespace std;

template<typename type>
void quantizeparameters(shared_ptr<Analysed> analysis,
                        unsigned precision) {
    analysis->qparameters.resize(analysis->order);

    // Determine the quantized coefficient bounds.
	int qmax = (1 << (precision - 1)) - 1;
	int qmin = -qmax - 1;

    // Determine the largest parameter (in norm).
    double cmax = 0.0;
	for (int i = 0; i < analysis->order; i++)
		cmax = max(cmax, fabs(analysis->parameters[i]));

    // Get the required shift.
    int log2cmax;
    frexp(cmax, &log2cmax);
    int shift = (int) precision - log2cmax - 1;

    // Quantize the parameters.
    double error = 0.0;
    int quantized;

    if (shift < 0) {
        cout << "Warning: negative shift!" << endl;
        cout << shift << endl;
        for (auto x : analysis->parameters)
            cout << x << " ";
        cout << endl;
    }

    for (int i = 0; i < analysis->order; i++) {
        error += analysis->parameters[i] * (1 << shift);
        quantized = lround(error);

        // Make sure we don't go out of bounds.
        if (quantized > qmax)
            quantized = qmax;
        else if (quantized < qmin)
            quantized = qmin;

        error -= quantized;
        analysis->qparameters[i] = quantized;
    }

    analysis->qprecision = precision;
    analysis->qlevel = shift;
}

template<typename type>
shared_ptr<Analysed> LPC<type>::analyse(RingBuffer<type>& input,
                                        int qlp_precision, LPCtype lpctype,
                                        int order_max, int Ndefault,
                                        int Nmin, int Nmax) {
    vector<double> a;
    int N = Nmin != -1 ? min(Ndefault, Nmax) : Ndefault;

    unsigned startfilled = input.filled();
    unsigned startpos = input.position();

    // Find the optimal a_k coefficients.
    switch (lpctype) {
        case LevinsonDurbinFixed:
            a = levinsondurbinmodified(input, order_max, N);
            break;
        case LatticeBurgFixed:
            a = burgfixed(input, order_max, N);
            break;
        case LatticeBurgAdaptive:
            if (Nmin != -1 && Nmax != -1)
                // Adaptive Burg determines a optimal block size, so
                // the size of `input` will be changed.
                a = burgadaptive(input, order_max, Nmin, Nmax);
            else
                a = burgfixed(input, order_max, N);
            break;
        case LatticeBurgAdaptiveAutocorrelation:
            if (Nmin != -1 && Nmax != -1)
                // Adaptive Burg determines a optimal block size, so
                // size of `input` will be changed.
                a = burgadaptive_autocor(input, order_max, Nmin, Nmax);
            else
                a = burgfixed(input, order_max, N);
            break;
    }


    // Update N if it was changed during an adapative method.
    N = startfilled - input.filled();
    input.reposition(startpos);

    shared_ptr<Analysed> analysis = make_shared<Analysed>(a.size(), N);
    analysis->parameters = a;
    int order = a.size();

    // Quantize parameters
    if (qlp_precision > 0)
        quantizeparameters<type>(analysis, qlp_precision);

    // Write the warmup values.
    for (int n = 0; n < order; n++)
        analysis->warmup[n] = input[n];

    // Calculate the residuals
    for (int n = order; n < N; n++) {
        long sum = 0;

        for (int k = 0; k < order; k++)
            sum += analysis->qparameters[k] * (int) input[n-k-1];

        analysis->residuals[n - order] = input[n] -
            (int) (sum >> analysis->qlevel);

        analysis->abs_sum += abs(analysis->residuals[n - order]);
    }

    input.reposition(startpos + N);

    return analysis;
}

template<typename type>
std::shared_ptr<vector<type>> LPC<type>::synthesize(Analysed& analysis) {
    int order = analysis.parameters.size();
    unsigned N = analysis.residuals.size();

    auto output = make_shared<vector<type>>(N, 0);

    for (int i = 0; i < order; i++)
        (*output)[i] = analysis.warmup[i];

    for (unsigned n = order; n < N; n++) {
        double res = analysis.residuals[n - order];

        for (int k = 0; k < order; k++)
            res -= (analysis.parameters[k] / (1 << analysis.qlevel))
                   * (double) (*output)[n-k-1];

        (*output)[n] = round(res);
    }

    return output;
}



template class LPC<int>;
template class LPC<int24>;
template class LPC<short>;
template class LPC<char>;
