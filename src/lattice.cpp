/* Lattice based methods for calculating the linear prediction
 * variables. Multiple methods are given.
 *
 * Author: Maxim van den Berg
 */

#include <memory>

#include "include/lattice.hpp"
#include "include/ringbuffer.hpp"
#include "include/levinson.hpp"

using namespace std;

vector<double> reflection2linear(vector<double>& k) {
    int order = k.size();
    vector<vector<double>> a(order, vector<double>());

    for (int m = 0; m < order; m++) {
        a[m].resize(m+1);
        a[m][m] = k[m];

        for (int j = 0; j < m; j++)
            a[m][j]= a[m-1][j] + k[m] * a[m-1][m-1-j];
    }

    return a[order-1];
}


/* Computes the autocorrelation matrix of N samples.
 * Since this matrix is Toeplitz, it is stored as a vector of
 * length `order+1`.
 */
template<typename type>
shared_ptr<vector<long>> autocor(RingBuffer<type>& input,
                                 int order, int N) {
    auto R = make_shared<vector<long>>(order+1);

    for (int l = 0; l <= order; l++)
        for (int n = l; n < N; n++)
            (*R)[l] += (long) input[n] * (long) input[n - l];

    return R;
}

/* Updates the given autocorrelation matrix by adding the specified sample.
 */
template<typename type>
void autocor_update(shared_ptr<vector<long>> R, RingBuffer<type>& input,
                    int order, unsigned n) {
    for (int l = 0; l <= order; l++)
        (*R)[l] += (long) input[n] * (long) input[n - l];
}


/* Updates the given autocorrelation matrix by adding the specified sample.
 */
template<typename type>
void autocor_update(long Rread[], long Rwrite[], RingBuffer<type>& input,
                    int order, unsigned n) {
    for (int l = 0; l <= order; l++)
        Rwrite[l] = Rread[l] + (long) input[n] * (long) input[n - l];
}

template<typename type>
shared_ptr<vector<long>> autocor2covar(RingBuffer<type>& input,
                                       shared_ptr<vector<long>> R,
                                       int N) {
    return autocor2covar(input, &(*R)[0], R->size(), N);
}


/* Creates the covariance matrix from the given autocorrelation matrix.
 */
template<typename type>
shared_ptr<vector<long>> autocor2covar(RingBuffer<type>& input,
                                       long* R, int size, int N) {
    auto phi = make_shared<vector<long>>(size * size);

    for (int l = 0; l < size; l++) {
        for (int r = 0; r <= l; r++) {
            if (l == 0 || r == 0)
                (*phi)[l*size + r] = R[max(r, l)];
            else {
                (*phi)[l*size + r] = (*phi)[(l-1)*size + r-1]
                    - (long) input[N - l] * (long) input[N - r];
            }
        }
    }

    return phi;
}

/* Computes and returns the covariance matrix.
 * This matrix is stored as a single vector, using that the size is
 * fixed to `order + 1`. Sets only the lower triangle of the matrix.
 */
template<typename type>
shared_ptr<vector<long>> covar(RingBuffer<type>& input, int order, int N) {
    auto R = autocor(input, order, N);
    auto phi = make_shared<vector<long>>((order+1) * (order+1));

    for (int l = 0; l <= order; l++) {
        for (int r = 0; r <= l; r++) {
            for (int n = max(l,r); n < N; n++)
                (*phi)[l*(order+1) + r] += (int) input[n - l] *
                    (int) input[n - r];
        }
    }

    return phi;
}

/* Updates the given covariance matrix by adding the specified sample.
 */
template<typename type>
void covar_update(shared_ptr<vector<long>> phi, RingBuffer<type>& input,
                  int order, unsigned n) {
    for (int l = 0; l <= order; l++) {
        for (int r = 0; r <= l; r++)
            (*phi)[l*(order+1) + r] += (long) input[n - l] * (long) input[n - r];
    }
}


/* Updates the given covariance matrix by adding the specified sample.
 */
template<typename type>
void covar_update(long phiread[], long phiwrite[], RingBuffer<type>& input,
                  int order, unsigned n) {
    for (int l = 0; l <= order; l++) {
        for (int r = 0; r <= l; r++)
            phiwrite[l*(order+1) + r] = phiread[l*(order+1) + r] +
                (long) input[n - l] * (long) input[n - r];
    }
}


/* Calculate either F[m], B[m] or C[m], based on wether the
 * iteration process should be inverted (specified by `backleft`
 * and `backright`).
 */
double sqrtav(shared_ptr<vector<long>> phi, vector<double> a,
              int m, int order, bool backward_l, bool backward_r) {
    double average = 0;

    for (int l = 0; l <= m; l++) {
        for (int r = 0; r <= m; r++) {
            int row = backward_l ? m + 1 - l : l;
            int column = backward_r ? m + 1 - r : r;

            average += a[l]*a[r] * (*phi)[max(row,column)*(order+1) + min(row,column)];
        }
    }

    return average;
}

/* Calculate either F[m], B[m] or C[m], based on wether the
 * iteration process should be inverted (specified by `backleft`
 * and `backright`).
 */
double sqrtav(long phi[], vector<double> a, int m, int order,
              bool backward_l, bool backward_r) {
    double average = 0;

    for (int l = 0; l <= m; l++) {
        for (int r = 0; r <= m; r++) {
            int row = backward_l ? m + 1 - l : l;
            int column = backward_r ? m + 1 - r : r;

            average += a[l]*a[r] * phi[max(row,column)*(order+1) + min(row,column)];
        }
    }

    return average;
}

/* Same as the normal `sqrtav` but using the autocorrelation matrix instead.
 */
double sqrtav_autocor(shared_ptr<vector<long>> R, vector<double> a,
                      int m, bool backward_l, bool backward_r) {
    double average = 0;

    for (int l = 0; l <= m; l++) {
        for (int r = 0; r <= m; r++) {
            int row = backward_l ? m + 1 - l : l;
            int column = backward_r ? m + 1 - r : r;

            average += a[l]*a[r] * (*R)[abs(row - column)];
        }
    }

    return average;
}

/* Same as the normal `sqrtav` but using the autocorrelation matrix instead.
 */
double sqrtav_autocor(long R[], vector<double> a, int m,
                      bool backward_l, bool backward_r) {
    double average = 0;

    for (int l = 0; l <= m; l++) {
        for (int r = 0; r <= m; r++) {
            int row = backward_l ? m + 1 - l : l;
            int column = backward_r ? m + 1 - r : r;

            average += a[l]*a[r] * R[abs(row - column)];
        }
    }

    return average;
}


/* Computes the Golomb-Rice expected codeword length from the
 * given squared erro Faverage = F_p / N = 1 / N * Σ u[n]^2.
 */
double expected_codeword_length(double Faverage) {
    // Estimate the distribution parameter.
    double p = sqrt(1 / (Faverage + 1));

    // Calculate the optimal guess for the Golomb-Rice parameter.
    double k = p >= 0.5 ? 0 : -log2(log2(1 / (1 - p)));

    // Calculate the exptected Golomb-Rice codeword length.
    return 1 + k + 1 / (1 - pow((1-p),pow(2,k)));
}


// static int blockcounter = 0;
vector<double> covarianceburg(shared_ptr<vector<long>> phi, int order) {
    vector<double> k(order+1);
    vector<vector<double>> a(order+1, vector<double>());
    a[0].insert(a[0].begin(), 1);

    // The forward, backward and center prediction squared
    // means. Used to compute the next reflection coefficient.
    double F, B, C;

    for (int m = 0; m < order; m++) {

        // Calculate the squared means.
        F = sqrtav(phi, a[m], m, order, false, false);
        B = sqrtav(phi, a[m], m, order, true, true);
        C = sqrtav(phi, a[m], m, order, false, true);

        // Calculate the next burg reflection coefficient.
        if (F + B != 0) {
            k[m+1] = -2 * C / (F + B);
        } else
            k[m+1] = 0;

        if (abs(k[m+1]) > 1) {
            cout << "Warning: Found instability, with k[" << m
                << "] = " << k[m+1] << endl;
            k[m+1] = k[m+1] > 0 ? 1 : -1;
        }

        // Update the LPC coefficients.
        a[m+1].resize(m+2);
        a[m+1][0] = 1;
        a[m+1][m+1] = k[m+1];

        for (int j = 1; j < m+1; j++)
            a[m+1][j] = a[m][j] + k[m+1] * a[m][m+1-j];
    }

    // Convert to the correct representation.
    a[order].erase(a[order].begin());
    for (unsigned i = 0; i < a[order].size(); i++)
        a[order][i] *= -1;

    return a[order];
}


template<typename type>
vector<double> burgfixed(RingBuffer<type>& input, int order, int N) {
    // The covariance matrix, which will be used repeatedly
    // to compute F_m, B_m and C_m in Burg's algorithm;
    auto phi = autocor2covar(input, autocor(input, order, N), N);
    input.reposition(input.position() + N);

    // Execute the covariance burg algorithm;
    return covarianceburg(phi, order);
}

// -- Required for plotting the residuals graphs -- //
// static int blockcounter = 0;
// ------------------------------------------------ //

template<typename type>
vector<double> burgadaptive_autocor(RingBuffer<type>& input, int order,
                                    int Nmin, int Nmax) {
    vector<double> k(order+1);
    vector<vector<double>> a(order+1, vector<double>());
    a[0].insert(a[0].begin(), 1);

    // The forward, backward and center prediction squared means.
    double F, B, C;

    // Initialisation of the covariance matrix.
    auto R0 = autocor(input, order, Nmin);
    auto R1 = make_shared<vector<long>>(*R0);
    int Rcur = 0;
    long* Rread = &(*R0)[0];
    long* Rwrite = &(*R0)[0];

    // Init the variables for detection of the best blocksize.
    double Fbest = sqrtav_autocor(Rwrite, a[0], 0, false, false) /
        (double) Nmin;
    int ordersmall = min(order, 4); // The order used in blocksize detection
    int Nbest = Nmin;   // The value of N for which F_p is lowest
    int counter = 0;    // n current - Nbest


    // -- Uncomment for plotting of the residual error ------- //
    // -- (uncomment together with the gnuplot calls below) -- //
    // Plot gnuplot("plots/data" + to_string(blockcounter++), true);
    // ------------------------------------------------------- //

    for (int n = Nmin; n < Nmax; n++) {
        int m = (n - Nmin) % ordersmall;

        autocor_update(Rread, Rwrite, input, order, n);
        Rread = Rwrite;

        // Calculate the squared means.
        F = sqrtav_autocor(Rwrite, a[m], m, false, false);
        B = sqrtav_autocor(Rwrite, a[m], m, true, true);
        C = sqrtav_autocor(Rwrite, a[m], m, false, true);

        // Calculate the next burg reflection coefficient.
        if (F + B != 0) {
            k[m+1] = -2 * C / (F + B);
        } else
            k[m+1] = 0;

        // Update the LPC coefficients.
        a[m+1].resize(m+2);
        a[m+1][0] = 1;
        a[m+1][m+1] = k[m+1];

        for (int j = 1; j < m+1; j++)
            a[m+1][j] = a[m][j] + k[m+1] * a[m][m+1-j];

        // Terminate iteration if the optimal N is found.
        int COUNTER_MAX = 2048;
        float RANGE_OPTIMAL = 1.01;
        if (m + 1 == ordersmall) {
            double Ferror = sqrtav_autocor(Rwrite, a[m+1], m+1, false, false) / (n+1);

            // -------------------- //
            // gnuplot.plot(n, Ferror);
            // -------------------- //

            if (Ferror <= Fbest) {
                // Set the new value for N which minimalizes F_p.
                Nbest = n + 1;
                Fbest = Ferror;

                // Swap the R buffer.
                Rcur = !Rcur;
                Rwrite = Rcur ? &(*R1)[0] : &(*R0)[0];

                counter = 0;
            } else if (Ferror <= Fbest * RANGE_OPTIMAL) {
                // Set the optimal value for N, whose error F_p is close to
                // the absolute minimum Fbest, but maximizes block length.
                Nbest = n + 1;

                // Swap the R buffer.
                Rcur = !Rcur;
                Rwrite = Rcur ? &(*R1)[0] : &(*R0)[0];

                counter = 0;
            } else {
                counter += ordersmall;
                if (counter > COUNTER_MAX) {
                    // -------------------- //
                    // gnuplot.plot(Nbest, 0);
                    // -------------------- //

                    break;
                }
            }
        }
    }

    // Calculate the remaining order using the created covariance matrix.
    auto phi = autocor2covar(input, Rcur ? R0 : R1, Nbest);

    auto a_n = covarianceburg(phi, order);
    input.reposition(input.position() + Nbest);

    return a_n;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

template<typename type>
vector<double> burgadaptive(RingBuffer<type>& input, int order,
                            int Nmin, int Nmax) {
    vector<double> k(order+1);
    vector<vector<double>> a(order+1, vector<double>());
    a[0].insert(a[0].begin(), 1);

    // The forward, backward and center prediction squared means.
    double F, B, C;

    // Initialisation of the covariance matrices.
    auto R = autocor(input, order, Nmin);
    auto phi0 = autocor2covar(input, R, Nmin);
    auto phi1 = make_shared<vector<long>>(*phi0);
    int phicur = 0;
    long* phiread = &(*phi0)[0];
    long* phiwrite = &(*phi0)[0];

    // Init the variables for detection of the best blocksize.
    double Fbest = sqrtav(phi1, a[0], 0, order, false, false) /
        (double) Nmin;

    // -- Alternative error measure. Assumes residuals ~ 2Geom -- //
    // double Ebest = expected_codeword_length(Fstart);

    int ordersmall = min(order, 4);  // The order used in blocksize detection
    int Nbest = Nmin;    // The highest value of N for which F_p is lowest
    int counter = 0;     // n current - Nbest

    // -- Uncomment for plotting of the residual error ------- //
    // -- (uncomment together with the gnuplot calls below) -- //
    // Plot gnuplot("plots/data" + to_string(blockcounter++), true);
    // ------------------------------------------------------- //

    for (int n = Nmin; n < Nmax; n++) {
        int m = (n - Nmin) % ordersmall;

        covar_update(phiread, phiwrite, input, order, n);
        phiread = phiwrite;

        // Calculate the squared means.
        F = sqrtav(phiwrite, a[m], m, order, false, false);
        B = sqrtav(phiwrite, a[m], m, order, true, true);
        C = sqrtav(phiwrite, a[m], m, order, false, true);

        // Calculate the next burg reflection coefficient.
        if (F + B != 0) {
            k[m+1] = -2 * C / (F + B);
        } else
            k[m+1] = 0;

        // Update the LPC coefficients.
        a[m+1].resize(m+2);
        a[m+1][0] = 1;
        a[m+1][m+1] = k[m+1];

        for (int j = 1; j < m+1; j++)
            a[m+1][j] = a[m][j] + k[m+1] * a[m][m+1-j];

        // Terminate iteration if the optimal N is found.
        int COUNTER_MAX = 2048;
        float RANGE_OPTIMAL = 1.01;
        if (m + 1 == ordersmall) {
            double Ferror = sqrtav(phiwrite, a[m+1], m+1, order, false, false)
                / (n+1);
            // -- Alternative error measure. Assumes residuals ~ 2Geom -- //
            // double E = expected_codeword_length(Ferror);

            // -------------------- //
            // gnuplot.plot(n, Ferror);
            // -------------------- //

            if (Ferror <= Fbest) {
                // Set the new value for N which minimalizes F_p.
                Nbest = n + 1;
                Fbest = Ferror;

                // Swap the phi buffer.
                phicur = !phicur;
                phiwrite = phicur ? &(*phi1)[0] : &(*phi0)[0];

                counter = 0;
            } else if (Ferror <= Fbest * RANGE_OPTIMAL) {
                // Set the optimal value for N, whose error F_p is close to
                // the absolute minimum Fbest, but maximizes block length.
                Nbest = n + 1;

                // Swap the phi buffer.
                phicur = !phicur;
                phiwrite = phicur ? &(*phi1)[0] : &(*phi0)[0];

                counter = 0;
            } else {
                counter += ordersmall;
                if (counter > COUNTER_MAX) {
                    // -------------------- //
                    // gnuplot.plot(Nbest, 0);
                    // -------------------- //

                    break;
                }
            }
        }
    }

    auto a_n = covarianceburg(phicur ? phi0 : phi1, order);
    input.reposition(input.position() + Nbest);

    return a_n;
}


template vector<double> burgfixed(RingBuffer<char>& input, int order, int N);
template vector<double> burgfixed(RingBuffer<short>& input, int order, int N);
template vector<double> burgfixed(RingBuffer<int24>& input, int order, int N);
template vector<double> burgfixed(RingBuffer<int>& input, int order, int N);

template vector<double> burgadaptive(RingBuffer<char>& input, int order,
                                     int Nmin, int Nmax);
template vector<double> burgadaptive(RingBuffer<short>& input, int order,
                                     int Nmin, int Nmax);
template vector<double> burgadaptive(RingBuffer<int24>& input, int order,
                                     int Nmin, int Nmax);
template vector<double> burgadaptive(RingBuffer<int>& input, int order,
                                     int Nmin, int Nmax);

template vector<double> burgadaptive_autocor(RingBuffer<char>& input, int order,
                                     int Nmin, int Nmax);
template vector<double> burgadaptive_autocor(RingBuffer<short>& input, int order,
                                     int Nmin, int Nmax);
template vector<double> burgadaptive_autocor(RingBuffer<int24>& input, int order,
                                     int Nmin, int Nmax);
template vector<double> burgadaptive_autocor(RingBuffer<int>& input, int order,
                                     int Nmin, int Nmax);
