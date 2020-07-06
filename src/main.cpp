/* A FLAC encoder using four possible linear predictive coding algorithms.
 *
 * Usage: encoder filename [order]
 *
 * The given file should be a .wav file.
 *
 * The order is optional, and specifies the orders all LPC algorithms should use.
 * The default is 8.
 *
 * The algorithms are:
 *  - Fixed-blocksize Levinson-Durbin algorithm
 *  - Fixed-blocksize covariance Burg algorithm
 *  - Adaptive-blocksize covariance Burg algorithm
 *  - Adaptive-blocksize autocorrelation Burg algorithm
 *
 * Compresses the given file using all four algorithms, and
 * prints the new file sizes, compression ratios and encoding times.
 *
 * Author: Maxim van den Berg
 */

#include <iostream>
#include <bitset>
#include <string>
#include <ctime>
#include <filesystem>

#include "include/bitarray.hpp"
#include "include/golomb.hpp"
#include "include/levinson.hpp"
#include "include/lpc.hpp"
#include "include/wav.hpp"
#include "include/flac.hpp"

using namespace std;


// -- Uncomment to save the statistics to data files -- //
// -- (Also see the function calls below -------------- //
// static Plot plotcompression("out/compression", false);
// static Plot plotexectime("out/exectime", false);
// ---------------------------------------------------- //

static int order = 8;

void run_flac(FLAC* flac, string method, string filename, string out) {
    cout << method << ":" << endl;

    time_t timer;
    double t;
    long filesize = filesystem::file_size(filename);

    timer = clock();
    flac->encode(filename, out, false);
    t = (clock() - timer) / (double) CLOCKS_PER_SEC;

    long outsize = filesystem::file_size(out);
    float outsize_ratio = outsize / (double) filesize * 100;

    cout << " done!" << endl;
    cout << " it took " << t << " seconds" << endl;
    cout << " and the filesize is " << outsize << " (" <<
        outsize_ratio << "%)" << endl;

    // -- Uncomment to save the statistics to data files -- //
    // -- (Also see the class inits above) ---------------- //
    // Plot plotcompression("out/compression-lpc" + to_string(flac->lpctype)
    //         + "-order" + to_string(order), false);
    // plotcompression.write(outsize_ratio);
    // Plot plotexectime("out/exectime-lpc" + to_string(flac->lpctype)
    //         + "-order" + to_string(order), false);
    // plotexectime.write(t);
    // ---------------------------------------------------- //
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Usage: encoder [.wav file] [order]" << endl;
        cout << "The order is optional (default is 8)" << endl;

        return 1;
    }

    if (argc == 3)
        order = stoi(argv[2]);

    FLAC* flacburgadaptive = new FLAC(LatticeBurgAdaptive, order);
    FLAC* flacburgadaptiveautocor = new FLAC(LatticeBurgAdaptiveAutocorrelation, order);
    FLAC* flacburgfixed = new FLAC(LatticeBurgFixed, order);
    FLAC* flaclevinsonfixed = new FLAC(LevinsonDurbinFixed, order);

    cout << "encoding " << string(argv[1]) << endl;
    cout << "with filesize " << filesystem::file_size(string(argv[1])) << endl << endl;

    run_flac(flaclevinsonfixed, "levinson-durbin fixed", string(argv[1]),
            "out/out.flac-levinson-fixed" + to_string(order));
    run_flac(flacburgfixed, "burg fixed", string(argv[1]),
            "out/out.flac-burg-fixed" + to_string(order));
    run_flac(flacburgadaptive, "burg adaptive", string(argv[1]),
            "out/out.flac-burg-adaptive" + to_string(order));
    run_flac(flacburgadaptiveautocor, "burg adaptive (using autocorrelation)",
             string(argv[1]), "out/out.flac-burg-adaptive-autocor" + to_string(order));

    delete(flacburgadaptive);
    delete(flacburgadaptiveautocor);
    delete(flacburgfixed);
    delete(flaclevinsonfixed);

    return 0;
}
