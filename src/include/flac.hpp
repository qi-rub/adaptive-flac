/* FLAC encoder. Support for 4 LPC-types:
 *  - Fixed-blocksize Levinson-Durbin algorithm
 *  - Fixed-blocksize covariance Burg algorithm
 *  - Adaptive-blocksize covariance Burg algorithm
 *  - Adaptive-blocksize autocorrelation Burg algorithm
 *
 * (Decoding is not yet implemented)
 *
 * Reference: https://xiph.org/flac/format.html
 * Author: Maxim van den Berg
 */

#include <string>

#include "bitarray.hpp"
#include "golomb.hpp"
#include "lpc.hpp"
#include "wav.hpp"
#include "ringbuffer.hpp"



/* FLAC decoder and encoder.
 * Decoding is not yet fully implemented.
 */
class FLAC {
public:
    FLAC(LPCtype type, int order = 8);

    /* Source should be a .wav file, destination a .flac file. */
    void encode(std::string source, std::string destination,
                bool showprogress = false);

    /* Source should be a .flac file, destination a .wav file. */
    void decode(std::string source, std::string destination);

    LPCtype lpctype;

private:
    void header();

    template<typename type>
    void writedata(WavFile& wav, bool showprogress = false);

    void frameheader(int blocksize, bool midside,
                     unsigned long samplestart, unsigned bitpos);
    void subframeheader(int lpc_order);

    /* File general properties */
    int samplerate;
    int channels;
    int bitdepth;
    unsigned size;

    /* LPC data. */
    int order_default;

    std::shared_ptr<BitStreamOut> out;

    /* For decoding, currently unused. */
    std::shared_ptr<BitStreamIn> in;
};
