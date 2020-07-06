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

#include <bitset>
#include <limits>
#include <algorithm>
#include <random>

#include "include/flac.hpp"

using namespace std;

/* Returns the bit code used by FLAC in every frame header for
 * the given block size. When returning 0b0111, a 16 bit sample
 * rate should be specified at the end of the header.
 */
inline unsigned header_blocksize(unsigned blocksize) {
    switch (blocksize) {
        case 192:     return 0b0001;
        case 576:     return 0b0010;
        case 1152:    return 0b0011;
        case 2304:    return 0b0100;
        case 4608:    return 0b0101;
        case 256:     return 0b1000;
        case 512:     return 0b1001;
        case 1024:    return 0b1010;
        case 2048:    return 0b1011;
        case 4096:    return 0b1100;
        case 8192:    return 0b1101;
        case 16384:   return 0b1110;
        case 32768:   return 0b1111;
        default:      return 0b0111;
    }
}

/* Returns the bit code used by FLAC in every frame header for
 * the given samplerate. When returning 0b1101, a 16 bit sample
 * rate should be specified at the end of the header.
 */
inline unsigned header_samplerate(unsigned samplerate) {
    switch (samplerate) {
        case 88200:   return 0b0001;
        case 176400:  return 0b0010;
        case 192000:  return 0b0011;
        case 8000:    return 0b0100;
        case 16000:   return 0b0101;
        case 22050:   return 0b0110;
        case 24000:   return 0b0111;
        case 32000:   return 0b1000;
        case 44100:   return 0b1001;
        case 48000:   return 0b1010;
        case 96000:   return 0b1011;
        default:      return 0b1101;
    }
}

/* Returns the bit code used by FLAC in every frame header for
 * the given bitrate. When returning 0b000, the value should
 * be read from the FLAC metadata header.
 */
inline unsigned header_bitdepth(unsigned bitdepth) {
    switch (bitdepth) {
        case 8:     return 0b001;
        case 12:    return 0b010;
        case 16:    return 0b100;
        case 20:    return 0b101;
        case 24:    return 0b110;
        default:    return 0b000;
    }
}

/* Computes the CRC-8 or CRC-16 of the given bitarray, up until the
 * last written byte, with the given generator polynomial `divisor`.
 * Uses a shift register for the computations.
 *
 * Since the bitarray is not required to be 32-bit aligned (the size
 * of an bitarray element), the bit position where reading has to
 * begin is specified as well.
 */
template <int T>
unsigned computeCRC(BitArray* array, const unsigned divisor,
                    unsigned bitpos) {
    // Initialise the reader and shift register.
    array->setreader(bitpos);
    unsigned reg = 0;

    try {
        while (true) {
            // Load next byte into the register.
            reg ^= array->read(8) << (T - 8);

            // Perform the shift.
            for (int i = 0; i < 8; i++) {
                reg <<= 1;

                // Division with repeated XOR's.
                if ((1 << T) & reg)
                    reg ^= divisor;
            }
        }
    } catch (const std::range_error& e) {
    }

    return reg;
}

/* Writes the UTF8 variable width character encoding of `value` to
 * the `out` bitarray.
 */
void writeUTF8(BitArray* out, unsigned long value) {
    if (value < 128)
        out->write(value, 8);
    else {
        // Determine amount of required byte fields.
        int bits = floor(log2(value)) + 1;
        int fields = bits / 6;
        fields += bits % 6 <= 6 - fields ? 0 : 1;

        // Write the first field.
        out->write(0xff, fields + 1);
        out->write(0, 1);
        out->write((value >> 6 * fields),  6 - fields);

        // Write remaining fields.
        for (int i = 0; i < fields; i++) {
            unsigned masked = 0b111111 & (value >> 6 * (fields - i - 1));
            out->write((0b10 << 6) | masked, 8);
        }
    }
}

/* Fills the given vector with randomly sampled values from the
 * two-sided geometric distribution with parameter p. Overwrite
 * the data in `fill`.
 */
template<typename type>
void write2geom(vector<type>& fill, double p) {
    default_random_engine gen;
    geometric_distribution<int> geom(p);
    bernoulli_distribution bernoulli(0.5);

    for (unsigned i = 0; i < fill.size(); i++)
        fill[i] = bernoulli(gen) ? geom(gen) : -geom(gen);
}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* The constructor sets the LPC type. */
FLAC::FLAC(LPCtype type, int order) {
    lpctype = type;
    order_default = order;
}


/* Write the FLAC-header, consisting of only the required
 * STREAMINFO block.
 */
void FLAC::header() {
    // Write "fLaC".
    out->write('f', 8);
    out->write('L', 8);
    out->write('a', 8);
    out->write('C', 8);

    // Write the STREAMINFO metadata block.
    out->write(1, 1);      // Specify this is the last metadata block
    out->write(0, 7);      // Specify this is the STREAMINFO block
    out->write(272/8, 24); // Specify length of the block

    out->write(16, 16);    // Minimum block size (min value)
    out->write(65535, 16); // Maximum block size (max value)

    out->write(0, 24);     // Minimum frame size (0 implies unknown)
    out->write(0, 24);     // Maximum frame size (0 implies unknown)

    out->write(samplerate, 20);  // Sample rate (Hz)
    out->write(channels - 1, 3); // Channel count
    out->write(bitdepth - 1, 5); // Bits per sample

    // Set the total samples in stream (36 bits).
    unsigned totalsize = (size / channels);
    out->write((totalsize >> 4), 32);
    out->write(totalsize & 0xf, 4);

    out->write(0, 32); // MD5 signature of the unencoded audio data
    out->write(0, 32); //  ^                            (skipped)
    out->write(0, 32); //  ^
    out->write(0, 32); //  ^
}

/* Write the FLAC frameheader. This signals the start of the processing of a
 * new block. One or more subframeheaders always follows (see below).
 */
void FLAC::frameheader(int blocksize, bool midside, unsigned long samplestart,
                       unsigned bitpos) {

    out->write(0b11111111111110, 14); // Sync code
    out->write(0, 1); // Mandatory bit
    out->write(1, 1); // Specify variable-blocksize stream

    out->write(header_blocksize(blocksize), 4);       // Block size
    out->write(header_samplerate(samplerate), 4);     // Samplerate code
    out->write(midside ? 0b1000 : (channels - 1), 4); // Write channel assignment
    out->write(header_bitdepth(bitdepth), 3);         // Sample size in bits

    out->write(0, 1); // Mandatory bit

    // Write the sample number, UTF-8 encoded.
    writeUTF8(out.get(), samplestart);

    // Write the explicit block size value if necessary.
    if (header_blocksize(blocksize) == 0b0111)
        out->write(blocksize - 1, 16);

    // Write samplerate as well if necessary.
    if (header_samplerate(samplerate) == 0b1101)
        out->write(samplerate, 16);

    // Write CRC of the header.
    out->write(computeCRC<8>(out.get(), 0b100000111, bitpos), 8);
}

/* Write the FLAC subframeheader. Encodes the encoding method.
 * Only LPC or constant is implemented.
 */
void FLAC::subframeheader(int lpc_order) {
    out->write(0, 1); // Zero-bit padding

    // Write the subframe type.
    if (lpc_order)
        out->write(0b100000 | (lpc_order - 1), 6); // LPC
    else
        out->write(0, 6); // Constant

    // "Wasted bits-per-sample" flag. Note that wav-files
    // don't support wasted bits, so we set this to zero.
    out->write(0, 1);
}



/* Main writing loop. Reads all data from the given wav file
 * and writes the FLAC encoded output to the specified
 * destination file.
 */
template<typename type>
void FLAC::writedata(WavFile& wav, bool showprogress) {
    unsigned blocksize_min = 256;
    unsigned blocksize_default = 4096;
    unsigned blocksize_max = 65535;
    unsigned blocksize = blocksize_default;

    // Init/set some of the parameters.
    int order_max = order_default;
    int qlp_precision = 14;
    bool midside = false; // channels == 2;

    LPC<type> lpc;
    GolombRice<type> golomb;

    auto ringbuffers = make_shared<vector<RingBuffer<type>>>(wav.channels,
            RingBuffer<type>(blocksize_max));
    wav.readchunk(ringbuffers, blocksize_max, midside);

    unsigned long samplenumber = 0;
    unsigned bitpos = 0;

    // -- Uncomment for plotting of general encoding info -- //
    // -- (Also see actual plot functions calls) ----------- //
    // Plot plotresiduals("out/residuals-lpc" + to_string(lpctype)
    //         + "-order" + to_string(order_max), false);
    // Plot plotblocksize("out/blocksize-lpc" + to_string(lpctype)
    //         + "-order" + to_string(order_max), false);
    // Plot plotbitspersample("out/bitspersample-lpc" + to_string(lpctype)
    //         + "-order" + to_string(order_max), false);
    // -------------- //

    while (true) {

        if (showprogress) {
            printf("\rtime: %5.1f sec", samplenumber / (double) samplerate);
            fflush(stdout);
        }

        for (int c = 0; c < channels; c++) {
            shared_ptr<typename LPC<type>::Analysed> analysis;

            // The first channel determines the blocksize.
            if (c == 0) {
                // Do the LPC analysis.
                analysis = lpc.analyse((*ringbuffers)[0], qlp_precision, lpctype,
                                       order_max, blocksize_default,
                                       blocksize_min, (*ringbuffers)[0].filled());
                blocksize = analysis->size;

                // Write the frameheader.
                bitpos = out->getbitpos();
                frameheader(blocksize, midside, samplenumber, bitpos);

                samplenumber += blocksize;
            } else {
                // Do the LPC analysis.
                analysis = lpc.analyse((*ringbuffers)[c], qlp_precision, lpctype,
                                       order_max, blocksize);
            }

            // Write the subframe header.
            subframeheader(analysis->order);

            // Write warmup residuals.
            for (int i = 0; i < analysis->order; i++)
                out->write(analysis->warmup[i], bitdepth);

            // Write quantized LPC parameters.
            out->write(analysis->qprecision - 1, 4);
            out->write(analysis->qlevel, 5);
            for (int i = 0; i < analysis->order; i++)
                out->write((type) analysis->qparameters[i], analysis->qprecision);

            // Calcuate and set the optimal Rice parameter.
            double average = analysis->abs_sum / (double) analysis->residuals.size();
            double p = 1 / (average + 1);
            double k_guess = p >= 0.5 ? 0 : -log2(log2(1 / (1 - p)));
            golomb.setk(round(clamp(k_guess, 0.0, 32.0)));

            // -- Uncomment for soft pops at the block edges -- //
            // analysis->residuals[analysis->residuals.size() - 1] = (type) 35000;
            // ------------------------------------------------ //

            // Write residuals
            out->write(0b01, 2); // Use 5-bit RICE parameter
            out->write(0x00, 4); // Rice partition order
            out->write(golomb.getk() + 1, 5); // Rice parameter
            golomb.encode(analysis->residuals, out);

            // -- Uncomment for plotting of general encoding info -- //
            // -- (Also see the initialisations above) ------------- //
            // plotresiduals.write(p, " ");
            // plotresiduals.write(analysis->residuals, end=" ");
            // plotblocksize.write(blocksize, " ");
            // plotbitspersample.write(out->getvector()->size() * sizeof(uint) *
            //         8 / (double) blocksize / channels, " ");
            // ----------------------------------------------------- //
        }

        // Pad to byte allignment.
        out->writebytepadding();
        // Write the CRC-16 of the whole frame.
        out->write(computeCRC<16>(out.get(), 0b11000000000000101, bitpos), 16);

        // Write to file.
        out->flush2file();

        // Read next chunk from the wav-file.
        wav.readchunk(ringbuffers, blocksize, midside);

        // End of file has been reached.
        if ((*ringbuffers)[0].filled() < blocksize_min)
            break;
    }
}

void FLAC::encode(string source, string destination,
                  bool showprogress) {
    WavFile wav(source);
    out = make_shared<BitStreamOut>();
    out->open(destination);

    try {
        wav.open();
    } catch (const char* msg) {
        return;
    }

    samplerate = wav.samplerate;
    channels = wav.channels;
    bitdepth = wav.bitdepth;
    size = wav.size;

    // Write the FLAC header to file.
    header();
    out->flush2file();

    // Write the complete file.
    switch (bitdepth) {
        case 8:  writedata<char>(wav, showprogress); break;
        case 16: writedata<short>(wav, showprogress); break;
        case 24: writedata<int24>(wav, showprogress); break;
        case 32: writedata<int>(wav, showprogress); break;
    }
}


void FLAC::decode(string source, string destination) {
    in = make_shared<BitStreamIn>();
    in->open(source);
}
