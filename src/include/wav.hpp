/* Class for .wav file decoding.
 * Support reading to fixed-size buffer or ringbuffer.
 *
 * Author: Maxim van den Berg
 */

#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <memory>

#include "general.hpp"
#include "ringbuffer.hpp"

/* Simple wav file decoder.
 * Source and breakdown: http://soundfile.sapp.org/doc/WaveFormat/
 *
 * Assumes little endian notation.
 */
class WavFile {
public:
    /* The constructor sets the filename internally. */
    WavFile(std::string filename);
    ~WavFile();

    /* Opens the filestream and calls `readheader`.
     * (see below). */
    void open();

    /* Tries to read a chunk of `length` audio samples, for every audio
     * channel. Fills the specified 2 dimensional vector with a list of
     * audio samples for all channels.
     *
     * The length of every channel vector is lower than `length` only if
     * the end of the audio file is reached, in which case it fill the
     * vector *only* until the last value.
     */
    template<typename T>
    std::shared_ptr<std::vector<std::vector<T>>> readchunk(int length);

    /* Same as `readchunk(int)`, but stores the data into `buffer`.
     * Begins writing in the buffer at the specified position `start`.
     */
    template<typename T>
    void readchunk(std::shared_ptr<std::vector<std::vector<T>>> buffer,
                   int start = 0);

    /* Same as `readchunk(int)`, but stores the data into the specified
     * ringbuffers. Writes for `length` amount of samples.
     */
    template<typename T>
    void readchunk(std::shared_ptr<std::vector<RingBuffer<T>>> ringbuffers,
                   int length, bool midside = false);

    /* File variables.
     * - size:       The total amount of samples
     * - channels:   Amount of audio channels
     * - samplerate: Amount of samples per second
     * - bitdepth:   The amount of bytes per sample (short, int, etc.)
     * - filename:   The name of the audio file.
     */
    unsigned size;
    unsigned short channels;
    unsigned samplerate;
    unsigned bitdepth;
    std::string filename;

private:
    /* Reads and processes the wav file header.
     * Sets ifstream to start of the `data` block, so
     * `readchunk` can be called immediatelly after.
     */
    void readheader();

    /* File stream. Is open as long as the WavFile object exists.  */
    std::ifstream ifs;

    /* Location (in bytes) where the `data` block begins.  */
    unsigned start;
};
