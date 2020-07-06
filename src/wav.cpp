/* Class for .wav file decoding.
 * Support reading to fixed-size buffer or ringbuffer.
 *
 * Author: Maxim van den Berg
 */

#include <algorithm>

#include "include/wav.hpp"

#define vec2d(type) vector<vector<type>>

using namespace std;

void WavFile::readheader() {
    // Reading buffers of both 2 and 4 bytes
    char buf2[5] = { 0 };
    char buf4[5] = { 0 };

    if (!ifs.is_open())
        throw "file might not exist";

    // The file format: asserts it equals 'RIFF'.
    ifs.read(buf4, 4);
    if (string(buf4) != "RIFF")
        throw "unsupported file type";

    // The size of the data in bytes (plus 36 bytes for general data info).
    ifs.read(buf4, 4);
    unsigned sizebytes = *((unsigned *) buf4) - 36;

    // The file type: asserts it equals 'WAVE'.
    ifs.read(buf4, 4);
    if (string(buf4) != "WAVE")
        throw "invalid .wav file";


    // START OF BLOCK 1 //

    // The sub chunk ID for the first sub chunk: assert is equals 'fmt '.
    ifs.read(buf4, 4);
    if (string(buf4) != "fmt ")
        throw "invalid .wav file";

    // Size of the first sub chunk (skipped).
    ifs.read(buf4, 4);

    // The audio format. Should equal 1.
    ifs.read(buf2, 2);
    unsigned short format = *((unsigned short *) buf2);
    if (format != 1)
        throw "unsupported format";

    // Channel count
    ifs.read(buf2, 2);
    channels = *((unsigned short *) buf2);

    // Samplerate
    ifs.read(buf4, 4);
    samplerate = *((unsigned *) buf4);

    // Byte rate (= samplerate * channels * bitsPerSample/8)
    ifs.read(buf4, 4);

    // Block (align = channels * bitsPerSample/8)
    ifs.read(buf2, 2);

    // Bits per sample. Equals 8, 16, etc.
    ifs.read(buf2, 2);
    bitdepth = *((unsigned *) buf2);


    // START OF BLOCK 2 //

    // Searches until encountering 'data'.
    unsigned index = 0;
    while (true) {
        ifs.read(buf2, 2);
        index += 2;

        if (string(buf2) == "da") {
            ifs.read(buf2, 2);
            index += 2;
            if (string(buf2) == "ta")
                break;
        }

        if (index >= sizebytes)
            throw "missing data chunk";
    }


    // START OF DATA BLOCK //

    // The size of the "data" sub chunk.
    ifs.read(buf4, 4);
    size = *((unsigned *) buf4) / (bitdepth / 8);

    // Set start location.
    start = ifs.tellg();
}


template<typename type>
void read(ifstream& ifs, int start, int end,
          int readsize, shared_ptr<vec2d(type)> data) {
    char buffer[4] = { 0 };
    int channels = (*data).size();

    for (int n = start; n < end; n++) {
        // When end is reached, shrink all channels to the correct size.
        if (ifs.eof()) {
            for (int c = 0; c < channels; c++)
                (*data)[c].resize(n);
            return;
        }

        for (int c = 0; c < channels; c++)  {

            ifs.read(buffer, readsize);
            (*data)[c][n] = *((type*) buffer);
        }
    }
}

template<typename type>
void read(ifstream& ifs, int length, int readsize,
          shared_ptr<vector<RingBuffer<type>>> ringbuffer,
          bool midside) {
    char buffer[4] = { 0 };
    int channels = (*ringbuffer).size();

    for (int n = 0; n < length; n++) {
        if (ifs.eof())
            return;

        if (!midside) {
            for (int c = 0; c < channels; c++)  {
                ifs.read(buffer, readsize);
                (*ringbuffer)[c].write(*((type*) buffer));
            }
        } else {
            ifs.read(buffer, readsize);
            type left = *((type*) buffer);
            ifs.read(buffer, readsize);
            type right = *((type*) buffer);

            (*ringbuffer)[0].write(left);
            (*ringbuffer)[1].write(left - right);
        }
    }
}

template<typename type>
shared_ptr<vec2d(type)> WavFile::readchunk(int length) {
    auto data = make_shared<vec2d(type)>(channels, vector<type>(length));
    readchunk(data);
    return data;
}

template<>
void WavFile::readchunk(shared_ptr<vec2d(int)> buffer, int start) {
    read<int>(ifs, start, (*buffer)[0].size(), 4, buffer);
}
template<>
void WavFile::readchunk(shared_ptr<vec2d(int24)> buffer, int start) {
    read<int24>(ifs, start, (*buffer)[0].size(), 3, buffer);
}
template<>
void WavFile::readchunk(shared_ptr<vec2d(short)> buffer, int start) {
    read<short>(ifs, start, (*buffer)[0].size(), 2, buffer);
}
template<>
void WavFile::readchunk(shared_ptr<vec2d(char)> buffer, int start) {
        read<char>(ifs, start, (*buffer)[0].size(), 1, buffer);
}


template<>
void WavFile::readchunk(shared_ptr<vector<RingBuffer<int>>> ringbuffer,
                        int length, bool midside) {
    read<int>(ifs, length, 4, ringbuffer, midside);
}
template<>
void WavFile::readchunk(shared_ptr<vector<RingBuffer<int24>>> ringbuffer,
                        int length, bool midside) {
    read<int24>(ifs, length, 3, ringbuffer, midside);
}
template<>
void WavFile::readchunk(shared_ptr<vector<RingBuffer<short>>> ringbuffer,
                        int length, bool midside) {
    read<short>(ifs, length, 2, ringbuffer, midside);
}
template<>
void WavFile::readchunk(shared_ptr<vector<RingBuffer<char>>> ringbuffer,
                        int length, bool midside) {
    read<char>(ifs, length, 1, ringbuffer, midside);
}


WavFile::WavFile(string file) {
    filename = file;
}

void WavFile::open() {
    ifs.open(filename);

    try {
        readheader();
    } catch (const char* msg) {
        ifs.close();
        cerr << "Error loading file " << filename << " "
            << msg << endl;
        throw "Could not load file";
    }
}

WavFile::~WavFile() {
    ifs.close();
}

template shared_ptr<vec2d(int)> WavFile::readchunk<int>(int length);
template shared_ptr<vec2d(int24)> WavFile::readchunk<int24>(int length);
template shared_ptr<vec2d(short)> WavFile::readchunk<short>(int length);
template shared_ptr<vec2d(char)> WavFile::readchunk<char>(int length);
