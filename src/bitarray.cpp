/* A simple datastructure for dealing with bit-level specific
 * encoding.
 *
 * Author: Maxim van den Berg
 */

#include <stdexcept>
#include <cmath>
#include <iostream>

#include "include/bitarray.hpp"

using namespace std;

BitArray::BitArray() {
    typesize = sizeof (uint) * 8;
    bitpos = 0;
    reader = 0;
    bitpos_r = 0;

    if (array.size() == 0)
        array.push_back(0);
}

BitArray::BitArray(vector<uint> data, int start) {
    BitArray();

    array = data;
    bitpos = start;
}

void BitArray::write(uint input, int length) {
    if (length > typesize)
        throw invalid_argument("Argument 'length' specifies larger "
            "size than possible for vector type.");

    // Ensure input is of length `length`.
    if (length != 32)
        input &= (1 << length) - 1;

    int shift = typesize - bitpos - length;

    if (shift > 0) {
        // Simply shifts the input when it fits.
        array[array.size() - 1] |= input << shift;
        bitpos += length;
    } else {
        // And creates a new element when necessary.
        array[array.size() - 1] |= input >> -shift;
        array.push_back(shift ? input << (typesize + shift) : 0);

        bitpos = -shift;
    }
}

void BitArray::writeunary(int length) {
    while (length > 0) {
        write(0, min(typesize, length));
        length -= typesize;
    }
}

void BitArray::writebytepadding() {
    if (bitpos % 8)
        write(0, 8 - bitpos % 8);
}

void BitArray::setreader(unsigned long start) {
    reader = start / typesize;
    bitpos_r = start % typesize;
}

uint BitArray::read(int length) {
    if (length > typesize)
        throw invalid_argument("Argument 'length' specifies larger "
            "size than possible for vector type.");

    if (reader >= array.size() - 1 && length + bitpos_r > bitpos)
        throw range_error("Reading reached end of data");

    uint output = 0;
    int shift = typesize - bitpos_r - length;

    if (shift > 0) {
        // Read the `length` bits in one go.
        output = array[reader] >> shift;
        bitpos_r += length;
    } else {
        // Go to the next element when necessary.
        output = array[reader] << -shift;
        output |= shift ? array[reader + 1] >> (typesize + shift) : 0;

        reader++;
        bitpos_r = -shift;
    }

    // Ensure only `length` bits are returned.
    return output & ((1 << length) - 1);
}

uint BitArray::readunary() {
    uint output = 0;

    // Iterate bits until encountering a 1.
    while (((array[reader] >> (typesize - bitpos_r - 1)) & 1) == 0) {
        output++;
        bitpos_r += 1;

        // Continue to the next element.
        if (bitpos_r >= 32) {
            reader++;
            bitpos_r = 0;
        }

        if (reader >= array.size() - 1 && bitpos_r >= bitpos)
            throw range_error("Reading reached end of data");
    }

    return output;
}

void BitArray::append(BitArray& bitarray) {
    int size = bitarray.array.size();
    cout << " : " << size << " " << bitarray.getbitpos()<<endl;

    for (int i = 0; i < size - 1; i++)
        write(bitarray.array[i], 32);

    write(bitarray.array[size - 1], bitarray.getbitpos());
}

vector<uint>* BitArray::getvector() {
    return &array;
}

int BitArray::getbitpos() {
    return bitpos;
}



BitStreamOut::BitStreamOut() : BitArray() { }

void BitStreamOut::open(string filename) {
    ofs.open(filename, ios_base::binary);
}
void BitStreamOut::close() {
    ofs.close();
}

inline uint64_t swapendian(uint64_t v) { return __builtin_bswap64(v); }
inline uint32_t swapendian(uint32_t v) { return __builtin_bswap32(v); }
inline uint16_t swapendian(uint16_t v) { return __builtin_bswap16(v); }
inline uint8_t swapendian(uint8_t v) { return v; }

void BitStreamOut::flush2file() {
    uint buf = 0;

    // Write all completed elements.
    for (unsigned n = 0; n < array.size() - 1; n++) {
        buf = swapendian(array[n]);
        ofs.write((const char*) &buf, sizeof (uint));
    }

    // Write the < 4 remaining bytes.
    int i = 0;
    for (; i < bitpos; i += 8) {
        buf = (array[array.size() - 1] >> (24 - i)) & 0xff;
        ofs.write((const char*) &buf, 1);
    }

    // Reset the array (and keep remaining bits).
    bitpos = i - bitpos;
    array[0] = array[array.size() - 1] << i;
    array.resize(1);
    reader = 0;
}



BitStreamIn::BitStreamIn() : BitArray() { }

void BitStreamIn::open(string filename) {
    ifs.open(filename, ios_base::binary);
}
void BitStreamIn::close() {
    ifs.close();
}

void BitStreamIn::getfromfile(int amount) {
    uint data = 0;

    for (int i = 0; i < amount; i++) {
        if (ifs.eof())
            throw range_error("Reading reached end of file");

        ifs.read((char*) &data, 32);
        write(swapendian(data), 32);
    }
}

uint BitStreamIn::read(int length) {
    if (reader >= array.size() - 1 && length + bitpos_r > bitpos)
        getfromfile(1);

    return BitArray::read(length);
}

uint BitStreamIn::readunary() {
    uint output = 0;

    // Iterate bits until encountering a 1.
    while (((array[reader] >> (typesize - bitpos_r - 1)) & 1) == 0) {
        output++;
        bitpos_r += 1;

        // Continue to the next element.
        if (bitpos_r >= 32) {
            reader++;
            bitpos_r = 0;
        }

        if (reader >= array.size() - 1 && bitpos_r >= bitpos)
            getfromfile(1);
    }

    return output;
}

void BitStreamIn::flush() {
    array[0] = array[array.size() - 1];
    array.resize(1);
    reader = 0;
}
