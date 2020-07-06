/* Implementation of Golomb-Rice and Exponential Golomb codess.
 *
 * Author: Maxim van den Berg
 */

#include <cmath>
#include <bitset>
#include <iostream>

#include "include/golomb.hpp"

using namespace std;

/*****************************************/
/*** Golomb-Rice Encoding and Decoding ***/
/*****************************************/

template<typename type>
shared_ptr<BitArray> GolombRice<type>::encode(vector<type>& datastream) {
    auto bitstream = make_shared<BitArray>();
    encode(datastream, bitstream);

    return bitstream;
}

template<typename type>
void GolombRice<type>::encode(vector<type>& datastream,
                              shared_ptr<BitArray> bitstream) {
    unsigned abs = 0;
    for (type x : datastream) {
        abs =  (x < 0 ? (unsigned) -(x+1) : (unsigned) x);

        // Write i, 1 terminated.
        bitstream->writeunary(abs >> k);
        bitstream->write(1, 1);

        // Write r.
        bitstream->write(abs, k);

        // Write sign bit
        bitstream->write(x < 0 ? 1 : 0, 1);
    }
}

template<typename type>
shared_ptr<vector<type>> GolombRice<type>::decode(BitArray& bitstream) {
    auto datastream = make_shared<vector<type>>();
    int abs;

    try {
        while (true) {
            int i = bitstream.readunary();
            bitstream.read(1);
            unsigned int r = bitstream.read(k);
            int sign = bitstream.read(1) ? -1 : 1;

            abs = r | (i << k);
            datastream->push_back(sign ? -abs - 1 : abs);
        }
    } catch (const range_error& e) {

        // Reached the end of the datastream.
        return datastream;
    }
}

template<typename type>
void GolombRice<type>::setk(int parameter) { k = parameter; }
template<typename type>
int GolombRice<type>::getk() { return k; }

/* Instantiate the template classes of integer types. */
template class GolombRice<int>;
template class GolombRice<int24>;
template class GolombRice<short>;
template class GolombRice<char>;



/************************************************/
/*** Exponential Golomb Encoding and Decoding ***/
/************************************************/

template<typename type>
ExpGolomb<type>::ExpGolomb(int order) {
    k = order;
}

template<typename type>
shared_ptr<BitArray> ExpGolomb<type>::encode(vector<type>* datastream) {
    auto bitstream = make_shared<BitArray>();

    for (type x : (*datastream)) {
        // Write sign bit
        bitstream->write(x < 0 ? 1 : 0, 1);

        // Write i.
        int i = floor(log2(abs((double) x))); // Efficiency could be improved
        bitstream->writeunary(i);

        // Write x in i bits.
        bitstream->write(abs(x), i + 1);
    }

    return bitstream;
}

template<typename type>
shared_ptr<vector<type>> ExpGolomb<type>::decode(BitArray* bitstream) {
    auto datastream = make_shared<vector<type>>();

    try {
        while (true) {
            int sign = bitstream->read(1) ? -1 : 1;
            int i = bitstream->readunary();
            datastream->push_back(sign * bitstream->read(i + 1));
        }
    } catch (const range_error& e) {

        // Reached the end of the datastream.
        return datastream;
    }
}

template<typename type>
void ExpGolomb<type>::setk(int order) {
    k = order;
}

/* Instantiate the template classes of integer types. */
template class ExpGolomb<int>;
template class ExpGolomb<int24>;
template class ExpGolomb<short>;
template class ExpGolomb<char>;
