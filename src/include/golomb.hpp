/* Implementation of Golomb-Rice and Exponential Golomb codess.
 *
 * Author: Maxim van den Berg
 */

#pragma once

#include <vector>
#include <iostream>
#include <string>
#include <memory>

#include "bitarray.hpp"
#include "general.hpp"

using byte = unsigned char;


/* Golomb-Rice encoding and decoding, for a given Rice-Parameter k.
 * Uses the BitArray class for storing the encoded results.
 *
 * Template class T is assumed to be an integer type, such as `char`,
 * `short` or `int`.
 */
template<typename T>
class GolombRice {
public:
    /* Encode the given data vector using Golomb-Rice coding.
     * Returns a pointer to a allocated BitArray object, which
     * the user should `delete` when necessary.
     */
    std::shared_ptr<BitArray> encode(std::vector<T>& datastream);
    void encode(std::vector<T>& datastream,
                std::shared_ptr<BitArray> bitstream);

    /* Decodes the given bitarray encoded with Golomb-Rice
     * coding. Returns a pointer to generated array object,
     * which the user should `delete` when necessary.
     */
    std::shared_ptr<std::vector<T>> decode(BitArray& bitstream);

    /* Updates the exponential Rice parameter k.
     */
    void setk(int parameter);
    int getk();

private:
    /* The exponential Rice parameter. */
    int k;
};



/* Exponential Golomob encoding and decoding, for a given order k.
 * Uses the BitArray class for storing the encoded results.
 *
 * NOTE: Currently only implements order 0, and fails in encoding `0`.
 *
 * Template class T is assumed to be an integer type, such as `char`,
 * `short` or `int`.
 */
template<typename T>
class ExpGolomb {
public:
    /* The constructor sets the order k. */
    ExpGolomb(int order);

    /* Encode the given data vector using Exponential Golomb coding.
     * Returns a pointer to a allocated BitArray object, which
     * the user should `delete` when necessary.
     */
    std::shared_ptr<BitArray> encode(std::vector<T>* datastream);

    /* Decodes the given bitarray encoded with Exponential Golomb
     * coding. Returns a pointer to generated array object,
     * which the user should `delete` when necessary.
     */
    std::shared_ptr<std::vector<T>> decode(BitArray* bitstream);

    /* Updates the order k.
     */
    void setk(int order);

private:
    /* The order of the exponential code. */
    int k;
};

