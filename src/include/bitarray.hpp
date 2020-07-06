/* A simple datastructure for dealing with bit-level specific
 * encoding.
 *
 * Author: Maxim van den Berg
 */

#pragma once

#include <string>
#include <vector>
#include <fstream>

typedef uint32_t uint;

/* Wrapper around std::vector which allows for writing and
 * reading on bit-level. Writes from left to right.
 */
class BitArray {

public:
    BitArray();

    /* Initialises the internal vector with the specified value, with
     * a starting bit position. Useful for reading binary content.
     */
    BitArray(std::vector<uint> data, int start);

    /* Writes the first `length` bits of the specified input
     * to the array. `writeunary` writes `length` amount of 0's.
     * `writebytepadding` fills the current byte with zeros.
     */
    void write(uint input, int length);
    void writeunary(int length);
    void writebytepadding();

    /* Reading is done in two steps:
     * - Set the reading pointer to the n'th bit using `setreader(n)`
     * - Read n bits using `read(n)`, which automatically updates
     *   the reader pointer.
     *
     * Calling `readunary` will read until encountering a 1 (or the
     * end of the vector), and return the amount of zero's found.
     * Note that this does *not* read past the 1.
     */
    void setreader(unsigned long start);
    uint read(int length);
    uint readunary();

    /* Returns the content of the internal vector. Note that the last
     * element does not necessarily have to be filled completely.
     */
    std::vector<uint>* getvector();

    /* Returns the current writing bit position.  */
    int getbitpos();

    /* Appends the given BitArray. */
    void append(BitArray& bitarray);

protected:
    std::vector<uint> array;

    /* Position of the bit which will be the next to which is written,
     * in the last vector element.
     */
    int bitpos;

    /* Position of the reader pointer.
     * The variables refer to the index in the vector and the bit
     * position  respectively.
     */
    unsigned long reader;
    int bitpos_r;

    /* The length of the array data type (unsigned int) in bits.
     * Equal to the maximum value of `bitpos` and `bitpos_r`
     */
    int typesize;
};


/* Similar to the BitArray class, but with additional functions
 * for writing vector elements bit.
 */
class BitStreamOut : public BitArray {
public:
    BitStreamOut();

    /* Open and close the file for writing. */
    void open(std::string filename);
    void close();

    /* Writes all fully filled bytes to the file and
     * removes them immediately thereafter.
     */
    void flush2file();

private:

    std::ofstream ofs;
};


/* Similar to the BitArray class, but with additional functions
 * for reading elements directly from a file.
 */
class BitStreamIn : public BitArray {
public:
    BitStreamIn();

    /* Open and close the file for writing. */
    void open(std::string filename);
    void close();

    /* Redefine the reader functions to read from the file
     * when necessary.
     */
    uint read(int length);
    uint readunary();

    /* Read the specified amount of (unsigned) ints from the
     * file and appends them to the data vector.
     */
    void getfromfile(int amount);

    /* Removes all fully filled vector elements ands resets
     * the reader position.
     */
    void flush();
private:

    std::ifstream ifs;
};
