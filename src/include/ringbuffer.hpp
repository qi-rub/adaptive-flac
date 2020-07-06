/* Simple ring buffer implementation for arbitrary types.
 *
 * Header only class.
 */

#pragma once

#include <vector>
#include <stdexcept>

#include "general.hpp"

/*
 * The ringbuffer offers basically no projection, so
 * it should be handled with care. For instance, reading and
 * repositioning of the reading pointer can be done irregardless
 * of whether data has been written to these sections.
 */
template<typename T>
class RingBuffer {
public:
    RingBuffer(unsigned size);


    void reposition(unsigned n);
    unsigned position();

    void write(T item);

    unsigned size();
    unsigned filled();

    T& operator[](int n);

private:
    std::vector<T> buffer;

    int sizemax;
    int start;
    int end;
};

using namespace std;

template<typename type>
RingBuffer<type>::RingBuffer(unsigned size) {
    sizemax = size + 1;
    buffer = vector<type>(sizemax);

    start = 0;
    end = 0;
}


template<typename type>
void RingBuffer<type>::reposition(unsigned n) {
    start = n % sizemax;
}

template<typename type>
unsigned RingBuffer<type>::position() {
    return start;
}

template<typename type>
void RingBuffer<type>::write(type item) {
    if ((end + 1) % sizemax == start) {
        throw length_error("Ringbuffer overflow");
    } else {
        buffer[end] = item;
        end = (end + 1) % sizemax;
    }
}

template<typename type>
unsigned RingBuffer<type>::size() {
    return sizemax;
}

template<typename type>
unsigned RingBuffer<type>::filled() {
    return (end - start + sizemax) % sizemax;
}

template<typename type>
type& RingBuffer<type>::operator[](int n) {
    return buffer[(start + n) % sizemax];
}
