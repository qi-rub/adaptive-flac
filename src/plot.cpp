/* Simple plot-to-file of the given data.
 *
 * Author: Maxim van den Berg
 */

#include "include/plot.hpp"

using namespace std;

Plot::Plot(string filename, bool gnuplot) {
    file.open(filename);
    if (gnuplot)
        file << "plot '-' using 1:2 with lines" << endl;
}

Plot::~Plot() {
    file.close();
}

void Plot::plot(int n, double data, bool newline) {
    file << n << " " << data;
    if (newline)
        file << endl;
}

void Plot::plot(vector<double>& data, bool newline) {
    for (unsigned n = 0; n < data.size(); n++) {
        file << n << " " << data[n];
        if (newline)
            file << endl;
    }
}


void Plot::write(double data, string end) {
    file << data << end;
}

template<typename type>
void Plot::write(vector<type>& data, string end) {
    for (type d : data)
        file << (int) d << end;
    file << endl;
}

template<typename type>
void Plot::array(vector<type>& data) {
    file << "[";
    for (unsigned n = 0; n < data.size() - 1; n++)
        file << (int) data[n] << ", ";
    file << (int) data[data.size() - 1] << "]" << endl;
}

template void Plot::array(vector<int>& data);
template void Plot::array(vector<int24>& data);
template void Plot::array(vector<short>& data);
template void Plot::array(vector<char>& data);

template void Plot::write(vector<int>& data, string end);
template void Plot::write(vector<int24>& data, string end);
template void Plot::write(vector<short>& data, string end);
template void Plot::write(vector<char>& data, string end);
