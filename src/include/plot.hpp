/* Simple plot-to-file of the given data. Automatically spaces data
 * points with a set interval of 1.
 *
 * When `gnuplot = true` is specified to the constructor, one can use
 *     `gnuplot -p < [filename]`
 * in bash to create a graph of the data.
 *
 * Author: Maxim van den Berg
 */

#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "general.hpp"

/* Class for writing general data to a given file,
 * specialized for creating plots.
 */
class Plot {
public:
    Plot(std::string filename, bool gnuplot);
    ~Plot();

    /* Writes the data as a gnuplot graph with integer x-axis,
     * specified as n. The vector input automatically assigns
     * x-axis values based on index.
     */
    void plot(int n, double data, bool newline=true);
    void plot(std::vector<double>& data, bool newline=false);

    /* Writes single value, or vector of values. Optional `end`
     * parameter adds given string after each written value.
     */
    void write(double data, std::string end="");
    template<typename T>
    void write(std::vector<T>& data, std::string end="");

    /* Write the given vector as a python array, adding commas
     * and parentheses [ and ], e.g. [0,-3,4].
     */
    template<typename T>
    void array(std::vector<T>& data);

private:
    std::ofstream file;
};
