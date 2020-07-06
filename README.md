## Linear prediction algorithms with adaptive blocksizes in the FLAC format.

This repository contains the source code used in my bachelor thesis, where I
describe an alternative approach to linear predictive coding in which the standard
Burg algorithm is modified to efficiently determine suitable blocksizes considering
the analysed signal in question.

Two adaptive algorithms are presented, the second a modification of the one described
above, which improves on compression speed.

These algorithms are implemented in the FLAC format specification, together with
the covariance Burg algorithm and the standard Levinson-Durbin method used by FLAC.
All FLAC-enabled audio players should be able to play the resulting files.

Note that only a subset of the FLAC feature set is implemented.

For more information, see `report/thesis.pdf`.

### Building
Simply run `make` to build the encoder.


### Running
One can simply the encoder in bash with

```
./encoder [filename] [order]
```

This will compress the file with all four algorithms, and print some statistics.

The order is the linear prediction order, which is fixed across the encoder
and used by all algorithms. Specifying the order is optimal, and the default is 8.

Output is saved to the `out/` folder, as the following four files:

```
out.flac-burg-adaptive-autocor[order]
out.flac-burg-adaptive[order]
out.flac-burg-fixed[order]
out.flac-levinson-fixed[order]
```

### Experiments
Apart from analysing the resulting file sizes and encoding times, more statistics can
be gathered into the `out/` folder, as well as to the `plots/` folder, if certain sections
in `main.cpp` (for file sizes and compression times), `flac.cpp` (for blocksizes and residuals) or `lattice.cpp` (for graphs of the error function) are uncommented.
This includes information about blocksizes, residuals, error plotting (for the adaptive algorithms), etc.

Note that the FLAC reference encoder has an useful build in analysis tool, which can be run with the following command.

```
flac -a --residual-text [filename]
```


__Author__

Maxim van den Berg

contact: `maximvdberg@gmail.com`


__Supervisors__
- Michael Walter
- Freek Witteveen
- Taco Walstra
