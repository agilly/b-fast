# b-fast

A basic C++ library to parse (and hopefully one day write) BGEN files. This was started mainly out of curiosity and is currently in an unfinished state.

`bgen_parser.h` currently exposes:
* `file_info parse_header(ifstream & bgenfile, bool VERBOSE)`, which reads the header of a BGEN file open in an `ifstream` passed as argument. This is returned as a `file_info` struct.
* `BGEN_RETURN_CODE write_header(ofstream & bgenfile, file_info f_info)`, which writes a header contained in a `file_info` struct to an ofstream. Currently unfinished.

`bgenfun.cpp` contains a few lines of code to read and parse genotype information.

`helper.h` contains the `info(...)` and `error(...)` methods, as well as the combinatorics `choose(k,n)` function.
