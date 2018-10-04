Raw data format change (from binary files to root files):

Compilation: g++ binary_readout.c `root-config --cflags --glibs` -o binary_readout
Usage: ./binary_readout <filename>

Estimation of channel time resolution

Compilation:  g++ Analyzer.cc `root-config --cflags --glibs` -o Analyzer -I .
Usage: ./Analyzer <filename>
