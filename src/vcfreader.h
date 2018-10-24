#ifndef VCF_READER_H
#define VCF_READER_H

#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <cstring>
#include <sstream>
#include <iostream>
#include "util.h"
#include "vcfreader.h"
#include "options.h"

using namespace std;

class vcfreader{
public:
    vcfreader(Options *opt);
    ~vcfreader();

public:
    void getvcfinfo();
    void getrepeatinfo();
    void fileExtract();

    map<string,vector<long*>> vcfPosChr;
    map<string,vector<string>> vcfLineSplit;
    map<string,vector<long*>> humanGenomeRepeat;

private:
    ifstream vcffp;
    ifstream repeatfp;
    Options * mOptions;

};
#endif // VCF_READER_H
