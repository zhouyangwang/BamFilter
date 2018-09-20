#ifndef OPTIONS_H
#define OPTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include "htslib/sam.h"

using namespace std;


class Options{
public:
    Options();
    bool validate();

public:
    string QueryVcf;
    string BamCfdna;
    string BamGdna;
    string output;
    string HumanRepeatFile;

    int LineSkip;
    int MinQualDrop;
    bool OutputSimpleOverlap;
    bool FastInferReadSize;
    bool DropInconsist;
    bool ReadDropXA;
    bool ReadMaxMismatch;
    bool DuplexUMIStat;
    bool snp;
    bool indel;
    bool debug;

};

#endif
