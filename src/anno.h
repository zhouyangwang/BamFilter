#ifndef ANNO_H
#define ANNO_H
#include <stdio.h>
#include <stdlib.h>
#include "htslib/sam.h"
#include "options.h"
#include "vcfreader.h"
#include "bamutil.h"

class Anno{
public:
    Anno(Options *opt);
    ~Anno();

public:
    void annotation();
    void readcollect(bam1_t* b,long *reference_start, long *reference_end);
    void vcfinforoutput(string chromosome, long pos);
    void readsFilter();
    void CountDifferentType();
    void extraInfoStat();

public:
    vcfreader *vcfInfo;
    map<string, vector<BamUtil*>> readcache;

private:
    Options *mOptions;
    bam_hdr_t *mBamHeader;
    ofstream vcfOUT;
    string lastVCFChr;
    long lastVCFPos;
};


#endif // ANNO_H

