#ifndef ANNO_H
#define ANNO_H
#include <stdio.h>
#include <stdlib.h>
#include "htslib/sam.h"
#include "options.h"
#include "vcfreader.h"
#include "bamutil.h"
#include "count.h"
#include "extrcat.h"

class Anno{
public:
    Anno(Options *opt);
    ~Anno();

public:
    void annotation();
    void readcollect(bam1_t* b,long *reference_start, long *reference_end);
    void vcfinforoutput(string chromosome,int ncontinuous);
    bool readsFilter(string ID);
    void CountDifferentType();
    void extraInfoStat();
    void readCombine(string chromosome, string ID, string continouID,long pos);
    void finaloutput(string newID, int j);
    void readsFree(string id, string nextid="");
    void readTempVariableFree();

    void singleSNPMutation(string idInLine,string idInAlt, bool altFromContinuous = true);
    void singleINDELMutation(string idInLine,string idInAlt);

public:
    vcfreader *vcfInfo;
    map<string, vector<BamUtil*>> readcache;
    map<string,vector<BamUtil*>> name_dict;
    map<vector<long>,vector<string>> unique_pairs;
    map<vector<long>,vector<string>> unique_single;
    vector<long> posid;
    string posbase;
    map<string,vector<char>> snv;
    vector<string> mutation;

    string refBase;
    string altBase;

private:
    Options *mOptions;
    bam_hdr_t *mBamHeader;
    ofstream vcfOUT;
    string lastVCFChr;
    long lastVCFPos;
    bool continuous;
    bool cfdna;

    string idInAlt;
    string idInLine;
    string nextidInAlt;
    string nextidInLine;
    string newidInLine;

};


#endif // ANNO_H

