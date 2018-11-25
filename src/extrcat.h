#ifndef EXTRA_INFO
#define EXTRA_INFO
#include <map>
#include <vector>
#include "options.h"
#include "util.h"
#include "bamutil.h"

using namespace std;


class otherinfo{

public:
    otherinfo(Options *opt,string refbase,string altbase);
    ~otherinfo();

public:
    void umi_dedup();
    void good_molecule_stat(map<string,vector<BamUtil*>>* name_dict);
    void mismatches_2_dedup();
    void stringTransfer();

public:
    Options * moptions;
    string mrefBase;
    string maltBase;

    int nq10;
    int nterminal;
    int good_molecule;
    int good_read;
    int nmulti;
    int noverlap_pe;
    int noverlap_se;
    int good_noverlap_pe;
    int good_noverlap_se;
    int numi;
    int nCN1;
    int nCN2;
    int ntotal;
    int nmapqual;
    int ave_mapqual;
    int totalInsertSize;
    int insertSize;
    int ave_insertSize;
    string mismatchCigar;
    string extrainfo;

    map<string,vector<vector<int>>> location2mismatch;
    map<string,vector<vector<int>>> location2UMI;
    vector<int> segment;
    // 5L < <1,2>,<1,3> >
    map<string,vector<vector<int>>> UMI;

};


#endif // EXTRA_INFO
