#ifndef count_different_type
#define count_diffenent_type
#include <map>
#include <vector>
#include "options.h"
#include "util.h"
using namespace std;

class Counter{
public:
    Counter(Options *opt,string refbase,string altbase);
    ~Counter();

public:
    void unique_pair_count(map<vector<long>,vector<string>> *pairs);
    void unique_single_count(map<vector<long>,vector<string>> *single);
    void stringTransfer();

public:
    int mor, mnr, msr, oor, onr, osr, moa, mna, msa, ooa, ona, osa, inconsis;
    string mrefBase;
    string maltbase;
    Options * mopt;
    string overlap;

};

#endif // count_different_type
