#ifndef BAM_UTIL_H
#define BAM_UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "htslib/sam.h"
#include "options.h"
#include <map>

using namespace std;

class BamUtil {
public:
    BamUtil(bam1_t *b, long pos, Options *opt);
    BamUtil(BamUtil *readinfo, long pos);
    ~BamUtil();

public:
    static string getQName(bam1_t *b);
    static string getChr(bam1_t *b);

    string getSeq(bam1_t *b);
    string getQual(bam1_t *b);
    void getCigar(bam1_t *b);
    char fourbits2base(uint8_t val);
    short nmismatch(bam1_t *b);
    int tagInfoCalc(bam1_t *b, const char *tag);
    void getquerypos(bam1_t *b);

    bool continuous_bases();
    void getaligned_pairs();
    void updatealigned_pairs();
    bool ifmismatch();
    void AddingMoreInforforRead(bam1_t *b, char refBase);
    void updateReadInfo(bam1_t *b, char refBase);
    int PrintAllReadInfo();

    string getUMI(string qname, const string& prefix);
    string getUMI(bam1_t *b, const string& prefix);

    uint8_t base2fourbits(char base);
    void dump(bam1_t *b);
    bool isPartOf(bam1_t *part, bam1_t *whole, bool isLeft);
    void dumpHeader(bam_hdr_t* hdr);
    int getRefOffset(bam1_t *b, int bampos);
    void copyQName(bam1_t *from, bam1_t *to);
    bool isPrimary(bam1_t *b);
    bool isProperPair(bam1_t *b);
    int getRightRefPos(bam1_t *b);
    void getMOffsetAndLen(bam1_t *b, int& MOffset, int& MLen);

    bool test();

public:
    bam1_t *mb;
    Options *mopt;
    string chromosome;
    string query_name;
    long reference_start;
    long reference_end;

    int mismatchBaseQuality;
    char mismatchbase;
    int nmismatchbases;
    int query_alignment_start;
    int query_length;

    bool mismatched;
    bool has_tag_XA;
    int template_length;
    long next_reference_start;
    bool is_reverse;
    bool is_pair_mate_mapped;
    bool averageQualOver10;
    bool mismatchAtTerminal;
    // 1 true, 0 false; -1 for indel
    int copynumber;
    int nm;
    string MDString;
    int mapping_quality;
    int visualSeqlength;
    char varianttype;
    vector<int> cigar;
    string cigarstring;

    string sequence;
    bool copy;

private:
    long PosFromVCF;
    char referenceBase;
    int querypos;
    string quality;
    map <int,vector<int>> aligned_pairs;
    int idInAligned_pairs;
    // id,alt,ref,type

};

#endif
