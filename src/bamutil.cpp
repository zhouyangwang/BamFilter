#include "bamutil.h"
#include <sstream>
#include <vector>
#include <memory.h>

BamUtil::BamUtil(bam1_t *b, long pos, Options *opt){
    mopt = opt;
    PosFromVCF = pos;
    chromosome = getChr(b);
    reference_start =  b->core.pos;
    reference_end = bam_endpos(b)-1;
    // get 0-based end;
    query_name = getQName(b);
    mb = b;
}

BamUtil::BamUtil(BamUtil *readinfo, long pos){
    PosFromVCF = pos;
    mopt = readinfo->mopt;
    query_name = readinfo->query_name;

    reference_start = readinfo->reference_start;
    query_alignment_start = readinfo->query_alignment_start;
    query_length = readinfo->query_length;

    has_tag_XA = readinfo->has_tag_XA;
    template_length = readinfo->template_length;
    next_reference_start = readinfo->next_reference_start;
    is_reverse = readinfo->is_reverse;
    is_pair_mate_mapped = readinfo->is_pair_mate_mapped;

    averageQualOver10 = readinfo->averageQualOver10;
    copynumber = readinfo->copynumber;
    mapping_quality = readinfo->mapping_quality;

}

BamUtil::~BamUtil(){

}

void BamUtil::AddingMoreInforforRead(bam1_t *b, char refbase){

    char tagCN[3] = "CN";
    char tagXA[3] = "XA";
    char tagMD[3] = "MD";
    char tagNM[3] = "NM";

    // unchanged part for read with multiple mutation
    sequence = getSeq(b);
    quality = getQual(b);
    template_length = abs(b->core.isize);
    next_reference_start = b->core.mpos;
    is_reverse = b->core.flag & BAM_FREVERSE;
    is_pair_mate_mapped = (b->core.flag & BAM_FPAIRED) & ( !(b->core.flag & BAM_FUNMAP));
    mapping_quality = b->core.qual;
    query_length = sequence.size();
    //query length from the CIGAR alignment but does not include hard-clipped bases.

    getCigar(b);
    visualSeqlength = 0;
    getaligned_pairs();

    copynumber = tagInfoCalc(b,tagCN);
    nm = tagInfoCalc(b,tagNM);
    has_tag_XA = (tagInfoCalc(b,tagXA) == 1)? true:false;
    tagInfoCalc(b,tagMD);

    // changed part for read with multiple mutation
    referenceBase = refbase;
    varianttype = 'E';
    getquerypos(b);
    mismatched = ifmismatch();
    mismatchBaseQuality = (varianttype == 'D' || varianttype == 'I')? -1: int(quality[querypos]);
    mismatchbase = (varianttype == 'D' || varianttype == 'I')? varianttype:sequence[querypos];
}

void BamUtil::updateReadInfo(bam1_t *b,char refbase){
    referenceBase = refbase;
    varianttype = 'E';
    getquerypos(b);
    mismatched = ifmismatch();
    mismatchBaseQuality = (varianttype == 'D' || varianttype == 'I')? -1: int(quality[querypos]);
    mismatchbase = (varianttype == 'D' || varianttype == 'I')? varianttype:sequence[querypos];
}

int BamUtil::PrintAllReadInfo(){

    cout << query_name << "\n" << "base: " << mismatchbase <<"\tqual: "<< mismatchBaseQuality <<"\tstart:" << reference_start - query_alignment_start;
    cout << "\tquerylength: " << query_length << "\tmismatced: " << mismatched <<"\tXA: " << has_tag_XA << "\tTemplateLen: "<< template_length;
    cout << "\tnext ref start: " << next_reference_start << "\treverse: " << is_reverse << "\tpair and mapped: " << is_pair_mate_mapped;
    cout << "\tq10: " << averageQualOver10 << "\tterminal: " << mismatchAtTerminal << "\tcopyNumber: " << copynumber <<endl;
    if(mopt->debug){
        int i = 0;
        cerr << "cigar: " << cigarstring << "\tMD: " << MDString  << endl;
        while(i<aligned_pairs.size()){
            cerr <<i <<"\t" << aligned_pairs[i][0] <<"\t" << aligned_pairs[i][1]<<"\t"<< aligned_pairs[i][2] << endl;
            i++;
        }
    }
    return 1;
}

void BamUtil::getquerypos(bam1_t *b){
    int id  = 0;
    int nBaseBeforMut = 0;

    while(id < aligned_pairs.size()){
        if(aligned_pairs[id][2] != 4 && aligned_pairs[id][2] != 2){
            nBaseBeforMut++;
        }
        if(aligned_pairs[id][1] == PosFromVCF){
            querypos = aligned_pairs[id][0];
            if(id + 1 < aligned_pairs.size() && aligned_pairs[id+1][2] == 1){
                varianttype = 'I';
            }else if( sequence[querypos] == referenceBase){
                varianttype = 'M';
            }else {
                // Mismatch base
                varianttype = 'X';
            }
            break;

        }else if(aligned_pairs[id][1] > PosFromVCF){
            break;
        }
        id++;
    }
    idInAligned_pairs = id;
    varianttype = (varianttype == 'E')?'D':varianttype;
    varianttype = (varianttype == 'D' || varianttype == 'I')?varianttype: sequence[querypos];
    // removing deletions and softclipped bases before mutation site while keeping insertions left
    // mutation base wasn't included when calculating
    float percent = float(nBaseBeforMut) / visualSeqlength;
    mismatchAtTerminal = (varianttype == 'M')? -1:((percent > 0.1 && percent < 0.9)?false:true);
}

int BamUtil::tagInfoCalc(bam1_t *b, const char *tag){
    uint8_t* v= bam_aux_get(b,tag);
    stringstream ss;
    if(!strcmp(tag,"XA")){
        if(!v){
            return 0;
        }else{
            return 1;
        }
    }
    // copy from: https://github.com/pysam-developers/pysam/blob/master/pysam/libcalignedsegment.pyx
    if(v[0] == 'B'){
        ss << v[0];ss << v[1];
    }else{
        ss << v[0];
    }
    string auxtype = ss.str();
    if(auxtype == "I" || auxtype == "i" || auxtype == "C" || auxtype == "c" || auxtype == "S" || auxtype == "s"){
        int value = bam_aux2i(v);
        return value;

    }else if(auxtype == "f" || auxtype == "F"){
        double value = bam_aux2f(v);
        cout << auxtype << " value2 " << value << "\t" << b->core.qual << "\t"<< b->core.pos <<  query_name << endl;
        error_exit("unknown auxiliary type");

    }else if(auxtype == "d" || auxtype == "D"){
        double value = bam_aux2f(v);
        cout << auxtype << " value3 " << value << "\t" << b->core.qual << "\t"<< b->core.pos <<  query_name << endl;
        error_exit("unknown auxiliary type");

    }else if(auxtype == "A" || auxtype == "a"){
            // force A to a
            v[0] = 'A';
            // there might a more efficient way
            // to convert a char into a string
            char value = bam_aux2A(v);
            cout << auxtype << " value4 " << value << "\t" << b->core.qual << "\t"<< b->core.pos <<  query_name << endl;
            error_exit("unknown auxiliary type");

    }else if(auxtype == "Z" || auxtype == "H"){
            // Z and H are treated equally as strings in htslib
            char * value = bam_aux2Z(v);
            MDString = value;
            return 1;

    }else if(auxtype[0] == 'B'){
        //bytesize, nvalues, values = convert_binary_tag(v + 1)
        int value = v[1];
        cout << auxtype << " value6 " << value << "\t" << b->core.qual << "\t"<< b->core.pos <<  query_name << endl;
        error_exit("unknown auxiliary type");

    }else{
        error_exit("unknown auxiliary type");
    }
}

string BamUtil::getChr(bam1_t *b){
    string chr = "chr";
    stringstream ss;
    int chrid = b->core.tid;

    if(chrid == 22){
        ss << 'X';
    }else if(chrid == 23){
        ss << 'Y';
    }else{
        chrid++;
        ss << chrid;
    }
    chr += ss.str();
    return chr;
}

string BamUtil::getQName(bam1_t *b) {
    return string(bam_get_qname(b));
}

string BamUtil::getUMI(bam1_t *b, const string& prefix) {
    return getUMI(string(bam_get_qname(b)), prefix);
}

string BamUtil::getUMI(string qname, const string& prefix) {
    int len = qname.length();
    int sep = len-1;
    int prefixLen = prefix.length();

    bool foundSep = false;
    bool found= false;
    for(sep = len-1; sep>=0; sep--) {
        char c = qname[sep];
        if(c == ':') {
            foundSep = true;
            break;
        }
    }

    if(!foundSep || sep + prefixLen >=len-1)
        return "";


    // check prefix
    bool goodPrefix = true;
    for(int p=0; p<prefixLen; p++) {
        if(prefix[p] != qname[sep + 1 + p]) {
            goodPrefix = false;
            break;
        }
    }

    if(!goodPrefix)
        return "";

    int start = sep + 1 + prefixLen;
    if(start < len-1 && qname[start] == '_')
        start++;

    int numOfUnderscore = 0;
    for(int i=start; i<len; i++) {
        char c = qname[i];
        // UMI can be only A/T/C/G/N/_
        if(c != 'A' && c != 'T' && c != 'C' && c != 'G' && c != '_') {
            return "";
        }
        if(c == '_') {
            numOfUnderscore++;
            if(numOfUnderscore > 1)
                return "";
        }
    }
    return qname.substr(start, len-start);
}

string BamUtil::getQual(bam1_t *b) {
    uint8_t *data = bam_get_qual(b);
    int len = b->core.l_qseq;
    string s(len, '\0');
    int totalquality = 0;
    for(int i=0; i<len; i++) {
//        s[i] = (char)(data[i] + 33);
        s[i] = (char)(data[i]);
        totalquality += int(data[i]);
    }
    averageQualOver10 = (float(totalquality) / len >= 10)?true:false;
    return s;
}

string BamUtil::getSeq(bam1_t *b) {
    uint8_t *data = bam_get_seq(b);
    int len = b->core.l_qseq;
    string s(len, '\0');
    for(int i=0; i<len; i++) {
        char base;
        if(i%2 == 1)
            base = fourbits2base(data[i/2] & 0xF);
        else
            base = fourbits2base((data[i/2]>>4) & 0xF);
        s[i] = base;
    }
    return s;
}

//Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G, 8 for T and 15 for N.
char BamUtil::fourbits2base(uint8_t val) {
    switch(val) {
        case 1:
            return 'A';
        case 2:
            return 'C';
        case 4:
            return 'G';
        case 8:
            return 'T';
        case 15:
            return 'N';
        default:
            cerr << "ERROR: Wrong base with value "<< (int)val << endl ;
            return 'N';
    }
}
/*
@discussion In the CIGAR array, each element is a 32-bit integer. The
 lower 4 bits gives a CIGAR operation and the higher 28 bits keep the
 length of a CIGAR.
*/

void BamUtil::getCigar(bam1_t *b) {
    uint32_t *data = (uint32_t *)bam_get_cigar(b);
    int cigarNum = b->core.n_cigar;
    stringstream ss;
    for(int i=0; i<cigarNum; i++) {
        uint32_t val = data[i];
        char op = bam_cigar_opchr(val);
        uint32_t len = bam_cigar_oplen(val);
        ss << len << op;
        cigar.push_back(op);
        cigar.push_back(len);
        if(i==0){
            if(op == 'S'){
                query_alignment_start = len;
            }else{
                query_alignment_start = 0;
            }
        }
    }
    cigarstring = ss.str();
    /*
    cout << "cigar " << ss.str() <<endl;
    for(vector<int>::iterator iter1=cigar.begin();iter1!=cigar.end();iter1++){
        cout << *iter1 << "\t"; iter1++; cout << *iter1 << "\t";
    }
    cout <<endl;
    */
}

void BamUtil::getaligned_pairs(){
    int i=0;
    int nalt = 0;
    int nref = reference_start;
    char op; int len; int location = 0;
    nmismatchbases = 0; visualSeqlength = 0;

    while(location < cigar.size()){
        op = cigar[location];
        location++;
        len = cigar[location];
        location++;

        switch(op){
            case 'I':
            case 'D':
                nmismatchbases += len;
                nmismatchbases --;
        }
    // aligned_pairs [pos in seq, pos in ref, basetype]
    // 0:match 1:insertion 2:deletion 4:softclipped bases 10:mismatch
        while(len >0){
            len--;
            switch (op){
                case 'M':
                    aligned_pairs[i].push_back(nalt);
                    aligned_pairs[i].push_back(nref);
                    aligned_pairs[i].push_back(0);
                    i++; nalt++; nref++; visualSeqlength++;
                    break;
                case 'I':
                    aligned_pairs[i].push_back(nalt);
                    aligned_pairs[i].push_back(-1);
                    aligned_pairs[i].push_back(1);
                    i++; nalt++; visualSeqlength++;
                    break;
                case 'S':
                    aligned_pairs[i].push_back(nalt);
                    aligned_pairs[i].push_back(-1);
                    aligned_pairs[i].push_back(4);
                    i++; nalt++;
                    break;
                case 'D':
                    aligned_pairs[i].push_back(-1);
                    aligned_pairs[i].push_back(nref);
                    aligned_pairs[i].push_back(2);
                    i++; nref++;
                    break;
                case 'H':
                    break;
                default:
                    cout << "unseen cigar alphabet: " << op <<"\t"<< len << "\t" <<query_name << endl;
                    error_exit("unseen cigar alphabet");
            }
        }
    }
}

bool BamUtil::ifmismatch(){
    // 4S19M3D14M6D36M3D64M1I5M
    // MD:Z: 8G10^GAG9G4^TGGTGA36^CGG0 T11T8G4A5G3T1G2A1A17G7
    if(mopt->debug)
        updatealigned_pairs();

    if(nm - nmismatchbases > mopt->ReadMaxMismatch){
        return true;
    }else if(nm - nmismatchbases <= 1){
        return false;
    }

    updatealigned_pairs();
    // mismatch bases need to located around mutation, distance should be equal or less than 2;
    if(aligned_pairs.count(idInAligned_pairs) != 0){
        if(aligned_pairs[idInAligned_pairs][2] == 0){
            for(int i=idInAligned_pairs-2;i<idInAligned_pairs+3;i++){
                if(aligned_pairs.count(i) != 0){
                    if(aligned_pairs[i][2] != 0){
                        return false;
                    }
                }
            }
            return true;
        }
        return true;
    }
    return true;
}


void BamUtil::updatealigned_pairs(){
    bool deletion = false;
    map<int,char> baseCondition;
    int j =0; string id;
    for(int i=0;i<MDString.size();i++){
        switch(MDString[i]){
        case '^':
            if(id.size()!=0){
                j += str2int(id);
                id = "";
            }
            deletion = true;
            break;
        case 'A':
        case 'T':
        case 'G':
        case 'C':
        case 'N':
            if(id.size() !=0){
                j += str2int(id);
                id = "";
            }
            if(deletion){
                baseCondition[j] = 'D';
                deletion = true;
                j++;
                break;
            }
            deletion = false;
            baseCondition[j] = MDString[i];
            j++;
            break;
        case '0':
            if(id.size()!=0){
                id += MDString[i];
            }
            deletion = false;
            break;
        default:
            id += MDString[i];
            deletion = false;
        }
    }

    int number (0);
    map<int,char>::iterator itr1 = baseCondition.begin();
    for(int i=0;i<aligned_pairs.size();i++){
        if(aligned_pairs[i][1] == -1){
            // SI
            continue;
        }
        if(itr1 != baseCondition.end()){
            if(number == itr1->first && (itr1->second == 'A' || itr1->second == 'G' || itr1->second == 'C' || itr1->second == 'T' || itr1->second == 'N')){
                aligned_pairs[i][2] = 10;
                itr1++;
            }
        }
        number++;
    }
}

bool BamUtil::isPrimary(bam1_t *b) {
    if(b->core.flag & BAM_FSECONDARY || b->core.flag & BAM_FSUPPLEMENTARY)
        return false;
    else
        return true;
}

bool BamUtil::isProperPair(bam1_t *b) {
    return b->core.flag & BAM_FPROPER_PAIR;
}

int BamUtil::getRightRefPos(bam1_t *b) {
    if(b->core.pos<0)
        return -1;
    return b->core.pos + bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
}

bool BamUtil::test() {
    vector<string> qnames;
    qnames.push_back("NB551106:8:H5Y57BGX2:1:13304:3538:1404");
    qnames.push_back("NB551106:8:H5Y57BGX2:1:13304:3538:1404:UMI_GAGCATAC");
    qnames.push_back("NB551106:8:H5Y57BGX2:1:13304:3538:1404:UMI_GAGC_ATAC");
    qnames.push_back("NB551106:8:H5Y57BGX2:1:13304:3538:1404:GAGC_ATAC");
    qnames.push_back("NB551106:8:H5Y57BGX2:1:13304:3538:1404:UMI_X");

    vector<string> prefixes;
    prefixes.push_back("");
    prefixes.push_back("UMI");
    prefixes.push_back("UMI");
    prefixes.push_back("");
    prefixes.push_back("UMI");

    vector<string> umis;
    umis.push_back("");
    umis.push_back("GAGCATAC");
    umis.push_back("GAGC_ATAC");
    umis.push_back("GAGC_ATAC");
    umis.push_back("");

    for(int i=0; i<qnames.size(); i++) {
        string umi = getUMI(qnames[i], prefixes[i]);
        if(umi != umis[i]) {
            cerr << "get UMI from " << qnames[i] << ", expect " << umis[i] << ", but got " << umi << endl;
            return false;
        }
    }

    return true;

}
