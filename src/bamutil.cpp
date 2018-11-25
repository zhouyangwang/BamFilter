#include "bamutil.h"
#include <sstream>
#include <memory.h>

BamUtil::BamUtil(bam1_t *b, long pos, Options *opt){
    //InforCopy = false;
    mopt = opt;
    PosFromVCF = pos;
    chromosome = getChr(b);
    reference_start =  b->core.pos;
    reference_end = bam_endpos(b)-1;
    // get 0-based end;
    query_name = getQName(b);
    AddingMoreInforforRead(b);
}

BamUtil::BamUtil(BamUtil *readinfo, long pos){
    //InforCopy = true;
    PosFromVCF = pos;
    mopt = readinfo->mopt;

    query_name = readinfo->query_name;
    reference_start = readinfo->reference_start;
    reference_end = readinfo->reference_end;
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
    nmismatchbases = readinfo->nmismatchbases;
    aligned_pairs = readinfo->aligned_pairs;

    sequence = readinfo->sequence;
    quality = readinfo->quality;

    // debug
    if(mopt->debug){
        cigarstring = readinfo->cigarstring;
        MDString = readinfo->MDString;
        visualSeqlength = readinfo->visualSeqlength;
        nm = readinfo->nm;
    }
}

BamUtil::~BamUtil(){
    aligned_pairs.clear();
    // free memory;
}

void BamUtil::AddingMoreInforforRead(bam1_t *b){

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
    //read with secondary alignment is also hiden in GenomeBrowse
    query_length = sequence.size();
    //query length from the CIGAR alignment but does not include hard-clipped bases.
    isSecondary(b);
    if (mopt->DuplexUMIStat)
        getUMI();
    else
        umi = "";

    getCigar(b);
    visualSeqlength = 0;
    // visual bases: match, mismatched,and inserted bases, which meads deletion and softclipped bases are removed, intended for geneBrowse
    getaligned_pairs();

    copynumber = tagInfoCalc(b,tagCN);
    nm = tagInfoCalc(b,tagNM);
    has_tag_XA = (tagInfoCalc(b,tagXA) == 1)? true:false;
    tagInfoCalc(b,tagMD);
    getMismatchNumber();
    if(nmismatchbases == 2)
        updatealigned_pairs();
}

int BamUtil::PrintAllReadInfo(){
    cout << query_name << "\n" << "base: " << mismatchbase <<"\tqual: "<< mismatchBaseQuality <<"\tstart:" << reference_start;
    cout << "\tquerylength: " << query_length << "\n:mismatchbases " << nmismatchbases <<"\tXA: " << has_tag_XA << "\tTemplateLen: "<< template_length;
    cout << "\tnext ref start: " << next_reference_start << "\treverse: " << is_reverse << "\tpair and mapped: " << is_pair_mate_mapped;
    cout << "\tq10: " << averageQualOver10 << "\tterminal: " << mismatchAtTerminal << "\tcopyNumber: " << copynumber <<"\tmapping qual ";
    cout << mapping_quality <<endl;
    return 1;
}

void BamUtil::getquerypos(string refbase, bool continuous){
    int id  = 0;
    int nBaseBeforeMut = 0;
    referenceBase = refbase;
    varianttype = "D";
    // nBaseBeforeMut refers to length of first bases to mismatched bases,while softclipped bases are removed and mismatched bases is included
    // removing deletions and softclipped bases before mutation site while keeping insertions left
    // mutation base was included when calculating

    while(id < aligned_pairs.size()){
        if(aligned_pairs[id][2] != 4 && aligned_pairs[id][2] != 2){
            nBaseBeforeMut++;
        }
        if(aligned_pairs[id][1] == PosFromVCF){
            querypos = aligned_pairs[id][0];
            if(querypos != -1){
                if(id + 1 < aligned_pairs.size() && aligned_pairs[id+1][2] == 1){
                    varianttype = "I";

                }else if(continuous){
                    mismatchbase = sequence.substr(querypos,2);
                    if( mismatchbase == referenceBase){
                        varianttype = "M";
                    }else {
                        // Mismatch base
                        varianttype = "X";
                    }

                }else{
                    mismatchbase = sequence.substr(querypos,1);
                    if( mismatchbase == referenceBase){
                        varianttype = "M";
                    }else {
                        // Mismatch base
                        varianttype = "X";
                    }
                }
            }else{
                varianttype = "D";
            }
            break;
        }
        id++;
        if(id == aligned_pairs.size()){
            cout << "error in location " << query_name << "\t" << PosFromVCF << "\t" <<varianttype <<endl;
            exit(-1);
        }
    }
    idInAligned_pairs = id;
    mismatchAtTerminal = (varianttype == "M")? -1:((visualSeqlength*0.1<nBaseBeforeMut &&
                                                    visualSeqlength - nBaseBeforeMut + 1 > visualSeqlength*0.1)?false:true);
    mismatchBaseQuality = (varianttype == "D" || varianttype == "I")? -1:quality[querypos];
    mismatchbase = (varianttype == "D" || varianttype == "I") ? varianttype : mismatchbase;

    if( varianttype != "M" && nmismatchbases == 2 && !continuous){
        inferAnotherInfo();
    }else
        anotherMismatchLocation = "";

    if(mopt->debug && varianttype != "M" && nmismatchbases <= 2){
        cerr << "read name: " << query_name << "\tRange: " << reference_start <<"-" <<reference_end << "\tcigar: ";
        cerr << cigarstring << "\tMD: " << MDString << "\tNM: " << nm << "\tmutation base: "<< mismatchbase << endl;
        cerr << "FrontBases: " << nBaseBeforeMut << "\tBackBases: " << visualSeqlength - nBaseBeforeMut + 1;
        cerr << "\tvisualSeqlength: " << visualSeqlength << "\tterminal calculation: " << visualSeqlength*0.1 << "<" << nBaseBeforeMut << " && ";
        cerr << visualSeqlength - nBaseBeforeMut + 1 << ">" << visualSeqlength*0.1;
        if(mismatchAtTerminal)
            cerr << "\tmutations is located at terminal(10%) part of read\n";
        else
            cerr << "\tmutations is located at middle(11%-90%) part of read\n";

        int x(0);
        while(x<aligned_pairs.size()){
            if(aligned_pairs[x][2] != 0)
                cerr << aligned_pairs[x][0] <<"-"<< aligned_pairs[x][1]<<"-"<< aligned_pairs[x][2] << "\t";
            x++;
        }
        cerr << endl;
        if(anotherMismatchLocation.size() !=0)
            cerr << "Another mismatch happened, " << anotherMismatchLocation <<endl << endl;
    }
}


void BamUtil::inferAnotherInfo(){
    int id(0);
    while(id < aligned_pairs.size()){
        if(id != idInAligned_pairs && aligned_pairs[id][2] != 0 && aligned_pairs[id][2] != 4){
            if(id > idInAligned_pairs && aligned_pairs[id-1][2] != 0){
            // in case of a snp located after the indel,while the indel is query target
                id++;
                continue;
            }
            if(id < idInAligned_pairs)
                anotherMismatchLocation = int2string(idInAligned_pairs - id) + "L";
            else
                anotherMismatchLocation = int2string(id - idInAligned_pairs) + "R";
            break;
        }
        id++;
    }
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
            // Z and H are treated equally as strings in htslib, MD string
            char * value = bam_aux2Z(v);
            MDString = value;
            return 1;

    }else if(auxtype[0] == 'B'){
        //bytesize, nvalues, values = convert_binary_tag(v + 1)
        int value = v[1];
        cerr << auxtype << " value6 " << value << "\t" << b->core.qual << "\t"<< b->core.pos <<  query_name << endl;
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

void BamUtil::getUMI(){
    size_t st = query_name.find("UMI_");
    if(st == string::npos){
        cerr << "--UMI is offered, while unable to extract UMI from read's name: " << query_name << endl;
        cerr << "example: A00283:41:H5TLKDSXX:4:2174:6849:30452:UMI_GGA_GTA";
        exit(-1);
    }else{
        umi = query_name.substr(st+4,query_name.size()-st-4);
    }

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
    for(int i=0; i<cigarNum; i++){
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
}

void BamUtil::getaligned_pairs(){
    int i=0;
    int nalt = 0;
    int nref = reference_start;
    char op; int len; int location = 0;
    nmismatchbases = 0;
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

void BamUtil::getMismatchNumber(){
    // 4S19M3D14M6D36M3D64M1I5M
    // MD:Z: 8G10^GAG9G4^TGGTGA36^CGG0 T11T8G4A5G3T1G2A1A17G7
    nmismatchbases = nm - nmismatchbases;
}


void BamUtil::updatealigned_pairs(){
    // split "M" into mismatched base and matched base according to cigar and MD string
    bool deletion = false;
    map<int,char> baseCondition;
    int j =0; string id="";
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
                deletion = true;
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
        if(aligned_pairs[i][2] == 4 || aligned_pairs[i][2] == 2 || aligned_pairs[i][2] == 1 ){
            // softclipped bases, indel
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
    baseCondition.clear();
}

bool BamUtil::isPrimary(bam1_t *b) {
    if(b->core.flag & BAM_FSECONDARY || b->core.flag & BAM_FSUPPLEMENTARY)
        return false;
    else
        return true;
}

void BamUtil::isSecondary(bam1_t *b){
    if(b->core.flag & BAM_FSECONDARY)
        secodary_alignment = true;
    else
        secodary_alignment = false;
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
