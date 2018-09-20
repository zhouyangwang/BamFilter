#include "options.h"
#include "util.h"

Options::Options(){
    QueryVcf = "";
    BamCfdna = "";
    BamGdna = "";
    output = "";
    HumanRepeatFile = "";

    LineSkip = 0;
    MinQualDrop = 25;
    OutputSimpleOverlap =  false;
    FastInferReadSize = false;
    DropInconsist = false;
    ReadDropXA = false;
    ReadMaxMismatch = -1;
    DuplexUMIStat = false;
    snp = false;
    indel = false;
    debug = false;
}

bool Options::validate() {

    if(QueryVcf.empty()){
        error_exit("vcf file must be provided");
    }else{
        check_file_valid(QueryVcf);
    }

    if(BamCfdna.empty() && BamGdna.empty()){
        error_exit("At least one of --cfdna and --gdna should be specified");
    }else{
        if(!BamCfdna.empty()){
            string BaiCfdna = BamCfdna + ".bai";
            check_file_valid(BaiCfdna);
            check_file_valid(BamCfdna);
        }
        if(!BamGdna.empty()){
            string BaiGdna = BamGdna + ".bai";
            check_file_valid(BaiGdna);
            check_file_valid(BamGdna);
        }
    }

    if(snp && indel){
        error_exit("Only one kind of mutation can be dealed with at a time");
    }
    return true;
}
