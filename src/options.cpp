#include "options.h"
#include "util.h"
#include "helper.h"

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
    explain = true;
}

bool Options::validate() {

    if(explain){
        instruction::explain();
        return false;
    }

    if(QueryVcf.empty()){
        error_exit("vcf file must be provided");
    }else{
        check_file_valid(QueryVcf);
    }

    if(BamCfdna.empty() || BamGdna.empty()){
        error_exit("both --cfdna and --gdna should be specified");
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
        error_exit("Only one kind of mutation can be handled at a time");
    }else if(!snp && !indel){
        error_exit("--snp or --indel must be provided");
    }

    if (output.empty()){
        size_t location = QueryVcf.find(".txt");
        if(location != string::npos){
            output = QueryVcf.substr(0,location) + "_MrBam.txt";
        }else
            output = QueryVcf + "_MrBam.txt";
    }

    return true;
}
