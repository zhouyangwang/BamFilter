#include "vcfreader.h"
using namespace std;

vcfreader::vcfreader(Options *opt){
    mOptions = opt;
    vcffp.open(mOptions->QueryVcf);
    if(!vcffp.is_open()){
        error_exit("unable to open vcf file");
    }
    if(!(mOptions->HumanRepeatFile).empty()){
        repeatfp.open(mOptions->HumanRepeatFile);
        if(!repeatfp.is_open()){
            error_exit("repeat file is offered, but can't be opened");
        }
    }
}

vcfreader::~vcfreader(){
    vcffp.close();
    if(repeatfp.is_open())
        repeatfp.close();
}

void vcfreader::getvcfinfo(){

    string line;
    string LastChr = "";
    string lineId = "";
    string mutationId = "";
    long LastPos = -1;

    if(mOptions->LineSkip != 0){
        for(int i=0;i< (mOptions->LineSkip);i++){
            getline(vcffp,line,'\n');
        }
    }
    line = "";

   while(getline(vcffp, line,'\n')){

        vector<string> linesplit;
        if(split(line,linesplit,"\t")){
            long *pos = new long (str2long(linesplit[1]) - 1);
            // pos in vcf is 1-based

            string chr = linesplit[0];
            string alt = linesplit[4];
            lineId = chr + long2string(*pos) + alt;
            mutationId = chr + long2string(*pos);
            vcfAlt[mutationId].push_back(alt);

            // check if the vcf file is sorted
            if(chr != LastChr){
                LastChr = chr;
                LastPos = *pos;
                vcfPosChr[chr].push_back(pos);
            }else{
                if(*pos < LastPos){
                    error_exit("vcf file is unsorted\n");
                }else if(*pos != LastPos){
                    vcfPosChr[chr].push_back(pos);
                }
                LastPos = *pos;
            }
            vcfLineSplit[lineId] = linesplit;

        }else{
            cerr << line << " is unable to split\n";
            exit(-1);
       }
    }
}

void vcfreader::getrepeatinfo(){
    string line = "";
    vector<string> linesplit;
    while(getline(repeatfp,line,'\n')){
        if(split(line,linesplit,"\t"))
        {
            long * refstart = new long (str2long(linesplit[1]));
            long * refend = new long (str2long(linesplit[2]));

            humanGenomeRepeat[linesplit[0]].push_back(refstart);
            humanGenomeRepeat[linesplit[0]].push_back(refend);

            while(linesplit.size()!=0){
                linesplit.pop_back();
            }
        }else{
            error_exit("Repeat file is unable to split");
        }
    }
}

void vcfreader::fileExtract(){
    //vcfreader::getrepeatinfo();
    vcfreader::getvcfinfo();
}
