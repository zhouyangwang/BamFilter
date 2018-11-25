#include "count.h"

Counter::Counter(Options *opt,string refbase,string altbase){
    mopt = opt;
    if(mopt->indel){
        if(refbase == "-" || altbase.size() >= 2){
            mrefBase = "ATGC";
            maltbase = "I";
        }else if(altbase == "-" || refbase.size() >= 2){
            mrefBase = "ATGC";
            maltbase = "D";
        }
    }else{
        mrefBase = refbase;
        maltbase = altbase;
    }
    mor = 0, mnr = 0, msr = 0, oor = 0, onr = 0, osr = 0, moa = 0, mna = 0, msa = 0, ooa = 0, ona = 0, osa = 0, inconsis =0;
}

Counter::~Counter(){
}

void Counter::unique_pair_count(map<vector<long>,vector<string>> *pairs){
    map<vector<long>,vector<string>>::iterator itr1 = pairs->begin();
    map<string, vector<long>> base;
    map<string, vector<long>>::iterator itr2;

    for(;itr1!=pairs->end();itr1++){
        for(long i=0;i<itr1->second.size();i++)
            base[itr1->second[i]].push_back(1);

        for(itr2=base.begin(); itr2!=base.end(); itr2++){
            if(itr2->second.size() == 1){
                if(stFind(itr2->first,maltbase)){
                    if(itr1->first[2] == 1)
                        ooa++;
                    else
                        ona++;

                }else if(stFind(itr2->first, mrefBase)){
                    if(itr1->first[2] == 1)
                        oor++;
                    else
                        onr++;
                }

            }else{
                if(stFind(itr2->first,maltbase)){
                    if(itr1->first[2] == 1)
                        moa++;
                    else
                        mna++;

                }else if(stFind(itr2->first, mrefBase)){
                    if(itr1->first[2] == 1)
                        mor++;
                    else
                        mnr++;
                }
            }
        }
        base.clear();
    }
}

void Counter::unique_single_count(map<vector<long>,vector<string>> *single){
    map<vector<long>,vector<string>>::iterator itr1 = single->begin();
    map<string, vector<long>> base;
    map<string, vector<long>>::iterator itr2;

    for(;itr1!=single->end();itr1++){
        for(long i=0;i<itr1->second.size();i++)
            base[itr1->second[i]].push_back(1);

        for(itr2=base.begin(); itr2!=base.end(); itr2++){
            if(itr2->second.size() == 1){
                if(stFind(itr2->first,maltbase)){
                    osa++;
                }else if(stFind(itr2->first, mrefBase)){
                    osr++;
                }

            }else{
                if(stFind(itr2->first,maltbase)){
                    msa++;
                }else if(stFind(itr2->first, mrefBase)){
                    msr++;
                }
            }
        }
        base.clear();
    }
    stringTransfer();
}

void Counter::stringTransfer(){
    int arr[13] ={mor, mnr, msr, oor, onr, osr, moa, mna, msa, ooa, ona, osa, inconsis};
    overlap = "";
    for(int i=0;i<11;i++){
        overlap += int2string(arr[i]) + ",";
    }
    overlap += int2string(arr[12]);
}




