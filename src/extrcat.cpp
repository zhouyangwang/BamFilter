#include "extrcat.h"

otherinfo::otherinfo(Options *opt,string refbase,string altbase){
    moptions = opt;
    if(moptions->indel){
        if(refbase == "-" || altbase.size() >= 2){
            mrefBase = "ATGC";
            maltBase = "I";
        }else if(altbase == "-" || refbase.size() >= 2){
            mrefBase = "ATGC";
            maltBase = "D";
        }
    }else{
        mrefBase = refbase;
        maltBase = altbase;
    }
    nq10 = 0;
    nterminal = 0;
    good_molecule = 0;
    good_read = 0;
    nmulti = 0;
    noverlap_pe = 0;
    noverlap_se = 0;
    good_noverlap_pe = 0;
    good_noverlap_se = 0;
    numi = 0;
    nCN1 = 0;
    nCN2 = 0;
    ntotal = 0;
    nmapqual = 0;
    ave_mapqual = 0;
    totalInsertSize = 0;
    insertSize = 0;
    ave_insertSize = 0;
    mismatchCigar = "";
}

otherinfo::~otherinfo(){
    location2mismatch.clear();
    location2UMI.clear();
    vector<int>().swap(segment);
}


void otherinfo::good_molecule_stat(map<string,vector<BamUtil*>>* name_dict){
    const int MaxInsertSize = 1000;
    map<string,vector<BamUtil*>>::iterator itr1 = name_dict->begin();

    if(moptions->debug)
        cerr << "list of good molecules: "<< endl;

    for(;itr1!=name_dict->end();itr1++){
        if(itr1->second.size() == 1){
            if(itr1->second[0]->mismatchbase == maltBase){
                int i = 0;
                if(itr1->second[i]->mismatchAtTerminal)
                    nterminal++;
                if(itr1->second[i]->averageQualOver10){
                    nq10++;
                    noverlap_se++;
                }
                if(itr1->second[i]->has_tag_XA && itr1->second[i]->mapping_quality ==0)
                    nmulti++;

                if(itr1->second[i]->averageQualOver10 && !itr1->second[0]->mismatchAtTerminal &&
                    (!itr1->second[i]->has_tag_XA || (itr1->second[i]->has_tag_XA &&
                    (itr1->second[i]->mapping_quality >0 && !itr1->second[i]->secodary_alignment)))){
                        good_molecule++;
                        good_read++;
                        good_noverlap_se++;
                        if(moptions->debug){
                            cerr << "Good Molecule: SE: " << itr1->second[0]->query_name << endl;
                        }

                        segment.push_back(itr1->second[i]->reference_start);
                        segment.push_back(itr1->second[i]->reference_end);
                        if(itr1->second[i]->anotherMismatchLocation.size() !=0)
                            location2mismatch[itr1->second[i]->anotherMismatchLocation].push_back(segment);
                        segment.clear();

                       if(itr1->second[i]->copynumber >1)
                            nCN2++;
                        else
                            nCN1++;

                        nmapqual += itr1->second[i]->mapping_quality;
                        if (itr1->second[i]->template_length< MaxInsertSize)
                            totalInsertSize += itr1->second[i]->template_length;
                    }
                }
        }else{
            for(int i=0;i<2;i++){
                if(itr1->second[i]->mismatchbase == maltBase){

                    if(itr1->second[i]->mismatchAtTerminal)
                        nterminal++;
                    if(itr1->second[i]->averageQualOver10){
                        nq10++;
                        if(i == 0)
                            noverlap_pe++;
                    }
                    if(itr1->second[i]->has_tag_XA && itr1->second[i]->mapping_quality ==0)
                        nmulti++;

                    if(i == 0 && itr1->second[0]->averageQualOver10 && itr1->second[1]->averageQualOver10){
                        if(!itr1->second[0]->mismatchAtTerminal || !itr1->second[1]->mismatchAtTerminal){
                            if( !(itr1->second[0]->has_tag_XA && itr1->second[0]->mapping_quality == 0 &&
                               itr1->second[1]->has_tag_XA && itr1->second[1]->mapping_quality == 0 &&
                               !itr1->second[0]->secodary_alignment && !itr1->second[1]->secodary_alignment)){
                                good_molecule++;
                                good_read += 2;
                                good_noverlap_pe++;
                                if(moptions->debug){
                                    cerr << "Good Molecule: PE: " << itr1->second[0]->query_name << endl;
                                }
                                for(int j=0;j<2;j++){
                                    segment.push_back(itr1->second[i]->reference_start);
                                    segment.push_back(itr1->second[i]->reference_end);
                                    if(itr1->second[j]->anotherMismatchLocation.size() !=0){
                                        location2mismatch[itr1->second[j]->anotherMismatchLocation].push_back(segment);
                                    }
                                    nmapqual += itr1->second[j]->mapping_quality ;
                                    if (itr1->second[j]->template_length< MaxInsertSize)
                                        totalInsertSize += itr1->second[j]->template_length;
                                    segment.clear();
                                }
                                if(itr1->second[i]->copynumber >1)
                                    nCN2++;
                                else
                                    nCN1++;
                            }
                        }
                    }
                }
            }
        }
    }

    if(good_read !=0){
        ave_mapqual = nmapqual / good_read;
        ave_insertSize = totalInsertSize / good_read;
    }else{
        ave_insertSize = 0;
        ave_mapqual = 0;
    }
}

void otherinfo::mismatches_2_dedup(){
    // once two read owes two mismatches with same distance and owing the same start position are usually reliable
    map<string,vector<vector<int>>>::iterator itr1 = location2mismatch.begin();
    map<int,vector<int>> position;
    // start<end>

    if(moptions->debug){
        if(location2mismatch.size()){
            cerr <<"location of reads with two mismatches:\n";
            for(;itr1!=location2mismatch.end();itr1++){
                for(int j=0;j<itr1->second.size();j++){
                    cerr << itr1->first <<"\t" <<itr1->second[j][0] << "\t" << itr1->second[j][1] <<endl;
                }
            }
            cerr << "------------------------------------------------------------\n\n";
            itr1 = location2mismatch.begin();
        }else{
            cerr << "there are no reads that covered the position with 2 mismatches\n";
            cerr << "------------------------------------------------------------\n\n";
        }
    }

    for(;itr1!=location2mismatch.end();itr1++){
        mismatchCigar += itr1->first + "(";
        for(int j=0;j<itr1->second.size();j++){
            position[itr1->second[j][0]].push_back(itr1->second[j][1]);
        }
        map<int,vector<int>>::iterator itr2 = position.begin();
        int singleRead(0);
        for(;itr2 != position.end();itr2++){
            if(itr2->second.size()>1){
                mismatchCigar += int2string(itr2->second.size()) + "+";
            }else
                singleRead++;
        }

        if(singleRead == 0)
            mismatchCigar = mismatchCigar.substr(0,mismatchCigar.size()-1);
        else
            mismatchCigar += int2string(singleRead);

        mismatchCigar += ")";
        position.clear();
    }
    stringTransfer();
}


void otherinfo::umi_dedup(){
    map<string,vector<vector<int>>>::iterator itr1 = location2UMI.end();
    if(moptions->debug){
        cerr <<"UMI && position\n";
        for(; itr1!=location2UMI.end();itr1++){
            cerr << itr1->first << "\t";
            for(int i=1;i<itr1->second.size()-1;i++){
                cerr <<itr1->second[i][0] <<"-" << itr1->second[i][1] << "\t";
            }
            cerr << "\n";
        }
        itr1 = location2UMI.begin();
    }

    for(; itr1!=location2UMI.end();itr1++){
        for(int i=1;i<itr1->second.size()-1;i++){
            for(int j=i+1;i<itr1->second.size()-1;j++){
            }
        }
    }
}

void otherinfo::stringTransfer(){
    int arr[13] = {
     good_molecule,
     good_read,
     good_noverlap_pe,
     good_noverlap_se,
     noverlap_pe,
     noverlap_se,
     ave_mapqual,
     ave_insertSize,
     nq10,
     nterminal,
     nmulti,
     nCN1,
     nCN2,
    };

    for(int i=0;i<13;i++){
        extrainfo += int2string(arr[i]) + ",";
    }
    mismatchCigar = (mismatchCigar.size() == 0)?"N":mismatchCigar;
    extrainfo += mismatchCigar;
}
