#include "anno.h"

Anno::Anno(Options *opt){
    mOptions = opt;
    vcfInfo = new vcfreader(mOptions);
    vcfInfo->fileExtract();
    vcfOUT.open(mOptions->output);
    if(mOptions->debug){
        cerr << "vcf file contains " << vcfInfo->vcfPosChr.size() << " chromosomes, " <<  vcfInfo->vcfLineSplit.size() <<" mutations\n\n";
    }
}

Anno::~Anno(){
    delete vcfInfo;
}

void Anno::annotation(){
    string bamfile[2] ={mOptions->BamCfdna, mOptions->BamGdna};
    long *reference_start = new long;
    long *reference_end = new long;

    for(int i=0;i<1;i++){
        lastVCFChr = "chr1";
        lastVCFPos = 0;
        samFile *in;
        cfdna = (i==0)?true:false;
        in = sam_open(bamfile[i].c_str(),"r");
        if(!in){
            cerr << "ERROR: failed to open " << mOptions->BamCfdna << endl;
            exit(-1);
        }
        mBamHeader = sam_hdr_read(in);
        mOptions->bamHeader = mBamHeader;
        if (mBamHeader == NULL || mBamHeader->n_targets == 0) {
            cerr << "ERROR: this SAM file has no header " << mOptions->BamCfdna << endl;
            exit(-1);
        }
        //BamUtil::dumpHeader(mBamHeader);

        bam1_t *b = NULL;
        b = bam_init1();
        int r;
        int lastid = -1;
        int lastPos = -1;
        int linemax = 1000;
        int linecurrent = 0;

        while((r = sam_read1(in, mBamHeader,b)) >= 0){
            //check whether the BAM is sorted
            if(b->core.tid <lastid || (b->core.tid == lastid && b->core.pos <lastPos)) {

                // skip the -1:-1, which means unmapped
                if(b->core.tid >=0 && b->core.pos >= 0) {
                    cerr << "ERROR: the input is unsorted. Found unsorted read in " << b->core.tid << ":" << b->core.pos << endl;
                    cerr << "Please sort the input first." << endl << endl;
                    exit(-1);
                }
            }
            *reference_end = bam_endpos(b)-1;
            *reference_start = b->core.pos;
            readcollect(b, reference_start, reference_end);
            lastid = b->core.tid;
            lastPos = b->core.pos;
            bam_destroy1(b);
            b = NULL;
            b = bam_init1();
         }

        // The last read in last chromosome aligned with mutation site
        if(vcfInfo->vcfPosChr[lastVCFChr].size() ){
            long pos = * vcfInfo->vcfPosChr[lastVCFChr][0];
            vcfinforoutput(lastVCFChr,vcfInfo->vcfPosChr[lastVCFChr].size()-1);
        }
    }
    delete reference_start;
    delete reference_end;
    reference_start = NULL;
    reference_end = NULL;
}

void Anno::readcollect(bam1_t* b, long * reference_start, long * reference_end){
    string chromosome = BamUtil::getChr(b);

    if(chromosome != lastVCFChr){
        if(vcfInfo->vcfPosChr[lastVCFChr].size()){
        // the last read in one chromosome rightly mapped with mutation site,therefore have to deal when switch to next chromosome;
            vcfinforoutput(lastVCFChr,vcfInfo->vcfPosChr[lastVCFChr].size()-1);
        }
        readcache.clear();
        lastVCFChr = chromosome;
    }

    if (vcfInfo->vcfPosChr[chromosome].size()){
        long pos = * vcfInfo->vcfPosChr[chromosome][0];
        string idInAlt = chromosome + long2string(pos);

        if(pos < *reference_start){
            // | -------------
            if(vcfInfo->vcfPosChr[chromosome].size() == 1){
                vcfinforoutput(chromosome,0);
                return;
            }
            BamUtil *readinfo = new BamUtil(b, pos, mOptions);

            for(int i=1;i<vcfInfo->vcfPosChr[chromosome].size();i++){
                if( *vcfInfo->vcfPosChr[chromosome][i] >= *reference_start){

                    for(int j=i;j<vcfInfo->vcfPosChr[chromosome].size();j++){
                        long pos2 = *vcfInfo->vcfPosChr[chromosome][j];
                        idInAlt = chromosome + long2string(pos2);
                        if(*reference_start <= pos2 && pos2 <= *reference_end){
                            BamUtil * readinfoCopy = new BamUtil(readinfo,pos2);
                            readcache[idInAlt].push_back(readinfoCopy);
                        }else
                            break;
                    }
                    if(*vcfInfo->vcfPosChr[chromosome][i] - *vcfInfo->vcfPosChr[chromosome][i-1] != 1)
                        vcfinforoutput(chromosome,i-1);
                        // |+--+--- continuous mutations might be combined
                    break;
                }
            }
            delete readinfo;
            readinfo= NULL;

        }else if( pos <= *reference_end){
            //read aligned with first vcf point -----|----
            BamUtil *readinfo = new BamUtil(b, pos, mOptions);
            readcache[idInAlt].push_back(readinfo);

            for(int i=1;i < vcfInfo->vcfPosChr[chromosome].size();i++){
                pos = *vcfInfo->vcfPosChr[chromosome][i];
                idInAlt = chromosome + long2string(pos);
                if (pos <= *reference_end){
                    // ------------|-------|------
                    BamUtil * readinfoCopy = new BamUtil(readinfo,pos);
                    readcache[idInAlt].push_back(readinfoCopy);
                }else
                    break;
            }
        }
    }
}

void Anno::vcfinforoutput(string chromosome, int ncontinuous){

    for(int i=0;i<=ncontinuous;i++){
        long pos = *vcfInfo->vcfPosChr[chromosome][i];

        idInAlt = chromosome + long2string(pos);
        idInLine = idInAlt + vcfInfo->vcfAlt[idInAlt][0];

        if(mOptions->debug)
            cerr <<chromosome << "\t" << pos<< "\tdepth: " <<readcache[idInAlt].size() << "\tneighboring: " << ncontinuous <<endl;

        if(cfdna){
            if(mOptions->snp){
                if(i!= ncontinuous && pos - *vcfInfo->vcfPosChr[chromosome][i+1] == -1){

                    continuous = true;
                    string continuousID = "continuous";
                    nextidInAlt = chromosome + long2string(pos+1);
                    nextidInLine = nextidInAlt + vcfInfo->vcfAlt[nextidInAlt][0];
                    readCombine(chromosome,idInAlt,continuousID,pos);

                    if (readsFilter(continuousID)){
                        // continuous bases could be combined
                        if(mOptions->debug)
                            cerr << "More than 3 piece of reads supporting continuous mutations:" << pos <<"-"<< pos+1 << endl;

                        vcfInfo->vcfAlt[idInAlt].clear();
                        int j(0);
                        while(j<mutation.size()){
                            vcfInfo->vcfAlt[idInAlt].push_back(mutation[j]);
                            newidInLine = chromosome + long2string(pos) + mutation[j];

                            if(mutation[j][0] != refBase[0] && mutation[j][1] != refBase[1]){
                                // eg: ref: AT, all position mutations: GC AC GT, if GC exist first, dealing GC first, other requires recalculation
                                vcfInfo->vcfLineSplit[newidInLine] = vcfInfo->vcfLineSplit[idInLine];
                                vcfInfo->vcfLineSplit[newidInLine][2] = long2string(str2long(vcfInfo->vcfLineSplit[newidInLine][2])+1);
                                vcfInfo->vcfLineSplit[newidInLine][3] = refBase;
                                vcfInfo->vcfLineSplit[newidInLine][4] = mutation[j];

                            }else if(vcfInfo->vcfAlt[idInAlt][j][0] != refBase[0] && vcfInfo->vcfAlt[idInAlt][j][1] == refBase[1]){
                                // AT->GT AT->CT AT->TT
                                vcfInfo->vcfLineSplit[newidInLine] = vcfInfo->vcfLineSplit[idInLine];
                                vcfInfo->vcfLineSplit[newidInLine][4] = vcfInfo->vcfAlt[idInAlt][j][0];

                            }else if(vcfInfo->vcfAlt[idInAlt][j][0] == refBase[0] && vcfInfo->vcfAlt[idInAlt][j][1] != refBase[1]){
                                // AT->AC AT->AG AT->AA
                                vcfInfo->vcfLineSplit[newidInLine] = vcfInfo->vcfLineSplit[nextidInLine];
                                vcfInfo->vcfLineSplit[newidInLine][4] = vcfInfo->vcfAlt[idInAlt][j][1];
                            }
                            finaloutput(newidInLine,j);
                            j++;
                        }
                        readTempVariableFree();
                        readsFree(idInAlt,nextidInAlt);
                        i++;
                        // skip to next next position;

                    }else{
                        // continuous mutation couldn't be combined;
                        if(mOptions->debug)
                            cerr << "failing to combine continuous mutations, position: " << pos << "-" << pos+1 << endl;
                        singleSNPMutation(idInLine,idInAlt);
                    }
                }else{
                    // there are no continuous mutations
                    singleSNPMutation(idInLine,idInAlt);
                }
            }else
                singleINDELMutation(idInLine,idInAlt);

        }else{ // gdna
            if(mOptions->snp){
                if( vcfInfo->vcfAlt[idInAlt][0].size() == 2  && i != ncontinuous && pos - *vcfInfo->vcfPosChr[chromosome][i+1] == -1 ){

                    continuous = true;
                    string continuousID = "continuous";
                    nextidInAlt = chromosome + long2string(pos+1);
                    nextidInLine = nextidInAlt + vcfInfo->vcfAlt[nextidInAlt][0];
                    readCombine(chromosome,idInAlt,continuousID,pos);
                    readsFilter(continuousID);

                    for(int j=0;j<vcfInfo->vcfAlt[idInAlt].size();j++){
                        newidInLine = chromosome + long2string(pos) + vcfInfo->vcfAlt[idInAlt][j];
                        finaloutput(newidInLine,j);
                    }
                    readTempVariableFree();
                    readsFree(idInAlt,nextidInAlt);
                    i++;
                    // skip to next next position;
                }else
                    singleSNPMutation(idInLine,idInAlt);
            }else
                singleINDELMutation(idInLine,idInAlt);
        }
    }

    for(int i=0;i<=ncontinuous;i++){
        delete vcfInfo->vcfPosChr[chromosome][0];
        vcfInfo->vcfPosChr[chromosome][0] = NULL;
        vcfInfo->vcfPosChr[chromosome].erase(vcfInfo->vcfPosChr[chromosome].begin());
    }
}


bool Anno::readsFilter(string idInAlt){
    int i(0);long start(0);
    while(i< readcache[idInAlt].size()){
        name_dict[readcache[idInAlt][i]->query_name].push_back(readcache[idInAlt][i]);
        // split into PE and SE mode
        i++;
    }

    map<string,vector<BamUtil*>>::iterator itr1 = name_dict.begin();
    for(;itr1!=name_dict.end();){
        if(itr1->second.size() == 1){
        // non-overlap or single
            if (0 <= itr1->second[0]->mismatchBaseQuality && itr1->second[0]->mismatchBaseQuality <=  mOptions->MinQualDrop){
                if(mOptions->debug)
                    cerr << "low quality " << itr1->first << endl;
                itr1 = name_dict.erase(itr1);
                continue;
            }
            if(mOptions->ReadMaxMismatch != -1 && itr1->second[0]->nmismatchbases > mOptions->ReadMaxMismatch){
                if(mOptions->debug){
                    cerr << itr1->first << " has mismatch bases more than " << mOptions->ReadMaxMismatch;
                    cerr << ",it's abandoned, mismatched bases: " << itr1->second[0]->nmismatchbases <<endl;
                }
                itr1 = name_dict.erase(itr1);
                continue;
            }
            if(mOptions->ReadDropXA && itr1->second[0]->has_tag_XA){
                if(mOptions->debug)
                    cerr << "multiple alignment with mapping quality equal to 0: " << itr1->first << endl;
                itr1 = name_dict.erase(itr1);
                continue;
            }

            start = min(itr1->second[0]->reference_start - itr1->second[0]->query_alignment_start, itr1->second[0]->next_reference_start);
            posbase = itr1->second[0]->mismatchbase;
            if(itr1->second[0]->is_pair_mate_mapped){
                //if(mOptions->FastInferReadSize){
                //    start = min(itr1->second[0]->reference_start, itr1->second[0]->next_reference_start);
                //}
                posid.push_back(start); posid.push_back(itr1->second[0]->template_length); posid.push_back(0);
                unique_pairs[posid].push_back(posbase);
            }else{
                posid.push_back(start); posid.push_back(itr1->second[0]->query_length); posid.push_back(1);
                unique_single[posid].push_back(posbase);
            }
            posid.clear();

            if(cfdna && mOptions->snp && posbase != refBase)
                snv[posbase].push_back('1');

        }else if(itr1->second.size() == 2){

            if(0 <= itr1->second[0]->mismatchBaseQuality && itr1->second[0]->mismatchBaseQuality <=  mOptions->MinQualDrop &&
               0 <= itr1->second[1]->mismatchBaseQuality && itr1->second[1]->mismatchBaseQuality <=  mOptions->MinQualDrop){
                if(mOptions->debug)
                    cerr << "low quality " << itr1->first << endl;
                itr1 = name_dict.erase(itr1);
                continue;
               }

            if(itr1->second[0]->mismatchbase != itr1->second[1]->mismatchbase){
                if(mOptions->debug)
                    cerr << "pair inconsistent: " << itr1->first << endl;
                itr1 = name_dict.erase(itr1);
                continue;
            }

            if(mOptions->ReadMaxMismatch != -1 && (itr1->second[0]->nmismatchbases > mOptions->ReadMaxMismatch ||
                                                   itr1->second[1]->nmismatchbases > mOptions->ReadMaxMismatch)){
                if(mOptions->debug){
                    cerr << itr1->first << " has mismatch bases more than " << mOptions->ReadMaxMismatch;
                    cerr << ",it's abandoned, mismatched bases: " << itr1->second[0]->nmismatchbases <<endl;
                }
                itr1 = name_dict.erase(itr1);
                continue;
            }

            if(mOptions->ReadDropXA && (itr1->second[0]->has_tag_XA || itr1->second[1]->has_tag_XA)){
                if(mOptions->debug)
                    cerr << "multiple alignment: " << itr1->first << endl;
                itr1 = name_dict.erase(itr1);
                continue;
            }

            start = min(itr1->second[0]->reference_start - itr1->second[0]->query_alignment_start,
                        itr1->second[1]->reference_start - itr1->second[1]->query_alignment_start);
            long tlen = max(itr1->second[0]->reference_start - itr1->second[0]->query_alignment_start + itr1->second[0]->query_length ,
                            itr1->second[1]->reference_start - itr1->second[1]->query_alignment_start + itr1->second[1]->query_length) - start;

            posbase = itr1->second[0]->mismatchbase;
            posid.push_back(start); posid.push_back(tlen); posid.push_back(1);
            unique_pairs[posid].push_back(posbase);
            posid.clear();

            if(cfdna && mOptions->snp && posbase != refBase){
                snv[posbase].push_back('1');
                snv[posbase].push_back('1');
            }

        }else if(mOptions->debug){
            cerr << itr1->first << ": more than 2 reads share the same name; all droped\n";
            cerr << "total size: " << itr1->second.size() << endl;
        }
        itr1++;
    }

    for(int i=0;i<vcfInfo->vcfAlt[idInAlt].size();i++){
        mutation.push_back(vcfInfo->vcfAlt[idInAlt][i]);
    }

    // check if there is MNV
    // MNV: mutation with supporting reads >= 3 will be regarded as novel snv
    if(cfdna){
        if(mOptions->snp && !continuous){
            if(mOptions->debug){
                cerr << "mutation from vcf file: ";
                for(int i=0;i<vcfInfo->vcfAlt[idInAlt].size();i++)
                   cerr << "\"" << vcfInfo->vcfAlt[idInAlt][i] << "\"" << "\t";
                cerr << endl << "mutations after filtering: ";
            }

            map<string,vector<char>>::iterator itr1 = snv.begin();
            for(;itr1 != snv.end();itr1++){
                if(itr1->first != vcfInfo->vcfAlt[idInAlt][0]){
                    if(vcfInfo->vcfAlt[idInAlt].size()==2 && itr1->first == vcfInfo->vcfAlt[idInAlt][1])
                    // one position owing MNV
                        continue;
                        if(itr1->second.size()> 2 && (itr1->first == "A" || itr1->first == "G" || itr1->first == "C" || itr1->first == "T")){
                            mutation.push_back(itr1->first);
                    }
                }
            }
            if(mOptions->debug){
                int j(0);
                while(j <mutation.size()){
                    cerr << "\"" << mutation[j] << "\"" << "\tsupporting reads: ";
                    if(snv.count(mutation[j]))
                        cerr << snv[mutation[j]].size() << "\t";
                    else
                        cerr << "0\t";
                    j++;
                }
                cerr << endl;
            }
            snv.clear();
            return false;

        }else if(mOptions->snp && continuous){
            mutation.clear();
            map<string,vector<char>>::iterator itr2 = snv.begin();
            for(;itr2!=snv.end();itr2++){
                if(itr2->second.size() > 2){
                    mutation.push_back(itr2->first);
                }
            }
            snv.clear();

            for(int i=0;i<mutation.size();i++){
                if(mutation[i][0] != refBase[0] && mutation[i][1] != refBase[1]){
                    if(mOptions->debug){
                        cerr << "trying to combine continuous mutations, mutations after filtering: ";
                        for(int i=0;i<mutation.size();i++)
                            cerr << "\"" << mutation[i] << "\"" << "\tsupporting reads: " << snv[mutation[i]].size() << "\t starting identifying good molecule:";
                        cerr << endl << "comtinuous mutations could be combined\n";
                    }
                    return true;
                }
            }

            if(mOptions->debug){
                cerr << "trying to combine continuous mutations, mutations after filtering: ";
                for(int i=0;i<mutation.size();i++)
                    cerr << "\"" << mutation[i] << "\"" << "\tsupporting reads: " << snv[mutation[i]].size() << "\t starting identifying good molecule:";
                cerr << endl << "failing to combine continuous mutations\n";
            }
            return false;
        }
        // indel
        return false;

    }else{
        return false;
    }
}

void Anno::readCombine(string chromosome, string idInAlt, string continuousID,long pos){
    int i(0);
    string idInAltnext = chromosome + long2string(pos+1);
    refBase = vcfInfo->vcfLineSplit[chromosome+long2string(pos)+vcfInfo->vcfAlt[idInAlt][0]][3] +
                    vcfInfo->vcfLineSplit[idInAltnext + vcfInfo->vcfAlt[idInAltnext][0]][3];

    while(i < readcache[idInAlt].size()){
        if(readcache[idInAlt][i]->reference_start <= pos && pos+1 <= readcache[idInAlt][i]->reference_end){
        // extract reads that cover both mutation site
            readcache[idInAlt][i]->getquerypos(refBase,true);
            readcache[continuousID].push_back(readcache[idInAlt][i]);
        }
        i++;
    }
}

void Anno::finaloutput(string idInLine, int j){

    Counter overlap(mOptions, refBase, vcfInfo->vcfAlt[idInAlt][j]);
    overlap.unique_pair_count(&unique_pairs);
    overlap.unique_single_count(&unique_single);
    otherinfo allinfo(mOptions, refBase, vcfInfo->vcfAlt[idInAlt][j]);
    allinfo.good_molecule_stat(&name_dict);
    allinfo.mismatches_2_dedup();

    if(cfdna){
        vcfInfo->vcfLineSplit[idInLine][60] += ":UDP";
        vcfInfo->vcfLineSplit[idInLine][61] += ":gdna";
        vcfInfo->vcfLineSplit[idInLine][62] += ":" + overlap.overlap;
        vcfInfo->vcfLineSplit[idInLine].push_back("gdna");
        vcfInfo->vcfLineSplit[idInLine].push_back(allinfo.extrainfo);
    }else{
        vcfInfo->vcfLineSplit[idInLine][63] = overlap.overlap;
        vcfInfo->vcfLineSplit[idInLine][65] = allinfo.extrainfo;
        int k(0);
        while(k<vcfInfo->vcfLineSplit[idInLine].size()-1){
            vcfOUT << vcfInfo->vcfLineSplit[idInLine][k] << "\t";
            k++;
        }
        vcfOUT <<vcfInfo->vcfLineSplit[idInLine][j] << endl;
        vcfInfo->vcfLineSplit[idInLine].resize(0);
        vcfInfo->vcfLineSplit[idInLine].shrink_to_fit();
    }

    // test output for cfdna only
    if(mOptions->debug){
        int k(0);
        while(k<vcfInfo->vcfLineSplit[idInLine].size()-1){
            vcfOUT << vcfInfo->vcfLineSplit[idInLine][k] << "\t";
            k++;
        }
        vcfOUT << vcfInfo->vcfLineSplit[idInLine][k] << endl;
        vcfInfo->vcfLineSplit[idInLine].resize(0);
        vcfInfo->vcfLineSplit[idInLine].shrink_to_fit();
    }
}

void Anno::readTempVariableFree(){
    unique_pairs.clear();
    unique_single.clear();
    name_dict.clear();
    vector<string>().swap(mutation);
    string continuousId = "continuous";
    if(readcache[continuousId].size() != 0)
        vector<BamUtil*>().swap(readcache[continuousId]);
}

void Anno::readsFree(string id, string nextid){

    int i(0);
    if(nextid.size()){
        while(i < readcache[nextid].size()){
            delete readcache[nextid][i];
            readcache[nextid][i] = NULL;
            i++;
        }
        readcache[nextid].resize(0);
        readcache[nextid].shrink_to_fit();
    }

    i =0;
    while(i < readcache[id].size()){
            delete readcache[id][i];
            readcache[id][i] = NULL;
            i++;
    }
    readcache[id].resize(0);
    readcache[id].shrink_to_fit();
}

void Anno::singleSNPMutation(string idInLine,string idInAlt, bool altFromContinuous){
    continuous = false;
    refBase = vcfInfo->vcfLineSplit[idInLine][3];

    // in case of failing to combine continuous bases
    readTempVariableFree();

    int k(0);
    while(k < readcache[idInAlt].size()){
        readcache[idInAlt][k]->getquerypos(refBase);
        k++;
    }
    readsFilter(idInAlt);

    //  cfdna owing novel mutation while gnda not will be output
    //  gdna owing novel muatation while cfdna not won't be output

    if(cfdna){
        k = 0;
        vcfInfo->vcfAlt[idInAlt].clear();
        while( k<mutation.size()){
            // in case of information deleted
            newidInLine = idInAlt + mutation[k];
            vcfInfo->vcfAlt[idInAlt].push_back(mutation[k]);

            if( ! vcfInfo->vcfLineSplit.count(newidInLine)){
                vcfInfo->vcfLineSplit[newidInLine] = vcfInfo->vcfLineSplit[idInLine];
                vcfInfo->vcfLineSplit[newidInLine][4] = mutation[k];
            }
            k++;
        }
    }

    k = 0;
    while(k<vcfInfo->vcfAlt[idInAlt].size()){
        newidInLine = idInAlt + vcfInfo->vcfAlt[idInAlt][k];
        finaloutput(newidInLine,k);
        k++;
    }

    readTempVariableFree();
    readsFree(idInAlt);
}

void Anno::singleINDELMutation(string idInLine, string idInAlt){
    refBase = vcfInfo->vcfLineSplit[idInLine][3];
    int k=0;
    while(k < readcache[idInAlt].size()){
        readcache[idInAlt][k]->getquerypos(refBase);
        k++;
    }
    readsFilter(idInAlt);
    int j(0);
    while(j<mutation.size()){
        string newidInLine = idInAlt + mutation[j];
        finaloutput(newidInLine,j);
        j++;
    }
    readTempVariableFree();
    readsFree(idInAlt);
}

