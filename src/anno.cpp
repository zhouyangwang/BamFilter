#include "anno.h"

Anno::Anno(Options *opt){
    mOptions = opt;
    vcfInfo = new vcfreader(mOptions);
    vcfInfo->fileExtract();
    //vcfOUT.open(mOptions->output);
    lastVCFChr = "chr1";
    lastVCFPos = 0;
    if(mOptions->debug){
        cerr << "vcf file contains " << vcfInfo->vcfPosChr.size() << " chromosomes, " <<  vcfInfo->vcfLineSplit.size() <<" mutations\n";
    }
}

Anno::~Anno(){
    delete vcfInfo;
}

void Anno::annotation(){

    samFile *in;
    in = sam_open(mOptions->BamCfdna.c_str(),"r");
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
    int lastTid = -1;
    int lastPos = -1;
    int linemax = 1000;
    int linecurrent = 0;
    long *reference_start = new long;
    long *reference_end = new long;

    while((r = sam_read1(in, mBamHeader,b)) >= 0){

        if(linecurrent > linemax){
            //break;
            linecurrent ++;
        }else{
            linecurrent++;
        }
        //check whether the BAM is sorted
        if(b->core.tid <lastTid || (b->core.tid == lastTid && b->core.pos <lastPos)) {

            // skip the -1:-1, which means unmapped
            if(b->core.tid >=0 && b->core.pos >= 0) {
                cerr << "ERROR: the input is unsorted. Found unsorted read in " << b->core.tid << ":" << b->core.pos << endl;
                cerr << "Please sort the input first." << endl << endl;
                exit(-1);
            }
        }
        *reference_end = bam_endpos(b)-1;
        *reference_start = lastPos;
        readcollect(b, reference_start, reference_end);
        lastTid = b->core.tid;
        lastPos = b->core.pos;
        bam_destroy1(b);
        b = NULL;
        b = bam_init1();
     }

    // The last read in last chromosome aligned with mutation site, leaving mutation info undealed
    if(vcfInfo->vcfPosChr.count(lastVCFChr)){
        vector<long*>::iterator itr1 = vcfInfo->vcfPosChr[lastVCFChr].begin();
        for(;itr1 != vcfInfo->vcfPosChr[lastVCFChr].end();){
            vcfinforoutput(lastVCFChr,**itr1);
            delete vcfInfo->vcfPosChr[lastVCFChr][0];
            vcfInfo->vcfPosChr[lastVCFChr][0] = NULL;
            vcfInfo->vcfPosChr[lastVCFChr].erase(itr1);
        }
    }

    delete reference_start;
    delete reference_end;
    reference_start = NULL;
    reference_end = NULL;
    if(mOptions->debug){
        cerr <<"bam file contains " << linecurrent << " reads\n";
    }
}

void Anno::readcollect(bam1_t* b, long * reference_start, long * reference_end){
    string chromosome = BamUtil::getChr(b);

    if(chromosome != lastVCFChr && vcfInfo->vcfPosChr.count(lastVCFChr)){
//  the last read in one chromosome rightly mapped with mutation site,therefore have to deal when switch to next chromosome;
        vector<long*>::iterator itr1 = vcfInfo->vcfPosChr[lastVCFChr].begin();
        for(;itr1 != vcfInfo->vcfPosChr[lastVCFChr].end();){
            vcfinforoutput(lastVCFChr,**itr1);
            delete vcfInfo->vcfPosChr[lastVCFChr][0];
            vcfInfo->vcfPosChr[lastVCFChr][0] = NULL;
            vcfInfo->vcfPosChr[lastVCFChr].erase(itr1);
        }
        lastVCFChr = chromosome;
    }
    lastVCFChr = chromosome;

    if ( vcfInfo->vcfPosChr.count(chromosome) && ! vcfInfo->vcfPosChr[chromosome].empty()){
        long pos = * vcfInfo->vcfPosChr[chromosome][0];
        string ID = chromosome + val2string(pos);

        if( *reference_start <= pos && pos <= *reference_end){
            //read aligned with first vcf point -----|----
            BamUtil *readinfo = new BamUtil(b, pos, mOptions);
//            readinfo->AddingMoreInforforRead(b,'a');
            readcache[ID].push_back(readinfo);

            for(int i=1;i < vcfInfo->vcfPosChr[chromosome].size();i++){
                pos = *vcfInfo->vcfPosChr[chromosome][i];
                ID = chromosome + val2string(pos);
                if (pos <= *reference_end){
                    // ------------|-------|------
                    BamUtil * readinfoCopy = new BamUtil(readinfo,pos);
                    //readinfoCopy->updateReadInfo(b,pos);
                    readcache[ID].push_back(readinfoCopy);
                }else{
                    break;
                }
            }

        }else if(pos < *reference_start){
            // | -------------
            vcfinforoutput(chromosome,pos);
            vector<long*>::iterator itr1 = vcfInfo->vcfPosChr[chromosome].begin();
            delete vcfInfo->vcfPosChr[chromosome][0];
            vcfInfo->vcfPosChr[chromosome][0] = NULL;
            vcfInfo->vcfPosChr[chromosome].erase(itr1);

            for(;itr1 != vcfInfo->vcfPosChr[chromosome].end();){
                pos = **itr1;
                ID = chromosome + val2string(pos);

                if ( *reference_start <= pos && pos <= *reference_end){
                    // ------------|-------|------
                    BamUtil *readinfo = new BamUtil(b, pos, mOptions);
                    readcache[ID].push_back(readinfo);

                    for(int i=1;i < vcfInfo->vcfPosChr[chromosome].size();i++){
                        pos = *vcfInfo->vcfPosChr[chromosome][i];
                        ID = chromosome + val2string(pos);
                        if (pos <= *reference_end){
                            // ------------|-------|------
                            BamUtil * readinfoCopy = new BamUtil(readinfo,pos);
                            //readinfoCopy->updateReadInfo(b,pos);
                            readcache[ID].push_back(readinfoCopy);
                        }else{
                            break;
                        }
                    }
                    break;

                }else if (pos < *reference_start) {
                    vcfinforoutput(chromosome,pos);
                    delete vcfInfo->vcfPosChr[chromosome][0];
                    vcfInfo->vcfPosChr[chromosome][0] = NULL;
                    vcfInfo->vcfPosChr[chromosome].erase(itr1);

                }else{
                    break;
                }
            }
        }
    }else if(vcfInfo->vcfPosChr.find(chromosome) != vcfInfo->vcfPosChr.end()){
        vcfInfo->vcfPosChr.erase(chromosome);
        if(mOptions->debug){
            cerr << "all mutaions in " << chromosome << " were dealed\n";
        }
    }
}

void Anno::vcfinforoutput(string chromosome, long pos){
    string ID = chromosome + val2string(pos);
    long posSt = ++pos;
    cout << chromosome << "\t" << posSt << "\t" << readcache[ID].size()<< endl;
    readsFilter();
    while(readcache[ID].size()!=0){
        delete readcache[ID][0];
        readcache[ID][0] = NULL;
        readcache[ID].erase(readcache[ID].begin());
    }
}

void Anno::readsFilter(){

}
