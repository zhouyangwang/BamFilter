#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "cmdline.h"
#include "common.h"
#include <sstream>
#include "options.h"
#include "anno.h"

using namespace std;

string command;

int main(int argc, char* argv[]){

    cmdline::parser cmd;
    // input/output
    cmd.add<string>("query", 0, "vcf file contains mutations to query", false);
    cmd.add<string>("cfdna",'c',"bam file contains cfdna reads info. There must be a corresponding .bai file in the same directory",false);
    cmd.add<string>("gdna",'g', "bam file contains gdna reads info. There must be a corresponding .bai file in the same directory",false);
    cmd.add<string>("output",'o',"output vcf file. Will be overwritten if already exists",false);

    cmd.add<int>("skip", 0,"skip the first N lines",false, 0);
    cmd.add<int>("qual", 'q', "drop bases whose quality is less than this", false, 25);
    cmd.add<int>("mismatch-limit", 'm', "if set, drop reads that has more mismatches than the limit. requires a 'MD' or a 'NM' tag to be present.", false,-1);

    cmd.add("simple",'s',"annotate less informations into vcf output");
    cmd.add("fast",'f',"do not infer origin read size by CIGAR, it can be faster and consume less memory.");
    //cmd.add("drop-inconsist", 0,"drop different reads stack at the same position. This decreases sensitivity.");
    cmd.add("dropXA",0, "drop reads that has XA tag (multiple alignment)");
    cmd.add("UMI", 'u', "count umi sequences when sample sequenced by duplex UMI");
    cmd.add("indel",0, "only indel exists in vcf file");
    cmd.add("snp", 0, "only snp exists in vcf file");

    cmd.add("debug", 0, "output some debug information to STDERR.");
    cmd.add("instruction", 'h', "detail explanation for result");

    cmd.parse_check(argc, argv);

    Options opt;
    opt.QueryVcf = cmd.get<string>("query");
    opt.BamCfdna = cmd.get<string>("cfdna");
    opt.BamGdna  = cmd.get<string>("gdna");
    opt.output = cmd.get<string>("output");
    opt.LineSkip = cmd.get<int>("skip");
    opt.MinQualDrop = cmd.get<int>("qual");
    opt.ReadMaxMismatch = cmd.get<int>("mismatch-limit");
    opt.FastInferReadSize = cmd.exist("fast");
    opt.OutputSimpleOverlap = cmd.exist("simple");
    //opt.DropInconsist = cmd.exist("drop-inconsist");
    opt.ReadDropXA = cmd.exist("dropXA");
    opt.DuplexUMIStat = cmd.exist("UMI");
    opt.indel = cmd.exist("indel");
    opt.snp = cmd.exist("snp");
    opt.debug = cmd.exist("debug");
    opt.explain = cmd.exist("instruction");

    if(! opt.validate()){
        return(0);
    }
    time_t t1 = time(NULL);

    stringstream ss;
    for(int i=0;i<argc;i++){
        ss << argv[i] << " ";
    }
    command = ss.str();

    Anno anno(&opt);
    anno.annotation();

    time_t t2 = time(NULL);
    cerr << endl << command << endl;
    cerr << "MrBam v:" << VERSION_NUMBER << ", time used: " << (t2)-t1 << " seconds" << endl;

    return 0;
}
