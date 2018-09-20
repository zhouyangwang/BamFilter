#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "cmdline.h"
#include "common.h"
#include <sstream>
#include "options.h"

using namespace std;

string command;

int main(int argc, char* argv[]){

    cmdline::parser cmd;
    // input/output
    cmd.add<string>("query", 0, "vcf file contains mutations to query", true);
    cmd.add<string>("cfdna",'c',"bam file contains cfdna reads info. There must be a corresponding .bai file in the same directory",true);
    cmd.add<string>("gdna",'g', "bam file contains gdna reads info. There must be a corresponding .bai file in the same directory",true);
    cmd.add<string>("output",'o',"output vcf file. Will be overwritten if already exists",false);
    cmd.add<string>("repeat", 'r', "file contains repeat region in huam genome", false);

    cmd.add<int>("skip", 0,"skip the first N lines",false, 0);
    cmd.add<int>("qual", 'q', "drop bases whose quality is less than this", false, 25);
    cmd.add<int>("mismatch-limit", 'm', "if set, drop reads that has more mismatches than the limit. requires a 'MD' or a 'NM' tag to be present.", false,-1);

    cmd.add<bool>("simple",'s',"annotate less informations into vcf output",false, false);
    cmd.add<bool>("fast",'f',"do not infer origin read size by CIGAR, it can be faster and consume less memory.",false, false);
    cmd.add<bool>("drop-inconsist", 0,"drop different reads stack at the same position. This decreases sensitivity.",false, false);
    cmd.add<bool>("dropXA",0, "drop reads that has XA tag (multiple alignment)",false, false);
    cmd.add<bool>("UMI", 'u', "count umi sequences when sample sequenced by duplex UMI", false, false);
    cmd.add<bool>("indel",0, "only indel exists in vcf file",false, false);
    cmd.add<bool>("snp", 0, "only snp exists in vcf file", false, false);

    cmd.add("debug", 0, "output some debug information to STDERR.");

    cmd.parse_check(argc, argv);

    Options opt;
    opt.QueryVcf = cmd.get<string>("query");
    opt.BamCfdna = cmd.get<string>("cfdna");
    opt.BamGdna  = cmd.get<string>("gdna");
    opt.output = cmd.get<string>("output");
    if (opt.output.empty()){
        opt.output = opt.QueryVcf + "MrBam" + ".txt";

    }
    opt.HumanRepeatFile = cmd.get<string>("repeat");

    opt.LineSkip = cmd.get<int>("skip");
    opt.MinQualDrop = cmd.get<int>("qual");
    opt.ReadMaxMismatch = cmd.get<int>("mismatch-limit");

    opt.FastInferReadSize = cmd.get<bool>("fast");
    opt.OutputSimpleOverlap = cmd.get<bool>("simple");
    opt.DropInconsist = cmd.get<bool>("drop-inconsist");
    opt.ReadDropXA = cmd.get<bool>("dropXA");
    opt.DuplexUMIStat = cmd.get<bool>("UMI");
    opt.indel = cmd.get<bool>("indel");
    opt.snp = cmd.get<bool>("snp");
    opt.debug = cmd.exist("debug");

    opt.validate();
    time_t t1 = time(NULL);

    stringstream ss;
    for(int i=0;i<argc;i++){
        ss << argv[i] << " ";
    }
    command = ss.str();

    time_t t2 = time(NULL);
    cerr << endl << command << endl;
    cerr << "MrBam " << VERSION_NUMBER << ", time used: " << (t2)-t1 << " seconds" << endl;

    return 0;
}
