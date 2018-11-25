#include "helper.h"

instruction::instruction(){
}

instruction::~instruction(){
}

void instruction::explain(){
    cout << "instruction for result\n";
    cout << "The first added column: mor, mnr, msr, oor, onr, osr, moa, mna, msa, ooa, ona, osa\n";
    cout << "The second added column: nq10, nterminal,nmulti,noverlap_pe,noverlap_se,good_noverlap_pe,";
    cout << "good_noverlap_se,good_molecule,good_read";
    cout << "nCN1,nCN2,ave_mapqual,ave_insertSize\n\n";

    cout << "nq10\t" << "the number of reads support mutation\n";
    cout << "nterminal\t" <<  "the number of reads owing mutation located at the 5' or 3' of the read, distance cutoff: 10%\n";
    cout << "nmulti\t" <<  "the number of reads owing mutation: has multiple alignments by BWA mapping quality == 0 or is secondary alignment\n";
    cout <<"noverlap_pe\t" <<  "the number of molecules owing mutation located the overlap area of pair-end reads\n";
    cout <<"noverlap_se\t" <<  "the number of molecules that are single end reads or one of read owing mutation while the other not\n";
    cout << "good_noverlap_pe\t" << "after filtering of noverlap_pe\n";
    cout << "good_noverlap_se\t" << "aftering filtering of noverlap_se\n";
    cout <<"good_molecule\t" <<  "the number of molecules without XA, mapping quality more than 0 and mutation site located at the middle of read\n";
    cout <<"good_read\t" <<  "the number of reads from good_molecule\n";
    cout <<"nCN1\t" <<  "the number of molecules that got copy number equal to 1\n";
    cout <<"nCN2\t" <<  "the number of molecules that got copy number more than 1\n";

    cout <<"ave_mapqual\t" <<  "the average map quality of reads support mutation\n";
    cout <<"ave_insertSize\t" <<  "the average insert size of reads support mutation\n";
    cout << "last:\t" << "location of reads support two mismatches\n";
    exit(-1);
}
