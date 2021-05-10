#ifndef __LREZ_BARCODES_LOADING__
#define __LREZ_BARCODES_LOADING__

#include "utils.h"

using namespace std;

const string barcodes_Haplotagging_A[] = {
	"ACGGAA", "CCAACA", "AGATCG", "TTCTCC", "TTCCTG", "TTCGGT", "TTGTGG", "TTGCCT", "TTGGTC", "TTACGC", 
	"TTAGCG", "TCTTCG", "TCTCTC", "TCTGGA", "TCCACT", "TCGTAC", "TCGATG", "TCACAG", "TGTTGC", "TGTCCA", 
	"TGTGTG", "TGCTAG", "TGCATC", "TGGAGT", "TGAGAC", "TATCGG", "TATGCC", "TACCAC", "TAGGAG", "CTTCGT", 
	"CTTGCA", "CTCTGA", "CTCAAC", "CTGCTA", "CTGGAT", "CTAAGG", "CCTCAA", "CCTGTT", "CCATTC", "CGTTCT", 
	"CGTAGA", "CGGTAA", "CGACTT", "CATACG", "CACTTG", "CACGAA", "CACAGT", "CAGATC", "CAACGA", "CAAGCT", 
	"GTTCAC", "GTCGTA", "GTGTCA", "GTGAAG", "GTAACC", "GCTTGT", "GCCTAA", "GCACTA", "GCAGAT", "GGTGAA", 
	"GGCAAT", "GGATGA", "GGAATG", "GATCCT", "GATAGC", "GACACA", "GAGCAA", "GAGGTT", "ATTCCG", "ATTGGC", 
	"ATCGAG", "ACTACC", "ACCAGA", "ACGTCT", "ACACGT", "ACAGTG", "AGCTGT", "AGCCTA", "AGGTTC", "AGGCAT", 
	"AGGACA", "AGAAGC", "AACGTC", "AAGCTG", "CGAGTA", "GAATCC", "GAATGG", "AAGTGC", "AAGAGG", "TACAGG", 
	"CTGACT", "CTAGTC", "CCTAAG", "CCATAG", "CGTAAC", "CAATGC"
};

const string barcodes_Haplotagging_B[] = {
	"AACGGA", "ACCAAC", "GAGATC", "CTTCTC", "GTTCCT", "TTTCGG", "GTTGTG", "TTTGCC", "CTTGGT", "CTTACG", 
	"GTTAGC", "GTCTTC", "CTCTCT", "ATCTGG", "TTCCAC", "CTCGTA", "GTCGAT", "GTCACA", "CTGTTG", "ATGTCC", 
	"GTGTGT", "GTGCTA", "CTGCAT", "TTGGAG", "CTGAGA", "GTATCG", "CTATGC", "CTACCA", "GTAGGA", "TCTTCG", 
	"ACTTGC", "ACTCTG", "CCTCAA", "ACTGCT", "TCTGGA", "GCTAAG", "ACCTCA", "TCCTGT", "CCCATT", "TCGTTC", 
	"ACGTAG", "ACGGTA", "TCGACT", "GCATAC", "GCACTT", "ACACGA", "TCACAG", "CCAGAT", "ACAACG", "TCAAGC", 
	"CGTTCA", "AGTCGT", "AGTGTC", "GGTGAA", "CGTAAC", "TGCTTG", "AGCCTA", "AGCACT", "TGCAGA", "AGGTGA", 
	"TGGCAA", "AGGATG", "GGGAAT", "TGATCC", "CGATAG", "AGACAC", "AGAGCA", "TGAGGT", "GATTCC", "CATTGG", 
	"GATCGA", "CACTAC", "AACCAG", "TACGTC", "TACACG", "GACAGT", "TAGCTG", "AAGCCT", "CAGGTT", "TAGGCA", 
	"AAGGAC", "CAGAAG", "CAACGT", "GAAGCT", "ACGAGT", "CGAATC", "GGAATG", "CAAGTG", "GAAGAG", "GTACAG", 
	"TCTGAC", "CCTAGT", "GCCTAA", "GCCATA", "CCGTAA", "CCAATG"
};

const string barcodes_Haplotagging_C[] = {
	"GAAACG", "ACACCA", "TCGAGA", "TCCTTC", "CTGTTC", "GGTTTC", "TGGTTG", "CCTTTG", "GTCTTG", "CGCTTA", 
	"GCGTTA", "TCGTCT", "CTCTCT", "GGATCT", "ACTTCC", "TACTCG", "ATGTCG", "CAGTCA", "TGCTGT", "CCATGT", 
	"GTGTGT", "TAGTGC", "ATCTGC", "AGTTGG", "GACTGA", "CGGTAT", "GCCTAT", "CACTAC", "GAGTAG", "CGTCTT", 
	"GCACTT", "TGACTC", "AACCTC", "CTACTG", "GATCTG", "AGGCTA", "CAACCT", "GTTCCT", "TTCCCA", "TCTCGT", 
	"AGACGT", "TAACGG", "CTTCGA", "ACGCAT", "TTGCAC", "GAACAC", "AGTCAC", "ATCCAG", "CGACAA", "GCTCAA", 
	"CACGTT", "GTAGTC", "TCAGTG", "AAGGTG", "ACCGTA", "TGTGCT", "TAAGCC", "CTAGCA", "GATGCA", "GAAGGT", 
	"AATGGC", "TGAGGA", "ATGGGA", "CCTGAT", "AGCGAT", "ACAGAC", "CAAGAG", "GTTGAG", "CCGATT", "GGCATT", 
	"GAGATC", "ACCACT", "AGAACC", "TCTACG", "CGTACA", "GTGACA", "TGTAGC", "CTAAGC", "TTCAGG", "CATAGG", 
	"ACAAGG", "AGCAGA", "GTCAAC", "CTGAAG", "GTACGA", "TCCGAA", "TGGGAA", "TGCAAG", "AGGAAG", "AGGTAC", 
	"ACTCTG", "GTCCTA", "AAGCCT", "TAGCCA", "AACCGT", "TGCCAA"
};

const string barcodes_Haplotagging_D[] = {
	"GGAAAC", "AACACC", "ATCGAG", "CTCCTT", "CCTGTT", "CGGTTT", "GTGGTT", "GCCTTT", "GGTCTT", "ACGCTT", 
	"AGCGTT", "TTCGTC", "TCTCTC", "TGGATC", "CACTTC", "GTACTC", "GATGTC", "ACAGTC", "TTGCTG", "TCCATG", 
	"TGTGTG", "CTAGTG", "CATCTG", "GAGTTG", "AGACTG", "TCGGTA", "TGCCTA", "CCACTA", "GGAGTA", "TCGTCT", 
	"TGCACT", "CTGACT", "CAACCT", "GCTACT", "GGATCT", "AAGGCT", "TCAACC", "TGTTCC", "ATTCCC", "TTCTCG", 
	"TAGACG", "GTAACG", "ACTTCG", "TACGCA", "CTTGCA", "CGAACA", "CAGTCA", "GATCCA", "ACGACA", "AGCTCA", 
	"TCACGT", "CGTAGT", "GTCAGT", "GAAGGT", "AACCGT", "TTGTGC", "CTAAGC", "ACTAGC", "AGATGC", "TGAAGG", 
	"CAATGG", "ATGAGG", "AATGGG", "TCCTGA", "TAGCGA", "CACAGA", "GCAAGA", "GGTTGA", "TCCGAT", "TGGCAT", 
	"CGAGAT", "TACCAC", "CAGAAC", "GTCTAC", "ACGTAC", "AGTGAC", "CTGTAG", "CCTAAG", "GTTCAG", "GCATAG", 
	"GACAAG", "AAGCAG", "CGTCAA", "GCTGAA", "AGTACG", "ATCCGA", "ATGGGA", "GTGCAA", "GAGGAA", "CAGGTA", 
	"GACTCT", "AGTCCT", "TAAGCC", "ATAGCC", "TAACCG", "ATGCCA"
};

const string barcodes_stLFR[] {
	"TAACAGCCAA", "CTAAGAGTCC", "TTACTGCCTT", "CGCTGAATTC", "TGACGTCCTT", "TGTGTGTAAC", "AGATTGAGAG", "GATAATGATG", "TTGAAGAAGC", "CTGATACTTC", 
	"GGCTCCAAGC", "AATCTACAAG", "TGAGGTGGTT", "TAGAAGACTC", "CACCGAGCAT", "CGAATATAGG", "CGACGGACGA", "TACGCCGATT", "AAGATTAACC", "GTCGGACCTC", 
	"ACGGTACCTC", "CATGGACTCC", "TTCTCCTGTA", "CTGTGAGCAG", "TAACATAGCA", "ACGACGTCCA", "ACTACGGTCG", "TATCATCTCT", "AATGCCAGTT", "CATAGGACTG", 
	"CGATCACGAT", "CCGGTGGACA", "AGGTTAGGAG", "AAGATCAATT", "TAGGCTCAGG", "CTTCGATGGT", "CCGGTCACTA", "GCAACAGAAC", "AATCAGTGGA", "AGACTTACGT", 
	"TCCTCCTGGC", "GTCAGTACCA", "GGAACGTGCC", "ACTGCACAGG", "TTCCACTACA", "TCCGTGAAGC", "GATATAGCGG", "ATTATTGTGC", "ACGCCTCAGT", "CTAGGCCACG", 
	"CTCACCTAAG", "GTCGCGATTC", "TATTACTCTC", "GTTCGCGCGC", "GATTGATGGC", "GCTTCTGTAC", "GCCGCAATAG", "AAGTAAGCGT", "TATCGCGCCA", "CTCTCCATAG", 
	"ACTGGTAAGA", "TCTAATGACG", "AGTCTATCGT", "GTACATAGGC", "TCTCGCGACC", "AGCCTCCTAT", "TATATGCGTT", "CATCAGTTGG", "ACGGCGAGGC", "TATCGGAATA", 
	"AATACTCGAT", "CGGACCACCT", "GTTAAGGTGA", "CAACTCGTCG", "GACATCGGTG", "TCATCGGATT", "ACGTTGAATA", "TCCTGGTCGG", "TGGTATTGTC", "ATCTCGCTGC", 
	"AGTTCCTGTC", "CCTGCTTAGT", "TACGGACTCT", "TCCAGGTACT", "TTGATCTTAA", "TGACTTCGTA", "CCTTACCACA", "CTCACAGCGG", "ACACACCAGA", "TGCGGTTACG", 
	"TACTCAGGAA", "CTGACCTGGT", "TAGTTAGACG", "CCTGCGAACA", "TTCAGTGCGG", "GACTCACCTT", "CATCGGCATG", "TGGACGTGAC", "AGTCGACCAT", "GTGTCACGAG", 
	"GTGTGAAGTT", "ATGGCGTGCG", "ATGTCTGCCA", "CGTTCTACGA", "GATTCAGAGA", "GATGAGAGGT", "CGTCCGCTCT", "GGTGCAACCT", "CTTCGGCCGG", "GGTTGGTAAC", 
	"GACGAAGACT", "GTAGTTCTAG", "AGTGGTGGTA", "ATCCACTAAG", "ACAGTCCACC", "TCAATGGTGG", "GTAGACGGAG", "TCAGATAGAT", "CTAGACACTT", "ATACTGAGGA", 
	"GCTAGTAGAA", "GTTGTTCTCT", "GAGGTGCCGC", "CTTCCTACAG", "ACATAGGTGG", "CTTGTAGCAT", "AGTTGTCGAA", "AAGGAATGTA", "CTTCCGATAA", "AGAGAGTCGC", 
	"AGATTGCTCA", "TCGACATCGA", "CTAAGATATT", "TCCTGATCAA", "TCGATGCGCG", "TACGAAGTGC", "CTTCAATAAT", "TGTATTCAGA", "GTAATACGTC", "ATTATCTAGC", 
	"TACAGTACGC", "GAAGACAGAC", "TGTTCACTAG", "CGTTCGATGG", "ATTCCGCTCA", "TGGAGCTCAG", "GTAATTGACA", "TACGTTATGG", "TACCATGACT", "GCAGGTATCC", 
	"TGTGTAGCCG", "CTATCACTAG", "GTGGCAGGTG", "CTCGATTGGT", "TTAATTGCTC", "ACGAGGCCAT", "ACGACTATTA", "GCCGGTTGTA", "GTCTACCTTC", "TCTTCATCTC", 
	"AGCAGGCAGC", "GCGTTAGATG", "TGCTCAGAAG", "GTTGCGCCGG", "GTATTGGCTA", "CGCGAGTCGA", "CACATGTTAC", "AACCTTGCGG", "CTTGCCTACA", "TATGGCATAG", 
	"AACGAGGTCC", "GCTGTGAGAC", "GTAACGACGC", "GGTATCCTTA", "GACTCTCCAA", "TAATCGGTGC", "TGCTGCAGTA", "ACCAAGCCAG", "TACTAATGGC", "GTTGGATCCT", 
	"ATTAACCGGC", "AACGAGTATT", "GCTGCTCTTC", "TCTTGGAGAC", "TAGTTGTCAC", "AGAGCGGCAA", "GCTCTCGGTT", "GCCTCTGGCT", "AGCACAGTTG", "GCACATAGTA", 
	"CGCCTGGAGA", "CGTATCGCCT", "TCCGACCGTA", "GCCTCCGGTC", "GCGAACACGG", "GGAGCTATGG", "GCGGTGCAGA", "GTCTCACAGT", "AAGGATTGAT", "ATTAAGATAA", 
	"ACAGTTCATT", "GGAGTACGTA", "AGACCGTCGT", "GTTACGCGTC", "CTTGTCTAGT", "AGAATTCAGT", "GCAACTATTG", "TGTCAGGCAC", "CGATTCCGGC", "AGCAAGCGAC", 
	"CTAGTGTAAC", "GATGGAGTGC", "AATTAATCTG", "AACTTGCACC", "GAGTAGATTG", "CATATTAGGA", "AAGTCGATCC", "TGGATTAGAC", "CGAGGTTAAC", "TTCTCACTCC", 
	"AGGCACACCT", "TTGAGCCTGA", "TGTTCCAGAC", "CGTTCAATAA", "TAAGCGTGGC", "ATACCGAGCT", "CCTGTGAAGT", "CAAGCCTTAG", "GACCGTCTCT", "ACCGACATTA", 
	"TGCAGTGGTA", "ATCAGAGTCA", "TCGCATCTTG", "ATAACATACG", "AAGTTGGAGG", "ATCGGACAAG", "AGTAGTAATA", "ACTTCAACAC", "ACCTCTCCTG", "AGTGTTGCTT", 
	"AATTCGAAGG", "CGTATTGCTC", "TTACCAAGGT", "TTCACTGGCG", "TTGATATTCC", "AATAATGCAC", "ATAGAGAGCC", "TCACGAATCC", "ATAGCTGGAT", "TCGTAGTCGT", 
	"TTGTGGCATC", "GCACTTCACG", "TGCAAGATAC", "ACATCATTGC", "GTGGACAGGA", "CAAGAGGTCA", "AGACCGATCG", "CCGAATCAGT", "GCCAGCTTCG", "CACTCAGTCA", 
	"CTTCCATCCT", "TCCTCGTGCG", "ATGGAAGGCC", "GAACGAATTC", "CGGACTGAGG", "TTCGTAGCAC", "CGGACGGATT", "CGATTGATAA", "ACCAATCTAA", "TAGTCGGAAT", 
	"ATGTTCTATA", "CCACTGGATT", "CGCCTCATCT", "GCTTAGTCGG", "CCAGATTCGG", "CGATAAGCAC", "CATGCTGACT", "CGCACTGGAT", "CCGTGCAACC", "TGATGACGTA", 
	"CTCTAAGCAT", "CAGTCCAGGA", "GACTATCACA", "AGTCGTCCTA", "GCTTAGATCT", "ATGCGAGCGT", "GCTCGGTATG", "GCGTACAACT", "CTCTGTCTCC", "TTCGCATAAT", 
	"CGCGACGAGG", "CAGGAACGAG", "CAAGCTGAAG", "ACATATGCGA", "ATTACCGCGT", "ATCAGTTAAC", "CGATATGTGC", "GCCAGCACGT", "ACGACAATAT", "CCATAAGGAG", 
	"TCTGTTCATA", "TCCGAAGATT", "GCTCAGCTGG", "CACCGCTAGT", "CGCCTCTCGG", "GTCACGTGCT", "ATCGAGCGAA", "GTGAACGAGA", "TTGACATCTC", "CGTACCGCGA", 
	"ATATTCTCAC", "GCGCGCCGAG", "GTTACCACAT", "TATAGAAGCT", "TTAACCGCGA", "CTCAAGTCGC", "TATTCATATA", "AAGTCACGGC", "CTATCGCCTT", "GACTTACTCT", 
	"GGAGGTTGGT", "CAATCTCATG", "CTCGCAGGTC", "CATCAATTAA", "CTTCTTCGTA", "GAACGCCGTG", "GACCAGAATA", "TATGAACAGC", "AGTAGTTACG", "CGATACTAGC", 
	"GGCTATATAA", "CGAGGAGTCA", "AGGAGCTCCA", "GGACACTCTT", "GTGCGAACAT", "ATGTCCGTGC", "AGCGATCCAT", "CAGTGTATGC", "GGAGCAACAG", "GGCCTGCACA", 
	"GGCACTGAGC", "AAGTTCGACC", "GCGAATCGAC", "CTGCAGGATA", "TTCTATCAGC", "CAGCAAGCAG", "ACATGCCTGT", "ATACTCCTAC", "TAGCTGCTAC", "TCCAATAGAA", 
	"CGGCAGAGTC", "ATGACGGATA", "CAGCGAAGTC", "AAGCGCCGTA", "AATTCCAACC", "AGTGGACACC", "TGTAAGAGCT", "TTAACGTAGG", "CTTGCCGCGC", "GTTGTCTGAC", 
	"CTTACCGGTG", "CTCAACATGT", "GACTGTTCCT", "ACGGCCTGTA", "GAGACAACGG", "CGGCGATACC", "ATGATCATTA", "GTCCAGACGA", "CATTCCGGCT", "TAGCGCACGT", 
	"GGAGGTCCAA", "GTCGTCGATC", "GACCAGTACG", "CCGAGCGATT", "GTGTCCATAC", "CTCACTGTAG", "GGCTTCGTCG", "AACAGCATAG", "CCTTAAGGCT", "CTTGACTCAA", 
	"CCATCGTGAC", "GATATGGCAA", "GCGCAGCACC", "GCGACAGTGT", "TACGTAACAG", "GCATTGTACA", "ATGGAACAAT", "AGATCTCTTG", "CACTTAGCTA", "TCAATCGTCC", 
	"AACGGTTGAG", "CCGATGAGCA", "TAGGCGCATT", "AAGGTCGGTT", "CGGACCTTGG", "AGAAGATATA", "AGGCATACTC", "TGCTCTATTC", "GCGTAGGCTT", "CGACATCACA", 
	"AGGAGCACTG", "ATTCCGAGAG", "TAATATTCAG", "TACCGTCAGC", "GACGTCTACA", "AATGGTACGC", "TGGAGAACTC", "TGGTACTGCT", "CGGTAAGATA", "GACTCGCTAG", 
	"TTAGGCTATG", "GCGAACGATC", "GGCTCGTATA", "TCATCCACCT", "GCATTGGCGC", "CGGTTGGTGG", "AGTGCCTACT", "ACTTCTATGC", "ACGACGACTG", "TGATCGAGGT", 
	"TTCGACCGGC", "ACGTCTTATA", "GCTCAACTAA", "CGATACGCCA", "GCTTCCTAAC", "GTTGACGGCT", "AACCTGGTGA", "CTTGAAGAGA", "GAACTGCTCC", "AATGCGCCAC", 
	"ACTGATCTCG", "TCTGAATGAT", "TGATGATTAA", "GTGTTGGAAC", "AGTCCATCCA", "AGCGCTGGAC", "CTTCTACGAT", "TACAAGCACC", "CCGAGTGACC", "CCTGTTGCGC", 
	"CCACGCTGTG", "CCTGATCGCA", "TCGTGACAAG", "CGAAGGCAGA", "CTCCTGAGGC", "TCGCGTTAGG", "TGATAGCATG", "GTGTTAATGT", "TTCCTCATAG", "ATGTTGAAGC", 
	"GGCACGACGT", "GCGCACCAGG", "CATAGAATAT", "GAGATGCGTG", "AGGCACTTGG", "GTGTTAACTA", "TGTCTGTCCT", "AGTTGCGATA", "GGTGGAAGGT", "GACCATGTGC", 
	"ATATGATGCT", "CGAGCAGGTT", "CTACTCATCC", "GCATATACAA", "TACGATGCAC", "GACCTATCCT", "AGACGGTGCT", "ATATTAGAGC", "AGGATCTGCT", "GACCTGCGCG", 
	"AACCTAGTAG", "GTCACTAGGA", "TCTGGATAGT", "GCAAGACCGT", "GCCGGACCGA", "GAGAACGTGT", "TTGCAGCCGC", "TTGGCGCTTA", "TTACTCAGAC", "CTTAAGCTCA", 
	"ATTAGCGGCT", "CTCGAGCTAG", "AATACATTAA", "GTAGGATCAG", "GTCGCCTCTT", "ACCTGGCGCT", "CAAGATGCCG", "GCGTATAATC", "ACCTGATCCG", "GGAATGCCAC", 
	"CTCGGTTAAT", "GCTTAATCAA", "TTCTGCACGG", "TGACGACCAA", "GTTCGCAATG", "AGGCCAGTCG", "ACGAGGTGGA", "GTCAGGCGAA", "TTGTCTGCAG", "ATCACTGGAA", 
	"TGCTGGTGGC", "CCTGTTAATG", "GCGATCGTAC", "ATGATCTTCG", "AAGGTTGGCC", "TATACATTCG", "AATACCTTCC", "GGAGTATTAA", "TATAATTAGG", "ACATATTACC", 
	"GGCTTGGTGC", "AGCTTAGAGT", "TATCACTGAC", "GATCATCAGC", "GATGCGCATG", "GATTACTACA", "ATTGGAGTAG", "CCACAGGTAT", "GTAGATCATG", "CAGTAAGGTG", 
	"AGATTAAGGA", "TTACTGCTGA", "GACGTATAAC", "ATCGAGATCG", "GCAATTACCG", "TCTTCCACAG", "GTGACATACA", "TTATGTCTAG", "GCAACCTCGC", "CATCCTGCAC", 
	"CTGATATGAC", "CTCGAACTGA", "GTCCTAGTGC", "TCTGCGCACG", "CACGGAGACC", "CGATGCTGAC", "ACCATCGGAA", "GGTACTCTGC", "TATTATATGG", "GCGCGTGATG", 
	"CCGTCCGGTG", "AAGACCTAGG", "CACGGATCGA", "CACGGTCGTA", "GCCACCAGCT", "TGGCTCAAGG", "CCTTGACGGC", "ACGATGATCG", "GAGAATTAGT", "TCCTAACAAT", 
	"CCGATTCCAA", "CGACTACTCT", "TTCGCCGTTA", "ATAGAACTGC", "TTATAGTAAT", "ATACCACCTA", "CATTATCCTT", "TGGCATTCAC", "CATAGTATTA", "AGACCGCGAA", 
	"ATTCACGGTG", "ACGTATACAG", "AGCTATGCTT", "CTGGTTCGAG", "TAACCGCACA", "ACGATGCGAA", "GATGAACTCT", "TATAAGTATT", "GTAGGCGAGG", "CAGGCTAAGA", 
	"GATTCCTCAA", "GGTAGGATAC", "AGAGACTCCG", "CATACGCCAT", "AAGACTGTGG", "TGAGTTGTGT", "CATTGTTATA", "ACCGATGACG", "AACATATGGC", "ATTAGTGGTC", 
	"TTGGTAAGTT", "TTATCAATAA", "AGGTGTAAGC", "GTATTCTATG", "AGTTGACGTT", "GCAGCCTGTG", "TTCGCCGCGT", "GTTATTATAC", "GCCATCAGGA", "CTTCCAATGG", 
	"GAGGCAAGTC", "ATATTAACTG", "AATTGCGGTG", "GTCTGGTACG", "CCATGAGAGG", "ATGTTAGTAT", "CGGTAGACAC", "GTGGAGCTAC", "ACCATAGGCC", "TTCTGTCGAC", 
	"AGATCTTGAG", "AGTCGATGGA", "GGCACCTTGC", "TATGCCTGAT", "ATGCGCTAAT", "GAGTATTCCG", "GTTCCTGGTT", "AGGAACCACT", "TGGCTTGCCT", "GAGACATTCT", 
	"ATACATGGAC", "GTTATAACGC", "GCGCGTACGC", "CTGTACCTTG", "CTTAGCAACC", "GTCTACTGAC", "CTTGGTGGTT", "ACCAGCGCAT", "TTCTCGAGGC", "CTGCTAGTGG", 
	"GACTCACTGA", "TAATTCGCTG", "CTTCACGTTA", "ATACCTCCAT", "GCCGTGTCGG", "AGTTATCAGA", "GGACGCCATA", "AGATCGCCTA", "TAAGGTATTC", "TAGGAACTCG", 
	"TATACAATTA", "ACCTCGTGAA", "AGGCTAGCTG", "CACGTACCAA", "GGAGTGACCT", "GCCAGTGACG", "TACCTCGTAC", "AACTTCCAGG", "CTGCGAGGTG", "CCTCCAGATT", 
	"GGCTAGACAG", "TATCTACACA", "CGACCACTGA", "GGCATTATCG", "AACGTGCGTT", "GCCTGTGCGT", "GTAACGGATG", "GGTGAGGTGC", "ATTAATACAG", "TTCCGTCCTC", 
	"CCGAGTTCGA", "TTCGTTGTGC", "CAACATATAA", "AATTACACAC", "ACAGAACCGT", "TCCTCGAGTA", "GTGAGTCATA", "TTCTGCTTCT", "TGAAGTAATT", "TGCCTGCCTC", 
	"AAGGAGCTGG", "ATATCCGAAT", "GTACCGCAGC", "CTCCATGGAA", "CGAGTTCAGC", "ACATTCTTAT", "ACGCAACTCT", "ATAAGCGGAG", "ATGGCGAGTA", "CGACCTTGGT", 
	"ACCTGAACTA", "TATATCACAC", "GCAACACGCT", "AAGCGATTCC", "CGCAATCCAC", "ATCATTCAGC", "GCGAATGACT", "TCGACACGAT", "CAGGTTAACT", "GCCACTAGTC", 
	"TCCTAAGGCC", "TACATGTGCC", "GTGACCTAAC", "CTATGAACGA", "GGCATACGTT", "CGTTGGTGGT", "GCAATGCGAG", "ATACTCTGTC", "TTATGACCGG", "CACTCTTAAC", 
	"GCCTGGAAGC", "AGATCAAGCT", "AGTGAAGAGT", "AATCAACTTG", "ACCATTCATA", "GCCTACCTGA", "TAGTCGCGCC", "GTGAAGTCGG", "GGTCCGGTGT", "AGCCACCATT", 
	"CCTAGGTACC", "TTAGGCAACA", "TATTGTCACA", "GACAATTGAG", "GAGCGGTACT", "TTCCAGGCCG", "CAGAGTGCCA", "GTGTTATCCG", "ATCGTCGCAG", "GAATGTAGGC", 
	"TTCGAGATAA", "TGTCTCATCC", "GCCGCGACTT", "GCTTCAGCGC", "GACTAAGGTC", "GTATTCAACA", "ACGTCCGCGT", "CATTGTAACG", "ATTATGAATA", "ATGGTTCTAA", 
	"GAGCCACAGG", "CTGGATTAAG", "GCTAAGGTTC", "CATTAACCAA", "AGCATGTAAC", "GCAATGATCA", "CACAGGCTGC", "TCCGACTTAA", "CATCATCGAT", "ATGAGGCTAC", 
	"GTATACTTAG", "GCTAGTCTCG", "GTTAAGGCTT", "GGCATTCGAA", "GCTCGGAACA", "TGTCTACGCG", "CATCGTCAGT", "CTGGCTACTA", "GTGATATAGT", "GCTAAGTACT", 
	"CGAAGACAAG", "TACGTGACGA", "GACTCGTGTG", "CATCATTCGA", "AGTGTTGTGA", "ATATCACGAA", "AGGCAGTTCC", "AGAGCTGTAG", "GGTAGCTCAT", "AAGTCCTCCT", 
	"CCGATACTGA", "AATACGTTGG", "TGACCGCGCG", "GCGTATGCGG", "ATTCATCAAG", "TCAGGTAAGT", "CATCACTTCC", "GCGGAGTGGA", "TACAGCCGAG", "CCTGAATTCT", 
	"TGAGCATAAC", "TTGGACGGCG", "CGATTGCGCG", "CGCTGCCGTG", "GGCATCGATA", "TTATTGTCGA", "TTACGATTGG", "TAAGGCTCGG", "CTTGCCAATG", "CTGCACTCTG", 
	"AACTGAGCGG", "AGTGTCTAGA", "GGAGAATATA", "CGGCGTGTAA", "GGTTACAGTA", "ACTACCGTGC", "CTCGGCGTAT", "TTCAGCAACT", "ATTGGCAACT", "ACTAACTTAA", 
	"GCAGTCAGGT", "CCGCAGACTG", "GTGCCACCTG", "TTAACCTACC", "CAGAGCGTGC", "TGAGACCGGT", "AGTGCCGTTC", "ATCCTCATCA", "AATGGAATAC", "GTCCGTCACA", 
	"GTTCGTGTCA", "CAGCGCCTTG", "ACCTTAAGTT", "CAGCCAGACG", "CTGGCTATGT", "ACGATGTCGT", "GGCGGCATGA", "TGCAACCGGA", "GGTATATGCC", "GCTCTTCAAT", 
	"CTATCAAGCA", "GCTCTACATA", "GCTGGCCGTA", "GCTTGCGGAA", "CGACGCTCTC", "TAACATCTAG", "CTGGTGTTGT", "AACGTATCTG", "CCTAGAGCAG", "AGAACAGGTC", 
	"GAGATTCGGT", "TGCTGCTGCG", "TCGTAGATCG", "GACTCGAGCA", "CCTGCCGCTA", "TTCCAATAAC", "GAGAAGAACA", "CAGGACATAC", "GCTTGACAAT", "ATGGTACCGA", 
	"CGTATGAATT", "TGCACTTATC", "TGTAATCCAT", "GCTAGAAGTT", "CCTCGTTGAC", "CGCGAATCAG", "GGCACTTCCA", "TTGACTCGGC", "GTGATCGCAT", "CATAAGCACT", 
	"GCCATGTGTC", "TATTGAGGTC", "ACGTGCGGCT", "ATGTACTTAA", "AGAGGAGGCG", "AGCTCATCGC", "TCGAGCCTTC", "CAAGACAAGT", "CTACCGCGAT", "ACAGAACTTA", 
	"GCAAGTCTAT", "TGTCTTACGA", "TTGTACACGT", "ATTACTAACG", "CTCCTCCTAA", "TCTACAGTGG", "TGTAGAGTCC", "CTATCGCTGA", "GCCGCGATGA", "CGCGATACCT", 
	"TCGAACGCCA", "TACGCTTCGA", "GATCGTGACT", "TAGATAGCGT", "TAGTAGGCCT", "GAGAGCCTCC", "ATTGCTGGCG", "CATCTATATA", "AGTTCTTGCT", "TGATCTAGTG", 
	"TGTACTGGAC", "GGTGCGCGTA", "CGTCTAAGCA", "CTTGTGGTTG", "GACGCCTAGT", "GAGATATTGA", "TGAAGGAAGG", "CTTCATAACG", "ATTAAGCGCG", "ACGTTAGCAC", 
	"CTATCTCCGG", "TGAGGCGGCC", "ATTACATAAT", "TAAGGAGAAG", "GCGTCCACAT", "TTCGGCGGCT", "GAATTAGAAG", "CATACATGAG", "GATATTGTAG", "ACTCTGCCGC", 
	"GCTGGATTCC", "CTCAGAGGCG", "GGTTAGTGGC", "ATGACCTCTG", "AAGATGGCCT", "GAGAACGCTA", "AACACTCTCC", "AACGCAATTA", "CACAAGGTCT", "GTCAGGTCGT", 
	"CTTCGACCAA", "CTCCTCTGTA", "CCTAACTGAG", "CGACGTTCCT", "TGCACTGTCT", "TCACGCCGCG", "TGTCCGCGAT", "ATCTGAGATG", "ACAGAATGAA", "CGATGAAGTG", 
	"CGATAGTACG", "CAATCCGGAG", "GTAACTGAGT", "AGTCATTATT", "TAGAACTCGA", "GTCAAGCAGA", "TAGCGCTTCG", "GCCAAGCATC", "CGTGGCAACA", "GTCCTTGCAC", 
	"AGTGACCGGA", "CTGGTCACGC", "AGCCGCACCG", "ACGTAGCGCG", "GACATTCAAG", "CTGTGCTAGG", "ATATCTCGTT", "CGCTGCTTAG", "GTTATAGATG", "AATAGTTGCA", 
	"GCGGAACTTG", "AATGGCCGAG", "CTTCCGCGCG", "CCTCGCCTCT", "CTCTTCCGGT", "AGTCCAGAGC", "TTGAGGAGAC", "TTCGTGAATA", "GAATCGGACT", "TGTCCTACCT", 
	"GTGCTTAGTT", "GCAATCGAGT", "GAAGGCAAGC", "CACAGAAGCC", "CGTCTCAGAC", "CAATAGAGTG", "TAATTGAACG", "TGACTGTCAG", "GATTCGGAAG", "CGGTCCAAGG", 
	"CATCTCGTAT", "GGCGTCAGTA", "AATCACAGTC", "CCTATCCAAG", "ATGGACGGAA", "TTGAATGCGT", "GAGACAGATC", "GCCTGTGTTA", "GAACGCTTAG", "GTAACCACCG", 
	"GACTTGAGGT", "GACCTAATGG", "ATCGACATGC", "TGTTCTCCGG", "GTGCTAAGAA", "TCAGCGCAAT", "TCTGATAGCG", "CCTCTATCTC", "CAACATCGCG", "TGTCGTAGGT", 
	"GTTGCATGGT", "TGTGCAGTTG", "ACATCTTCAC", "GTGCGGTTCA", "GACTGGACGA", "TTATTGCGAT", "ATATTATCCA", "TGCGCACATG", "TAGATCTAAT", "TTATTCACTC", 
	"TTATGGCCAA", "GAACCGCCTC", "AACGCATTCG", "TTGTACATTA", "AAGGACCTCC", "TGACCATCCT", "ATGATGTTGC", "ACACGCCGAA", "TATGGAGACT", "AATCAATGAG", 
	"AGCTGCTGAA", "CTCTATGCTA", "ATCGTATAGG", "GTAGGTACCT", "GGTTGAGCCG", "AGACTAATAT", "ACGATTACCA", "AATAGTAGTG", "GCCGGTAGCG", "CGTTCACGCG", 
	"TGATTGCTAG", "AGCTTGATAC", "AGTCGGTGAG", "ACCAATTGTA", "CTCCGGATTC", "AAGCGCTTAA", "ATCTCAAGCC", "GAAGAGCCTT", "CCTTCCTGGT", "CGGCGAGCGA", 
	"TTCCTACGAC", "TCGACAGACC", "CAGATTCAAC", "GCCAGGTTGC", "CAATTGCACA", "TCATCTACTC", "CATCAACGTA", "CATGGCTGTA", "GAGAGTTGAT", "GCCGCACGCA", 
	"GGTGATTACG", "TCCTTACCGA", "ATATGTCTCA", "TGAGCGGCCG", "GGCACTACTG", "CTACTACGCG", "TCCTGCGAGA", "ACCATACAAT", "CCTGTCGTCA", "TTCGGACACA", 
	"GGAGCCTCTC", "TTATTCGAGG", "ACTGTCGGCG", "TTACGAACCT", "TAATCTTACG", "CTGAGCATGG", "GCCGTGATCT", "CGATCTCGTA", "TGCGAGCTCT", "GTTGCTCCTT", 
	"TCTACGGTAA", "TAATAGTTAA", "ACGCAGAGGT", "ATATATGCTC", "GGCATGACCA", "TCTCGCCGAT", "ACGTTCTAGC", "CACCGTGTGT", "TGATGGACCT", "AAGTCCATGG", 
	"CCAATTCATC", "CGATCGTCAG", "TACACTAGCC", "CAACAGACAG", "CACGTGCCGG", "ACTTCGGATA", "ATGATCACGT", "GAGATGTTAG", "GATTACAATG", "GCAACCGACA", 
	"TGAGCCTACA", "CGTGAGCTCC", "GATAACATCA", "CACGCATGCA", "AGACTCGACG", "TGTGGCCACG", "TCACGGCGGC", "ACTTGGTGTC", "TGAGCCGCGC", "AATCGCGCAG", 
	"TACGGTTGCA", "GGCTATCGCG", "ACAGGATAGA", "GACCGCTGAC", "AGTAAGCTCT", "GTCGAAGTCA", "GAGAGAAGCG", "ATTGGTAATC", "CGAATCGCAG", "GCCTCTCAAC", 
	"TCGACGTCAG", "ACTAACCGTA", "TAACATTGTG", "TTGTCGGTAA", "TGGCTGAACC", "AACATTCCTC", "TCTTATGTTA", "ACGTTCGCCA", "ACGTTGTACG", "CATTAATGGT", 
	"GTATATGAAG", "TGACCTACAG", "CCGGTCTCCG", "ACATTCCGTT", "CGCGAGGACC", "TAGTACAATT", "CGCGTGTAAT", "AGCTGCCTTA", "TGGTTCGGAC", "AGCAGCCACG", 
	"CGGTCGAACC", "TCCTGGATCT", "CAGAATCTTC", "CAGCCGTCAC", "GACCTGATAA", "ATGCTCCAGT", "CATGGCAGCG", "CGTACGTAGG", "TTCAGTAATC", "TGGAGAGAGG", 
	"GCGCTGTGCC", "ACACCACCGC", "TGAGTGGTTG", "CACCTAGGAA", "ACGACGGAGC", "AGTGTGAATC", "ATCCGGTGCC", "GTAGTACCGG", "CAGCTAATGC", "CTGGCAACAT", 
	"GAGCATAGAA", "TCCGTCAACG", "TATGCTCTCC", "TCGTGGCAGA", "CGTCGGTTAA", "ACACCTTGTC", "TATGTATGGC", "CTCCACCATA", "CACTAGATTC", "TCATCGTTCC", 
	"ATAATGTAAT", "TAGAAGGAGG", "CTGAGACGGC", "TCTTACAACG", "TATCAATGCA", "TGACAATAAT", "AGAAGTTAAT", "TACGCTGACC", "ATGTAATTCC", "TGAGCCAATG", 
	"TTATGGTGGT", "AACGGTCTTG", "ACCTACGGCG", "GCAACGACTA", "TCGCTCGGTG", "CACAACAAGA", "TCGAGAAGTG", "TTATTGGACC", "TTGATAGATT", "TGCTAAGCCG", 
	"CTTGACCGGT", "GCCGCTGATC", "TGGTTCCACT", "CTAGCGGCAT", "CGCGAGCGAT", "ACGGATCAGC", "GTGCGTATGT", "TACTGCAATG", "ATACGTTCCA", "CGACTTAGTG", 
	"AGCGAGAGCT", "ATTATCGCCA", "CCGACACTCT", "AAGCCTTGAC", "ATGTAAGATT", "CGTACTTATT", "CTCTTCTCAA", "TCATCCTTGG", "GTGTCAATCA", "AGCGAACTGT", 
	"GGATATCTCC", "CGATAGAATA", "CGTCTGCTGA", "GTGGTGCATC", "GTTAATAATC", "GGTGTTATAA", "AATACACGTA", "GACGTCAATG", "CGTAGCGGCA", "ACGCATTGCA", 
	"AATACCGATT", "GTGGAGTGTC", "TGTCGGAGTG", "CACGCTCTCT", "TTCTCCAGCG", "ATATTGATAT", "AGGACCTGGA", "TTACTATGTG", "GAAGAATGTG", "CATAGCGACG", 
	"CAGGTGTAGA", "CGTAGACACT", "CCTCAGGTCG", "TACTGTGTCA", "ATAGGAGTCT", "TCCGTCTATA", "GCTTATTCTT", "AGAGCCTAAG", "AAGGATCCGA", "TGATGTCGAT", 
	"TCTTCTCGGC", "TACGCATTAA", "ATACGTACTG", "CATCATGACC", "GTGAGCGGAA", "TCGCGGTATT", "GTACGAGTTG", "CGTCAACATG", "TGAGTCTAGT", "TCTGGCGCAT", 
	"ACACGAGAAT", "AGGATTCTAC", "CCTTGGATCC", "TACTGCTACA", "CCGCATACGT", "TTCTGTGACT", "CACGCAAGTG", "CCTGATATAG", "AACCATATAC", "GTCTGAGCAC", 
	"CACTATGAGA", "AGGAACGGAC", "GATAGGCACC", "TCAACCGCTC", "CTTGAGGAAG", "ACAGTCGGAT", "AGCGAACCTA", "AACGTCACAC", "CATGATGCAT", "CTATGTATAA", 
	"AGCCGGACGC", "TAATTGTATA", "GATAATACGC", "GCTACGTCAT", "CTATGTCGCG", "CTTCATTATA", "GTAATATTAC", "GCTATGCGCT", "GCTACCATAC", "ATATACGCCT", 
	"TGGTTAGGCA", "GATAAGGAGT", "ATACGGATTA", "TCACTGCTTC", "GCCGACGCAT", "TATCCACAGT", "TTACGTTCAG", "ACATTAGATA", "CGCACTCACC", "AATCCGCACT", 
	"CACCTTGGTT", "GTCTAATGCA", "TCTGTAGGCC", "TGGTGAGCCT", "CTACCGGACC", "TTCGAATCCT", "AAGTACTAAT", "TCCTTATGAT", "AGTTCACTCA", "TAGTCTACAC", 
	"ATTATAGCAC", "TCTGATCTAA", "ATATCAATCG", "TGGCCTTACC", "TCCAATCTCG", "AATATGACGA", "TAGGACAGCC", "GTCAGGATCG", "TTATTATCAG", "CATGAAGTGT", 
	"AATTGTCAAG", "ACGAGTAGCT", "CGCACAGGTA", "CCGGCAATAC", "GAGTATATGT", "CGGCCTCATC", "ACTACTGCCA", "TACCACACGG", "CCATGCTCAG", "TTCGTGTACG", 
	"GTAGTGCCAA", "CTAAGCGTAA", "CCGGCTGATG", "CCTCATGCCA", "AGAGCAGCGG", "TCCGTGGCTG", "TGCCTCAGAT", "CTGGCTTCCG", "GATCCGAGGC", "TTATCCATCC", 
	"ACAGATTGTT", "TGAGGCCAAT", "ATGCACTGGT", "ACGAGACTTG", "ATGACTATGC", "AATACGACCT", "ACGCATAGTG", "TGACGTCTGA", "ACGTAATCCT", "GCAAGATGAA", 
	"GTGTCGCGGA", "GGAAGCATGT", "GACAGTTAGG", "CTTGGCGGCC", "GGCCGATACG", "GCAGCAAGAC", "ATCCTGGAGT", "AGGTCGCAAC", "AACAGACGAC", "AAGGTTCAAT", 
	"GATGATTGCA", "GAGACTCGCA", "TCCTTGCCAG", "GCTATTACAT", "AGTTGTATCG", "ATCTCCAGAA", "ACCAATAGCG", "TGAGTAGCAT", "GTGACTGTAC", "TAGCGCATTA", 
	"CACTACTTGA", "ACCGACTTCG", "GAGGCCTGAG", "TCGGCGCTGC", "GGACATGAGA", "GGCAGACTGT", "CTACGTAGGT", "CTCAAGGACA", "TAACTTCATG", "GAGTCCGTAG", 
	"GCGGAGCCAT", "ATAGCAGGTA", "CAGCGCTGAG", "ACCGCGGTAA", "GCCGTCGAGA", "TATCGCTAGC", "GTGCCGTGTT", "ACATAAGTAA", "AACAGAGACT", "TGAAGCAACC", 
	"GGAGATTAAT", "AGCCGCTCTA", "CTCTCACGAC", "GAGTAGTTCA", "GCGATAGTCA", "TGACACGTTA", "CTGAACCACA", "TGTTGTCGCG", "ATCGATTCTG", "GTTGTATGCA", 
	"TTATCGATGG", "TAGCGGTTGC", "GGTGCATTGG", "ATGGTATGAT", "CGCTCGCCTC", "GTAACCTTGT", "TCAACGAACC", "TCGTTATAGG", "ACACTGAGTC", "CTTGTTGTGT", 
	"GCCAGCATTA", "GAATCCACTT", "CCGAGGTTGG", "GTTAACAACT", "GACGTGGCCG", "CTGAGGATCC", "GACGCGGTTG", "GTCAGCATGC", "CCTATTGGTG", "AGTTGTTCGT", 
	"ACGGCACTCC", "AGCAACATGA", "CACCTCCAAT", "TATCTCGGCT", "AAGAAGTCAC", "CAGATCGGTC", "CGCGCATACG", "TGACACGCGT", "GACTGTATGG", "GTCGCCTTGA", 
	"CATGTCTTGA", "TCGAAGTACG", "CTCGAACCTT", "TCGATATCCT", "TCTGAACCGA", "ACGACCGACG", "GTGACAGCGC", "GCAAGCTGCC", "CCGCAGGAGC", "CAAGCAATTC", 
	"CAGTTCAGCT", "CACACTATAC", "CTACAATACG", "GCGGTACAAG", "CAACCGTATA", "GTTGATGGTC", "TAGAATGATT", "AGTCCAACTG", "GATCCACTCC", "GGCAGGAGCT", 
	"TGACTATCGA", "GCCGTACGGT", "TCGCAGCTGT", "GGAGCGACGA", "GTTCGATAAC", "GACTGCTCTC", "ACAGATCTAT", "GACGGACAGC", "CTCGTGCATG", "AGGTATAGAC", 
	"ACGTAGATAA", "CCAGAGTTGA", "CCAATAGGCA", "AGGCGCCATG", "AGAGAGGACA", "TGATTACCTT", "CGGACGTTCC", "GAGAGTCCGA", "AGCTTCCGGA", "CTCTATGTGT", 
	"TGTAGCGTAA", "TGGTAACCGC", "TTCGTCGCCA", "CCGCGGCACT", "CTACCGTCGA", "ACATACGTCC", "TTACGGTTAA", "CATACACTTG", "AGATCGCTGT", "GGATACTGAT", 
	"GGCTCCGCTG", "GCTCACCTCC", "AGAAGCGTAT", "TGGCAATTGC", "AAGGTACATA", "GATGTAGGTC", "CTTGTTGCTA", "GAACGTGAAG", "ATAAGTCATG", "ACGGCCAGCG", 
	"CCGAGATTAA", "TACTATCCTC", "GTCGTCACGG", "GGTCTGGTCA", "AGTGCTGTCT", "GCACCGCATA", "TATTCTTAAT", "CGGTCCGCTC", "GCTCGTTAGT", "TGAGGAGGAA", 
	"CCGGTCATGT", "CCAGAATTAG", "CCGGCTACGC", "TAATCCGTCG", "CTACTGATGG", "GGATGTGCTA", "AAGACGTACC", "TGAACCGGTG", "TATAGGCTGT", "GGCAGACCTA", 
	"TCCGCGGTCG", "CCGTAGCTCA", "CGTTGTAGCA", "CCAGATATCT", "TATATGTTAT", "CGAGCTGGAA", "CATACCAGTC", "AGTCTACGAA", "AATTCTAATT", "CATTCTGGTC", 
	"GAATCCATGA", "GCGTAAGTAG", "GATCGTCGAC", "GTGCTACTCG", "CCTAGTGTGG", "GGAAGCACTA", "AATATTATGG", "TCTAACACGT", "CAGGTATAAG", "AATAGACTCT", 
	"TACAGCATCA", "CTACTAATAA", "GCGTTAACGC", "GCGCCTAGCC", "CCAATCGGAC", "AGTCGGCTTG", "CACTGTGGAA", "GAGAGATGTA", "TGACGGCCGG", "GAATTCTCGG", 
	"TCAGTACACG", "TAGGATAGTT", "GTAATCTTCA", "AATCGATAGG", "AACACCTGAT", "GTTGGTACAG", "GTGAAGATCT", "GTCCATATGG", "CCGACTTGCA", "AAGTGCTGGT", 
	"TCGCTTCAAG", "ACTATGGCTG", "TATAGGCCTA", "GTGTAGTACT", "GCGGCGCACT", "CAGGCGTACT", "GCCGAATAGT", "CATAGGTCCA", "ACCGAGTTGC", "GTCAGTTCTG", 
	"AAGTAAGTTA", "GCATCCGTTG", "TGGCCAGTAA", "GGTACACCAC", "GGAGTTCGAT", "CTTCTATCGA", "CATTCACACA", "TATCTTGGTC", "TTCCATGTCA", "AAGGAACTAA", 
	"TAATCTAATA", "GAGACTATAG", "TCTACTGCAG", "GGCACGTTCG", "GCTACAATCA", "TTCGTCTAGC", "GGACGTGGAA", "TTGTCAGTGG", "GACGCTGCTA", "GCACATTGCG", 
	"GCCTGAGTAT", "ATCCTAATAC", "CGTCGTTCAG", "GGAGTTGACC", "GATTACGCGC", "GCTTGTCATA", "CTTGGCCAAT", "TCTCCATGAC", "CTCGAGAGCA", "GGATTCGGCC", 
	"GACAAGCCTA", "GGAATCAGTT", "AGACTGGAGC", "AGATGTATAT", "CCGAGACGTA", "GTCTTCGGCT"
};

#endif