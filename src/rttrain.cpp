#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>

#include "genes.h"
#include "fasta.h"
#include "nt.h"
#include "codon.h"



using namespace std;



#define ALPHA  1000



static void
usage()
{
  cerr << "Usage: rttrain [-options] <KEGG GENES> <KEGG GENOME>";
  cerr << endl << endl;
  cerr << "Options:" << endl;;
  cerr << "  -w <x> : total amount of pseudocount; default " << ALPHA;
  cerr << endl << endl;
  exit(1);
};



static void
failedToOpen(char *filename)
{
  cerr << "ERROR: failed to open '" << filename << "'" << endl;
  exit(1);
}



int
main(int argc, char **argv)
{

  vector<Codon> codons;
  double alpha = ALPHA;


  if(argc >= 3 and strcmp(argv[1], "-w") == 0) {
    alpha = atof(argv[2]);
    argv += 2;
    argc -= 2;
  }
  if(argc != 3) usage();


  // get NTSEQ of CDS from GENES
  ifstream f_genes(argv[1]);
  string seq;
  bool is_CDS;
  CodonFreq cf_cds;

  if(f_genes.fail()) failedToOpen(argv[1]);

  while(getNtseqFromGenes(seq, is_CDS, f_genes)) {
    if(!is_CDS) continue;
    Nt nt(seq);
    codons = nt.encode(1, 1);
    cf_cds.add(codons);
  }


  // get genome sequences from FASTA files
  ifstream f_genome(argv[2]);
  Fasta fas;
  CodonFreq cf_genome;

  if(f_genome.fail()) failedToOpen(argv[2]);

  while(getFasta(fas, f_genome)) {
    Nt nt(fas.seq);
    for(int strand = -1; strand <= 1; strand += 2) {
      for(int frame = 1; frame <= 3; frame++) {
	codons = nt.encode(strand, frame);
	cf_genome.add(codons);
      }
    }
  }

  cf_genome.sub(cf_cds);

  

  // print log-odds
  cf_cds.print_lod(cout, cf_genome, alpha);



  return 0;
}



