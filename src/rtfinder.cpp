#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cstring>

#include "fasta.h"
#include "nt.h"
#include "codon.h"
#include "orf.h"



using namespace std;



#define MIN_TOTAL_LEN    100



void
usage()
{
  cerr << "Usage: rtfinder [-options] <FASTA>";
  cerr << endl << endl;
  cerr << "Options:" << endl;;
  cerr << "  -t     : only scan the top strand" << endl;
  cerr << "  -e     : use eukaryotic start codons" << endl;
  cerr << "  -m     : use Mycoplasma stop codons" << endl;
  cerr << "  -l <n> : minimal ORF length; default " << MIN_TOTAL_LEN << endl;
  cerr << "  -c <f> : score ORF by codon transition model <f>" << endl;
  cerr << "  -s <n> : minimal length of sense codon run (or min. score)";
  cerr << endl << endl;
  exit(1);
};



int
main(int argc, char **argv)
{

  bool both_strand = true;
  bool eukaryote = false;
  bool mycoplasma = false;
  int min_total_len = MIN_TOTAL_LEN;
  bool use_lod = false;
  char *lod_file = 0;
  int use_threshold = false;
  int threshold;
  char *genome_file = 0;

  
  // parse command line options
  argc--; argv++;
  while(argc > 0) {
    if(strcmp(argv[0], "-t") == 0) {
      both_strand = false;
      argc--; argv++;
    }
    else if(strcmp(argv[0], "-e") == 0) {
      eukaryote = true;
      argc--; argv++;
    }
    else if(strcmp(argv[0], "-m") == 0) {
      mycoplasma = true;
      argc--; argv++;
    }
    else if(strcmp(argv[0], "-l") == 0) {
      if(argc == 1) usage();
      min_total_len = atoi(argv[1]);
      argc -= 2; argv += 2;
    }
    else if(strcmp(argv[0], "-c") == 0) {
      if(argc == 1) usage();
      use_lod = true;
      lod_file = new char[strlen(argv[1]) + 1];
      strcpy(lod_file, argv[1]);
      argc -= 2; argv += 2;
    }
    else if(strcmp(argv[0], "-s") == 0) {
      if(argc == 1) usage();
      use_threshold = true;
      threshold = atoi(argv[1]);
      argc -= 2; argv += 2;
    }
    else if(argv[0][0] == '-')
      usage();
    else
      break;
  }
  if(argc != 1) usage();
  genome_file = new char[strlen(argv[0]) + 1];
  strcpy(genome_file, argv[0]);
  

  
  // read log-odds (if required)
  int lod[N_CODON][N_CODON];
  if(use_lod) {
    ifstream f_lod(lod_file);
    loadLodScores(lod, f_lod);
  }



  // read genome
  ifstream f_genome(genome_file);
  Fasta fas;

  while(getFasta(fas, f_genome)) {
    Nt nt(fas.seq);
    string entry_id = parse_entry_id(fas.definition);

    for(int strand = (both_strand ? -1 : 1); strand <= 1; strand += 2) {
      for(int frame = 1; frame <= 3; frame++) {
	vector<Codon> codons = nt.encode(strand, frame);

	// find ORFs
	int offset = (strand == 1) ? frame : (nt.length - frame + 1);
	vector<Orf> orfs = findOrfs(codons, entry_id, offset,
				    min_total_len, eukaryote, mycoplasma);

	// print ORFs
	vector<Orf>::iterator it = orfs.begin();
	while(it != orfs.end()){
	  if(use_lod)
	    it->evaluate(lod);
	  else
	    it->evaluate();

	  if(!use_threshold or it->score > threshold)
	      it->print(cout, use_lod);
	  it++;
	}

      }
    }
  }


  if(lod_file) delete [] genome_file;
  if(genome_file) delete [] lod_file;

  
  return 0;
}

