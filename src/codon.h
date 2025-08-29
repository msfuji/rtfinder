#ifndef _CODON_H_
#define _CODON_H_


#include <vector>
#include <ostream>
#include <fstream>


#define N_CODON  64
#define CODON_X  -1


// start codons
#define CODON_ATG     6  // ATG
#define CODON_TTG    22  // TTG
#define CODON_GTG    38  // GTG
#define CODON_CTG    54  // TTG

// stop codons
#define CODON_AMBER  18  // TAG
#define CODON_OCHRE  16  // TAA
#define CODON_OPAL   24  // TGA


typedef int Codon;


class CodonFreq
{

 public:
  CodonFreq();
  void add(std::vector<Codon> &codons);
  void sub(CodonFreq &cf);
  void print_lod(std::ostream &of, CodonFreq &cf, double alpha);

 private:
  int total;
  int visit[N_CODON];
  int trans[N_CODON][N_CODON];
  double prob(int cij, int ci, int cj, int c_tot, double alpha);

};



void loadLodScores(int lod[][N_CODON], std::ifstream &f_lod);
std::string toCodonStr(int i);


#endif // _CODON_H_
