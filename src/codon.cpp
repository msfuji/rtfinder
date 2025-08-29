#include <iostream>
#include <cmath>
#include <string>
#include <sstream>

#include "codon.h"
#include "nt.h"



using namespace std;



#define SIGNATURE "#RTTRAIN"



static string codonStr[N_CODON] = {
  "AAA", "AAT", "AAG", "AAC", 
  "ATA", "ATT", "ATG", "ATC", 
  "AGA", "AGT", "AGG", "AGC", 
  "ACA", "ACT", "ACG", "ACC", 
  "TAA", "TAT", "TAG", "TAC", 
  "TTA", "TTT", "TTG", "TTC", 
  "TGA", "TGT", "TGG", "TGC", 
  "TCA", "TCT", "TCG", "TCC", 
  "GAA", "GAT", "GAG", "GAC", 
  "GTA", "GTT", "GTG", "GTC", 
  "GGA", "GGT", "GGG", "GGC", 
  "GCA", "GCT", "GCG", "GCC", 
  "CAA", "CAT", "CAG", "CAC", 
  "CTA", "CTT", "CTG", "CTC", 
  "CGA", "CGT", "CGG", "CGC", 
  "CCA", "CCT", "CCG", "CCC"
};



CodonFreq::CodonFreq()
{
  total = 0;
  for(int i = 0; i < N_CODON; i++) {
    visit[i] = 0;
    for(int j = 0; j < N_CODON; j++)
      trans[i][j] = 0;
  }
}



void
CodonFreq::add(std::vector<Codon> &codons)
{
  total += codons.size();

  visit[codons[0]]++;
    
  for(int i = 1; i < codons.size(); i++) {
    visit[codons[i]]++;
    trans[codons[i - 1]][codons[i]]++;
  }

}



void
CodonFreq::sub(CodonFreq &cf)
{
  total -= cf.total;

  for(int i = 0; i < N_CODON; i++){
    visit[i] -= cf.visit[i];
    for(int j = 0; j < N_CODON; j++){
      trans[i][j] -= cf.trans[i][j];
    }
  }

}



void
CodonFreq::print_lod(ostream &of, CodonFreq &cf, double alpha)
{
  string str1, str2;
  double p1, p2;
  double lod;

  
  of << SIGNATURE << endl;

  for(int i = 0; i < N_CODON; i++){
    str1 = toCodonStr(i);

    for(int j = 0; j < N_CODON; j++){
      str2 = toCodonStr(j);

      p1 = prob(trans[i][j], visit[i], visit[j], total, alpha);
      p2 = prob(cf.trans[i][j], cf.visit[i], cf.visit[j], cf.total, alpha);
      lod = 1000 * (log(p1) - log(p2)) / log(2.0);

      of << str1 << "->" << str2 << "\t";
      of << i << "\t" << j << "\t";
      of << (int)lod << endl;
    }
  }

}



double
CodonFreq::prob(int cij, int ci, int cj, int c_tot, double alpha)
{
  double x = cij + (alpha * cj) / c_tot + 1;
  double y = ci + alpha + N_CODON;

  return x / y;
}



string
toCodonStr(int i)
{
  if(0 <= i and i < N_CODON)
    return codonStr[i];
  else
    return "NNN";
}



void
loadLodScores(int lod[][N_CODON], ifstream &f_lod)
{
  string line;
  string transition;
  int i, j, score;
  
  // check header
  getline(f_lod, line);
  if(line != SIGNATURE) {
    cerr << "ERROR: illegal file format" << endl;
    exit(1);
  }

  while(getline(f_lod, line)) {
    stringstream ss(line);
    ss >> transition >> i >> j >> score;
    lod[i][j] = score;
  }

}
