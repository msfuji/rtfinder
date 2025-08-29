#include <ostream>

#include "orf.h"


using namespace std;



Orf::Orf(string &ent_id, int strand, int offset,
	 vector<Codon> &genome_codons,
	 int start_pos, int rt_pos, int stop_pos)
{
  
  entry_id = ent_id;
  stop_pos--;

  if(strand == 1) {
    from = 3 * start_pos + offset;
    to   = 3 * stop_pos + offset + 2;
  }
  else {
    from = offset - 3 * start_pos;
    to   = offset - 3 * stop_pos - 2;
  }

  rt_idx = rt_pos - start_pos;

  for(int i = start_pos; i <= stop_pos; i++)
    codons.push_back(genome_codons[i]);

}



void
Orf::evaluate(int lod[][N_CODON])
{

  int score_N = 0;
  for(int i = 1; i < rt_idx; i++) {
    if(codons[i - 1] != CODON_X and codons[i] != CODON_X)
      score_N += lod[codons[i - 1]][codons[i]];
  }

  int score_C = 0;
  for(int i = rt_idx + 2; i < codons.size(); i++) {
    if(codons[i - 1] != CODON_X and codons[i] != CODON_X)
      score_C += lod[codons[i - 1]][codons[i]];
  }
  
  score = (score_N < score_C) ? score_N : score_C;
  
}


void
Orf::evaluate(void)
{

  int score_N = rt_idx;
  int score_C = codons.size() - rt_idx - 1;
  score = (score_N < score_C) ? score_N : score_C;
  
}



void
Orf::print(ostream &of, bool use_lod)
{

  of << ">" << entry_id << ":" << from << ".." << to << " ";
  of << "stop_codon=" << rt_idx + 1 << "";
  of << "(" << toCodonStr(codons[rt_idx]) << ")";
  if(use_lod)
    of << " score=" << score;
  of << endl;
  
  vector<Codon>::iterator it = codons.begin();
  while(it != codons.end()) {
    of << toCodonStr(*it);
    it++;
  }

  of << endl;
    
}



vector<Orf>
findOrfs(vector<Codon> &codons, string entry_id, int offset,
	 int min_total_len, bool eukaryote, bool mycoplasma)
{

  int strand;
  if(1 <= offset and offset <= 3)
    strand = 1;
  else
    strand = -1;


  // scan stop positions
  vector<int> stop_pos;
  stop_pos.push_back(-1);
  for(int i = 0; i < codons.size(); i++) {
    int c = codons[i];
    if(c == CODON_AMBER or c == CODON_OCHRE or (c == CODON_OPAL and !mycoplasma))
      stop_pos.push_back(i);
  }
  stop_pos.push_back(codons.size());


  // find long ORFs
  vector <Orf> orfs;
  for(int i = 1; i <= stop_pos.size() - 2; i++) {
    int len_tot = stop_pos[i + 1] - stop_pos[i - 1] - 1;
    if(len_tot < min_total_len) continue;

    int start_pos;
    for(start_pos = stop_pos[i - 1] + 1; start_pos <= stop_pos[i]; start_pos++) {
      int c = codons[start_pos];
      if(eukaryote) {
	if(c == CODON_ATG) break;
      }
      else {
	if(c == CODON_ATG or c == CODON_TTG or c == CODON_GTG) break;	
      }	  
    }
    if(start_pos >= stop_pos[i]) continue;
    len_tot = stop_pos[i + 1] - start_pos;
    if(len_tot < min_total_len) continue;

    Orf orf(entry_id, strand, offset, codons,
	    start_pos, stop_pos[i], stop_pos[i + 1]);
    orfs.push_back(orf);
  }


  return orfs;
}


