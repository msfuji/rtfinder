#ifndef _ORF_H_
#define _ORF_H_


#include <vector>
#include <ostream>
#include <string>

#include "codon.h"



class Orf
{

 public:
  Orf(std::string &ent_id, int strand, int offset,
      std::vector<Codon> &genome_codons,
      int start_pos, int rt_pos, int stop_pos);

  void evaluate(int lod[][N_CODON]);
  void evaluate(void);
  void print(std::ostream &of, bool use_lod);
  int score;

 private:
  std::string entry_id;
  int from;
  int to;
  int rt_idx;
  std::vector<Codon> codons;

};



std::vector<Orf> findOrfs(std::vector<Codon> &codons, std::string entry_id,
			  int offset, int min_total_len,
			  bool eukaryote, bool mycoplasma);




#endif // _ORF_H_
