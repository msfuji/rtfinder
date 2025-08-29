#ifndef _FASTA_H_
#define _FASTA_H_


#include <fstream>
#include <string>


class Fasta
{
 public:
  std::string definition;
  std::string seq;
};



bool getFasta(Fasta &fas, std::ifstream &f);
std::string parse_entry_id(std::string definition);



#endif // _FASTA_H_
