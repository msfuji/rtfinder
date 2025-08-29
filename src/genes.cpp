#include "genes.h"


using namespace std;


bool
getNtseqFromGenes(string &seq, bool &is_CDS, ifstream &f)
{
  string line;
  string::size_type pos;

  
  getline(f, line);
  if(line.compare(0, 5, "ENTRY") != 0)
    return false;
  if(line.find(" CDS ") != string::npos)
    is_CDS = true;
  else
    is_CDS = false;

  
  while(getline(f, line)) {
    if(line.compare(0, 5, "NTSEQ") == 0) break;
  }
  seq = "";
  while(getline(f, line)) {
    if(line == "///") break;
    pos = line.find_first_not_of(" ");
    seq += line.substr(pos, line.length() - pos);
  }

  
  // remove stop codon
  if(seq.length() >= 3) 
    seq.erase(seq.length() - 3);


  return true;
}
