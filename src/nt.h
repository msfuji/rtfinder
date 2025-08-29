#ifndef _NT_H_
#define _NT_H_


#include <string>
#include <vector>

#include "codon.h"


#define NT_A   0
#define NT_T   1
#define NT_G   2
#define NT_C   3
#define NT_N  -1


class Nt
{
 public:
  Nt(std::string seq);
  std::vector<Codon> encode(int strand, int frame);
  int length;

 private:
  std::vector<int> forward;
  std::vector<int> reverse;
};


#endif // _NT_H_
