#include <iostream>
#include "nt.h"



using namespace std;



Nt::Nt(string seq)
{

  length = seq.length();


  // forward strand
  for(int i = 0; i < length; i++) {

    switch(seq[i]) {

    case 'A': case 'a':
      forward.push_back(NT_A);
      break;

    case 'T': case 't':
      forward.push_back(NT_T);
      break;

    case 'G': case 'g':
      forward.push_back(NT_G);
      break;

    case 'C': case 'c':
      forward.push_back(NT_C);
      break;

    default:
      forward.push_back(NT_N);

    }

  }


  // reverse strand
  for(int i = length - 1; i >= 0; i--) {

    switch(seq[i]) {

    case 'A': case 'a':
      reverse.push_back(NT_T);
      break;

    case 'T': case 't':
      reverse.push_back(NT_A);
      break;

    case 'G': case 'g':
      reverse.push_back(NT_C);
      break;

    case 'C': case 'c':
      reverse.push_back(NT_G);
      break;

    default:
      reverse.push_back(NT_N);

    }

  }

}



vector<Codon> Nt::encode(int strand, int frame)
{
  vector<Codon> codons;
  vector<int> *s;


  if(strand == 1) {
    s = &forward;
  }
  else if(strand == -1) {
    s = &reverse;
  }
  else {
    cerr << "ERROR: strand should be +1 or -1" << endl;
    exit(1);
  }

  if(1 <= frame && frame <= 3)
    frame--;
  else {
    cerr << "ERROR: frame should be 1, 2 or 3" << endl;
    exit(1);
  }


  for(int i = frame; i <= s->size() - 3; i += 3) {
    int n0 = s->at(i);
    int n1 = s->at(i + 1);
    int n2 = s->at(i + 2);

    if(n0 == NT_N || n1 == NT_N || n2 == NT_N)
      codons.push_back(CODON_X);
    else
      codons.push_back(16 * n0 + 4 * n1 + n2);
  }


  return codons;
}
