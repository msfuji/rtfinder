#include "fasta.h"
#include <iostream>



using namespace std;



bool
getFasta(Fasta &fas, ifstream &f)
{
  string line;

  getline(f, line);
  if(line[0] != '>')
    return false;
  fas.definition = line.substr(1, line.length() - 1);


  fas.seq = "";
  while(getline(f, line)) {
    if(line[0] == '>') {
      f.seekg(-line.length() - 1, ios::cur);
      break;
    }
    fas.seq += line;
  }
  
  return true;
}



string
parse_entry_id(string definition)
{
  if(definition.compare(0, 3, "gn:") == 0) 
    definition = definition.substr(3, definition.length() - 3);

  string::size_type to = definition.find_first_of(" \t\r\n");
  return definition.substr(0, to);
}


