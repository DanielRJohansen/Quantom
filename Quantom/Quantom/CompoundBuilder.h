#pragma once


#include "Bodies.cuh"
#include <string.h>
#include <fstream>



using namespace std;


enum LineType { atom, pairbond, anglebond, torsionbond };

struct ParsedLine {
	LineType linetype;
	int pdb_index;
	string id;
	string monomer;
	Float3 position;
	string atom;
};





using namespace std;
class CompoundBuilder
{
public:
	CompoundBuilder() {}
	Molecule1 buildMolecule(string pdb_path);

private:
	ParsedLine parseLine(int line_index);
	ParsedLine parseAtom(string line);
	ParsedLine parseConnection(string line);
	int countMonomers(string path);
};

