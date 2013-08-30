//
//  main.cpp
//  harvest
//
//  Created by Brian Ondov on 6/19/13.
//
//

#include <iostream>
#include "HarvestIO.h"

using namespace::std;

int main(int argc, const char * argv[])
{
	const char * output = 0;
	const char * fasta = 0;
	const char * genbank = 0;
	const char * newick = 0;
	const char * vcf = 0;
	const char * xmfa = 0;
	
	for ( int i = 0; i < argc; i++ )
	{
		if ( argv[i][0] == '-' )
		{
			switch ( argv[i][1] )
			{
				case 'b': genbank = argv[++i]; break;
				case 'f': fasta = argv[++i]; break;
				case 'n': newick = argv[++i]; break;
				case 'o': output = argv[++i]; break;
				case 'v': vcf = argv[++i]; break;
				case 'x': xmfa = argv[++i]; break;
			}
		}
	}
	
	cout << "Output:" << output << '\n';
	cout << "Fasta:" << fasta << '\n';
	cout << "GenBank:" << genbank << "\n";
	cout << "Newick:" << newick << "\n";
	cout << "VCF:" << vcf << "\n";
	cout << "XMFA:" << xmfa << "\n";
	
	HarvestIO hio;
	
	hio.loadFasta(fasta);
	
	if ( genbank )
	{
		hio.loadGenbank(genbank);
	}
	
	hio.loadXmfa(xmfa);
	hio.loadNewick(newick);
	hio.loadVcf(vcf);
	
	hio.writeHarvest(output);
	
    return 0;
}

