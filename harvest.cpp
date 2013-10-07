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
	const char * mfa = 0;
	const char * newick = 0;
	const char * vcf = 0;
	const char * xmfa = 0;
	const char * outSnp = 0;
	
	for ( int i = 0; i < argc; i++ )
	{
		if ( argv[i][0] == '-' )
		{
			switch ( argv[i][1] )
			{
				case 'b': genbank = argv[++i]; break;
				case 'f': fasta = argv[++i]; break;
				case 'm': mfa = argv[++i]; break;
				case 'n': newick = argv[++i]; break;
				case 'o': output = argv[++i]; break;
				case 'S': outSnp = argv[++i]; break;
				case 'v': vcf = argv[++i]; break;
				case 'x': xmfa = argv[++i]; break;
			}
		}
	}
	
	HarvestIO hio;
	
	if ( mfa )
	{
		hio.loadMFA(mfa);
		hio.loadNewick(newick);
		hio.writeHarvest(output);
		return 0;
	}
	
	if ( fasta )
	{
		hio.loadFasta(fasta);
	}
	
	if ( genbank )
	{
		hio.loadGenbank(genbank);
	}
	
	if ( xmfa )
	{
		printf("Loading %s...\n", xmfa);
		hio.loadXmfa(xmfa, vcf == 0);
	}
	
	if ( newick )
	{
		hio.loadNewick(newick);
	}
	
	if ( vcf )
	{
		hio.loadVcf(vcf);
	}
	
	if ( output )
	{
		hio.writeHarvest(output);
	}
	
	if ( outSnp )
	{
		printf("Writing %s...\n", outSnp);
		hio.writeSnp(outSnp);
	}
	
    return 0;
}

