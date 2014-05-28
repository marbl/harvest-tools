//
//  main.cpp
//  harvest
//
//  Created by Brian Ondov on 6/19/13.
//
//

#include <iostream>
#include <fstream>
#include "HarvestIO.h"
#include <string.h>

using namespace::std;

int main(int argc, const char * argv[])
{
	const char * input = 0;
	const char * output = 0;
	vector<const char *> bed;
	const char * fasta = 0;
	vector<const char *> genbank;
	const char * mfa = 0;
	const char * newick = 0;
	const char * vcf = 0;
	const char * xmfa = 0;
	const char * outNewick = 0;
	const char * outSnp = 0;
	const char * outVcf = 0;
	const char * outBB = 0;
	const char * outXmfa = 0;
	bool help = false;
	bool quiet = false;
	
	//stdout flag
	string out1("-");
	
	for ( int i = 0; i < argc; i++ )
	{
		if ( argv[i][0] == '-' )
		{
			switch ( argv[i][1] )
			{
				case 'b': bed.push_back(argv[++i]); break;
				case 'B': outBB = argv[++i]; break;
				case 'f': fasta = argv[++i]; break;
				case 'g': genbank.push_back(argv[++i]); break;
				case 'h': help = true; break;
				case 'i': input = argv[++i]; break;
				case 'm': mfa = argv[++i]; break;
				case 'n': newick = argv[++i]; break;
				case 'N': outNewick = argv[++i]; break;
				case 'o': output = argv[++i]; break;
				case 'q': quiet = true; break;
				case 'S': outSnp = argv[++i]; break;
				case 'V': outVcf = argv[++i]; break;
				case 'v': vcf = argv[++i]; break;
				case 'x': xmfa = argv[++i]; break;
				case 'X': outXmfa = argv[++i]; break;
				
				default:
					printf("ERROR: Unrecognized option ('%s').\n", argv[i]);
					help = true;
			}
		}
	}
	
	if (help || argc < 2)
	{
		cout << "harVest usage: harvest " << endl;
		cout << "   -i <harvest input>" << endl;
		cout << "   -b <bed filter intervals>,<filter name>,\"<description>\"" << endl;
		cout << "   -B <output backbone intervals>" << endl;
		cout << "   -f <reference fasta>" << endl;
		cout << "   -g <reference genbank>" << endl;
		cout << "   -m <multi-fasta alignment input>" << endl;
		cout << "   -n <Newick tree input>" << endl;
		cout << "   -N <Newick tree output>" << endl;
		cout << "   -o <hvt output>" << endl;
		cout << "   -S <output for multi-fasta SNPs>" << endl;
		cout << "   -v <VCF input>" << endl;
		cout << "   -V <VCF output>" << endl;
		cout << "   -x <xmfa alignment file>" << endl;
		cout << "   -X <output xmfa alignment file>" << endl;
		cout << "   -h (show this help)" << endl;
		cout << "   -q (quiet mode)" << endl;
		exit(0);
	}
	
	HarvestIO hio;
	
	if ( input )
	{
		hio.loadHarvest(input);
	}
	
	if ( mfa )
	{
		hio.loadMFA(mfa);
	}
	
	if ( fasta )
	{
		hio.loadFasta(fasta);
	}
	
	for ( int i = 0; i < genbank.size(); i++ )
	{
		hio.loadGenbank(genbank[i]);
	}
	
	if ( xmfa )
	{
		if (!quiet)
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
	
	for ( int i = 0; i < bed.size(); i++ )
	{
		char * arg = new char[strlen(bed[i]) + 1];
		
		strcpy(arg, bed[i]);
		
		const char * file = strtok(arg, ",");
		const char * name = strtok(0, ",");
		const char * desc = strtok(0, "");
		
		if ( name == 0 )
		{
			printf("ERROR: no filter name for bed file %s\n", file);
			return 1;
		}
		
		if ( desc == 0 )
		{
			printf("ERROR: no filter description for bed file %s\n", file);
			return 1;
		}
		
		hio.loadBed(file, name, desc);
		delete [] arg;
	}
	
	if ( output )
	{
		hio.writeHarvest(output);
	}
	
	if ( outNewick )
	{
		if (!quiet) printf("Writing %s...\n", outNewick);
		
		std::ostream* fp = &cout;
		std::ofstream fout;
		
		if (out1.compare(outNewick) != 0) 
		{
			fout.open(outNewick);
			fp = &fout;
		}
		
		hio.writeNewick(*fp);
	}
	
	if ( outSnp )
	{
		if (!quiet) printf("Writing %s...\n", outSnp);
		
		std::ostream* fp = &cout;
		std::ofstream fout;
		
		if (out1.compare(outSnp) != 0) 
		{
			fout.open(outSnp);
			fp = &fout;
		}
		
		hio.writeSnp(*fp, true);
	}

	if ( outBB )
	{
		if (!quiet) printf("Writing %s...\n", outBB);
		
		std::ostream* fp = &cout;
		std::ofstream fout;
		
		if (out1.compare(outBB) != 0) 
		{
			fout.open(outBB);
			fp = &fout;
		}
		
		hio.writeBackbone(*fp);
	}
	if ( outXmfa )
	{
		if (!quiet) printf("Writing %s...\n", outXmfa);
		
		std::ostream* fp = &cout;
		std::ofstream fout;
		
		if (out1.compare(outXmfa) != 0) 
		{
			fout.open(outXmfa);
			fp = &fout;
		}
		
		hio.writeXmfa(*fp);
	}

	if ( outVcf )
	{
		if (!quiet) printf("Writing %s...\n", outVcf);
		
		std::ostream* fp = &cout;
		std::ofstream fout;
		
		if (out1.compare(outVcf) != 0) 
		{
			fout.open(outVcf);
			fp = &fout;
		}
		
		hio.writeVcf(*fp, true);

	}
	
    return 0;
}

