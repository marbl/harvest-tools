// Copyright © 2014, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "parse.h"
#include <string.h>

using namespace::std;

char complement(char base)
{
	switch ( base )
	{
		case 'A': return 'T'; break;
		case 'C': return 'G'; break;
		case 'G': return 'C'; break;
		case 'T': return 'A'; break;
		default: return base;
	}
}

char * removePrefix(char * string, const char * substring)
{
	size_t len = strlen(substring);
	
	if ( strncmp(string, substring, len) == 0 )
	{
		return string + len;
	}
	else
	{
		return 0;
	}
}

void reverseComplement(string & sequence)
{
	string copy = sequence;
	
	for ( int i = 0; i < sequence.length(); i++ )
	{
		sequence[sequence.length() - i - 1] = complement(copy[i]);
	}
}

void ungap(string & gapped)
{
	int pos = 0;
	
	for ( int i = 0; i < gapped.length(); i++ )
	{
		if ( gapped[i] != '-' )
		{
			gapped[pos] = gapped[i];
			pos++;
		}
	}
	
	gapped.resize(pos);
}
