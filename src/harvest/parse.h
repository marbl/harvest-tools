// Copyright © 2014, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef parse_h
#define parse_h

#include <string>

char complement(char base);
char * removePrefix(char * string, const char * substring);
void reverseComplement(std::string & sequence);
void ungap(std::string & gapped);

#endif
