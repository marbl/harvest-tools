// Copyright © 2014, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "AnnotationList.h"
#include <fstream>
#include "parse.h"

using namespace std;

void AnnotationList::clear()
{
	annotations.clear();
}

void AnnotationList::initFromGenbank(const char * file, ReferenceList & referenceList, bool useSeq)
{
	ifstream in(file);
	char * line = new char[1 << 20];
	Annotation * annotation = 0;
	int offset;
	
	while ( ! in.eof() )
	{
		string locus;
		string definition;
		
		if ( useSeq )
		{
			offset = 0;
			
			for ( int i = 0; i < referenceList.getReferenceCount(); i++ )
			{
				offset += referenceList.getReference(i).sequence.length();
			}
		}
		
		// header
		
		while ( in.getline(line, (1 << 20) - 1) )
		{
			char * token;
			
			if ( useSeq && (token = removePrefix(line, "LOCUS")) )
			{
				locus = strtok(token, " \t");
			}
			else if ( useSeq && (token = removePrefix(line, "DEFINITION")) )
			{
				do
				{
					while ( *token == ' ' )
					{
						token++;
					}
				
					if ( definition.length() )
					{
						definition.append(" ");
					}
				
					definition.append(token);
				
					in.getline(token, (1 << 20) - 1);
				}
				while ( *token == ' ' );
			}
			else if ( ! useSeq && removePrefix(line, "VERSION") )
			{
				if ( const char * giToken = strstr(line, " GI:") )
				{
					long int gi = atol(giToken + 4);
					int sequence = referenceList.getReferenceSequenceFromGi(gi);
					
					offset = 0;
				
					for ( int i = 0; i < sequence; i++ )
					{
						offset += referenceList.getReference(i).sequence.length();
					}
				}
				else
				{
					throw NoGiException(file);
				}
			}
			else if ( removePrefix(line, "FEATURES") )
			{
				break;
			}
		}
	
		// annotations
	
		while ( ! in.eof() )
		{
			in.getline(line, (1 << 20) - 1);
		
			if ( in.eof() || strcmp(line, "//") == 0 || removePrefix(line, "ORIGIN") )
			{
				break;
			}
		
			char * token = line;
			char * suffix;
		
			while ( *token == ' ' )
			{
				token++;
			}
		
			if ( token == line + 5 )
			{
				int start;
				int end;
				bool reverse;
			
				string feature = strtok(token, " ");
			
				token = strtok(0, " ");
			
				suffix = removePrefix(token, "complement(");
				reverse = suffix;
			
				if ( suffix )
				{
					token = suffix;
				}
			
				suffix = removePrefix(token, "join(");
			
				if ( suffix )
				{
					token = suffix;
				}
			
				suffix = removePrefix(token, "order(");
			
				if ( suffix )
				{
					token = suffix;
				}
			
				suffix = removePrefix(token, "complement(");
			
				if ( suffix )
				{
					reverse = true;
					token = suffix;
				}
			
				if ( *token == '<' || *token == '>' )
				{
					token++;
				}
			
				start = atoi(strtok(token, ".")) + offset - 1;
				end = atoi(strtok(0, ".,)<>")) + offset - 1;
			
				if ( ! annotation || start != annotation->start || end != annotation->end || reverse != annotation->reverse )
				{
					if ( feature != "source" && feature != "misc_feature" )
					{
						annotations.resize(annotations.size() + 1);
						annotation = &annotations.at(annotations.size() - 1);
					
						annotation->start = start;
						annotation->end = end;
						annotation->reverse = reverse;
					}
					else
					{
						annotation = 0;
					}
				}
			
				if ( annotation )
				{
					annotation->feature = feature;
				}
			}
			else if ( annotation )
			{
				char * suffix;
		
				if ( (suffix = removePrefix(token, "/locus_tag=\"")) )
				{
					annotation->locus = strtok(suffix, "\"");
				}
				else if ( (suffix = removePrefix(token, "/gene=\"")) && annotation->feature == "gene" )
				{
					annotation->name = strtok(suffix, "\"");
				}
				else if ( (suffix = removePrefix(token, "/product=\"")) )
				{
					annotation->description = suffix;
			
					if ( annotation->description[annotation->description.length() - 1] == '"' )
					{
						annotation->description.resize(annotation->description.length() - 1);
					}
				
					while ( suffix[strlen(suffix) - 1] != '"' )
					{
						in.getline(line, (1 << 20) - 1);
						suffix = line;
				
						while ( *suffix == ' ' )
						{
							suffix++;
						}
					
						annotation->description.append(suffix - 1, strlen(suffix));
					}
				}
			}
		}
	
		// sequence
	
		string sequence;
	
		while ( ! in.eof() && strcmp(line, "//") != 0 )
		{
			in.getline(line, (1 << 20) - 1);
			
			if ( useSeq )
			{
				if ( in.eof() || strcmp(line, "//") == 0 )
				{
					break;
				}
		
				strtok(line, " "); // eat number
				const char * token;
		
				while ( (token = strtok(0, " ")) )
				{
					sequence.append(token);
				}
			}
		}
		
		if ( in.eof() )
		{
			break;
		}
		
		if ( useSeq )
		{
			if ( sequence.length() )
			{
				for ( int i = 0; i < sequence.length(); i++ )
				{
					sequence[i] = toupper(sequence.at(i));
				}
			
				referenceList.addReference(locus, definition, sequence);
			}
			else
			{
				clear();
				throw NoSequenceException(file);
			}
		}
	}
	
	for ( int i = 0; false && i < annotations.size(); i++ )
	{
		annotation = &annotations.at(i);
		printf("%s\t%d\t%d\t%c\t%s\t%s\n", annotation->locus.c_str(), annotation->start, annotation->end, annotation->reverse ? '-' : '+', annotation->name.c_str(), annotation->description.c_str());
	}
	
	delete [] line;
	in.close();
}

void AnnotationList::initFromProtocolBuffer(const Harvest::AnnotationList & msg, const ReferenceList & referenceList)
{
	int sequence = 0;
	int offset = 0;
	
	annotations.resize(0);
	annotations.resize(msg.annotations_size());
	
	for ( int i = 0; i < msg.annotations_size(); i++ )
	{
		Annotation & annotation = annotations.at(i);
		
		const Harvest::AnnotationList::Annotation & msgAnn = msg.annotations(i);
		
		while ( sequence < msgAnn.sequence() && sequence < referenceList.getReferenceCount() )
		{
			offset += referenceList.getReference(sequence).sequence.length();
			sequence++;
		}
		
		if ( sequence == referenceList.getReferenceCount() )
		{
			offset = 0;
			//printf("ERROR: sequence %d not found in reference or annotation out of order in protobuf.\n", msgAnn.sequence());
			//exit(1);
		}
		
		annotation.start = msgAnn.regions(0).start() + offset;
		annotation.end = msgAnn.regions(0).end() + offset;
		annotation.reverse = msgAnn.reverse();
		annotation.name = msgAnn.name();
		annotation.locus = msgAnn.locus();
		annotation.description = msgAnn.description();
		annotation.feature = msgAnn.feature();
		
		//printf("%s\t%d\t%d\t%d\t%c\t%s\t%s\n", annotation.locus.c_str(), msgAnn.sequence(), annotation.start, annotation.end, annotation.reverse ? '-' : '+', annotation.name.c_str(), annotation.description.c_str());
	}
}

void AnnotationList::writeToProtocolBuffer(Harvest * harvest, const ReferenceList & referenceList) const
{
	int sequence = 0;
	int offset = 0;
	int offsetEnd = referenceList.getReference(0).sequence.length();
	
	for ( int i = 0; i < annotations.size(); i++ )
	{
		const Annotation & annotation = annotations.at(i);
		
		Harvest::AnnotationList::Annotation * msgAnn = harvest->mutable_annotations()->add_annotations();
		
		while ( annotation.start >= offsetEnd && sequence < referenceList.getReferenceCount() - 1 )
		{
			sequence++;
			offset = offsetEnd;
			offsetEnd += referenceList.getReference(sequence).sequence.length();
		}
		
		msgAnn->set_sequence(sequence);
	
		msgAnn->add_regions();
		msgAnn->mutable_regions(0)->set_start(annotation.start - offset);
		msgAnn->mutable_regions(0)->set_end(annotation.end - offset);
		msgAnn->set_reverse(annotation.reverse);
		msgAnn->set_name(annotation.name);
		msgAnn->set_locus(annotation.locus);
		msgAnn->set_description(annotation.description);
		msgAnn->set_feature(annotation.feature);
	}
}
