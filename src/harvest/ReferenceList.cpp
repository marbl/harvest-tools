// Copyright © 2014, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include <fstream>
#include "ReferenceList.h"

using namespace::std;

void ReferenceList::addReference(string name, string desc, string sequence)
{
	references.resize(references.size() + 1);
	references[references.size() - 1].name = name;
	references[references.size() - 1].description = desc;
	references[references.size() - 1].sequence = sequence;
}

long int ReferenceList::getConcatenatedPosition(int sequence, long int position) const
{
	for ( int i = 0; i < sequence; i++ )
	{
		position += references.at(i).sequence.length();
	}
	
	return position;
}

int ReferenceList::getPositionFromConcatenated(int sequence, long int position) const
{
	long int sum = 0;
	
	for ( int i = 0; i < sequence; i++ )
	{
		sum += references.at(i).sequence.length();
	}
	
	return position - sum;
}

int ReferenceList::getReferenceSequenceFromConcatenated(long int position) const
{
	long int sum = 0;
	
	for ( int i = 0; i < references.size(); i++ )
	{
		sum += references.at(i).sequence.length();
		
		if ( sum > position )
		{
			return i;
		}
	}
	
	return undef;
}

int ReferenceList::getReferenceSequenceFromGi(long int gi) const
{
	for ( int i = 0; i < references.size(); i++ )
	{
		size_t giToken = references.at(i).name.find("gi|");
		
		if ( giToken != string::npos )
		{
			if ( gi == atol(references.at(i).name.c_str() + giToken + 3))
			{
				return i;
			}
		}
	}
	
	char giString[1024];
	sprintf(giString, "%ld", gi);
	
	throw GiNotFoundException(giString);
	
	return undef;
}

int ReferenceList::getReferenceSequenceFromName(string name) const
{
	for ( int i = 0; i < references.size(); i++ )
	{
		if ( name == references.at(i).name )
		{
			return i;
		}
	}
	
	throw NameNotFoundException(name);
	
	return undef;
}

void ReferenceList::initFromFasta(const char * file)
{
	ifstream in(file);
	string line;
	Reference * reference;
	
	while ( getline(in, line) )
	{
		if ( line[0] == '>' )
		{
			references.resize(references.size() + 1);
			reference = &references.at(references.size() - 1);
			
			string tag = line.substr(1);
			
			reference->name = parseNameFromTag(tag);
			reference->description = parseDescriptionFromTag(tag);
		}
		else if ( line[0] != '#' )
		{
			int length = line.length();
			
			if ( line[length - 1] == '\r' )
			{
				line.resize(length - 1);
			}
			
			reference->sequence.append(line);
		}
	}
	
	in.close();
}

void ReferenceList::initFromProtocolBuffer(const Harvest::Reference & msg)
{
	references.resize(0);
	references.resize(msg.references_size());
	
	for ( int i = 0; i < msg.references_size(); i++ )
	{
		references[i].name = parseNameFromTag(msg.references(i).tag());
		references[i].description = parseDescriptionFromTag(msg.references(i).tag());
		references[i].sequence = msg.references(i).sequence();
	}
}

void ReferenceList::writeToFasta(ostream & out) const
{
	for ( int i = 0; i < references.size(); i++ )
	{
		out << '>' << references[i].name;
		
		if ( references[i].description.length() )
		{
			out << ' ' << references[i].description;
		}
		
		out << endl;
		
		const string & sequence = references[i].sequence;
		int width = 0;
		
		for ( int j = 0; j < sequence.length(); j++ )
		{
			if ( width == 70 )
			{
				out << endl;
				width = 0;
			}
			
			out << sequence.at(j);
			width++;
		}
		
		out << endl;
	}
}

void ReferenceList::writeToProtocolBuffer(Harvest * msg) const
{
	for ( int i = 0; i < references.size(); i++ )
	{
		Harvest::Reference::Sequence * msgRef = msg->mutable_reference()->add_references();
		const Reference & reference = references.at(i);
		
		msgRef->mutable_tag()->append(reference.name);
		
		if ( reference.description.length() )
		{
			msgRef->mutable_tag()->append(" ");
			msgRef->mutable_tag()->append(reference.description);
		}
		
		msgRef->set_sequence(reference.sequence);
	}
}

string parseNameFromTag(string tag)
{
	for ( int i = 0; i < tag.length(); i++ )
	{
		if ( tag[i] == ' ' )
		{
			return tag.substr(0, i);
			break;
		}
	}
	
	return tag;
}

string parseDescriptionFromTag(string tag)
{
	for ( int i = 0; i < tag.length(); i++ )
	{
		if ( tag[i] == ' ' )
		{
			return tag.substr(i);
			break;
		}
	}
	
	return "";
}