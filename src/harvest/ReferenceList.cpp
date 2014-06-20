
#include <fstream>
#include "ReferenceList.h"

using namespace::std;

void ReferenceList::addReference(string name, string sequence)
{
	references.resize(references.size() + 1);
	references[references.size() - 1].name = name;
	references[references.size() - 1].sequence = sequence;
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
	
	return undef;
}

void ReferenceList::initFromFasta(const char * file)
{
	ifstream in(file);
	char line[1 << 20];
	Reference * reference;
	
	while ( ! in.eof() )
	{
		in.getline(line, (1 << 20) - 1);
		
		if ( *line == '>' )
		{
			references.resize(references.size() + 1);
			reference = &references.at(references.size() - 1);
			
			reference->name = line + 1;
		}
		else if ( *line != '#' )
		{
			int length = strlen(line);
			
			if ( line[length - 1] == '\r' )
			{
				line[length - 1] = 0;
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
		references[i].name = msg.references(i).tag();
		references[i].sequence = msg.references(i).sequence();
	}
}

void ReferenceList::writeToFasta(ostream & out) const
{
	for ( int i = 0; i < references.size(); i++ )
	{
		out << '>' << references[i].name << endl;
		
		const string & sequence = references[i].sequence;
		int width = 0;
		
		for ( int j = 0; j < sequence.length(); j++ )
		{
			if ( width == 80 )
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
		
		msgRef->set_tag(references.at(i).name);
		msgRef->set_sequence(references.at(i).sequence);
	}
}
