// Copyright © 2014, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef ReferenceList_h
#define ReferenceList_h

#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>

#include "harvest/capnp/harvest.capnp.h"
#include "harvest/pb/harvest.pb.h"

static const int undef = -1;
	
struct Reference
{
	std::string name;
	std::string description;
	std::string sequence;
};

class ReferenceList
{
public:
	
	class GiNotFoundException : public std::exception
	{
	public:
		
		GiNotFoundException(const std::string & giNew)
		{
			gi = giNew;
		}
		
		virtual ~GiNotFoundException() throw() {}
		
		std::string gi;
	};
	
	class NameNotFoundException : public std::exception
	{
	public:
		
		NameNotFoundException(const std::string & nameNew)
		{
			name = nameNew;
		}
		
		virtual ~NameNotFoundException() throw() {}
		
		std::string name;
	};
	
	void addReference(std::string name, std::string desc, std::string sequence);
	void clear();
	long int getConcatenatedPosition(int sequence, long int position) const;
	int getPositionFromConcatenated(int sequence, long int position) const;
	const Reference & getReference(int index) const;
	int getReferenceCount() const;
	int getReferenceSequenceFromConcatenated(long int position) const;
	int getReferenceSequenceFromGi(long int gi) const;
	int getReferenceSequenceFromName(std::string name) const;
	void initFromCapnp(const capnp::Harvest::Reader & harvestReader);
	void initFromFasta(const char * file);
	void initFromProtocolBuffer(const Harvest::Reference & msg);
	void writeToCapnp(capnp::Harvest::Builder & harvestBuilder) const;
	void writeToFasta(std::ostream & out) const;
	void writeToProtocolBuffer(Harvest * msg) const;
	
private:
	
	std::vector<Reference> references;
};

std::string parseNameFromTag(std::string tag);
std::string parseDescriptionFromTag(std::string tag);

inline void ReferenceList::clear() { references.resize(0); }
inline int ReferenceList::getReferenceCount() const { return references.size(); }
inline const Reference & ReferenceList::getReference(int index) const { return references.at(index); }

#endif
