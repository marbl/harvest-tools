// Copyright © 2014, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef AnnotationList_h
#define AnnotationList_h

#include <string>
#include <vector>
#include <stdexcept>

#include "harvest/capnp/harvest.capnp.h"
#include "harvest/pb/harvest.pb.h"
#include "harvest/ReferenceList.h"

struct Annotation
{
	int start;
	int end;
	bool reverse;
	std::string name;
	std::string locus;
	std::string description;
	std::string feature;
};

class AnnotationList
{
public:
	
	class NoSequenceException : public std::exception
	{
	public:
		
		NoSequenceException(const std::string & fileNew)
		{
			file = fileNew;
		}
		
		virtual ~NoSequenceException() throw() {}
		
		std::string file;
	};
	
	class NoAccException : public std::exception
	{
	public:
		
		NoAccException(const std::string & fileNew)
		{
			file = fileNew;
		}
		
		virtual ~NoAccException() throw() {}
		
		std::string file;
	};
	
	void clear();
	int getAnnotationCount() const;
	const Annotation & getAnnotation(int index) const;
	void initFromCapnp(const capnp::Harvest::Reader & harvestReader, const ReferenceList & referenceList);
	void initFromGenbank(const char * file, ReferenceList & referenceList, bool useSeq);
	void initFromProtocolBuffer(const Harvest::AnnotationList & msg, const ReferenceList & referenceList);
	void writeToCapnp(capnp::Harvest::Builder & harvestBuilder, const ReferenceList & referenceList) const;
	void writeToProtocolBuffer(Harvest * msg, const ReferenceList & referenceList) const;
	
private:
	
	std::vector<Annotation> annotations;
};

inline int AnnotationList::getAnnotationCount() const { return annotations.size(); }
inline const Annotation & AnnotationList::getAnnotation(int index) const { return annotations.at(index); }

#endif
