// Copyright © 2014, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef AnnotationList_h
#define AnnotationList_h

#include <string>
#include <vector>
#include <map>
#include <stdexcept>

#include "harvest/capnp/harvest.capnp.h"
#include "harvest/pb/harvest.pb.h"
#include "harvest/ReferenceList.h"

static const std::map<std::string, std::string> translations =
{
	{"TTT", "F"},
	{"TCT", "S"},
	{"TAT", "Y"},
	{"TGT", "C"},
	{"TTC", "F"},
	{"TCC", "S"},
	{"TAC", "Y"},
	{"TGC", "C"},
	{"TTA", "L"},
	{"TCA", "S"},
	{"TAA", "*"},
	{"TGA", "*"},
	{"TTG", "L"},
	{"TCG", "S"},
	{"TAG", "*"},
	{"TGG", "W"},
	{"CTT", "L"},
	{"CCT", "P"},
	{"CAT", "H"},
	{"CGT", "R"},
	{"CTC", "L"},
	{"CCC", "P"},
	{"CAC", "H"},
	{"CGC", "R"},
	{"CTA", "L"},
	{"CCA", "P"},
	{"CAA", "Q"},
	{"CGA", "R"},
	{"CTG", "L"},
	{"CCG", "P"},
	{"CAG", "Q"},
	{"CGG", "R"},
	{"ATT", "I"},
	{"ACT", "T"},
	{"AAT", "N"},
	{"AGT", "S"},
	{"ATC", "I"},
	{"ACC", "T"},
	{"AAC", "N"},
	{"AGC", "S"},
	{"ATA", "I"},
	{"ACA", "T"},
	{"AAA", "K"},
	{"AGA", "R"},
	{"ATG", "M"},
	{"ACG", "T"},
	{"AAG", "K"},
	{"AGG", "R"},
	{"GTT", "V"},
	{"GCT", "A"},
	{"GAT", "D"},
	{"GGT", "G"},
	{"GTC", "V"},
	{"GCC", "A"},
	{"GAC", "D"},
	{"GGC", "G"},
	{"GTA", "V"},
	{"GCA", "A"},
	{"GAA", "E"},
	{"GGA", "G"},
	{"GTG", "V"},
	{"GCG", "A"},
	{"GAG", "E"},
	{"GGG", "G"},
};

static const std::map<std::string, std::string> translationsRc =
{
	{"TTT", "K"},
	{"GTT", "N"},
	{"CTT", "K"},
	{"ATT", "N"},
	{"TGT", "T"},
	{"GGT", "T"},
	{"CGT", "T"},
	{"AGT", "T"},
	{"TCT", "R"},
	{"GCT", "S"},
	{"CCT", "R"},
	{"ACT", "S"},
	{"TAT", "I"},
	{"GAT", "I"},
	{"CAT", "M"},
	{"AAT", "I"},
	{"TTG", "Q"},
	{"GTG", "H"},
	{"CTG", "Q"},
	{"ATG", "H"},
	{"TGG", "P"},
	{"GGG", "P"},
	{"CGG", "P"},
	{"AGG", "P"},
	{"TCG", "R"},
	{"GCG", "R"},
	{"CCG", "R"},
	{"ACG", "R"},
	{"TAG", "L"},
	{"GAG", "L"},
	{"CAG", "L"},
	{"AAG", "L"},
	{"TTC", "E"},
	{"GTC", "D"},
	{"CTC", "E"},
	{"ATC", "D"},
	{"TGC", "A"},
	{"GGC", "A"},
	{"CGC", "A"},
	{"AGC", "A"},
	{"TCC", "G"},
	{"GCC", "G"},
	{"CCC", "G"},
	{"ACC", "G"},
	{"TAC", "V"},
	{"GAC", "V"},
	{"CAC", "V"},
	{"AAC", "V"},
	{"TTA", "*"},
	{"GTA", "Y"},
	{"CTA", "*"},
	{"ATA", "Y"},
	{"TGA", "S"},
	{"GGA", "S"},
	{"CGA", "S"},
	{"AGA", "S"},
	{"TCA", "*"},
	{"GCA", "C"},
	{"CCA", "W"},
	{"ACA", "C"},
	{"TAA", "L"},
	{"GAA", "F"},
	{"CAA", "L"},
	{"AAA", "F"},
};

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
