// Copyright © 2014, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef VariantList_h
#define VariantList_h

#include <vector>
#include "harvest/pb/harvest.pb.h"
#include "harvest/LcbList.h"
#include "harvest/PhylogenyTree.h"
#include "harvest/ReferenceList.h"
#include "harvest/TrackList.h"

typedef long long unsigned int uint64;

class VariantList
{
public:
	
	struct Filter
	{
		uint64 flag;
		std::string name;
		std::string description;
	};
	
	struct Variant
	{
		int sequence;
		int position;
		int offset;
		std::string alleles;
		long long int filters;
		int quality;
	};
	
	void addFilterFromBed(const char * file, const char * name, const char * desc);
	void addVariantsFromAlignment(const std::vector<std::string> & seqs, const ReferenceList & referenceList, int sequence, int position, bool lcbfilt);
	void clear();
	const Filter & getFilter(int index) const;
	int getFilterCount() const;
	const Variant & getVariant(int index) const;
	int getVariantCount() const;
	void init();
	void initFromProtocolBuffer(const Harvest::Variation & msgVariation);
	void initFromVcf(const char * file, const ReferenceList & referenceList, TrackList * trackList, LcbList * lcbList, PhylogenyTree * phylogenyTree);
	void sortVariants();
	void writeToMfa(std::ostream &out, bool indels, const TrackList & trackList) const;
	void writeToProtocolBuffer(Harvest * harvest) const;
	void writeToVcf(std::ostream &out, bool indels, const ReferenceList & referenceList, const TrackList & trackList) const;
	
	static bool variantLessThan(const Variant & a, const Variant & b)
	{
		if ( a.sequence == b.sequence )
		{
			if ( a.position == b.position )
			{
				return a.offset < b.offset;
			}
			else
			{
				return a.position < b.position;
			}
		}
		else
		{
			return a.sequence < b.sequence;
		}
	}
	
private:
	
	enum FilterFlag
	{
		FILTER_indel = 1,
		FILTER_n = 2,
		FILTER_lcb = 4,
		FILTER_conservation = 8,
		FILTER_gaps = 16,
	};
	
	void addFilter(long long int flag, std::string name, std::string description);
	
	std::vector<Filter> filters;
	std::vector<Variant> variants;
};

inline const VariantList::Filter & VariantList::getFilter(int index) const { return filters.at(index); }
inline int VariantList::getFilterCount() const { return filters.size(); }
inline const VariantList::Variant & VariantList::getVariant(int index) const { return variants.at(index); }
inline int VariantList::getVariantCount() const { return variants.size(); }

#endif
