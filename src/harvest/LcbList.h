// Copyright © 2014, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef LcbList_h
#define LcbList_h

#include <vector>
#include "harvest/ReferenceList.h"
#include "harvest/PhylogenyTree.h"
#include "harvest/TrackList.h"
#include <stdexcept>

#include "harvest/capnp/harvest.capnp.h"
#include "harvest/pb/harvest.pb.h"

class VariantList;

class LcbList
{
public:
	
	class NoCoreException : public std::exception
	{
	public:
		
		NoCoreException(int queryCountNew)
		{
			queryCount = queryCountNew;
		}
		
		virtual ~NoCoreException() throw() {}
		
		int queryCount;
	};
	
	struct Interval
	{
		int sequence;
		int start;
		int end;
		
		Interval(int sequenceNew, int startNew, int endNew)
			: sequence(sequenceNew), start(startNew), end(endNew) {}
	};
	
	struct Region
	{
		int position;
		int length;
		bool reverse;
	};
	
	struct Lcb
	{
		std::vector<Region> regions;
		int sequence;
		int position;
		int length;
		float concordance;
	};
	
	void addLcbByReference(int startSeq, int startPos, int endSeq, int endPos, const ReferenceList & referenceList, const TrackList & trackList);
	void clear();
	const Lcb & getLcb(int index) const;
        double getCoreSize() const;
	int getLcbCount() const;
	void initFromCapnp(const capnp::Harvest::Reader & harvestReader);
	void initFromMaf(const char * file, ReferenceList * referenceList, TrackList * trackList, PhylogenyTree * phylogenyTree, VariantList * variantList, const char * referenceFileName);
	void initFromMfa(const char * file, ReferenceList * referenceList, TrackList * trackList, PhylogenyTree * phylogenyTree, VariantList * variantList);
	void initFromProtocolBuffer(const Harvest::Alignment & msgAlignment);
	void initFromXmfa(const char * file, ReferenceList * referenceList, TrackList * trackList, PhylogenyTree * phylogenyTree, VariantList * variantList);
	void initWithSingleLcb(const ReferenceList & referenceList, const TrackList & trackList);
	void writeToCapnp(capnp::Harvest::Builder & harvestBuilder) const;
	void writeToMfa(std::ostream & out, const ReferenceList & referenceList, const TrackList & trackList, const VariantList & variantList) const;
	void writeFilteredToMfa(std::ostream & out, std::ostream & out2, const ReferenceList & referenceList, const TrackList & trackList, const VariantList & variantList) const;
	void writeToProtocolBuffer(Harvest * msg) const;
	void writeToXmfa(std::ostream & out, const ReferenceList & referenceList, const TrackList & trackList, const VariantList & variantList) const;
	
private:
	
	std::vector<Lcb> lcbs;
};

inline const LcbList::Lcb & LcbList::getLcb(int index) const { return lcbs.at(index); }
inline int LcbList::getLcbCount() const { return lcbs.size(); }

#endif
