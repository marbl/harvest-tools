#ifndef LcbList_h
#define LcbList_h

#include <vector>
#include "harvest/ReferenceList.h"
#include "harvest/PhylogenyTree.h"
#include "harvest/TrackList.h"

#include "harvest/pb/harvest.pb.h"

class VariantList;

class LcbList
{
public:
	
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
	
	const Lcb & getLcb(int index) const;
	int getLcbCount() const;
	void initFromMfa(const char * file, ReferenceList * referenceList, TrackList * trackList, PhylogenyTree * phylogenyTree, VariantList * variantList);
	void initFromProtocolBuffer(const Harvest::Alignment & msgAlignment);
	void initFromXmfa(const char * file, const ReferenceList & referenceList, TrackList * trackList, PhylogenyTree * phylogenyTree, VariantList * variantList);
	void initWithSingleLcb(const ReferenceList & referenceList, const TrackList & trackList);
	void writeToProtocolBuffer(Harvest * msg) const;
	void writeToXmfa(std::ostream & out, const ReferenceList & referenceList, const TrackList & trackList, const VariantList & variantList) const;
	
private:
	
	std::vector<Lcb> lcbs;
};

inline const LcbList::Lcb & LcbList::getLcb(int index) const { return lcbs.at(index); }
inline int LcbList::getLcbCount() const { return lcbs.size(); }

#endif