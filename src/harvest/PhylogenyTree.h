// Copyright © 2014, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef harvest_PhylogenyTree
#define harvest_PhylogenyTree

#include <vector>
#include <iostream>

#include "harvest/pb/harvest.pb.h"
#include "harvest/PhylogenyTreeNode.h"
#include "harvest/TrackList.h"

class PhylogenyTree
{
public:
	
	PhylogenyTree();
	~PhylogenyTree();
	
	void clear();
	const PhylogenyTreeNode * getLeaf(int id) const;
	void getLeafIds(std::vector<int> & ids) const;
	int getNodeCount() const;
	void initFromNewick(const char * file, TrackList * trackList);
	void initFromProtocolBuffer(const Harvest::Tree & msg);
	float leafDistance(int leaf1, int leaf2) const;
	void midpointReroot();
	void setOutgroup(const PhylogenyTreeNode * node);
	void setTrackIndeces(int * trackIndecesNew);
	void writeToNewick(std::ostream &out, const TrackList & trackList) const;
	void writeToProtocolBuffer(Harvest * msg) const;
	
	PhylogenyTreeNode * getRoot() const;
	double mult;
private:
	
	void init();
	void reroot(const PhylogenyTreeNode * rootNew, float distance, bool reorder = false);
	std::vector<PhylogenyTreeNode *> leaves;
	PhylogenyTreeNode * root;
	int nodeCount;

};

inline const PhylogenyTreeNode * PhylogenyTree::getLeaf(int id) const {return leaves[id];}
inline int PhylogenyTree::getNodeCount() const {return nodeCount;}
inline PhylogenyTreeNode * PhylogenyTree::getRoot() const {return this->root;}

#endif
