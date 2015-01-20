// Copyright © 2014, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef harvest_PhylogenyTree
#define harvest_PhylogenyTree

#include <vector>
#include <iostream>

#include "harvest/capnp/harvest.capnp.h"
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
	double getMult() const;
	int getNodeCount() const;
	void initFromCapnp(const capnp::Harvest::Reader & harvestReader);
	void initFromNewick(const char * file, TrackList * trackList);
	void initFromProtocolBuffer(const Harvest::Tree & msg);
	float leafDistance(int leaf1, int leaf2) const;
	void midpointReroot();
	void setMult(double multNew);
	void setOutgroup(const PhylogenyTreeNode * node);
	void setTrackIndeces(int * trackIndecesNew);
	void writeToCapnp(capnp::Harvest::Builder & harvestBuilder) const;
	void writeToNewick(std::ostream &out, const TrackList & trackList, bool useMult) const;
	void writeToProtocolBuffer(Harvest * msg) const;
	
	PhylogenyTreeNode * getRoot() const;
private:
	
	void init();
	void reroot(const PhylogenyTreeNode * rootNew, float distance, bool reorder = false);
	std::vector<PhylogenyTreeNode *> leaves;
	PhylogenyTreeNode * root;
	int nodeCount;
	double mult;
};

inline const PhylogenyTreeNode * PhylogenyTree::getLeaf(int id) const {return leaves[id];}
inline int PhylogenyTree::getNodeCount() const {return nodeCount;}
inline double PhylogenyTree::getMult() const { return mult; }
inline PhylogenyTreeNode * PhylogenyTree::getRoot() const {return this->root;}
inline void PhylogenyTree::setMult(double multNew) { mult = multNew; }

#endif
