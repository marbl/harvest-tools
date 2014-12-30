// Copyright © 2014, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef PhylogenyTreeNode_h
#define PhylogenyTreeNode_h

#include <iostream>
#include <vector>
#include "harvest/pb/harvest.pb.h"
#include "harvest/TrackList.h"

class PhylogenyTreeNode
{
public:
	
	PhylogenyTreeNode(const Harvest::Tree::Node & msgNode, PhylogenyTreeNode * parent = 0);
	PhylogenyTreeNode(char *& token, TrackList * trackList, bool useNames, PhylogenyTreeNode * parent = 0);
	PhylogenyTreeNode(PhylogenyTreeNode * parent, PhylogenyTreeNode * child); // for edge bisection
	~PhylogenyTreeNode();
	
	PhylogenyTreeNode * bisectEdge(float distanceLower);
	PhylogenyTreeNode * collapse();
	float getBootstrap() const;
	PhylogenyTreeNode * getChild(unsigned int index) const;
	int getChildrenCount() const;
	float getDepth() const;
	double getDistance() const;
	int getId() const;
	int getTrackId() const;
	int getLeafCount() const;
	int getLeafMax() const;
	int getLeafMin() const;
	void getLeaves(std::vector<PhylogenyTreeNode *> & leaves);
	void getLeafIds(std::vector<int> & ids) const;
	void getPairwiseDistances(float ** matrix, int size);
	const PhylogenyTreeNode * getParent() const;
	void initialize(int & newId, int & leaf, float depthParent = 0);
	void invert(PhylogenyTreeNode * fromChild = 0);
	void setAlignDist(float dist, float dep);
	void setParent(PhylogenyTreeNode * parentNew, float distanceNew);
	void setTrackId(int trackIdNew);
	void swapSiblings();
	void writeToNewick(std::ostream &out, const TrackList & trackList, const double mult = 1.0) const;
	void writeToProtocolBuffer(Harvest::Tree::Node * msgNode) const;
	
private:
	
	enum ParseState
	{
		STATE_start,
		STATE_children,
		STATE_nameLeaf,
		STATE_nameInternal,
		STATE_length,
		STATE_end,
	};
	
	std::vector<PhylogenyTreeNode *> children;
	PhylogenyTreeNode * parent;
	int id;
	int trackId;
	double distance;
	float depth;
	float depthAlign;
	float distanceAlign;
	int leafMin;
	int leafMax;
	float bootstrap;
};

inline float PhylogenyTreeNode::getBootstrap() const {return bootstrap;}
inline PhylogenyTreeNode * PhylogenyTreeNode::getChild(unsigned int index) const {return children[index];};
inline int PhylogenyTreeNode::getChildrenCount() const {return children.size();}
inline float PhylogenyTreeNode::getDepth() const {return depth;}
inline double PhylogenyTreeNode::getDistance() const {return distance;}
inline int PhylogenyTreeNode::getId() const {return id;}
inline int PhylogenyTreeNode::getTrackId() const {return trackId;}
inline int PhylogenyTreeNode::getLeafCount() const {return leafMax - leafMin + 1;}
inline int PhylogenyTreeNode::getLeafMax() const {return leafMax;}
inline int PhylogenyTreeNode::getLeafMin() const {return leafMin;}
inline const PhylogenyTreeNode * PhylogenyTreeNode::getParent() const {return parent;}
inline void PhylogenyTreeNode::setTrackId(int trackIdNew) { trackId = trackIdNew; }

#endif
