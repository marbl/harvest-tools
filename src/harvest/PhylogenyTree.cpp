
#include "PhylogenyTree.h"
#include <fstream>

using namespace::std;

PhylogenyTree::PhylogenyTree()
{
	root = 0;
}

PhylogenyTree::~PhylogenyTree()
{
	if ( root )
	{
		delete root;
	}
}

void PhylogenyTree::getLeafIds(vector<int> & ids) const
{
	ids.resize(0);
	root->getLeafIds(ids);
}

void PhylogenyTree::init()
{
	int leaf = 0;
	
	root->initialize(nodeCount, leaf);
	leaves.resize(0);
	root->getLeaves(leaves);
}

void PhylogenyTree::initFromNewick(const char * file, TrackList * trackList)
{
	if ( root )
	{
		delete root;
	}
	
	ifstream in(file);
	char line[1 << 20];
	
	bool useNames = trackList->getTrackCount() == 0;
	
	while ( in.getline(line, (1 << 20) - 1) )
	{
		char * token = line;
		root = new PhylogenyTreeNode(token, trackList, useNames);
	}
	
	in.close();
	init();
}

void PhylogenyTree::initFromProtocolBuffer(const Harvest::Tree & msg)
{
	if ( root )
	{
		delete root;
	}
	
	nodeCount = 0;
	int leaf = 0;
	root = new PhylogenyTreeNode(msg.root());
	init();
}

float PhylogenyTree::leafDistance(int leaf1, int leaf2) const
{
	const PhylogenyTreeNode * node1 = leaves[leaf1];
	const PhylogenyTreeNode * node2 = leaves[leaf2];
	
	float distance = 0;
	
	while ( node1 != node2 )
	{
		if ( node1->getDepth() > node2->getDepth() )
		{
			distance += node1->getDistance();
			node1 = node1->getParent();
		}
		else if ( node2->getDepth() > node1->getDepth() )
		{
			distance += node2->getDistance();
			node2 = node2->getParent();
		}
		else
		{
			distance += node1->getDistance();
			distance += node2->getDistance();
			node1 = node1->getParent();
			node2 = node2->getParent();
		}
	}
	
	return distance;
}

void PhylogenyTree::midpointReroot()
{
	// lower triangular matrix of pairwise distances between leaves
	//
	int leavesCount = leaves.size();
	float ** distance = new float*[leavesCount - 1];
	
	for ( int i = 0; i < leavesCount - 1; i++ )
	{
		distance[i] = new float[i + 1];
		memset(distance[i], 0, sizeof(float) * (i + 1));
	}
	
	root->getPairwiseDistances(distance, leavesCount);
	
	float max = 0;
	int maxLeaf1;
	int maxLeaf2;
	
	for ( int i = 0; i < leavesCount - 1; i++ )
	{
		for ( int j = 0; j < i + 1; j++ )
		{
			if ( distance[i][j] > max )
			{
				max = distance[i][j];
				maxLeaf1 = i + 1;
				maxLeaf2 = j;
			}
		}
	}
	
	float midDistance = distance[maxLeaf1 - 1][maxLeaf2] / 2;
	
	for ( int i = 0; i < leavesCount - 1; i++ )
	{
		delete [] distance[i];
	}
	
	delete [] distance;
	
	const PhylogenyTreeNode * node;
	
	if ( leaves[maxLeaf1]->getDepth() > leaves[maxLeaf2]->getDepth() )
	{
		node = leaves[maxLeaf1];
	}
	else
	{
		node = leaves[maxLeaf2];
	}
	
	float depth = 0;
	
	while ( depth + node->getDistance() < midDistance && node->getParent() )
	{
		depth += node->getDistance();
		node = node->getParent();
	}
	
	if ( node != root )
	{
		reroot(node, midDistance - depth);
	}
}

void PhylogenyTree::setOutgroup(const PhylogenyTreeNode * node)
{
	reroot(node, node->getParent() == root ? (root->getChild(0)->getDistance() + root->getChild(1)->getDistance()) / 2 : node->getDistance() / 2, true);
}

void PhylogenyTree::setTrackIndeces(int * trackIndecesNew)
{
	for ( int i = 0; i < leaves.size(); i++ )
	{
		leaves[i]->setTrackId(trackIndecesNew[leaves[i]->getTrackId()]);
	}
}

void PhylogenyTree::reroot(const PhylogenyTreeNode * rootNew, float distance, bool reorder)
{
	int leaf = 0;
	nodeCount = 0;
	
	if ( rootNew->getParent() == root )
	{
		PhylogenyTreeNode * rootNewMutable;
		PhylogenyTreeNode * sibling;
		
		if ( root->getChild(0) == rootNew )
		{
			rootNewMutable = root->getChild(0);
			sibling = root->getChild(1);
		}
		else
		{
			sibling = root->getChild(0);
			rootNewMutable = root->getChild(1);
			
			if ( reorder )
			{
				root->swapSiblings();
			}
		}
		
		sibling->setParent(root, rootNew->getDistance() + sibling->getDistance() - distance);
		rootNewMutable->setParent(root, distance);
	}
	else
	{
		root = const_cast<PhylogenyTreeNode *>(rootNew)->bisectEdge(distance);
	}
	
	root->initialize(nodeCount, leaf);
	leaves.resize(0);
	root->getLeaves(leaves);
	//root->setAlignDist(root->getDistanceMax(), 0);
}

void PhylogenyTree::writeToNewick(std::ostream &out, const TrackList & trackList) const
{
	root->writeToNewick(out, trackList);
	out << '\n';
}

void PhylogenyTree::writeToProtocolBuffer(Harvest * msg) const
{
	root->writeToProtocolBuffer(msg->mutable_tree()->mutable_root());
}
