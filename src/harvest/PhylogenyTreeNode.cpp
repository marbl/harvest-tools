// Copyright © 2014, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "PhylogenyTreeNode.h"

#include <stdlib.h>

using namespace::std;

PhylogenyTreeNode::PhylogenyTreeNode(const Harvest::Tree::Node & msgNode, PhylogenyTreeNode * parent)
{
	// load from protobuf
	
	this->parent = parent;
	children.resize(msgNode.children_size());
	
	trackId = msgNode.track();
	distance = msgNode.branchlength();
	bootstrap = msgNode.bootstrap();
	
	for ( int i = 0; i < children.size(); i++ )
	{
		children[i] = new PhylogenyTreeNode(msgNode.children(i), this);
	}
	
//	collapse = false;
}

PhylogenyTreeNode::PhylogenyTreeNode(char *& token, TrackList * trackList, bool useNames, PhylogenyTreeNode * parent)
{
	ParseState state = STATE_start;
	char * valueStart;
	
	this->parent = parent;
	bootstrap = 0;
	distance = 0;
	
	while ( state != STATE_end )
	{
		while ( *token == '\n' || *token == '\r' )
		{
			token++;
		}
		
		if ( state == STATE_start )
		{
			if ( *token == '(' )
			{
				state = STATE_children;
			}
			else
			{
				state = STATE_nameLeaf;
				valueStart = token;
			}
			
			token++;
		}
		else if ( state == STATE_children )
		{
			if ( *token == ')' )
			{
				state = STATE_nameInternal;
				valueStart = token + 1;
				token++;
			}
			else if ( *token == ',' )
			{
				token++;
			}
			else
			{
				children.push_back(new PhylogenyTreeNode(token, trackList, useNames, this));
			}
		}
		else if ( state == STATE_nameLeaf || state == STATE_nameInternal )
		{
			if ( *token == 0 )
			{
				state = STATE_end;
			}
			else if ( *token == ':' )
			{
				if ( valueStart != token )
				{
					*token = 0;
					
					if
					(
						(*valueStart == '"' && *(token - 1) == '"') ||
						(*valueStart == '\'' && *(token - 1) == '\'')
					)
					{
						// remove quotes
						
						valueStart++;
						*(token - 1) = 0;
					}
					
					if ( state == STATE_nameInternal )
					{
						bootstrap = atof(valueStart);
					}
					else
					{
						if ( useNames )
						{
							trackId = trackList->addTrack(valueStart);
						}
						else
						{
							trackId = trackList->getTrackIndexByFile(valueStart);
						}
					}
				}
				
				state = STATE_length;
				valueStart = token + 1;
			}
			
			token++;
		}
		else if ( state == STATE_length )
		{
			if ( *token == ',' || *token == ')' )
			{
				//*token = 0;
				distance = atof(valueStart);
				state = STATE_end;
			}
			else
			{
				token++;
			}
		}
	}
}

PhylogenyTreeNode::PhylogenyTreeNode(PhylogenyTreeNode * child1, PhylogenyTreeNode * child2)
{
	// edge bisection
	
	distance = 0;
	parent = 0;
	bootstrap = 1;
	
	children.resize(2);
	
	children[0] = child1;
	children[1] = child2;
}

PhylogenyTreeNode::~PhylogenyTreeNode()
{
	for ( int i = 0; i < children.size(); i++ )
	{
		delete children[i];
	}
}

PhylogenyTreeNode * PhylogenyTreeNode::bisectEdge(float distanceLower)
{
	PhylogenyTreeNode * parentNew = new PhylogenyTreeNode(this, parent);
	
	parent->invert(this);
	parent->setParent(parentNew, distance - distanceLower);
	distance = distanceLower;
	parent = parentNew;
	
	return parentNew;
}

PhylogenyTreeNode * PhylogenyTreeNode::collapse()
{
	PhylogenyTreeNode * child = children[0];
	
	children[0]->setParent(parent, distance + children[0]->getDistance());
	children.resize(0);
	
	return child;
}

void PhylogenyTreeNode::getLeafIds(vector<int> & ids) const
{
	if ( children.size() == 0 )
	{
		ids.push_back(trackId);
	}
	else
	{
		for ( int i = 0; i < children.size(); i++ )
		{
			children[i]->getLeafIds(ids);
		}
	}
}

void PhylogenyTreeNode::getLeaves(vector<PhylogenyTreeNode *> & leaves)
{
	if ( children.size() == 0 )
	{
		leaves.push_back(this);
	}
	else
	{
		for ( int i = 0; i < children.size(); i++ )
		{
			children[i]->getLeaves(leaves);
		}
	}
}

void PhylogenyTreeNode::getPairwiseDistances(float ** matrix, int size)
{
	for ( int i = 0; i < size; i++ )
	{
		if ( i < leafMin || i > leafMax )
		{
			for ( int j = leafMin; j <= leafMax; j++ )
			{
				int row;
				int col;
				
				if ( i > j )
				{
					row = i;
					col = j;
				}
				else
				{
					row = j;
					col = i;
				}
				
				matrix[row - 1][col] += distance;
			}
		}
	}
	
	for ( int i = 0; i < children.size(); i++ )
	{
		children[i]->getPairwiseDistances(matrix, size);
	}
}

void PhylogenyTreeNode::initialize(int & newId, int &leaf, float depthParent)
{
	id = newId;
	newId++;
	leafMin = leaf;
	depth = depthParent + distance;
	
	for ( int i = 0; i < children.size(); i++ )
	{
		children[i]->initialize(newId, leaf, depth);
	}
	
	if ( children.size() == 0 )
	{
		leaf++;
	}
	
	leafMax = leaf - 1;
}

void PhylogenyTreeNode::invert(PhylogenyTreeNode * fromChild)
{
	vector<PhylogenyTreeNode *> childrenNew;
	
	for ( int i = 0; i < children.size(); i++ )
	{
		if ( children[i] != fromChild )
		{
			childrenNew.push_back(children[i]);
		}
	}
	
	if ( parent )
	{
		childrenNew.push_back(parent);
	}
	
	children.resize(childrenNew.size());
	
	for ( int i = 0; i < children.size(); i++ )
	{
		children[i] = childrenNew.at(i);
	}
	
	if ( parent )
	{
		parent->invert(this);
		
		if ( parent->getChildrenCount() == 1 )
		{
			children[children.size() - 1] = parent->collapse();
			delete parent;
		}
	}
	
	if ( fromChild )
	{
		distance = fromChild->getDistance();
	}
	
	parent = fromChild;
}

void PhylogenyTreeNode::setParent(PhylogenyTreeNode *parentNew, float distanceNew)
{
	parent = parentNew;
	distance = distanceNew;
}

void PhylogenyTreeNode::swapSiblings()
{
	PhylogenyTreeNode * temp = children[0];
	children[0] = children[1];
	children[1] = temp;
}

void PhylogenyTreeNode::writeToNewick(std::ostream &out, const TrackList & trackList, const double mult) const
{
	if ( children.size() )
	{
		out << '(';
		
		for ( int i = 0; i < children.size(); i++ )
		{
		        children[i]->writeToNewick(out, trackList,mult);
		
			if ( i < children.size() - 1 )
			{
				out << ',';
			}
		}
		
		out << ')';
		
		if ( bootstrap != 0 )
		{
			out << bootstrap;
		}
	}
	else
	{
		out << '\'' << trackList.getTrack(trackId).file << '\'';
	}
        //by default, always use multiplier
        //can be 1.0, or an adjusted value
        //alternatively, this could be conditional based on the parameter, instead of using 1.0 vs non-1.0 values
	out << ':' << distance * mult;
}

void PhylogenyTreeNode::writeToProtocolBuffer(Harvest::Tree::Node * msgNode) const
{
	if ( children.size() )
	{
		for ( int i = 0; i < children.size(); i++ )
		{
			children[i]->writeToProtocolBuffer(msgNode->add_children());
		}
		
		if ( bootstrap != 0 )
		{
			msgNode->set_bootstrap(bootstrap);
		}
	}
	else
	{
		msgNode->set_track(trackId);
	}
	
	msgNode->set_branchlength(distance);
}

