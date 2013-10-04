#ifndef harvest_TreeNode
#define harvest_TreeNode

#include "harvest.pb.h"

class TreeNode
{
public:

	TreeNode
	(
		const Harvest::Tree::Node & msgNode,
		int & newId,
		int & leaf,
		const TreeNode * parent = 0, float depth = 0);
	~TreeNode();
	
	float getBootstrap() const;
	TreeNode * getChild(unsigned int index) const;
	int getChildrenCount() const;
	float getDistance() const;
	int getId() const;
	int getTrack() const;
	const TreeNode * getParent() const;
	
private:
	
	TreeNode** children;
	const TreeNode * parent;
	int childrenCount;
	int track;
	float distance;
	float bootstrap;
};

inline float TreeNode::getBootstrap() const {return bootstrap;}
inline TreeNode * TreeNode::getChild(unsigned int index) const {return children[index];};
inline int TreeNode::getChildrenCount() const {return childrenCount;}
inline double TreeNode::getDistance() const {return distance;}
inline int TreeNode::getId() const {return id;}
inline int TreeNode::getTrack() const {return trackId;}
inline const TreeNode * TreeNode::getParent() const {return parent;}

#endif
