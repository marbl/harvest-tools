#ifndef ReferenceList_h
#define ReferenceList_h

#include <string>
#include <vector>

#include "harvest/pb/harvest.pb.h"

static const int undef = -1;
	
struct Reference
{
	std::string name;
	std::string sequence;
};

class ReferenceList
{
public:
	
	void addReference(std::string name, std::string sequence);
	int getPositionFromConcatenated(int sequence, long int position) const;
	const Reference & getReference(int index) const;
	int getReferenceCount() const;
	int getReferenceSequenceFromConcatenated(long int position) const;
	int getReferenceSequenceFromGi(long int gi) const;
	void initFromFasta(const char * file);
	void initFromProtocolBuffer(const Harvest::Reference & msg);
	void writeToProtocolBuffer(Harvest * msg) const;
	
private:
	
	std::vector<Reference> references;
};

inline int ReferenceList::getReferenceCount() const { return references.size(); }
inline const Reference & ReferenceList::getReference(int index) const { return references.at(index); }

#endif
