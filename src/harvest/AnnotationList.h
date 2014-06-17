#ifndef AnnotationList_h
#define AnnotationList_h

#include <string>
#include <vector>

#include "harvest/pb/harvest.pb.h"
#include "harvest/ReferenceList.h"

struct Annotation
{
	int start;
	int end;
	bool reverse;
	std::string name;
	std::string locus;
	std::string description;
	std::string feature;
};

class AnnotationList
{
public:
	
	int getAnnotationCount() const;
	const Annotation & getAnnotation(int index) const;
	void initFromGenbank(const char * file, const ReferenceList & referenceList);
	void initFromProtocolBuffer(const Harvest::AnnotationList & msg, const ReferenceList & referenceList);
	void writeToProtocolBuffer(Harvest * msg, const ReferenceList & referenceList) const;
	
private:
	
	std::vector<Annotation> annotations;
};

inline int AnnotationList::getAnnotationCount() const { return annotations.size(); }
inline const Annotation & AnnotationList::getAnnotation(int index) const { return annotations.at(index); }

#endif
