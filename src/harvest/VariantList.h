
#ifndef VariantList_h
#define VariantList_h

#include <vector>
#include "harvest/pb/harvest.pb.h"
#include "harvest/ReferenceList.h"

class VariantList
{
public:
	
	struct Filter
	{
		long long int flag;
		std::string name;
		std::string description;
	};
	
	struct Variant
	{
		int sequence;
		int position;
		int offset;
		std::string alleles;
		long long int filters;
		int quality;
	};
	
	void addVariantsFromAlignment(const std::vector<std::string> & seqs, const ReferenceList & referenceList, int sequence, int position, bool lcbfilt);
	const Variant & getVariant(int index) const;
	int getVariantCount() const;
	void init();
	void initFromProtocolBuffer(const Harvest::Variation & msgVariation);
	void sortVariants();
	void writeToProtocolBuffer(Harvest * harvest) const;
	
	static bool variantLessThan(const Variant & a, const Variant & b)
	{
		if ( a.sequence == b.sequence )
		{
			if ( a.position == b.position )
			{
				return a.offset < b.offset;
			}
			else
			{
				return a.position < b.position;
			}
		}
		else
		{
			return a.sequence < b.sequence;
		}
	}
	
private:
	
	enum FilterFlag
	{
		FILTER_indel = 1,
		FILTER_n = 2,
		FILTER_lcb = 4,
	};
	
	void addFilter(long long int flag, std::string name, std::string description);
	
	std::vector<Filter> filters;
	std::vector<Variant> variants;
};

inline const VariantList::Variant & VariantList::getVariant(int index) const { return variants.at(index); }
inline int VariantList::getVariantCount() const { return variants.size(); }

#endif
