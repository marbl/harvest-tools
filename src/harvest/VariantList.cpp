
#include "harvest/VariantList.h"

using namespace::std;

void VariantList::addVariantsFromAlignment(const vector<string> & seqs, const ReferenceList & referenceList, int sequence, int position, bool lcbfilt)
{
//	Harvest::Variation * msg = harvest.mutable_variation();
	char col[seqs.size() + 1];
	int offset = 0;
	
	col[seqs.size()] = 0;
	
	// Since insertions to the reference take on the left-most reference
	// position, this allows the alignment to start with an insertion
	// (possibly at reference position -1).
	//
	position--;
	
	for ( int i = 0; i < seqs[0].length(); i++ )
	{
		bool variant = false;
		bool n = false;

		col[0] = seqs[0][i];
		
		bool indel = col[0] == '-';
		
		if ( indel )
		{
			// insertion relative to the reference
			offset++;
		}
		else
		{
			position++;
			offset = 0;
		}
		
		for ( int j = 1; j < seqs.size(); j++ )
		{
			col[j] = seqs[j][i];
			
			if ( ! variant && col[j] != col[0] )
			{
				variant = true;
			}
			
			if ( ! indel && col[j] == '-' )
			{
				indel = true;
			}
			
			if ( col[j] == 'N' || col[j] == 'n' )
			{
				n = true;
			}
		}
		
		if ( variant )
		{
			variants.resize(variants.size() + 1);
			Variant * varNew = &variants[variants.size() - 1];
			
			while ( position >= 0 && position >= referenceList.getReference(sequence).sequence.length() )
			{
				position -= referenceList.getReference(sequence).sequence.length();
				sequence++;
			}
			
			varNew->sequence = sequence;
			varNew->position = position;
			varNew->offset = offset;
			varNew->alleles = col;
			varNew->filters = 0;
			
			if ( indel )
			{
				varNew->filters |= FILTER_indel;
			}
			
			if ( n )
			{
				varNew->filters |= FILTER_n;
			}
			
			if (lcbfilt)
			{
				varNew->filters |= FILTER_lcb;
			}
			
			varNew->quality = 0;
		}
	}
}

void VariantList::init()
{
	filters.resize(0);
	addFilter(FILTER_indel, "IND", "Column contains indel");
	addFilter(FILTER_n, "N", "Column contains N");
	addFilter(FILTER_lcb, "LCB", "LCB smaller than 200bp");
	
	variants.resize(0);
}

void VariantList::initFromProtocolBuffer(const Harvest::Variation & msgVariation)
{
	filters.resize(msgVariation.filters_size());
	
	for ( int i = 0; i < msgVariation.filters_size(); i++ )
	{
		filters[i].flag = msgVariation.filters(i).flag();
		filters[i].name = msgVariation.filters(i).name();
		filters[i].description = msgVariation.filters(i).description();
	}
	
	variants.resize(msgVariation.variants_size());
	
	for ( int i = 0; i < msgVariation.variants_size(); i++ )
	{
		Variant & variant = variants[i];
		const Harvest::Variation::Variant & msgVariant = msgVariation.variants(i);
		
		variant.sequence = msgVariant.sequence();
		variant.position = msgVariant.position();
		variant.alleles = msgVariant.alleles();
		variant.filters = msgVariant.filters();
		variant.quality = msgVariant.quality();
	}
}

void VariantList::sortVariants()
{
	sort(variants.begin(), variants.end(), variantLessThan);
}

void VariantList::writeToProtocolBuffer(Harvest * msg) const
{
	Harvest::Variation * msgVar = msg->mutable_variation();
	
	for ( int i = 0; i < filters.size(); i++ )
	{
		Harvest::Variation::Filter * msgFilter = msgVar->add_filters();
		
		msgFilter->set_flag(filters[i].flag);
		msgFilter->set_name(filters[i].name);
		msgFilter->set_description(filters[i].description);
	}
	
	for ( int i = 0; i < variants.size(); i++ )
	{
		Harvest::Variation::Variant * variant = msgVar->add_variants();
		
		variant->set_sequence(variants[i].sequence);
		variant->set_position(variants[i].position);
		variant->set_alleles(variants[i].alleles);
		variant->set_filters(variants[i].filters);
	}
}

void VariantList::addFilter(long long int flag, string name, string description)
{
	filters.resize(filters.size() + 1);
	filters[filters.size() - 1].flag = flag;
	filters[filters.size() - 1].name = name;
	filters[filters.size() - 1].description = description;
}
