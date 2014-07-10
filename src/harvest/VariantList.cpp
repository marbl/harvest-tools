
#include "harvest/VariantList.h"
#include <fstream>

using namespace::std;

void VariantList::addFilterFromBed(const char * file, const char * name, const char * desc)
{
	ifstream in(file);
	char line[1 << 20];
	int i = 0;
	long long int flag = 1 << filters.size();
	
	addFilter(flag, name, desc);
	
	while ( ! in.eof() )
	{
		in.getline(line, (1 << 20) - 1);
		
		if ( in.eof() )
		{
			break;
		}
		
		int seq = atoi(strtok(line, "\t")) - 1;
		int start = atoi(strtok(0, "\t")) - 1;
		int end = atoi(strtok(0, "\t")) - 1;
		
		// seek to interval start
		//
		while
		(
			i < variants.size() &&
			(
				variants.at(i).sequence < seq ||
				(
					variants.at(i).sequence == seq &&
					variants.at(i).position < start
				)
			)
		)
		{
			i++;
		}
		
		// set flags through interval
		//
		while
		(
			i < variants.size() &&
			variants.at(i).sequence == seq &&
			variants.at(i).position <= end
		)
		{
			variants[i].filters = variants.at(i).filters | flag;
			i++;
		}
	}
}

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
			
			while ( referenceList.getReferenceCount() > 0 && position >= 0 && position >= referenceList.getReference(sequence).sequence.length() )
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

void VariantList::clear()
{
	filters.clear();
	variants.clear();
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
//		cout << "Filter " << i << '\n';
		filters[i].flag = msgVariation.filters(i).flag();
		filters[i].name = msgVariation.filters(i).name();
		filters[i].description = msgVariation.filters(i).description();
	}
	
	variants.resize(msgVariation.variants_size());
	
	for ( int i = 0; i < msgVariation.variants_size(); i++ )
	{
//		cout << "Variant " << i << '\n';
		Variant & variant = variants[i];
		const Harvest::Variation::Variant & msgVariant = msgVariation.variants(i);
		
		variant.sequence = msgVariant.sequence();
		variant.position = msgVariant.position();
		variant.alleles = msgVariant.alleles();
		variant.filters = msgVariant.filters();
		variant.quality = msgVariant.quality();
	}
}

void VariantList::initFromVcf(const char * file, const ReferenceList & referenceList, TrackList * trackList, LcbList * lcbList, PhylogenyTree * phylogenyTree)
{
	filters.resize(0);
	variants.resize(0);
	
	ifstream in(file);
	
	char line[1 << 20];
	map<string, long long int> flagsByFilter;
	map<string, int> refByTag;
	unsigned int alleleCount = 0;
	
	bool oldTags = trackList->getTrackCount();
	int * trackIndecesNew;
	
	if ( oldTags )
	{
		trackIndecesNew = new int[trackList->getTrackCount()];
	}
	
	for ( int i = 0; i < referenceList.getReferenceCount(); i++ )
	{
		refByTag[referenceList.getReference(i).name] = i;
	}
	
	while ( ! in.eof() )
	{
		if ( in.peek() == '#' )
		{
			in.getline(line, (1 << 20) - 1);
			
			if ( strncmp(line, "##FILTER", 8) == 0 )
			{
				char * token;
				
				filters.resize(filters.size() + 1);
				Filter & filter = filters[filters.size() - 1];
				
				token = strtok(line, "<");
				
				while ( (token = strtok(0, ",>\"")) )
				{
					if ( strncmp(token, "ID=", 3) == 0 )
					{
						filter.name = token + 3;
					}
					else if ( strcmp(token, "Description=") == 0 )
					{
						filter.description = strtok(0, "\"");
						strtok(0, ">"); // eat
					}
				}
				
				uint64 flag = 1 << flagsByFilter.size();
				flagsByFilter[filter.name] = flag;
				filter.flag = flag;
				//printf("FILTER:\t%d\t%s\t%s\n", filter->flag(), filter->name().c_str(), filter->description().c_str());
			}
			else if ( strncmp(line, "#CHROM", 6) == 0 )
			{
				TrackList::Track * track;
				char * token;
				int n = 0;
				
				strtok(line, "\t");
				
				for ( int i = 0; i < 8; i++ )
				{
					strtok(0, "\t"); // eat headers
				}
				
				while ( (token = strtok(0, "\t")) )
				{
					if ( oldTags )
					{
						track = &trackList->getTrackMutable(n); // TODO: clear track
						trackIndecesNew[trackList->getTrackIndexByFile(token)] = n;
						n++;
					}
					else
					{
						track = &trackList->getTrackMutable(trackList->addTrack(token));
					}
			
					track->file = token;
				}
			}
		}
		else
		{
			in.getline(line, (1 << 20) - 1);
			
			if ( in.eof() )
			{
				break;
			}
			
			variants.resize(variants.size() + 1);
			Variant & variant = variants[variants.size() - 1];
			
			variant.sequence = refByTag[strtok(line, "\t")];
			variant.position = atoi(strtok(0, "\t")) - 1;
			strtok(0, "\t"); // eat id
			char * alleles = strtok(0, "\t"); // ref allele
			strtok(0, "\t"); // eat alt alleles
			variant.quality = atoi(strtok(0, "\t"));
			
			uint64 filters = 0;
			char * filterString = strtok(0, "\t");
			
			if ( filterString[-1] == '\t' )
			{
				filterString--;
				*filterString = 0;
			}
			else
			{
				strtok(0, "\t"); // eat info
			}
			char * del = filterString;
			
			while ( del )
			{
				del = strchr(filterString, ':');
				
				if ( del )
				{
					*del = 0;
				}
				
				if ( *filterString )
				{
					filters |= flagsByFilter[filterString];
				}
				
				filterString = del + 1;
			}
			
			variant.filters = filters;
			
			strtok(0, "\t"); // eat format
			
			char * alleleString;
			variant.alleles.resize(alleleCount);
			alleleCount = 0;
			
			while ( (alleleString = strtok(0, "\t")) )
			{
				if ( variant.alleles.size() < alleleCount + 1 )
				{
					variant.alleles.resize(alleleCount + 1);
				}
				
				if ( *alleleString == '.' )
				{
					variant.alleles[alleleCount] = 'N';
				}
				else
				{
					variant.alleles[alleleCount] = alleles[atoi(alleleString) * 2];
				}
				
				alleleCount++;
			}
			
//			(*variant->mutable_alleles())[alleleCount] = '\0';
			
//			printf("VARIANT:\t%d\t%d\t%d\t%d\t%s\n", variant.sequence, variant.position, variant.filters, variant.quality, variant.alleles.c_str());
		}
	}
	
	if ( oldTags )
	{
		trackList->setTracksByFile();
		phylogenyTree->setTrackIndeces(trackIndecesNew);
		delete [] trackIndecesNew;
	}
	
	if ( lcbList->getLcbCount() == 0 )
	{
		lcbList->initWithSingleLcb(referenceList, *trackList);
	}
	
	in.close();
}

void VariantList::sortVariants()
{
	sort(variants.begin(), variants.end(), variantLessThan);
}

void VariantList::writeToMfa(std::ostream &out, bool indels, const TrackList & trackList) const
{
	int wrap = 80;
	int col;
	
	for ( int i = 0; i < trackList.getTrackCount(); i++ )
	{
		const TrackList::Track & track = trackList.getTrack(i);
		
		out << '>' << (track.file.length() ? track.file : track.name) << endl;
		col = 0;
		
		for ( int j = 0; j < variants.size(); j++ )
		{
			if ( ! indels && (variants.at(j).filters & FILTER_indel) )
			{
				continue;
			}
			
			col++;
			
			if ( wrap && col > wrap )
			{
				out << endl;
				col = 1;
			}
			
			out << variants.at(j).alleles[i];
		}
		
		out << endl;
	}
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

void VariantList::writeToVcf(std::ostream &out, bool indels, const ReferenceList & referenceList, const TrackList & trackList) const
{
	//tjt: Currently outputs SNPs, no indels
	//tjt: next pass will add standard VCF output for indels, plus an attempt at qual vals
	//tjt: also filters need to be added to findVariants to populate FILTer column

	//indel char, to skip columns with indels (for now)
	char indl = '-';
	//the VCF output file

	for ( int i = 0; i < filters.size(); i++ )
	{
		const Filter & filter = filters.at(i);
		
		out << "##FILTER=<ID=" << filter.name << ",Description=\"" << filter.description << "\">\n";
	}
	
	//the VCF header line (skipping previous lines for simplicity, can/will add in later)
	//#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  AA1 
	out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";


	//output the file name for each column
	for ( int i = 0; i < trackList.getTrackCount(); i++ )
	{
		const TrackList::Track & track = trackList.getTrack(i);
		//out << '\t' << (msgTrack.has_name() ? msgTrack.name() : msgTrack.file());
		out << '\t' << track.file;
	}
	
	out << '\n';
	
	//now iterate over variants and output
	for ( int j = 0; j < variants.size(); j++ )
	{
		const Variant & variant = variants.at(j);

		//no indels for now..
		if (variant.alleles[0] == indl)
			continue;
			
		if (find(variant.alleles.begin(), variant.alleles.end(),indl) != variant.alleles.end())
			continue;

		//capture the reference position of variant
		int pos = variant.position;

		//output first few columns, including context (+/- 7bp for now)
		int ws = 10;
		int lend = pos-ws;
		int rend = ws;
		
		const string & refseq = referenceList.getReference(variant.sequence).sequence;
		
		if (lend < 0)
			lend = 0;
		if (pos+ws >= refseq.size())
			rend = refseq.size()-pos;
		if (pos+rend >= refseq.size())
			rend = 0;
			
		out << variant.sequence + 1 << "\t" << pos + 1 << "\t" << refseq.substr(lend,ws) << "." << refseq.substr(pos,rend);

		//build non-redundant allele list from cur alleles
		vector<char> allele_list;
		//first allele is ref allele (0)
		allele_list.push_back(variant.alleles[0]);
		bool prev_var = false;
		for ( int i = 0; i < trackList.getTrackCount(); i++ )
		{
			if (find(allele_list.begin(), allele_list.end(), variant.alleles[i]) == allele_list.end())
			{
				if (variant.alleles[i] == indl) 
					continue;
					
				//to know if we need to output a preceding comma
				if (!prev_var)
				out << variant.alleles[i];
				else
				out << "," << variant.alleles[i];

				allele_list.push_back(variant.alleles[i]);
				prev_var = true;
			}
			//to see if we are in REF column
			else if (i == 0)
			{
				out << "\t" << variant.alleles[i] << "\t";
			}
		}
		
		//below values, punt for now, fill in with actual values later..
		//QUAL
		out << "\t40";

		//FILT
		//
		out << '\t';
		int filterCount = 0;
		//
		for ( int i = 0; i < filters.size(); i++ )
		{
			const Filter & filter = filters.at(i);
			
			if ( variant.filters & filter.flag )
			{
				if ( filterCount > 0 )
				{
					out << ':';
				}
				
				out << filter.name;
				filterCount++;
			}
		}
		//
		if ( filterCount == 0 )
		{
			out << "PASS";
		}
		
		//INFO
		out << "\tNA";
		//FORMAT
		out << "\tGT";

		//catch last one for newline
		int i = 0;
		
		map<char, int> indexByAllele;
		
		for ( int i = 0; i < allele_list.size(); i++ )
		{
			indexByAllele[allele_list[i]] = i;
		}
		
		for (i = 0; i < trackList.getTrackCount(); i++ )
		{
			out << "\t" << indexByAllele[variant.alleles[i]];
		}
		
		out << "\n";

	}
	//done! should be well-formated VCF (see above notes)
	//out.close();
}

void VariantList::addFilter(long long int flag, string name, string description)
{
	filters.resize(filters.size() + 1);
	filters[filters.size() - 1].flag = flag;
	filters[filters.size() - 1].name = name;
	filters[filters.size() - 1].description = description;
}
