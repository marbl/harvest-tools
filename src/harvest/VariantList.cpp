// Copyright © 2014, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "harvest/VariantList.h"
#include <fstream>
#include <sstream>
#include "harvest/parse.h"
#include <set>
#include <algorithm>

using namespace::std;
using namespace::flatbuffers;

bool operator<(const VariantList::VariantSortKey & a, const VariantList::VariantSortKey & b)
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

void VariantList::addFilterFromBed(const char * file, const char * name, const char * desc)
{
	ifstream in(file);
	char * line = new char[1 << 20];
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
	
	delete [] line;
}

void VariantList::addVariantsFromAlignment(const vector<string> & seqs, const ReferenceList & referenceList, int sequence, int position, int length, bool reverse)
{
//	Harvest::Variation * msg = harvest.mutable_variation();
	char col[seqs.size() + 1];
        //add arrays for tracking conserved,poorly aligned columns
	vector<bool> conserved(seqs[0].length()+1,true);
	vector<bool> gaps(seqs[0].length()+1,false);
	vector<bool> nns(seqs[0].length()+1,false);
	int offset = 0;
	
	col[seqs.size()] = 0; // null-terminate for use as a c-style string
	
        //simple loop to check for column conservation
        //this could be done via SP-score all-v-all pairs
        //but for now, simply use to flag SNPs that are within a window of 100bp
        //with less than 50% column conservation (w.r.t ref, not consensus)
	for ( int i = 0; i < seqs[0].length(); i++ )
	{
		bool variant = false;
		bool indel = false;
		vector<int> nt_cnt(5,0);
		//vector<int>::iterator maxval;
		
		for (int j = 0; j < seqs.size(); j++)
		{
			if ( reverse )
			{
				col[j] = seqs[j][seqs[0].length() - i - 1];
			}
			else
			{
				col[j] = seqs[j][i];
			}

  		    if (col[j] == 'a' || col[j] == 'A')
		      nt_cnt[0] =1;
  		    else if (col[j] == 't' || col[j] == 'T')
		      nt_cnt[1] =1;
  		    else if (col[j] == 'g' || col[j] == 'G')
		      nt_cnt[2] =1;
  		    else if (col[j] == 'c' || col[j] == 'C')
		      nt_cnt[3] =1;
  		    else if (col[j] == 'n' || col[j] == 'N')
		      nns[i] = true;
  		    else if (col[j] == '-')
		    {
		      gaps[i] = true;
                      nt_cnt[4] = 1;
		    }
                }
                //maxval = std::max_element(nt_cnt.begin(),nt_cnt.end());
                if ((nt_cnt[0] + nt_cnt[1] + nt_cnt[2] + nt_cnt[3] +nt_cnt[4]) > 1)
                  conserved[i] = false;
                /*multi-allelic
                if ((nt_cnt[0] + nt_cnt[1] + nt_cnt[2] + nt_cnt[3] +nt_cnt[4]) > 2)
		{
		  conserved[i] = false;
		}
                */
	}
	
	// Since insertions to the reference take on the left-most reference
	// position, this allows the alignment to start with an insertion
	// (possibly at reference position -1).
	//
	position--;
	
	for ( int i = 0; i < seqs[0].length(); i++ )
	{
		bool variant = false;
		bool n = false;
		
		if ( reverse )
		{
			col[0] = seqs[0][seqs[0].length() - i - 1];
		}
		else
		{
			col[0] = seqs[0][i];
		}
		
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
			if ( reverse )
			{
				col[j] = seqs[j][seqs[0].length() - i - 1];
			}
			else
			{
				col[j] = seqs[j][i];
			}
		
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
			
			int windowsize = 0;
                        int window = 50;
                        if (window > i)
			{
			  window = i -1;
			}
                        windowsize+=window;
                        int conserved_cnt = 0;
                        int gap_cnt = 0;
                        for (int z = 1; z<=window;z++)
			{
			  if (conserved.at(i-z))
			  {
			    conserved_cnt+=1;
			  }

			  if (gaps.at(i-z))
			  {
			    gap_cnt+=1;
			  }
            
			}
                        window = 50;
                        if (window+i > seqs[0].length())
			{
			  window = seqs[0].length() - i;
			}
                        windowsize+=window;
                        for (int z = 1; z<=window;z++)
			{
			  if (conserved.at(i+z))
			  {
			    conserved_cnt+=1;
			  }
			  if (gaps.at(i+z))
			  {
			    gap_cnt+=1;
			  }
			}
			
			if ( reverse )
			{
				for ( int j = 0; j < seqs.size(); j++ )
				{
					col[j] = complement(col[j]);
				}
			}
			
			varNew->sequence = sequence;
			varNew->position = position;
			varNew->offset = offset;
			
			if ( referenceList.getReferenceCount() )
			{
				varNew->reference = referenceList.getReference(sequence).sequence[position];
			}
			else
			{
				varNew->reference = col[0];
			}
			
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
			
			if ( length < 200 )
			{
				varNew->filters |= FILTER_lcb;
			}
			
                        if ( ((float)conserved_cnt/(float)windowsize) < 0.5 )
			{
				varNew->filters |= FILTER_conservation;
			}
                        if ( ((float)gap_cnt/(float)windowsize) > 0.2 )
			{
				varNew->filters |= FILTER_gaps;
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
	addFilter(FILTER_conservation, "CID", "SNP in aligned 100bp window with < 50% column % ID");
	addFilter(FILTER_gaps, "ALN", "SNP in aligned 100b window with > 20 indels");
	
	variants.resize(0);
}

void VariantList::initFromFlatBuffers(const fbHarvest::VariantList * pointerVariantList)
{
	filters.resize(pointerVariantList->filters()->Length());
	
	for ( int i = 0; i < filters.size(); i++ )
	{
		const fbHarvest::Filter * pointerFilter = pointerVariantList->filters()->Get(i);
		
		filters[i].flag = pointerFilter->flag();
		filters[i].name = pointerFilter->name()->c_str();
		filters[i].description = pointerFilter->description()->c_str();
		printf("FILTER:\t%d\t%s\t%s\n", filters[i].flag, filters[i].name.c_str(), filters[i].description.c_str());
	}
	
	variants.resize(pointerVariantList->variants()->Length());
	
	for ( int i = 0; i < variants.size(); i++ )
	{
		Variant & variant = variants[i];
		const fbHarvest::Variant * pointerVariant = pointerVariantList->variants()->Get(i);
		
		variant.sequence = pointerVariant->sequence();
		variant.position = pointerVariant->position();
		variant.alleles = pointerVariant->alleles()->c_str();
		variant.filters = pointerVariant->filters();
		variant.quality = pointerVariant->quality();
		
		printf("VARIANT: %d\t%d\t%s\t%ld\t%d\n", variant.sequence, variant.position, variant.alleles.c_str(), variant.filters, variant.quality);
	}
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
		
		if ( msgVariant.has_reference() )
		{
			variant.reference = msgVariant.reference();
		}
		else
		{
			variant.reference = variant.alleles[0];
		}
	}
}

void VariantList::initFromVcf(const char * file, const ReferenceList & referenceList, TrackList * trackList, LcbList * lcbList, PhylogenyTree * phylogenyTree)
{
	filters.resize(0);
	variants.resize(0);
	
	ifstream in(file);
	
	// Since we will be transposing multi-base alleles to our column-based
	// representation, we will refer to columns multiple times and will use a
	// map to look up existing columns efficiently.
	//
	map<VariantSortKey, int> variantIndecesBySortKey;
	
	// Insertions where any allele inserted more than one base are ambiguous and
	// will be replaced with an LCB boundary; also, insertions or deletions with
	// missing ('.') alleles are considered non-core and removed. This map keeps
	// track of such cases for reference during transposition and for creating
	// LCBs later. For these keys, offset is 0 for deletions and 1 for
	// insertions; this determines whether the base itself is included in an LCB
	//
	set<VariantSortKey> ambiguousIndels;
	
	string line;
	map<string, long long int> flagsByFilter;
	map<string, int> refByTag;
	//unsigned int alleleCount = 0;
	
	const bool oldTags = phylogenyTree->getRoot();
	int * trackIndecesNew;
	int lineIndex = 1;
	
	if ( oldTags )
	{
		trackIndecesNew = new int[trackList->getTrackCount()];
	}
	else
	{
		trackList->clear();
	}
	
	for ( int i = 0; i < referenceList.getReferenceCount(); i++ )
	{
		refByTag[referenceList.getReference(i).name] = i;
	}
	
	while ( getline(in, line) )
	{
		if ( line[0] == '#' )
		{
			if ( strncmp(line.c_str(), "##FILTER", 8) == 0 )
			{
				char * token;
				
				filters.resize(filters.size() + 1);
				Filter & filter = filters[filters.size() - 1];
				
				size_t pos = line.find("ID=", 10);
				
				if ( pos != string::npos )
				{
					pos += 3;
					size_t end = line.find_first_of(",>", pos);
					filter.name = line.substr(pos, end - pos);
				}
				
				pos = line.find("Description=", 10);
				
				if ( pos != string::npos )
				{
					pos += 13;
					size_t end = line.find_first_of("\"", pos);
					filter.description = line.substr(pos, end - pos);
				}
				
				uint64 flag = 1 << flagsByFilter.size();
				flagsByFilter[filter.name] = flag;
				filter.flag = flag;
				//printf("FILTER:\t%d\t%s\t%s\n", filter->flag(), filter->name().c_str(), filter->description().c_str());
			}
			else if ( strncmp(line.c_str(), "#CHROM", 6) == 0 )
			{
				TrackList::Track * track;
				int n = 0;
				stringstream lineStream(line);
				string field;
				
				// eat headers
				//
				for ( int i = 0; i < 9; i++ )
				{
					lineStream >> field;
				}
				
				// get names
				//
				while ( (lineStream >> field) )
				{
					if ( oldTags )
					{
						track = &trackList->getTrackMutable(n); // TODO: clear track
						
						try
						{
							trackIndecesNew[trackList->getTrackIndexByFile(field)] = n;
						}
						catch ( const TrackList::TrackNotFoundException & e )
						{
							delete [] trackIndecesNew;
							throw;
							return;
						}
						
						n++;
					}
					else
					{
						track = &trackList->getTrackMutable(trackList->addTrack(field));
					}
			
					track->file = field;
				}
			}
		}
		else
		{
			stringstream lineStream(line);
			
			string refName;
			int position;
			string eaten;
			string ref;
			string altAlleles;
			float quality;
			string info;
			string filterString;
			int offset = 0;
			
			lineStream >> refName >> position >> eaten >> ref >> altAlleles >> quality >> filterString >> eaten >> eaten;
			int sequence = refByTag[refName];
			position--;
			
			vector<string> alleleStrings;
			string alleleString;
			stringstream alleleStream(altAlleles);
			
			while ( getline(alleleStream, alleleString, ',') )
			{
				alleleStrings.push_back(alleleString);
			}
			
			uint64 filters = 0;
			stringstream filterStream(filterString);
			
			while ( getline(filterStream, filterString, ':') )
			{
				if ( filterString.compare(".") != 0 && filterString.compare("PASS") != 0 )
				{
					filters |= flagsByFilter[filterString];
				}
			}
			
			string alleleIndex;
			vector<int> alleleIndeces;
			bool missing = false;
			
			while ( lineStream >> alleleIndex )
			{
				if ( alleleIndex[0] == '.' )
				{
					missing = true;
					alleleIndeces.push_back(-1);
				}
				else
				{
					alleleIndeces.push_back(atoi(alleleIndex.c_str()));
				}
			}
			
			for ( int i = 0; i < alleleStrings.size(); i++ )
			{
				if ( alleleStrings[i].find_first_of("<>[]*X") != string::npos )
				{
					// we don't yet handle symbolic alleles, breakends, or other
					// weird stuff
					
					continue;
				}
				
				if ( alleleStrings[i].length() != ref.length() )
				{
					if ( alleleStrings[i][0] != ref[0] )
					{
						throw CompoundVariantException(lineIndex);
					}
				}
				
				int lengthVariant;
				
				if ( alleleStrings[i].length() > ref.length() )
				{
					lengthVariant = alleleStrings[i].length();
				}
				else
				{
					lengthVariant = ref.length();
				}
				
				for ( int j = 0; j < lengthVariant; j++ )
				{
					if ( j < ref.length() && j < alleleStrings[i].length() && alleleStrings[i].at(j) == ref.at(j) )
					{
						continue;
					}
					
					int positionVariant;
					int offset;
					
					if ( j >= ref.length() )
					{
						positionVariant = position + ref.length() - 1;
						offset = j - ref.length() + 1;
					}
					else
					{
						positionVariant = position + j;
						offset = 0;
					}
					
					if ( offset > 0 )
					{
						if ( ambiguousIndels.count(VariantSortKey(sequence, position + ref.length() - 1, 1)) )
						{
							// another variant tried to insert more than one
							// base here; this is now ambiguous
							
							break;
						}
					}
					
					if ( offset > 1 || ( offset == 1 && missing) )
					{
						// insertions of more than one base become ambiguous;
						// replace with LCB boundary, destroy any single base
						// insertions at this spot, and prevent more
						
						VariantSortKey key(sequence, position + ref.length() - 1, 1);
						
						if ( variantIndecesBySortKey.count(key) )
						{
							variants.erase(variants.begin() + variantIndecesBySortKey.at(key));
						}
						
						ambiguousIndels.insert(key);
						
						break;
					}
					
					VariantSortKey key(sequence, positionVariant, offset);
					Variant * variant;
					
					if ( ambiguousIndels.count(key) )
					{
						// ambiguous deletion here; no variants allowed
						
						continue;
					}
					
					if ( missing && j >= alleleStrings[i].length() )
					{
						// ambiguous deletion; destroy any variants at this base
						// (including insertions) and prevent more
						
						if ( variantIndecesBySortKey.count(key) )
						{
							variants.erase(variants.begin() + variantIndecesBySortKey.at(key));
						}
						
						ambiguousIndels.insert(key);
						
						VariantSortKey keyInsertion(sequence, positionVariant, 1);
						
						if ( variantIndecesBySortKey.count(keyInsertion) )
						{
							variants.erase(variants.begin() + variantIndecesBySortKey.at(keyInsertion));
						}
						
						ambiguousIndels.insert(keyInsertion);
						
						continue;
					}
					
					if ( variantIndecesBySortKey.count(key) )
					{
						// existing variant at this column
						
						variant = & variants[variantIndecesBySortKey.at(key)];
						
						// use the minimum quality to be conservative
						//
						if ( quality < variant->quality )
						{
							variant->quality = quality;
						}
						
						// use the union of the filters
						//
						variant->filters |= filters;
					}
					else
					{
						variantIndecesBySortKey[key] = variants.size();
						variants.resize(variants.size() + 1);
						variant = & variants[variants.size() - 1];
						
						if ( offset )
						{
							variant->reference = '-';
						}
						else
						{
							variant->reference = ref.at(j);
						}
						
						variant->sequence = sequence;
						variant->position = positionVariant;
						variant->offset = offset;
						variant->quality = quality;
						variant->filters = filters;
						variant->alleles.resize(trackList->getTrackCount(), 0);
					}
					
					char snp;
					
					if ( j < alleleStrings[i].length() )
					{
						snp = alleleStrings[i].at(j);
					}
					else
					{
						snp = '-';
					}
					
					for ( int k = 0; k < alleleIndeces.size(); k++ )
					{
						if ( alleleIndeces[k] - 1 == i || alleleIndeces[k] == -1 )
						{
							// we only set alternate bases, since reference alleles
							// might not reflect other variants
							
							char snpAllele = alleleIndeces[k] == -1 ? 'N' : snp;
							
							if ( variant->alleles[k] != 0 && variant->alleles[k] != snpAllele)
							{
								throw ConflictingVariantException
								(
									lineIndex,
									trackList->getTrack(k).file,
									variant->alleles[k],
									snpAllele
								);
							}
							
							variant->alleles[k] = snpAllele;
						}
					}
				}
			}
		}
		
		lineIndex++;
	}
	
	// since indel and snp changes can be cumulative in VCF, we only set
	// alternate alleles above and will now fill in any missing values with
	// their reference bases
	//
	for ( int i = 0; i < variants.size(); i++ )
	{
		for ( int j = 0; j < trackList->getTrackCount(); j++ )
		{
			if ( variants.at(i).alleles.at(j) == 0 )
			{
				variants[i].alleles[j] = variants.at(i).reference;
			}
		}
	}
	
	sortVariants();
	
	if ( oldTags )
	{
		trackList->setTracksByFile();
		phylogenyTree->setTrackIndeces(trackIndecesNew);
		delete [] trackIndecesNew;
	}
	
	if ( lcbList->getLcbCount() == 0 )
	{
		// use ambiguous indels as breakpoints for LCBs
		
		VariantSortKey keyLast(0, 0, 0);
		set<VariantSortKey>::iterator key = ambiguousIndels.begin();
		
		while ( true )
		{
			int endSeq;
			int endPos;
			
			if ( key == ambiguousIndels.end() )
			{
				endSeq = referenceList.getReferenceCount() - 1;
				endPos = referenceList.getReference(endSeq).sequence.length() - 1;
			}
			else
			{
				endSeq = key->sequence;
				
				if ( key->offset == 0 )
				{
					endPos = key->position - 1;
				}
				else
				{
					endPos = key->position;
				}
			}
			
			lcbList->addLcbByReference(keyLast.sequence, keyLast.position, endSeq, endPos, referenceList, *trackList);
			
			if ( key == ambiguousIndels.end() )
			{
				break;
			}
			
			// increment, skipping runs of adjacent deletions
			//
			do
			{
				keyLast.sequence = key->sequence;
				keyLast.position = key->position + 1; // next lcb should start at next base
				keyLast.offset = key->offset;
				
				if
				(
					keyLast.position == referenceList.getReference(key->sequence).sequence.length() &&
					key->sequence < referenceList.getReferenceCount() - 1
				)
				{
					// roll over to next sequence TODO: error?
					
					keyLast.sequence++;
					keyLast.position = 0;
				}
				
					
				key++;
			}
			while
			(
				key != ambiguousIndels.end() &&
				key->sequence == keyLast.sequence && 
				key->position <= keyLast.position &&
				
				// allow 1-base LCB for consecutive ambiguous insertions
				//
				(key->offset == 0 || keyLast.offset == 0)
			);
		}
	}
	
	in.close();
}

void VariantList::sortVariants()
{
	sort(variants.begin(), variants.end(), variantLessThan);
}

Offset<fbHarvest::VariantList> VariantList::writeToFlatBuffers(FlatBufferBuilder & fbb) const
{
	Offset<fbHarvest::Filter> * offsetFilters = new Offset<fbHarvest::Filter>[filters.size()];
	
	for ( int i = 0; i < filters.size(); i++ )
	{
		fbHarvest::FilterBuilder filterBuilder(fbb);
		
		filterBuilder.add_flag(filters[i].flag);
		filterBuilder.add_name(fbb.CreateString(filters[i].name));
		filterBuilder.add_description(fbb.CreateString(filters[i].description));
		
		offsetFilters[i] = filterBuilder.Finish();
	}
	
	delete [] offsetFilters;
	
	Offset<fbHarvest::Variant> * offsetVariants = new Offset<fbHarvest::Variant>[variants.size()];
	
	for ( int i = 0; i < variants.size(); i++ )
	{
		fbHarvest::VariantBuilder variantBuilder(fbb);
		
		variantBuilder.add_sequence(variants[i].sequence);
		variantBuilder.add_reference(variants[i].reference);
		variantBuilder.add_position(variants[i].position);
		variantBuilder.add_alleles(fbb.CreateString(variants[i].alleles));
		variantBuilder.add_filters(variants[i].filters);
		
		offsetVariants[i] = variantBuilder.Finish();
	}
	
	delete [] offsetVariants;
	
	return fbHarvest::CreateVariantList
	(
		fbb,
		fbb.CreateVector(offsetFilters, filters.size()),
		fbb.CreateVector(offsetVariants, variants.size())
	);
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
			if ( ! indels && variants.at(j).filters && variants.at(j).filters != FILTER_n )
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
		variant->set_reference(variants[i].reference);
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
			
		out << referenceList.getReference(variant.sequence).name << "\t" << pos + 1 << "\t" << refseq.substr(lend,ws) << "." << refseq.substr(pos,rend);

		//build non-redundant allele list from cur alleles
		vector<char> allele_list;
		//first allele is ref allele (0)
		out << "\t" << variant.reference << "\t";
		allele_list.push_back(variant.reference);
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
		}
		
		//below values, punt for now, fill in with actual values later..
		//QUAL
		if ( variant.quality != 0 )
		{
			out << '\t' << variant.quality; // currently only exists if imported from VCF
		}
		else
		{
			out << "\t40";
		}

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
