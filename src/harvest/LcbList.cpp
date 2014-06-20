
#include "harvest/LcbList.h"
#include <iostream>
#include <fstream>
#include "harvest/parse.h"

using namespace::std;

bool lcbLessThan(const LcbList::Lcb a, const LcbList::Lcb & b)
{
	if ( a.sequence == b.sequence )
	{
		return a.position < b.position;
	}
	else
	{
		return a.sequence < b.sequence;
	}
}

void LcbList::initFromMfa(const char * file, ReferenceList * referenceList, TrackList * trackList, PhylogenyTree * phylogenyTree, VariantList * variantList)
{
	lcbs.resize(0);
	
	ifstream in(file);
	char line[1 << 20];
	vector<string> seqs;
	bool oldTags = trackList->getTrackCount();
	string refTag;
	int * trackIndecesNew;
	
	if ( oldTags )
	{
		trackIndecesNew = new int[trackList->getTrackCount()];
	}
	
	while ( ! in.eof() )
	{
		if ( in.peek() == '#' )
		{
			continue;
		}
		
		if ( in.peek() == '>' )
		{
			in.getline(line, (1 << 20) - 1);
			string tag(strtok(line + 1, " "));
			string desc(strtok(0, "\n"));
			
			if ( seqs.size() == 0 )
			{
				refTag = line + 1;
			}
			
			TrackList::Track * track;
			
			if ( oldTags )
			{
				track = &trackList->getTrackMutable(seqs.size());
				trackIndecesNew[trackList->getTrackIndexByFile(tag.c_str())] = seqs.size();
			}
			else
			{
				track = &trackList->getTrackMutable(trackList->addTrack(tag.c_str()));
			}
			
			track->file = tag;
			seqs.resize(seqs.size() + 1);
		}
		
		in.getline(line, (1 << 20) - 1);
		seqs[seqs.size() - 1].append(line);
	}
	
	in.close();
	
	if ( oldTags )
	{
		trackList->setTracksByFile();
		phylogenyTree->setTrackIndeces(trackIndecesNew);
		delete [] trackIndecesNew;
	}
	
	string ref = seqs[0];
	ungap(ref);
	
	referenceList->clear();
	referenceList->addReference(refTag, ref);
	
	lcbs.resize(1);
	LcbList::Lcb * lcb = &lcbs[0];
	lcb->position = 0;
	
	for ( int i = 0; i < trackList->getTrackCount(); i++ )
	{
		lcb->regions.resize(lcb->regions.size() + 1);
		LcbList::Region * region = &lcb->regions[lcb->regions.size() - 1];
		
		region->position = 0;
		region->length = ref.length();
		region->reverse = false;
	}
	
	vector<const VariantList::Variant *> vars;
	
	if ( variantList )
	{
		variantList->init();
		variantList->addVariantsFromAlignment(seqs, *referenceList, lcb->sequence, lcb->position, false);
		variantList->sortVariants();
	}
}

void LcbList::initFromProtocolBuffer(const Harvest::Alignment & msgAlignment)
{
	lcbs.resize(msgAlignment.lcbs_size());
	
	for ( int i = 0; i < msgAlignment.lcbs_size(); i++ )
	{
		Lcb & lcb = lcbs[i];
		const Harvest::Alignment::Lcb & msgLcb = msgAlignment.lcbs(i);
		
		lcb.sequence = msgLcb.sequence();
		lcb.position = msgLcb.position();
		lcb.length = msgLcb.length();
		lcb.concordance = msgLcb.concordance();
		lcb.regions.resize(msgLcb.regions_size());
		
		for ( int j = 0; j < msgLcb.regions_size(); j++ )
		{
			Region & region = lcb.regions[j];
			const Harvest::Alignment::Lcb::Region & msgRegion = msgLcb.regions(j);
			
			// TODO: track id?
			
			region.position = msgRegion.position();
			region.length = msgRegion.length();
			region.reverse = msgRegion.reverse();
		}
	}
}

void LcbList::initFromXmfa(const char * file, const ReferenceList & referenceList, TrackList * trackList, PhylogenyTree * phylogenyTree, VariantList * variantList)
{
	lcbs.resize(0);
	
	ifstream in(file);
	char line[1 << 20];
	int trackIndex = 0;
	vector<string> seqs;
	bool oldTags = trackList->getTrackCount();
	int * trackIndecesNew;
	
	if ( variantList )
	{
		variantList->init();
	}
	
	if ( oldTags )
	{
		trackIndecesNew = new int[trackList->getTrackCount()];
	}
	
	LcbList::Lcb * lcb;
	
	int lcbLength = 0;
	
	bool lcbfilt = false; // if lcb < 200bp filter SNPs inside
	TrackList::Track * track;
	
	while ( ! in.eof() )
	{
		in.getline(line, (1 << 20) - 1);
		
		if ( *line == '#' )
		{
			const char * suffix;
			
			if ( (suffix = removePrefix(line, "##SequenceFile ")) )
			{
				if ( oldTags )
				{
					track = &trackList->getTrackMutable(trackIndex);
					trackIndecesNew[trackList->getTrackIndexByFile(suffix)] = trackIndex;
				}
				else
				{
					track = &trackList->getTrackMutable(trackList->addTrack(suffix));
				}
				
				track->file = suffix;
				trackIndex++;
			}
			else if ( (suffix = removePrefix(line, "##SequenceHeader ")) )
			{
				track->name = suffix;
			}
			else if ( (suffix = removePrefix(line, "##SequenceLength ")) )
			{
				string length(suffix);
				string length_t(length.begin(),length.end()-2);
				track->size = atoi(length_t.c_str());
			}
		}
		else if ( *line == '>' )
		{
			trackIndex = atoi(strtok(line, ":") + 1) - 1;
			
			if ( trackIndex == 0 )
			{
				if ( variantList && lcbs.size() == 0 )
				{
					seqs.resize(trackList->getTrackCount());
				}
				
				lcbs.resize(lcbs.size() + 1);
				lcb = &lcbs[lcbs.size() - 1];
			}
			
			if ( trackIndex >= lcb->regions.size() )
			{
				int sizeOld = lcb->regions.size();
				lcb->regions.resize(trackIndex + 1);
				
				for ( int i = sizeOld; i < trackIndex; i++ )
				{
					lcb->regions[i].position = 0;
					lcb->regions[i].length = 0;
					lcb->regions[i].reverse = false;
				}
			}
			
			LcbList::Region * region = &lcb->regions[trackIndex];
			
			region->position = atoi(strtok(0, "-"));
			int end = atoi(strtok(0, " "));
			region->length = end - region->position + 1;
			region->reverse = *strtok(0, " ") == '-';
			
			if ( trackIndex == 0 )
			{
				lcb->sequence = referenceList.getReferenceSequenceFromConcatenated(region->position);
				lcb->position = referenceList.getPositionFromConcatenated(lcb->sequence, region->position);
			}
		}
		else if ( variantList && *line == '=' )
		{
			bool all = true;
			
			for ( int i = 0; i < seqs.size(); i++ )
			{
				if ( seqs[i].length() == 0 )
				{
					all = false;
					break;
				}
			}
			
			if ( all )
			{
				variantList->addVariantsFromAlignment(seqs, referenceList, lcb->sequence, lcb->position, lcbLength < 200 ? true : false);
			}
			
			for ( int i = 0; i < seqs.size(); i++ )
			{
				seqs[i].clear();
			}
		}
		else if ( variantList )
		{
			seqs[trackIndex].append(line);
		}
		if ( *line != '=' && *line != '>' && *line != '#')
		{
			if ( trackIndex == 0 )
			{
				lcbLength += strlen(line);
			}
		}
		else if (*line == '=')
		{
			lcb->length = lcbLength;
			lcbLength = 0;
		}
	}
	
	if ( oldTags )
	{
		trackList->setTracksByFile();
		phylogenyTree->setTrackIndeces(trackIndecesNew);
		delete [] trackIndecesNew;
	}
	
	sort(lcbs.begin(), lcbs.end(), lcbLessThan);
	
	if ( variantList )
	{
		variantList->sortVariants();
	}
	
	in.close();
}

void LcbList::writeToProtocolBuffer(Harvest * msg) const
{
	Harvest::Alignment * msgAlignment = msg->mutable_alignment();
	
	for ( int i = 0; i < lcbs.size(); i++ )
	{
		Harvest::Alignment::Lcb * msgLcb = msgAlignment->add_lcbs();
		const LcbList::Lcb & lcb = lcbs.at(i);
		
		msgLcb->set_sequence(lcb.sequence);
		msgLcb->set_position(lcb.position);
		msgLcb->set_length(lcb.length);
		msgLcb->set_concordance(lcb.concordance);
		
		for ( int j = 0; j < lcb.regions.size(); j++ )
		{
			// TODO: empty tracks?
			
			Harvest::Alignment::Lcb::Region * msgRegion = msgLcb->add_regions();
			const LcbList::Region & region = lcb.regions.at(j);
			
			msgRegion->set_track(j);
			msgRegion->set_position(region.position);
			msgRegion->set_length(region.length);
			msgRegion->set_reverse(region.reverse);
		}
	}
}

void LcbList::writeToXmfa(ostream & out, const ReferenceList & referenceList, const TrackList & trackList, const VariantList & variantList) const
{
/* EXAMPLE header
#FormatVersion MultiSNiP
#SequenceCount 8
##SequenceIndex 1
##SequenceFile b1.fna
##SequenceHeader >gi|76577973|gb|CP000124.1| Burkholderia pseudomallei 1710b chromosome I, complete sequence
##SequenceLength 4126292bp
*/
	out << "#FormatVersion ParSNP v1.0" << endl;
	out << "#SequenceCount " << trackList.getTrackCount() << endl;
	
	for ( int i = 0; i < trackList.getTrackCount(); i++ )
	{
		const TrackList::Track & track = trackList.getTrack(i);
		out << "##SequenceIndex " << i + 1 << endl;
		out << "##SequenceFile " << track.file << endl;
		out << "##SequenceHeader " << track.name << endl;
		
		if ( track.size )
		{
			out << "##SequenceLength " << track.size << "bp" << endl;
		}
	}
	
	out << "#IntervalCount " << lcbs.size() << endl;
	
	// now iterate over alignments
	
	int totrefgaps = 0;
	int currvar = 0;
	
	for ( int j = 0; j < lcbs.size(); j++ )
	{
		const LcbList::Lcb & lcb = lcbs.at(j);
		int refstart = lcb.regions.at(0).position;
		int refend = 0;
		int blockVarStart;
		
		for ( int r = 0; r < lcb.regions.size(); r++)
		{
			// >1:8230-11010 + cluster174 s1:p8230
			
			const LcbList::Region & region = lcb.regions.at(r);
			int start = region.position + 1;
			int end = start + region.length - 1;
			
			if (r == 0)
			{
				blockVarStart = currvar;
			}
			else
			{
				currvar = blockVarStart;
			}
			
			out << ">" << r+1 << ":" << start << "-" << end << " ";
			
			if ( ! region.reverse )
			{
				out << "+ ";
			}
			else
			{
				out << "- ";
			}
			
			out << "cluster" << j + 1 << endl;
			int currpos = refstart;
			int width = 80;
			int col = 0;
			int variantsSize = variantList.getVariantCount();
			const VariantList::Variant * currvarref = &variantList.getVariant(currvar);
			
			//var =  harvest.variation().variants(currvar);//.alleles()[r];
			//string ref_slice(harvest.reference().references(0).sequence().substr(refstart,(end-start)+1));
			
			while
			(
				currpos - refstart < lcb.regions.at(0).length ||
				(
					currvar < variantsSize &&
					currvarref->position - refstart < lcb.regions.at(0).length
				)
			)
			{
				if ( currpos == referenceList.getReference(0).sequence.size() )
				{
					printf("ERROR: LCB %d extends beyond reference (position %d)\n", j, currpos);
				}
				
				if ( col == width )
				{
					out << endl;
					col = 0;
				}
				
				if
				(
					currvar == variantsSize ||
					currpos != currvarref->position ||
					(
						currvarref->alleles[0] == '-' &&
						currvar > 0 &&
						variantList.getVariant(currvar - 1).position != currpos
					)
				)
				{
					out << referenceList.getReference(0).sequence.substr(currpos, 1);
					col++;
				}
				
				if ( currvar < variantsSize && currpos == currvarref->position )
				{
					out << currvarref->alleles[r];
					currvar++;
					
					if ( currvar < variantsSize )
					{
						currvarref = &variantList.getVariant(currvar);
					}
					
					col++;
				}
				
				if ( currvar == variantsSize || currvarref->position > currpos )
				{
					currpos++;
				}
			}
			
			out << endl;
		}
		
		out << "=" << endl;
	}
}
