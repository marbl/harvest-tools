// Copyright © 2014, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "harvest/LcbList.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "harvest/parse.h"
#include <set>
#include <stdlib.h>
#include "harvest/exceptions.h"
#include "harvest/VariantList.h"
#include <algorithm>
#include <limits>

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

double LcbList::getCoreSize(void) const
{
    double coreSize = 0;
	
    //LcbList::Lcb * lcb = &lcbs[lcbs.size() - 1];
    
    for ( int i = 0; i < lcbs.size(); i++ )
    {
		const LcbList::Lcb & lcb = lcbs.at(i);
		
		coreSize += lcb.regions.at(0).length;
    }
    
    return coreSize;
}

bool operator<(const LcbList::Interval & a, const LcbList::Interval & b)
{
	if ( a.sequence == b.sequence )
	{
		if ( a.start == b.start )
		{
			return a.end < b.end;
		}
		else
		{
			return a.start < b.start;
		}
	}
	else
	{
		return a.sequence < b.sequence;
	}
}

void LcbList::addLcbByReference(int startSeq, int startPos, int endSeq, int endPos, const ReferenceList & referenceList, const TrackList & trackList)
{
	// add an LCB that has the same position in all queries as in the reference
	
	int startConcat = referenceList.getConcatenatedPosition(startSeq, startPos);
	int endConcat = referenceList.getConcatenatedPosition(endSeq, endPos);
	int length = endConcat - startConcat + 1;
	
	lcbs.resize(lcbs.size() + 1);
	LcbList::Lcb * lcb = &lcbs[lcbs.size() - 1];
	
	lcb->sequence = startSeq;
	lcb->position = startPos;
	lcb->length = length;
	lcb->regions.resize(trackList.getTrackCount());
	
	for ( int i = 0; i < trackList.getTrackCount(); i++ )
	{
		LcbList::Region * region = &lcb->regions[i];
		
		region->position = startConcat;
		region->length = length;
		region->reverse = false;
	}
}

void LcbList::clear()
{
	lcbs.clear();
}

void LcbList::initFromCapnp(const capnp::Harvest::Reader & harvestReader)
{
	auto lcbListReader = harvestReader.getLcbList();
	auto lcbsReader = lcbListReader.getLcbs();
	
	lcbs.resize(lcbsReader.size());
	
	for ( int i = 0; i < lcbsReader.size(); i++ )
	{
		Lcb & lcb = lcbs[i];
		auto lcbReader = lcbsReader[i];
		
		lcb.sequence = lcbReader.getSequence();
		lcb.position = lcbReader.getPosition();
		lcb.length = lcbReader.getLength();
		lcb.concordance = lcbReader.getConcordance();
		
		auto regionsReader = lcbReader.getRegions();
		lcb.regions.resize(regionsReader.size());
		
		for ( int j = 0; j < regionsReader.size(); j++ )
		{
			Region & region = lcb.regions[j];
			auto regionReader = regionsReader[j];
			
			// TODO: track id?
			
			region.position = regionReader.getPosition();
			region.length = regionReader.getLength();
			region.reverse = regionReader.getReverse();
		}
	}
}

void LcbList::initFromMaf(const char * file, ReferenceList * referenceList, TrackList * trackList, PhylogenyTree * phylogenyTree, VariantList * variantList, const char * referenceFileName)
{
	lcbs.resize(0);
	
	ifstream in(file);
	string line;
	int trackCount = 0;
	vector<string> seqs;
	const bool oldTags = phylogenyTree->getRoot();
	int * trackIndecesNew;
	bool mauve = false;
	vector<string> refs;
	vector<string> refNames;
	map<string, int> refIndexByName;
	vector<map<string, int> > regionOffsetBySeqNameByTrack;
	vector<int> totalOffsetByTrack;
	
	set<Interval> lcbIntervals;
	
	bool createReference = referenceList->getReferenceCount() == 0;
	
	if ( variantList )
	{
		variantList->init();
	}
	
	if ( oldTags )
	{
		trackIndecesNew = new int[trackList->getTrackCount()];
	}
	else
	{
		trackList->clear();
	}
	
	if ( referenceFileName )
	{
		char referenceBaseName[strlen(referenceFileName)];
		
		for ( const char * i = referenceFileName; *i != 0; i++ )
		{
			if ( *i == '/' )
			{
				referenceFileName = i + 1;
			}
		}
		
		strcpy(referenceBaseName, referenceFileName);
		strtok(referenceBaseName, ".");
		
		int trackIndex;
		
		try
		{
			trackIndex = trackList->getTrackIndexByFile(referenceBaseName);
		}
		catch ( const TrackList::TrackNotFoundException & e )
		{
			if ( oldTags )
			{
				delete [] trackIndecesNew;
				throw;
			}
			else
			{
				trackIndex = trackList->addTrack(referenceBaseName);
			}
		}
		
		if ( oldTags )
		{
			trackIndecesNew[trackIndex] = trackCount;
			trackCount++;
		}
		
	}
	
	LcbList::Lcb * lcb = 0;
	bool lcbReverse;
	
	TrackList::Track * track;
	char * suffix;
	int queryCount = 0;
	int lcbQueryCount = 0;
	vector<int> queryCountsByLcb;
	set<string> seqNames;
	
	// scan through once to determine number of query sequences...
	
	if ( ! in.is_open() )
	{
		cerr << "ERROR: " << file << " could not be opened.";
		return;
	}
	
	while ( getline(in, line) )
	{
		if ( line[0] == 's' )
		{
			stringstream lineStream(line);
			
			lineStream.ignore(2);
			
			string trackName;
			getline(lineStream, trackName, '.');
			seqNames.insert(trackName);
			
			lcbQueryCount++;
		}
		else if ( line[0] == 0 )
		{
			queryCountsByLcb.push_back(lcbQueryCount);
			lcbQueryCount = 0;
		}
	}
	
	queryCount = seqNames.size();
	regionOffsetBySeqNameByTrack.resize(queryCount);
	totalOffsetByTrack.resize(queryCount, 0);
	
	// ...now parse only those that are core
	
	in.clear();
	in.seekg(0);
	int lcbIndex = 0;
	
	while ( getline(in, line) )
	{
		if ( line[0] == 'a' )
		{
			// new alignment block; first check if it's core
			
			if ( queryCountsByLcb[lcbIndex] == queryCount )
			{
				// core; create a new Lcb
				
				if ( variantList && lcbs.size() == 0 )
				{
					seqs.resize(queryCount);
				}
			
				lcbs.resize(lcbs.size() + 1);
				lcb = &lcbs[lcbs.size() - 1];
			}
			else
			{
				// not core; eat it
				
				while ( in.peek() != '\n' )
				{
					in.ignore(numeric_limits<streamsize>::max(), '\n');
				}
				
				lcb = 0;
			}
			
			lcbIndex++;
		}
		else if ( line[0] == 's' )
		{
			stringstream lineStream(line);
			
			lineStream.ignore(2);
			
			string trackName;
			getline(lineStream, trackName, '.');
			
			int trackIndex;
			
			// create track if name is new
			
			try
			{
				trackIndex = trackList->getTrackIndexByFile(trackName);
			}
			catch ( const TrackList::TrackNotFoundException & e )
			{
				if ( oldTags )
				{
					delete [] trackIndecesNew;
					throw;
					return;
				}
				else
				{
					trackIndex = trackList->addTrack(trackName);
				}
			}
			
			if ( oldTags && trackCount < trackList->getTrackCount() )
			{
				trackIndecesNew[trackIndex] = trackCount;
				trackCount++;
			}
			
			track = &trackList->getTrackMutable(trackIndex);
			
			// basename parsing; not sure if this is needed for MAF
			//
/*					char * name = suffix;
			
			for ( int i = 0; i < strlen(suffix) - 1; i++ )
			{
				if ( suffix[i] == '/' )
				{
					name = suffix + i + 1;
				}
			}
*/					
			track->file = trackName;
			
			// parse positional info
			
			string seqName;
			lineStream >> seqName;
			
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
			
			int position;
			int length;
			char reverseChar;
			int seqLength;
			
			lineStream >> position >> length >> reverseChar >> seqLength;
			
			bool reverse = reverseChar == '-';
			
			if ( reverse )
			{
				position = seqLength - position - length; // MAF is stupid
			}
			
			if ( regionOffsetBySeqNameByTrack[trackIndex].count(seqName) == 0 )
			{
				// new sequence for this track; it's offset is the current total
				
				regionOffsetBySeqNameByTrack[trackIndex][seqName] = totalOffsetByTrack[trackIndex];
				totalOffsetByTrack[trackIndex] += seqLength;
			}
			
			int refIndex;
			
			if ( trackIndex == 0 )
			{
				// translate ref seq name to index and create if needed
				
				try
				{
					// check our local cache of names
					
					refIndex = refIndexByName.at(seqName);
				}
				catch ( const out_of_range & e )
				{
					// not in local cache of names
					
					if ( createReference )
					{
						// create ref and cache its index
						
						refIndex = refs.size();
						refs.resize(refs.size() + 1);
						refNames.resize(refNames.size() + 1);
						refNames[refIndex] = seqName;
						refs[refIndex].resize(seqLength, 'N');
						
						refIndexByName[seqName] = refIndex;
					}
					else
					{
						// reference should already exist
						
						try
						{
							refIndex = referenceList->getReferenceSequenceFromName(seqName);
						}
						catch ( ReferenceList::NameNotFoundException & e )
						{
							// reference doesn't exist; error
							throw;
						}
						
						// cache for faster lookup next time
						
						refIndexByName[seqName] = refIndex;
					}
				}
				
				set<Interval>::iterator lowerBound = lcbIntervals.lower_bound(Interval(refIndex, position, position + length - 1));
				bool overlap = false;
				
				if ( lowerBound != lcbIntervals.end() )
				{
					if ( refIndex == lowerBound->sequence && position + length - 1 >= lowerBound->start )
					{
						overlap = true;
					}
					
					lowerBound--;
					
					if ( ! overlap && lowerBound != lcbIntervals.begin() )
					{
						if ( refIndex == lowerBound->sequence && position <= lowerBound->end )
						{
							overlap = true;
						}
					}
				}
				
				if ( overlap )
				{
					// destroy lcb and eat alignment
					
					lcbs.resize(lcbs.size() - 1);
					
					while ( in.peek() != '\n' )
					{
						in.ignore(numeric_limits<streamsize>::max(), '\n');
					}
					
					lcb = 0;
					
					continue;
				}
				
				lcbIntervals.insert(Interval(refIndex, position, position + length - 1));
				
				lcb->sequence = refIndex;
				lcb->position = position;
				lcbReverse = reverse;
			}
			
			LcbList::Region * region = &lcb->regions[trackIndex];
			
			region->position = position + regionOffsetBySeqNameByTrack[trackIndex][seqName];
			region->length = length;
			region->reverse = reverse;
			
			// parse sequence
			
			if ( trackIndex == 0 || variantList )
			{
				string seq;
				lineStream >> seq;
				
				if ( createReference && trackIndex == 0 )
				{
					string ungapped = seq;
					ungap(ungapped);
					
					if ( lcbReverse )
					{
						reverseComplement(ungapped);
					}
					
					refs[refIndex].replace(lcb->position, ungapped.length(), ungapped);
				}
				
				lcb->length = seq.length();
				
				if ( variantList )
				{
					seqs[trackIndex] = seq;
				}
			}
		}
		else if ( variantList && line[0] == 0 && lcb != 0 )
		{
			for ( int i = 0; i < seqs.size(); i++ )
			{
				for ( int j = 0; j < seqs[i].length(); j++ )
				{
					seqs[i][j] = toupper(seqs[i][j]);
				}
			}
			
			variantList->addVariantsFromAlignment(seqs, *referenceList, lcb->sequence, lcb->position, lcb->length, lcbReverse);
			
			for ( int i = 0; i < seqs.size(); i++ )
			{
				seqs[i].clear();
			}
			
			lcb = 0;
		}
	}
	
	if ( queryCount && lcbs.size() == 0 )
	{
		throw NoCoreException(queryCount);
	}
	
	if ( oldTags )
	{
		trackList->setTracksByFile();
		phylogenyTree->setTrackIndeces(trackIndecesNew);
		delete [] trackIndecesNew;
	}
	
	if ( createReference )
	{
		for ( int i = 0; i < refNames.size(); i++ )
		{
			referenceList->addReference(refNames.at(i), "", refs.at(i));
		}
	}
	
	sort(lcbs.begin(), lcbs.end(), lcbLessThan);
	
	if ( variantList )
	{
		variantList->sortVariants();
	}
	
	in.close();
}

void LcbList::initFromMfa(const char * file, ReferenceList * referenceList, TrackList * trackList, PhylogenyTree * phylogenyTree, VariantList * variantList)
{
	lcbs.resize(0);
	
	ifstream in(file);
	string line;
	vector<string> seqs;
	const bool oldTags = phylogenyTree->getRoot();
	string refTag;
	string refDesc;
	int * trackIndecesNew;
	
	if ( oldTags )
	{
		trackIndecesNew = new int[trackList->getTrackCount()];
	}
	else
	{
		trackList->clear();
	}
	
	while ( getline(in, line) )
	{
		if ( line[0] == '>' )
		{
			string tag = line.substr(1);
			
			string name = parseNameFromTag(tag);
			string desc = parseDescriptionFromTag(tag);
			
			if ( seqs.size() == 0 )
			{
				refTag = tag;
				refDesc = desc;
			}
			
			TrackList::Track * track;
			
			if ( oldTags )
			{
				track = &trackList->getTrackMutable(seqs.size());
				
				try
				{
					trackIndecesNew[trackList->getTrackIndexByFile(name.c_str())] = seqs.size();
				}
				catch ( const TrackList::TrackNotFoundException & e )
				{
					delete [] trackIndecesNew;
					throw;
					return;
				}
			}
			else
			{
				track = &trackList->getTrackMutable(trackList->addTrack(name.c_str()));
			}
			
			track->file = name;
			seqs.resize(seqs.size() + 1);
		}
		else if ( line[0] != '#' )
		{
			if ( seqs.size() == 0 )
			{
				throw BadInputFileException();
			}
			
			seqs[seqs.size() - 1].append(line);
		}
	}
	
	if ( seqs.size() == 0 )
	{
		throw BadInputFileException();
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
	referenceList->addReference(refTag, refDesc, ref);
	
	initWithSingleLcb(*referenceList, *trackList);
	
	if ( variantList )
	{
		variantList->init();
		variantList->addVariantsFromAlignment(seqs, *referenceList, 0, 0, lcbs[0].length, false);
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

void LcbList::initFromXmfa(const char * file, ReferenceList * referenceList, TrackList * trackList, PhylogenyTree * phylogenyTree, VariantList * variantList)
{
	lcbs.resize(0);
	
	ifstream in(file);
	char * line = new char[1 << 20];
	int trackIndex = 0;
	vector<string> seqs;
	const bool oldTags = phylogenyTree->getRoot();
	int * trackIndecesNew;
	bool mauve = false;
	string ref;
	bool createReference = referenceList->getReferenceCount() == 0;
	
	if ( variantList )
	{
		variantList->init();
	}
	
	if ( oldTags )
	{
		trackIndecesNew = new int[trackList->getTrackCount()];
	}
	else
	{
		trackList->clear();
	}
	
	LcbList::Lcb * lcb = 0;
	
	int lcbLength = 0;
	
	TrackList::Track * track;
	
	while ( ! in.eof() )
	{
		in.getline(line, (1 << 20) - 1);
		
		if ( *line == '#' )
		{
			char * suffix;
			
			if ( (suffix = removePrefix(line, "#FormatVersion ")) )
			{
				if ( removePrefix(suffix, "Mauve") )
				{
					mauve = true;
				}
			}
			else if ( mauve && (suffix = removePrefix(line, "#Sequence")) )
			{
				while ( *suffix >= '0' && *suffix <= '9' )
				{
					suffix++;
				}
				
				if ( (suffix = removePrefix(suffix, "File\t")) )
				{
					if ( oldTags )
					{
						track = &trackList->getTrackMutable(trackIndex);
						
						try
						{
							trackIndecesNew[trackList->getTrackIndexByFile(suffix)] = trackIndex;
						}
						catch ( const TrackList::TrackNotFoundException & e )
						{
							delete [] trackIndecesNew;
							throw;
							return;
						}
					}
					else
					{
						track = &trackList->getTrackMutable(trackList->addTrack(suffix));
					}
					
					char * name = suffix;
					
					for ( int i = 0; i < strlen(suffix) - 1; i++ )
					{
						if ( suffix[i] == '/' )
						{
							name = suffix + i + 1;
						}
					}
					
					track->file = name;
					trackIndex++;
				}
			}
			else if ( (suffix = removePrefix(line, "##SequenceFile ")) )
			{
				if ( oldTags )
				{
					track = &trackList->getTrackMutable(trackIndex);
					
					try
					{
						trackIndecesNew[trackList->getTrackIndexByFile(suffix)] = trackIndex;
					}
					catch ( const TrackList::TrackNotFoundException & e )
					{
						delete [] trackIndecesNew;
						throw;
						return;
					}
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
			char * suffix = line + 1;
			
			while ( *suffix == ' ' )
			{
				suffix++;
			}
			
			trackIndex = atoi(strtok(suffix, ":")) - 1;
			
			if ( lcb == 0 )
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
			
			if ( mauve )
			{
				region->position--;
			}
			
			int end = atoi(strtok(0, " "));
			
			if ( mauve )
			{
				end--;
			}
			
			region->length = end - region->position + 1;
			region->reverse = *strtok(0, " ") == '-';
			
			if ( trackIndex == 0 )
			{
				if ( createReference )
				{
					lcb->sequence = 0;
					lcb->position = region->position;
					
					if ( end > ref.length() )
					{
						ref.resize(end, 'N');
					}
				}
				else
				{
					lcb->sequence = referenceList->getReferenceSequenceFromConcatenated(region->position);
					lcb->position = referenceList->getPositionFromConcatenated(lcb->sequence, region->position);
				}
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
				for ( int i = 0; i < seqs.size(); i++ )
				{
					for ( int j = 0; j < seqs[i].length(); j++ )
					{
						seqs[i][j] = toupper(seqs[i][j]);
					}
				}
				
				if ( createReference )
				{
					string ungapped = seqs[0];
					ungap(ungapped);
					ref.replace(lcb->position, ungapped.length(), ungapped);
				}
				
				variantList->addVariantsFromAlignment(seqs, *referenceList, lcb->sequence, lcb->position, lcbLength);
			}
			else
			{
				// not core; destroy
				
				lcbs.resize(lcbs.size() - 1);
			}
			
			lcb = 0;
			
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
			if ( lcb )
			{
				lcb->length = lcbLength;
			}
			
			lcbLength = 0;
		}
	}
	
	if ( oldTags )
	{
		trackList->setTracksByFile();
		phylogenyTree->setTrackIndeces(trackIndecesNew);
		delete [] trackIndecesNew;
	}
	
	if ( createReference )
	{
		referenceList->addReference(trackList->getTrack(0).file, "", ref);
	}
	
	sort(lcbs.begin(), lcbs.end(), lcbLessThan);
	
	if ( variantList )
	{
		variantList->sortVariants();
	}
	
	delete [] line;
	in.close();
}

void LcbList::initWithSingleLcb(const ReferenceList & referenceList, const TrackList & trackList)
{
	int totalLength = 0;
	
	for ( int i = 0; i < referenceList.getReferenceCount(); i++ )
	{
		totalLength += referenceList.getReference(i).sequence.length();
	}
	
	lcbs.resize(1);
	LcbList::Lcb * lcb = &lcbs[0];
	lcb->position = 0;
	lcb->regions.resize(trackList.getTrackCount());
	lcb->length = totalLength;
	
	for ( int i = 0; i < trackList.getTrackCount(); i++ )
	{
		LcbList::Region * region = &lcb->regions[i];
		
		region->position = 0;
		region->length = totalLength;
		region->reverse = false;
	}
}

void LcbList::writeToCapnp(capnp::Harvest::Builder & harvestBuilder) const
{
	auto lcbListBuilder = harvestBuilder.initLcbList();
	auto lcbsBuilder = lcbListBuilder.initLcbs(lcbs.size());
	
	for ( int i = 0; i < lcbs.size(); i++ )
	{
		auto lcbBuilder = lcbsBuilder[i];
		const LcbList::Lcb & lcb = lcbs.at(i);
		
		lcbBuilder.setSequence(lcb.sequence);
		lcbBuilder.setPosition(lcb.position);
		lcbBuilder.setLength(lcb.length);
		lcbBuilder.setConcordance(lcb.concordance);
		
		auto regionsBuilder = lcbBuilder.initRegions(lcb.regions.size());
		
		for ( int j = 0; j < lcb.regions.size(); j++ )
		{
			// TODO: empty tracks?
			
			auto regionBuilder = regionsBuilder[j];
			const LcbList::Region & region = lcb.regions.at(j);
			
			regionBuilder.setTrack(j);
			regionBuilder.setPosition(region.position);
			regionBuilder.setLength(region.length);
			regionBuilder.setReverse(region.reverse);
		}
	}
}

void LcbList::writeToMfa(ostream & out, const ReferenceList & referenceList, const TrackList & trackList, const VariantList & variantList) const
{
	// now iterate over alignments
	
	int totrefgaps = 0;
	
	for ( int i = 0; i < trackList.getTrackCount(); i++)
	{
		int refend = 0;
		int currvar = 0;
		int width = 80;
		int col = 0;
		
		out << '>' << trackList.getTrack(i).file << endl;
		
		for ( int j = 0; j < lcbs.size(); j++ )
		{
			const LcbList::Lcb & lcb = lcbs.at(j);
			int refIndex = lcb.sequence;
			int refstart = lcb.position;
			const LcbList::Region & region = lcb.regions.at(i);
			
			int currpos = refstart;
			int variantsSize = variantList.getVariantCount();
			const VariantList::Variant * currvarref;
			
			if ( currvar < variantsSize )
			{
				currvarref = &variantList.getVariant(currvar);
			
				if ( currvarref->alleles[0] == '-' )
				{
					currpos--;
				}
			}
			
			while
			(
				currpos - refstart < lcb.regions.at(0).length ||
				(
					currvar < variantsSize &&
					currvarref->sequence == refIndex &&
					currvarref->position - refstart < lcb.regions.at(0).length
				)
			)
			{
				if ( currpos == referenceList.getReference(refIndex).sequence.size() )
				{
					printf("ERROR: LCB %d extends beyond reference (position %d)\n", j,
					currpos);
					return;
				}
				
				if ( col == width )
				{
					out << endl;
					col = 0;
				}
				
				if
				(
					currvar == variantsSize ||
					(currpos != currvarref->position && currpos >= refstart) ||
					(
						currvarref->reference == '-' &&
						currvar > 0 &&
						variantList.getVariant(currvar - 1).position != currpos &&
						currpos >= refstart
					)
				)
				{
					out << referenceList.getReference(refIndex).sequence.at(currpos);
					col++;
					
					if ( col == width )
					{
						out << endl;
						col = 0;
					}
				}
				
				if ( currvar < variantsSize && currpos == currvarref->position )
				{
					out << currvarref->alleles[i];
					currvar++;
					
					if ( currvar < variantsSize )
					{
						currvarref = &variantList.getVariant(currvar);
					}
					
					col++;
				}
				
				if ( currvar == variantsSize || currvarref->position > currpos || currvarref->sequence != refIndex )
				{
					currpos++;
				}
			}
			
		}
		
		out << endl;
	}
}
void LcbList::writeFilteredToMfa(ostream & out, ostream & out2, const ReferenceList & referenceList, const TrackList & trackList, const VariantList & variantList) const
{
	// now iterate over alignments
	
	int totrefgaps = 0;
	
	for ( int i = 0; i < trackList.getTrackCount(); i++)
	{
		int refend = 0;
		int currvar = 0;
		int width = 80;
		int col = 0;
		
		out << '>' << trackList.getTrack(i).file << endl;
		
		for ( int j = 0; j < lcbs.size(); j++ )
		{
			const LcbList::Lcb & lcb = lcbs.at(j);
			int refIndex = lcb.sequence;
			int refstart = lcb.position;
			const LcbList::Region & region = lcb.regions.at(i);
			
			int currpos = refstart;
			int variantsSize = variantList.getVariantCount();
			const VariantList::Variant * currvarref;
			
			if ( currvar < variantsSize )
			{
				currvarref = &variantList.getVariant(currvar);
			
				if ( currvarref->alleles[0] == '-' )
				{
					currpos--;
				}
			}
			
			while
			(
				currpos - refstart < lcb.regions.at(0).length ||
				(
					currvar < variantsSize &&
					currvarref->sequence == refIndex &&
					currvarref->position - refstart < lcb.regions.at(0).length
				)
			)
			{
				if ( currpos == referenceList.getReference(refIndex).sequence.size() )
				{
					printf("ERROR: LCB %d extends beyond reference (position %d)\n", j,
					currpos);
					return;
				}
				
				if ( col == width )
				{
					out << endl;
					col = 0;
				}
				
				if
				(
					currvar == variantsSize ||
					(currpos != currvarref->position && currpos >= refstart) ||
					(
						currvarref->reference == '-' &&
						currvar > 0 &&
						variantList.getVariant(currvar - 1).position != currpos &&
						currpos >= refstart
					)
				)
				{
				  // ALB -- do not output if this SNP has been filtered for some reason
				  //if ( currvarref->filters == 0 ) {
					out << referenceList.getReference(refIndex).sequence.at(currpos);
					if(i == 0) {
					  // believe our internal genome sequence index is 0-based, so add one here to be compatible with GenBank 1-based system
					  out2 << currpos + 1 << ",";
					}
					//}
					col++;
					
					if ( col == width )
					{
						out << endl;
						col = 0;
					}
				}
				
				if ( currvar < variantsSize && currpos == currvarref->position )
				{
				  // ALB -- do not output if this SNP has been filtered for some reason
				  if ( currvarref->filters == 0 ) {
					out << currvarref->alleles[i];
					if(i == 0) {
					  out2 << currpos + 1 << ",";
					}
				  }
					currvar++;
					
					if ( currvar < variantsSize )
					{
						currvarref = &variantList.getVariant(currvar);
					}
					
					col++;
				}
				
				if ( currvar == variantsSize || currvarref->position > currpos || currvarref->sequence != refIndex )
				{
					currpos++;
				}
			}
			
		}
		
		out << endl;
		if(i == 0) {
		  out2 << endl;
		}
	}
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
		int refstart = lcb.position;
		int refIndex = lcb.sequence;
		int refend = 0;
		int blockVarStart;
		
		for ( int r = 0; r < lcb.regions.size(); r++)
		{
			// >1:8230-11010 + cluster174 s1:p8230
			
			const LcbList::Region & region = lcb.regions.at(r);
			int start = region.position;
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
			const VariantList::Variant * currvarref;
			
			if ( currvar < variantsSize )
			{
				currvarref = &variantList.getVariant(currvar);
			
				if ( currvarref->alleles[0] == '-' )
				{
					currpos--;
				}
			}
			
			//var =  harvest.variation().variants(currvar);//.alleles()[r];
			//string ref_slice(harvest.reference().references(0).sequence().substr(refstart,(end-start)+1));
			
			while
			(
				currpos - refstart < lcb.regions.at(0).length ||
				(
					currvar < variantsSize &&
					currvarref->sequence == refIndex &&
					currvarref->position - refstart < lcb.regions.at(0).length
				)
			)
			{
				if ( currpos == referenceList.getReference(refIndex).sequence.size() )
				{
					printf("ERROR: LCB %d extends beyond reference (position %d)\n", j, currpos);
					return;
				}
				
				if ( col == width )
				{
					out << endl;
					col = 0;
				}
				
				if
				(
					currvar == variantsSize ||
					(currpos != currvarref->position && currpos >= refstart) ||
					(
						currvarref->alleles[0] == '-' &&
						currvar > 0 &&
						variantList.getVariant(currvar - 1).position != currpos &&
						currpos >= refstart
					)
				)
				{
					out << referenceList.getReference(refIndex).sequence.at(currpos);
					col++;
				
					if ( col == width )
					{
						out << endl;
						col = 0;
					}
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
				
				if ( currvar == variantsSize || currvarref->position > currpos || currvarref->sequence != refIndex )
				{
					currpos++;
				}
			}
			
			out << endl;
		}
		
		out << "=" << endl;
	}
}
