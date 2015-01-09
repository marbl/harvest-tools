// Copyright © 2014, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "HarvestIO.h"

#include <google/protobuf/io/gzip_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/coded_stream.h>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include "parse.h"
#include <sys/stat.h>
#include "harvest/fb/harvest_generated.h"
#include <string.h>

using namespace::std;
using namespace::google::protobuf::io;

HarvestIO::HarvestIO()
{
	GOOGLE_PROTOBUF_VERIFY_VERSION;
}

void HarvestIO::clear()
{
	referenceList.clear();
	trackList.clear();
	lcbList.clear();
	variantList.clear();
	annotationList.clear();
	phylogenyTree.clear();
}

void HarvestIO::loadBed(const char * file, const char * name, const char * desc)
{
	variantList.addFilterFromBed(file, name, desc);
}

void HarvestIO::loadFasta(const char * file)
{
	referenceList.initFromFasta(file);
}

void HarvestIO::loadGenbank(const char * file, bool useSeq)
{
	annotationList.initFromGenbank(file, referenceList, useSeq);
}

bool HarvestIO::loadHarvest(const char * file)
{
	ifstream in(file);
	
	char header[flatBufferHeaderLength];
	
	in.read(header, flatBufferHeaderLength);
	in.close();
	
	if ( strncmp(header, flatBufferHeader, flatBufferHeaderLength) == 0 )
	{
		return loadHarvestFlatBuffer(file);
	}
	else
	{
		return loadHarvestProtocolBuffer(file);
	}
}

bool HarvestIO::loadHarvestFlatBuffer(const char * file)
{
	struct stat info;
	stat(file, &info);
	
	char * data = new char[info.st_size];
	
	ifstream in(file);
	in.read(data, info.st_size);
	in.close();
	
	auto harvest = fbHarvest::GetHarvest(data + flatBufferHeaderLength);
	
	variantList.initFromFlatBuffers(harvest->variantList());
	
	return true;
}

bool HarvestIO::loadHarvestProtocolBuffer(const char * file)
{
	Harvest harvest;
	
	int fd = open(file, O_RDONLY);
	
	if ( fd < 0 )
	{
		return false;
	}
	
	FileInputStream raw_input(fd);
	GzipInputStream gz(&raw_input);
	CodedInputStream coded_input(&gz);
	
	coded_input.SetTotalBytesLimit(INT_MAX, INT_MAX);
	
	if ( ! harvest.ParseFromCodedStream(&coded_input) )
	{
		printf("FAIL!\n");
		return false;
	}
	
	if ( harvest.has_reference() )
	{
		referenceList.initFromProtocolBuffer(harvest.reference());
	}
	
	if ( harvest.has_annotations() )
	{
		annotationList.initFromProtocolBuffer(harvest.annotations(), referenceList);
	}
	
	trackList.initFromProtocolBuffer(harvest.tracks());
	
	if ( harvest.has_tree() )
	{
		phylogenyTree.initFromProtocolBuffer(harvest.tree());
	}
	
	if ( harvest.has_alignment() )
	{
		lcbList.initFromProtocolBuffer(harvest.alignment());
	}
	
	if ( harvest.has_variation() )
	{
		variantList.initFromProtocolBuffer(harvest.variation());
	}
	
	harvest.Clear();
	
	close(fd);
	google::protobuf::ShutdownProtobufLibrary();
	return true;
}

void HarvestIO::loadMaf(const char * file, bool findVariants, const char * referenceFileName)
{
	lcbList.initFromMaf(file, &referenceList, &trackList, &phylogenyTree, findVariants ? &variantList : 0, referenceFileName);
}

void HarvestIO::loadMfa(const char * file, bool findVariants)
{
	lcbList.initFromMfa(file, &referenceList, &trackList, &phylogenyTree, findVariants ? &variantList : 0);
}

void HarvestIO::loadNewick(const char * file)
{
	if ( lcbList.getLcbCount() == 0 )
	{
		trackList.clear();
	}
	
	phylogenyTree.initFromNewick(file, &trackList);
}

void HarvestIO::loadVcf(const char * file)
{
	variantList.initFromVcf(file, referenceList, &trackList, &lcbList, &phylogenyTree);
}

void HarvestIO::loadXmfa(const char * file, bool findVariants)
{
	lcbList.initFromXmfa(file, &referenceList, &trackList, &phylogenyTree, findVariants ? &variantList : 0);
}

void HarvestIO::writeFasta(std::ostream &out) const
{
	referenceList.writeToFasta(out);
}

void HarvestIO::writeHarvest(const char * file)
{
	flatbuffers::FlatBufferBuilder fbb;
	fbHarvest::HarvestBuilder harvestBuilder(fbb);
	
	/*
	if ( referenceList.getReferenceCount() )
	{
		referenceList.writeToProtocolBuffer(&harvest);
	}
	
	if ( annotationList.getAnnotationCount() )
	{
		annotationList.writeToProtocolBuffer(&harvest, referenceList);
	}
	
	trackList.writeToProtocolBuffer(&harvest);
	
	if ( phylogenyTree.getRoot() )
	{
		phylogenyTree.writeToProtocolBuffer(&harvest);
	}
	
	if ( lcbList.getLcbCount() )
	{
		lcbList.writeToProtocolBuffer(&harvest);
	}
	*/
	
	if ( variantList.getVariantCount() || variantList.getFilterCount() )
	{
		harvestBuilder.add_variantList(variantList.writeToFlatBuffers(fbb));
	}
	
	FinishHarvestBuffer(fbb, harvestBuilder.Finish());
	
	std::ofstream fout(file);
	
	fout.write(flatBufferHeader, flatBufferHeaderLength);
	fout.write((char *)fbb.GetBufferPointer(), fbb.GetSize());
	fout.close();
	
//	int fd = open(file, O_CREAT | O_WRONLY | O_TRUNC, 0644);
	
//	FileOutputStream stream(fd);
//	GzipOutputStream zip_stream(&stream);
	
//	if ( ! harvest.SerializeToZeroCopyStream(&zip_stream) )
//	{
//		printf("Failed to write.\n");
//	}
	
//	zip_stream.Close();
//	stream.Close();
//	close(fd);
	
//	harvest.Clear();
//	google::protobuf::ShutdownProtobufLibrary();
}

void HarvestIO::writeMfa(std::ostream &out) const
{
	lcbList.writeToMfa(out, referenceList, trackList, variantList);
}

void HarvestIO::writeNewick(std::ostream &out, bool useMult) const
{
	if ( ! phylogenyTree.getRoot() )
	{
		printf("Cannot write Newick; no tree loaded.\n");
		exit(1);
	}
	
	phylogenyTree.writeToNewick(out, trackList, useMult);
}

void HarvestIO::writeXmfa(std::ostream &out, bool split) const
{
	lcbList.writeToXmfa(out, referenceList, trackList, variantList);
}

void HarvestIO::writeBackbone(std::ostream &out) const
{
	int i = 0;
	
	for ( i = 0; i < trackList.getTrackCount(); i++ )
	{
		const TrackList::Track & track = trackList.getTrack(i);
		out << (track.name.length() ? track.name : track.file) << "_start\t";
		out << (track.name.length() ? track.name : track.file) << "_end";
		
		if ( i == trackList.getTrackCount() - 1 )
		{
			out << endl;
		}
		else
		{
			out << '\t';
		}
	}
	
	//now iterate over alignments
	for ( int j = 0; j < lcbList.getLcbCount(); j++ )
	{
  	    const LcbList::Lcb & lcb = lcbList.getLcb(j);
            int r = 0;
            for (r= 0; r < lcb.regions.size(); r++)
	    {
	          const LcbList::Region & region = lcb.regions.at(r);
                  int start = region.position + 1;
                  int end = region.position + region.length;
                  if (!region.reverse)
		  {
		    out << start << "\t" << end;
		  }
                  else
		  {
		    out << -1*start << "\t" << -1*end;
		  }
		  
		  if ( r == lcb.regions.size() - 1 )
		  {
		  	out << endl;
		  }
		  else
		  {
		  	out << '\t';
		  }
	    }
	}

}

void HarvestIO::writeSnp(std::ostream &out, bool indels) const
{
	variantList.writeToMfa(out, indels, trackList);
}

void HarvestIO::writeVcf(std::ostream &out, bool indels) const
{
	variantList.writeToVcf(out, indels, referenceList, trackList);
}

