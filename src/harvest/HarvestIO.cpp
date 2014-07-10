

#include "HarvestIO.h"

#include <google/protobuf/io/gzip_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/coded_stream.h>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include "parse.h"

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

void HarvestIO::loadGenbank(const char * file)
{
	annotationList.initFromGenbank(file, referenceList);
}

bool HarvestIO::loadHarvest(const char * file)
{
	int fd = open(file, O_RDONLY);
	
	FileInputStream raw_input(fd);
	GzipInputStream gz(&raw_input);
	CodedInputStream coded_input(&gz);
	
	coded_input.SetTotalBytesLimit(1 << 30, 1 << 30);
	
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

void HarvestIO::loadMFA(const char * file, bool findVariants)
{
	lcbList.initFromMfa(file, &referenceList, &trackList, &phylogenyTree, findVariants ? &variantList : 0);
}

void HarvestIO::loadNewick(const char * file)
{
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
	
	if ( variantList.getVariantCount() || variantList.getFilterCount() )
	{
		variantList.writeToProtocolBuffer(&harvest);
	}
	
	int fd = open(file, O_CREAT | O_WRONLY | O_TRUNC, 0644);
	FileOutputStream stream(fd);
	GzipOutputStream zip_stream(&stream);
	
	if ( ! harvest.SerializeToZeroCopyStream(&zip_stream) )
	{
		printf("Failed to write.\n");
	}
	
	zip_stream.Close();
	stream.Close();
	close(fd);
	
	harvest.Clear();
	google::protobuf::ShutdownProtobufLibrary();
}

void HarvestIO::writeNewick(std::ostream &out) const
{
	if ( ! phylogenyTree.getRoot() )
	{
		printf("Cannot write Newick; no tree loaded.\n");
		exit(1);
	}
	
	phylogenyTree.writeToNewick(out, trackList);
}

void HarvestIO::writeXmfa(std::ostream &out, bool split) const
{
	lcbList.writeToXmfa(out, referenceList, trackList, variantList);
}

void HarvestIO::writeBackbone(std::ostream &out) const
{
	int i = 0;
	
	for ( i = 0; i < harvest.tracks().tracks_size()-1; i++ )
	{
		const Harvest::TrackList::Track & msgTrack = harvest.tracks().tracks(i);
		out << (msgTrack.has_name() ? msgTrack.name() : msgTrack.file()) << "_start\t";
		out << (msgTrack.has_name() ? msgTrack.name() : msgTrack.file()) << "_end\t";
	}
	
	const Harvest::TrackList::Track & msgTrack = harvest.tracks().tracks(i);
	out << (msgTrack.has_name() ? msgTrack.name() : msgTrack.file()) << "_start\t";
	out << (msgTrack.has_name() ? msgTrack.name() : msgTrack.file()) << "_end\n";
	
	//now iterate over alignments
	const Harvest::Alignment & msgAlignment = harvest.alignment();
	for ( int j = 0; j < msgAlignment.lcbs_size(); j++ )
	{
        
  	    const Harvest::Alignment::Lcb & msgLcb = msgAlignment.lcbs(j);
            int r = 0;
            for (r= 0; r < msgLcb.regions_size()-1; r++)
	    {
	          const Harvest::Alignment::Lcb::Region & msgRegion = msgLcb.regions(r);
                  int start = msgRegion.position();
                  int end = msgRegion.position()+msgRegion.length();
                  if (!msgRegion.reverse())
		  {
		    out << start << "\t" << end << "\t";
		  }
                  else
		  {
		    out << -1*start << "\t" << -1*end << "\t";
		  }
                  
	    }
	   const Harvest::Alignment::Lcb::Region & msgRegion = msgLcb.regions(r);
           int start = msgRegion.position();
           int  end = msgRegion.position()+msgRegion.length();
           if (!msgRegion.reverse())
	   {
		out << start << "\t" << end << "\n";
	   }
           else
	   {
	        out << -1*start << "\t" << -1*end << "\n";
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

