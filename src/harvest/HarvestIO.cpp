

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

void HarvestIO::loadBed(const char * file, const char * name, const char * desc)
{
	ifstream in(file);
	char line[1 << 20];
	int i = 0;
	google::protobuf::uint64 flag = 1 << harvest.variation().filters_size();
	
	Harvest::Variation::Filter * msgFilter = harvest.mutable_variation()->add_filters();
	
	msgFilter->set_flag(flag);
	msgFilter->set_name(name);
	msgFilter->set_description(desc);
	
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
			i < harvest.variation().variants_size() &&
			(
				harvest.variation().variants(i).sequence() < seq ||
				(
					harvest.variation().variants(i).sequence() == seq &&
					harvest.variation().variants(i).position() < start
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
			i < harvest.variation().variants_size() &&
			harvest.variation().variants(i).sequence() == seq &&
			harvest.variation().variants(i).position() <= end
		)
		{
			Harvest::Variation::Variant * msgVariant = harvest.mutable_variation()->mutable_variants(i);
			msgVariant->set_filters(msgVariant->filters() | flag);
			i++;
		}
	}
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
	
	referenceList.initFromProtocolBuffer(harvest.reference());
	annotationList.initFromProtocolBuffer(harvest.annotations(), referenceList);
	trackList.initFromProtocolBuffer(harvest.tracks());
	phylogenyTree.initFromProtocolBuffer(harvest.tree());
	lcbList.initFromProtocolBuffer(harvest.alignment());
	variantList.initFromProtocolBuffer(harvest.variation());
	
	close(fd);
	google::protobuf::ShutdownProtobufLibrary();
	return true;
}

void HarvestIO::loadMFA(const char * file, bool findVariants, int * trackIndecesNew)
{
	lcbList.initFromMfa(file, &referenceList, &trackList, &phylogenyTree, findVariants ? &variantList : 0);
}

void HarvestIO::loadNewick(const char * file)
{
	phylogenyTree.initFromNewick(file, &trackList);
}

void HarvestIO::loadVcf(const char * file)
{
	ifstream in(file);
	
	Harvest::Variation * msg = harvest.mutable_variation();
	char line[1 << 20];
	map<string, google::protobuf::uint64> flagsByFilter;
	unsigned int alleleCount = 0;
	
	while ( ! in.eof() )
	{
		if ( in.peek() == '#' )
		{
			in.getline(line, (1 << 20) - 1);
			
			if ( strncmp(line, "##FILTER", 8) == 0 )
			{
				char * token;
				
				Harvest::Variation::Filter * filter = msg->add_filters();
				
				token = strtok(line, "<");
				
				while ( (token = strtok(0, ",>\"")) )
				{
					if ( strncmp(token, "ID=", 3) == 0 )
					{
						filter->set_name(token + 3);
					}
					else if ( strcmp(token, "Description=") == 0 )
					{
						filter->set_description(strtok(0, "\""));
						strtok(0, ">"); // eat
					}
				}
				
				google::protobuf::uint64 flag = 1 << flagsByFilter.size();
				flagsByFilter[filter->name()] = flag;
				filter->set_flag(flag);
				//printf("FILTER:\t%d\t%s\t%s\n", filter->flag(), filter->name().c_str(), filter->description().c_str());
			}
			else if ( ! harvest.has_tracks() && strncmp(line, "#CHROM", 6) == 0 )
			{
				char * token;
				
				strtok(line, "\t");
				
				for ( int i = 0; i < 8; i++ )
				{
					strtok(0, "\t"); // eat headers
				}
				
				while ( (token = strtok(0, "\t")) )
				{
					Harvest::TrackList::Track * msgTrack = harvest.mutable_tracks()->add_tracks();
					msgTrack->set_name(token);
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
			
			Harvest::Variation::Variant * variant = msg->add_variants();
			
			variant->set_sequence(atoi(strtok(line, "\t")) - 1);
			variant->set_position(atoi(strtok(0, "\t")) - 1);
			strtok(0, "\t"); // eat id
			char * alleles = strtok(0, "\t"); // ref allele
			strtok(0, "\t"); // eat alt alleles
			variant->set_quality(atoi(strtok(0, "\t")));
			
			google::protobuf::uint64 filters = 0;
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
			
			variant->set_filters(filters);
			
			strtok(0, "\t"); // eat format
			
			char * alleleString;
			variant->mutable_alleles()->resize(alleleCount);
			alleleCount = 0;
			
			while ( (alleleString = strtok(0, "\t")) )
			{
				if ( variant->alleles().size() < alleleCount + 1 )
				{
					variant->mutable_alleles()->resize(alleleCount + 1);
				}
				
				(*variant->mutable_alleles())[alleleCount] = alleles[atoi(alleleString) * 2];
				alleleCount++;
			}
			
//			(*variant->mutable_alleles())[alleleCount] = '\0';
			
			//printf("VARIANT:\t%d\t%d\t%d\t%d\t%s\n", variant->sequence(), variant->position(), variant->filters(), variant->quality(), variant->alleles().c_str());
		}
	}
	
	in.close();
}

void HarvestIO::loadXmfa(const char * file, bool findVariants, int * trackIndecesNew)
{
	lcbList.initFromXmfa(file, referenceList, &trackList, &phylogenyTree, findVariants ? &variantList : 0);
}

void HarvestIO::writeHarvest(const char * file)
{
	referenceList.writeToProtocolBuffer(&harvest);
	annotationList.writeToProtocolBuffer(&harvest, referenceList);
	trackList.writeToProtocolBuffer(&harvest);
	phylogenyTree.writeToProtocolBuffer(&harvest);
	lcbList.writeToProtocolBuffer(&harvest);
	variantList.writeToProtocolBuffer(&harvest);
	
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
	int wrap = 80;
	int col;
	
	for ( int i = 0; i < harvest.tracks().tracks_size(); i++ )
	{
		const Harvest::TrackList::Track & msgTrack = harvest.tracks().tracks(i);
		out << '>' << (msgTrack.has_file() ? msgTrack.file() : msgTrack.name()) << '\n';
		col = 0;
		
		for ( int j = 0; j < harvest.variation().variants_size(); j++ )
		{
			if ( ! indels && harvest.variation().variants(j).filters() )
			{
				continue;
			}
			
			col++;
			
			if ( wrap && col > wrap )
			{
				out << '\n';
				col = 1;
			}
			
			out << harvest.variation().variants(j).alleles()[i];
		}
		
		out << '\n';
	}
	
}

void HarvestIO::writeVcf(std::ostream &out, bool indels) const
{
        //tjt: Currently outputs SNPs, no indels
        //tjt: next pass will add standard VCF output for indels, plus an attempt at qual vals
        //tjt: also filters need to be added to findVariants to populate FILTer column

        //indel char, to skip columns with indels (for now)
        char indl = '-';
        //the VCF output file
	
	for ( int i = 0; i < harvest.variation().filters_size(); i++ )
	{
		const Harvest::Variation::Filter & msgFilter = harvest.variation().filters(i);
		
		out << "##FILTER=<ID=" << msgFilter.name() << ",Description=\"" << msgFilter.description() << "\">\n";
	}
	
        //the reference sequence
	string refseq;
        refseq = harvest.reference().references(0).sequence();

        //the VCF header line (skipping previous lines for simplicity, can/will add in later)
        //#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  AA1 
        out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";


	//output the file name for each column
	for ( int i = 0; i < harvest.tracks().tracks_size(); i++ )
	{
		const Harvest::TrackList::Track & msgTrack = harvest.tracks().tracks(i);
		//out << '\t' << (msgTrack.has_name() ? msgTrack.name() : msgTrack.file());
		out << '\t' << msgTrack.file();
	}
	
	out << '\n';
	
        //now iterate over variants and output
	for ( int j = 0; j < harvest.variation().variants_size(); j++ )
	{
          
	  const Harvest::Variation::Variant & msgSnp = harvest.variation().variants(j);

          //no indels for now..
          if (harvest.variation().variants(j).alleles()[0] == indl)
	    continue;
          if (find(harvest.variation().variants(j).alleles().begin(), harvest.variation().variants(j).alleles().end(),indl) != harvest.variation().variants(j).alleles().end())
            continue;

          //capture the reference position of variant
	  int pos = msgSnp.position();
          
          //output first few columns, including context (+/- 7bp for now)
          int ws = 10;
          int lend = pos-ws;
          int rend = ws;
 
          if (lend < 0)
	    lend = 0;
          if (pos+ws >= refseq.size())
            rend = refseq.size()-pos;
          if (pos+rend >= refseq.size())
            rend = 0;
    	  out << msgSnp.sequence() + 1 << "\t" << pos + 1 << "\t" << refseq.substr(lend,ws) << "." << refseq.substr(pos,rend);
          
	  //build non-redundant allele list from cur alleles
          vector<char> allele_list;
          //first allele is ref allele (0)
          allele_list.push_back(harvest.variation().variants(j).alleles()[0]);
          bool prev_var = false;
	  for ( int i = 0; i < harvest.tracks().tracks_size(); i++ )
	  {
	    if (find(allele_list.begin(), allele_list.end(), harvest.variation().variants(j).alleles()[i]) == allele_list.end())
	    {
              if (harvest.variation().variants(j).alleles()[i] == indl) 
		continue;
              //to know if we need to output a preceding comma
              if (!prev_var)
                  out << harvest.variation().variants(j).alleles()[i];
              else
		out << "," << harvest.variation().variants(j).alleles()[i];

              allele_list.push_back(harvest.variation().variants(j).alleles()[i]);
              prev_var = true;
	    }
            //to see if we are in REF column
            else if (i == 0)
	    {
		out << "\t" << harvest.variation().variants(j).alleles()[i] << "\t";
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
		for ( int i = 0; i < harvest.variation().filters_size(); i++ )
		{
			const Harvest::Variation::Filter & msgFilter = harvest.variation().filters(i);
			
			if ( msgSnp.filters() & msgFilter.flag() )
			{
				if ( filterCount > 0 )
				{
					out << ':';
				}
				
				out << msgFilter.name();
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
          //now get allele index and output. must be an easier way than this? I've been spoiled by python it seems..
	  for (i = 0; i < harvest.tracks().tracks_size()-1; i++ )
	  {
	    int idx = distance(allele_list.begin(),(find(allele_list.begin(), allele_list.end(), harvest.variation().variants(j).alleles()[i])));
            out << "\t" << idx;
            
	  }
	  int idx = distance(allele_list.begin(),(find(allele_list.begin(), allele_list.end(), harvest.variation().variants(j).alleles()[i])));
          out << "\t" << idx << "\n";
		
	}
	//done! should be well-formated VCF (see above notes)
	//out.close();
}
