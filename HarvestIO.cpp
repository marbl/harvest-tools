

#include "HarvestIO.h"

#include <google/protobuf/io/gzip_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/coded_stream.h>
#include <fcntl.h>
#include <fstream>
#include <iostream>

using namespace::std;
using namespace::google::protobuf::io;

HarvestIO::HarvestIO()
{
	GOOGLE_PROTOBUF_VERIFY_VERSION;
}

void HarvestIO::loadBed(const char * file)
{
	ifstream in(file);
	char line[1 << 20];
	int i = 0;
	google::protobuf::uint64 flag = 2 << harvest.variation().filters_size();
	
	Harvest::Variation::Filter * msgFilter = harvest.mutable_variation()->add_filters();
	
	msgFilter->set_flag(flag);
	msgFilter->set_name("REC");
	msgFilter->set_description("Recombination filter");
	
	while ( ! in.eof() )
	{
		in.getline(line, (1 << 20) - 1);
		
		if ( in.eof() )
		{
			break;
		}
		
		strtok(line, "\t"); // eat chromosome
		
		int start = atoi(strtok(0, "\t")) - 1;
		int end = atoi(strtok(0, "\t")) - 1;
		
		// seek to interval start
		//
		while ( i < harvest.variation().variants_size() && harvest.variation().variants(i).position() < start )
		{
			i++;
		}
		
		// set flags through interval
		//
		while ( i < harvest.variation().variants_size() && harvest.variation().variants(i).position() <= end )
		{
			Harvest::Variation::Variant * msgVariant = harvest.mutable_variation()->mutable_variants(i);
			msgVariant->set_filters(msgVariant->filters() | flag);
			i++;
		}
	}
}

void HarvestIO::loadFasta(const char * file)
{
	ifstream in(file);
	char line[1 << 20];
	Harvest::Reference::Sequence * msgSequence;
	
	while ( ! in.eof() )
	{
		in.getline(line, (1 << 20) - 1);
		
		if ( *line == '>' )
		{
			msgSequence = harvest.mutable_reference()->add_references();
			msgSequence->set_tag(line + 1);
		}
		else if ( *line != '#' )
		{
			msgSequence->mutable_sequence()->append(line);
		}
	}
	
	in.close();
}

void HarvestIO::loadGenbank(const char * file)
{
	ifstream in(file);
	char line[1 << 20];
	Harvest::AnnotationList::Annotation * msgAnn = 0;
	google::protobuf::uint32 sequence;
	
	while ( ! in.eof() )
	{
		while ( in.getline(line, (1 << 20) - 1) )
		{
			if ( in.eof() )
			{
				//printf("bad genbank\n");
				return;
			}
		
			if ( removePrefix(line, "VERSION") )
			{
				if ( const char * giToken = strstr(line, " GI:") )
				{
					long int gi = atol(giToken + 4);
					sequence = getReferenceSequenceFromGi(gi);
				
					if ( sequence == undef )
					{
						printf("Couldn't find GI %ld in reference.\n", gi);
						return;
					}
				}
			}
			else if ( removePrefix(line, "FEATURES") )
			{
				break;
			}
		}
	
		while ( ! in.eof() )
		{
			in.getline(line, (1 << 20) - 1);
		
			if ( in.eof() || strcmp(line, "//") == 0 || removePrefix(line, "ORIGIN") )
			{
				break;
			}
			
			char * token = line;
			char * suffix;
			
			while ( *token == ' ' )
			{
				token++;
			}
			
			if ( token == line + 5 )
			{
				int start;
				int end;
				bool reverse;
				
				string feature = strtok(token, " ");
				
				token = strtok(0, " ");
				
				suffix = removePrefix(token, "complement(");
				reverse = suffix;
				
				if ( suffix )
				{
					token = suffix;
				}
				
				suffix = removePrefix(token, "join(");
				
				if ( suffix )
				{
					token = suffix;
				}
				
				start = atoi(strtok(token, ".")) - 1;
				end = atoi(strtok(0, ".,)")) - 1;
				
				if ( ! msgAnn || start != msgAnn->regions(0).start() || end != msgAnn->regions(0).end() || reverse != msgAnn->reverse() )
				{
					if ( feature != "source" )
					{
						msgAnn = harvest.mutable_annotations()->add_annotations();
						msgAnn->set_sequence(sequence);
					
						msgAnn->add_regions();
						msgAnn->mutable_regions(0)->set_start(start);
						msgAnn->mutable_regions(0)->set_end(end);
						msgAnn->set_reverse(reverse);
					}
				}
				
				if ( msgAnn )
				{
					msgAnn->set_feature(feature);
				}
			}
			else
			{
				char * suffix;
			
				if ( (suffix = removePrefix(token, "/locus_tag=\"")) )
				{
					msgAnn->set_locus(strtok(suffix, "\""));
				}
				else if ( (suffix = removePrefix(token, "/gene=\"")) && msgAnn->feature() == "gene" )
				{
					msgAnn->set_name(strtok(suffix, "\""));
				}
				else if ( (suffix = removePrefix(token, "/product=\"")) )
				{
					msgAnn->set_description(suffix);
				
					if ( msgAnn->description()[msgAnn->description().length() - 1] == '"' )
					{
						msgAnn->mutable_description()->resize(msgAnn->description().length() - 1);
					}
				
					while ( suffix[strlen(suffix) - 1] != '"' )
					{
						in.getline(line, (1 << 20) - 1);
						suffix = line;
					
						while ( *suffix == ' ' )
						{
							suffix++;
						}
					
						msgAnn->mutable_description()->append(suffix - 1, strlen(suffix));
					}
				}
			}
		}
	}
	
	for ( int i = 0; i < harvest.annotations().annotations_size(); i++ )
	{
		const Harvest::AnnotationList::Annotation msgAnn = harvest.annotations().annotations(i);
		//printf("%s\t%d\t%d\t%c\t%s\t%s\n", msgAnn.locus().c_str(), msgAnn.regions(0).start(), msgAnn.regions(0).end(), msgAnn.reverse() ? '-' : '+', msgAnn.name().c_str(), msgAnn.description().c_str());
	}
	
	in.close();
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
	
	close(fd);
	google::protobuf::ShutdownProtobufLibrary();
	return true;
}

void HarvestIO::loadMFA(const char * file)
{
	ifstream in(file);
	
	char line[1 << 20];
	vector<string> seqs;
	bool oldTags = tracksByFile.size();
	Harvest::Reference::Sequence * msgSequence;
	msgSequence = harvest.mutable_reference()->add_references();
	Harvest::Alignment * msgAlignment = harvest.mutable_alignment();
	Harvest::TrackList * msgTracks = harvest.mutable_tracks();
	
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
				msgSequence->set_tag(line + 1);
			}
			
			msgTracks->add_tracks()->set_file(tag);
			tracksByFile[tag] = seqs.size();
			seqs.resize(seqs.size() + 1);
		}
		
		in.getline(line, (1 << 20) - 1);
		seqs[seqs.size() - 1].append(line);
	}
	
	in.close();
	msgSequence->set_sequence(seqs[0]);
	ungap(*msgSequence->mutable_sequence());
	
	Harvest::Alignment::Lcb * msgLcb = msgAlignment->add_lcbs();
	msgLcb->set_position(0);
	
	for ( int i = 0; i < msgTracks->tracks_size(); i++ )
	{
		Harvest::Alignment::Lcb::Region * msgRegion = msgLcb->add_regions();
		
		msgRegion->set_track(i);
		msgRegion->set_position(0);
		msgRegion->set_length(msgSequence->sequence().length());
		msgRegion->set_reverse(false);
	}
	
	vector<const Variant *> vars;
	
	findVariants(seqs, vars, 0);
	
	Harvest::Variation * msgVar = harvest.mutable_variation();
	
	//sort(vars.begin(), vars.end(), variantLessThan);
	
	for ( int i = 0; i < vars.size(); i++ )
	{
		Harvest::Variation::Variant * variant = msgVar->add_variants();
		
		variant->set_sequence(vars[i]->sequence);
		variant->set_position(vars[i]->position);
		variant->set_alleles(vars[i]->alleles);
		variant->set_filters(vars[i]->filters);
		
		//printf("VARIANT:\t%d\t%d\t%llu\t%d\t%s\n", variant->sequence(), variant->position(), variant->filters(), variant->quality(), variant->alleles().c_str());
		delete vars[i];
	}
}

void HarvestIO::loadNewick(const char * file)
{
	ifstream in(file);
	char line[1 << 20];
	
	bool useNames = ! harvest.has_tracks() || harvest.tracks().tracks_size() == 0;
	
	while ( in.getline(line, (1 << 20) - 1) )
	{
		loadNewickNode(line, harvest.mutable_tree()->mutable_root(), useNames);
		harvest.mutable_tree()->set_newick(line);
	}
	
	in.close();
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

void HarvestIO::loadXmfa(const char * file, bool variants)
{
	ifstream in(file);
	char line[1 << 20];
	int track = 0;
	vector<string> seqs;
	vector<const Variant *> vars;
	
	if ( variants )
	{
		Harvest::Variation::Filter * msgFilter = harvest.mutable_variation()->add_filters();
	
		msgFilter->set_flag(1);
		msgFilter->set_name("IND");
		msgFilter->set_description("Indel");
	}
	
	Harvest::Alignment * msgAlignment = harvest.mutable_alignment();
	Harvest::TrackList * msgTracks = harvest.mutable_tracks();
	Harvest::Alignment::Lcb * msgLcb;
	
	google::protobuf::uint32 lcb_length = 0;
	
	while ( ! in.eof() )
	{
		in.getline(line, (1 << 20) - 1);
		
		if ( *line == '#' )
		{
			const char * suffix;
			
			if ( (suffix = removePrefix(line, "##SequenceFile ")) )
			{
				msgTracks->add_tracks()->set_file(suffix);
//				printf("track %d - %s\n", track, suffix);
				tracksByFile[suffix] = track;
				track++;
			}
			else if ( (suffix = removePrefix(line, "##SequenceHeader ")) )
			{
				msgTracks->mutable_tracks(track-1)->set_name(suffix);
			}
			else if ( (suffix = removePrefix(line, "##SequenceLength ")) )
			{
				string length (suffix);
				string length_t (length.begin(),length.end()-2);
				msgTracks->mutable_tracks(track-1)->set_size(atoi(length_t.c_str()));
			}
		}
		else if ( *line == '>' )
		{
			track = atoi(strtok(line, ":") + 1) - 1;
			
			if ( track == 0 )
			{
				if ( variants && msgAlignment->lcbs_size() == 0 )
				{
					seqs.resize(tracksByFile.size());
				}
				
				msgLcb = msgAlignment->add_lcbs();
			}
			
			Harvest::Alignment::Lcb::Region * msgRegion = msgLcb->add_regions();
			
			msgRegion->set_track(track);
			msgRegion->set_position(atoi(strtok(0, "-")));
			google::protobuf::uint32 end = atoi(strtok(0, " "));
			msgRegion->set_length(end - msgRegion->position() + 1);
			msgRegion->set_reverse(*strtok(0, " ") == '-');
			
			if ( track == 0 )
			{
				msgLcb->set_sequence(getReferenceSequenceFromConcatenated(msgRegion->position()));
				msgLcb->set_position(getPositionFromConcatenated(msgLcb->sequence(), msgRegion->position()));
			}
		}
		else if ( variants && *line == '=' )
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
				findVariants(seqs, vars, msgLcb->sequence(), msgLcb->position());
			}
			
			for ( int i = 0; i < seqs.size(); i++ )
			{
				seqs[i].clear();
			}
		}
		else if ( variants )
		{
			seqs[track].append(line);
		}
		if ( *line != '=' && *line != '>' && *line != '#')
		{
			if ( track == 0 )
			{
				lcb_length += strlen(line);
			}
		}
		else if (*line == '=')
		{
			msgLcb->set_length(lcb_length);
			lcb_length = 0;
		}
	}
	
	if ( variants )
	{
		Harvest::Variation * msgVar = harvest.mutable_variation();
		
		sort(vars.begin(), vars.end(), variantLessThan);
		
		for ( int i = 0; i < vars.size(); i++ )
		{
			Harvest::Variation::Variant * variant = msgVar->add_variants();
			
			variant->set_sequence(vars[i]->sequence);
			variant->set_position(vars[i]->position);
			variant->set_alleles(vars[i]->alleles);
			variant->set_filters(vars[i]->filters);
			
			//printf("VARIANT:\t%d\t%d\t%llu\t%d\t%s\n", variant->sequence(), variant->position(), variant->filters(), variant->quality(), variant->alleles().c_str());
			delete vars[i];
		}
	}
	
	in.close();
}

void HarvestIO::writeHarvest(const char * file) const
{
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
	if ( ! harvest.has_tree() )
	{
		printf("Cannot write Newick; no tree loaded.\n");
		exit(1);
	}
	
	writeNewickNode(out, harvest.tree().root());
	out << ";\n";
}

void HarvestIO::writeXmfa(std::ostream &out, bool split) const
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
	out << "#SequenceCount " << harvest.tracks().tracks_size() << endl;
	
	for ( int i = 0; i < harvest.tracks().tracks_size(); i++ )
	{
		const Harvest::TrackList::Track & msgTrack = harvest.tracks().tracks(i);
		out << "##SequenceIndex " << i+1 << endl;
		out << "##SequenceFile " << msgTrack.file() << endl;
		out << "##SequenceHeader " << msgTrack.name() << endl;
		
		if ( msgTrack.has_size() )
		{
			out << "##SequenceLength " << msgTrack.size() << "bp" << endl;
		}
	}
	
	out << "#IntervalCount " << harvest.alignment().lcbs_size() << endl;
	
	//now iterate over alignments
	const Harvest::Alignment & msgAlignment = harvest.alignment();
	int totrefgaps = 0;
	int currvar = 0;
	
	for ( int j = 0; j < msgAlignment.lcbs_size(); j++ )
	{
		
		const Harvest::Alignment::Lcb & msgLcb = msgAlignment.lcbs(j);
		int refstart = msgLcb.regions(0).position();
		int refend = 0;
		int blockVarStart;
		
		for ( int r= 0; r < msgLcb.regions_size(); r++)
		{
			//>1:8230-11010 + cluster174 s1:p8230
			const Harvest::TrackList::Track & msgTrack = harvest.tracks().tracks(r);
			const Harvest::Alignment::Lcb::Region & msgRegion = msgLcb.regions(r);
			int start = msgRegion.position() + 1;
			int end = start + msgRegion.length() - 1;
			
			if (r == 0)
			{
				blockVarStart = currvar;
			}
			else
			{
				currvar = blockVarStart;
			}

			out << ">" << r+1 << ":" << start << "-" << end << " ";
			
			if (!msgRegion.reverse())
			{
				out << "+ ";
			}
			else
			{
				out << "- ";
			}
			
			out << "cluster" << j+1 << endl;
			google::protobuf::uint32 currpos = refstart;
			int width = 80;
			int col = 0;
			int variantsSize = harvest.variation().variants_size();
			const Harvest::Variation::Variant * currvarref = & harvest.variation().variants(currvar);
			
			//var =  harvest.variation().variants(currvar);//.alleles()[r];
			//string ref_slice(harvest.reference().references(0).sequence().substr(refstart,(end-start)+1));
			
			while
			(
				currpos - refstart < msgLcb.regions(0).length() ||
				(
					currvar < variantsSize &&
					currvarref->position() - refstart < msgLcb.regions(0).length()
				)
			)
			{
				if ( currpos == harvest.reference().references(0).sequence().size() )
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
					currpos != currvarref->position() ||
					(
						currvarref->alleles()[0] == '-' &&
						currvar > 0 &&
						harvest.variation().variants(currvar - 1).position() != currpos
					)
				)
				{
					out << harvest.reference().references(0).sequence().substr(currpos,1);
					col++;
				}
				
				if ( currvar < variantsSize && currpos == currvarref->position() )
				{
					out << currvarref->alleles()[r];
					currvar++;
					
					if ( currvar < variantsSize )
					{
						currvarref = & harvest.variation().variants(currvar);
					}
					
					col++;
				}
				
				if ( currvar == variantsSize || currvarref->position() > currpos )
				{
					currpos++;
				}
			}
			out << endl;
		}
		out << "=" << endl;
	
	}
		
}

void HarvestIO::writeBackbone(std::ostream &out) const
{
        
        int i  =0;
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

        //the reference sequence
	string refseq;
        refseq = harvest.reference().references(0).sequence();

        //the VCF header line (skipping previous lines for simplicity, can/will add in later)
        //#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  AA1 
        out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";


        int i = 0;
        //output the file name for each column
	for ( i = 0; i < harvest.tracks().tracks_size()-1; i++ )
	{
            const Harvest::TrackList::Track & msgTrack = harvest.tracks().tracks(i);
	    out << (msgTrack.has_name() ? msgTrack.name() : msgTrack.file()) << '\t';
    	    
	}
        //grab the last one and add a new line to avoid extra tab at end of line
        const Harvest::TrackList::Track & msgTrack = harvest.tracks().tracks(i);
        out << (msgTrack.has_name() ? msgTrack.name() : msgTrack.file()) << '\n';
        
        
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
    	  out << msgSnp.sequence() + 1 << "\t" << pos << "\t" << refseq.substr(lend,ws) << "." << refseq.substr(pos,rend);
          
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
          out << "\tNA";
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

int HarvestIO::getPositionFromConcatenated(int sequence, long int position) const
{
	long int sum = 0;
	
	for ( int i = 0; i < sequence; i++ )
	{
		sum += harvest.reference().references(i).sequence().length();
	}
	
	return position - sum;
}

google::protobuf::uint32 HarvestIO::getReferenceSequenceFromConcatenated(long int position) const
{
	long int sum = 0;
	
	for ( int i = 0; i < harvest.reference().references_size(); i++ )
	{
		sum += harvest.reference().references(i).sequence().length();
		
		if ( sum > position )
		{
			return i;
		}
	}
	
	return undef;
}

google::protobuf::uint32 HarvestIO::getReferenceSequenceFromGi(long int gi) const
{
	if ( ! harvest.has_reference() )
	{
		return undef;
	}
	
	for ( int i = 0; i < harvest.reference().references_size(); i++ )
	{
		size_t giToken = harvest.reference().references(i).tag().find("gi|");
		
		if ( giToken != string::npos )
		{
			if ( gi == atol(harvest.reference().references(i).tag().c_str() + giToken + 3))
			{
				return i;
			}
		}
	}
	
	return undef;
}

void HarvestIO::findVariants(const vector<string> & seqs, vector<const Variant *> & vars, int sequence, int position)
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
		}
		
		if ( variant )
		{
			Variant * varNew = new Variant();
			
			while ( position >= harvest.reference().references(sequence).sequence().length() )
			{
				position -= harvest.reference().references(sequence).sequence().length();
				sequence++;
			}
			
			varNew->sequence = sequence;
			varNew->position = position;
			varNew->offset = offset;
			varNew->alleles = col;
			varNew->filters = indel ? 1 : 0;
			varNew->quality = 0;
			/*
			Harvest::Variation::Variant * variant = msg->add_variants();
			
			while ( position >= harvest.reference().references(sequence).sequence().length() )
			{
				position -= harvest.reference().references(sequence).sequence().length();
				sequence++;
			}
			
			variant->set_sequence(sequence);
			variant->set_position(position);
			variant->set_alleles(col);
			variant->set_filters(indel ? 1 : 0);
			*/
			vars.push_back(varNew);
		}
	}
}

char * HarvestIO::loadNewickNode(char * token, Harvest::Tree::Node * msg, bool useNames)
{
	ParseState state = STATE_start;
	char * valueStart;
	
	while ( state != STATE_end )
	{
		if ( state == STATE_start )
		{
			if ( *token == '(' )
			{
				state = STATE_children;
			}
			else
			{
				state = STATE_nameLeaf;
				valueStart = token;
			}
			
			token++;
		}
		else if ( state == STATE_children )
		{
			if ( *token == ')' )
			{
				state = STATE_nameInternal;
				valueStart = token + 1;
				token++;
			}
			else if ( *token == ',' )
			{
				token++;
			}
			else
			{
				token = loadNewickNode(token, msg->add_children(), useNames);
			}
		}
		else if ( state == STATE_nameLeaf || state == STATE_nameInternal )
		{
			if ( *token == ';' )
			{
				state = STATE_end;
			}
			else if ( *token == ':' )
			{
				if ( valueStart != token )
				{
					*token = 0;
					
					if
					(
						(*valueStart == '"' && *(token - 1) == '"') ||
						(*valueStart == '\'' && *(token - 1) == '\'')
					)
					{
						// remove quotes
						
						valueStart++;
						*(token - 1) = 0;
					}
					
					//printf("leaf: %s\n", valueStart);
					if ( state == STATE_nameInternal )
					{
						msg->set_bootstrap(atof(valueStart));
					}
					else
					{
						if ( useNames )
						{
							msg->set_track(harvest.tracks().tracks_size());
							tracksByFile[valueStart] = msg->track();
							harvest.mutable_tracks()->add_tracks()->set_file(valueStart);
						}
						else
						{
							msg->set_track(tracksByFile.at(valueStart));
						}
					}
				}
				
				state = STATE_length;
				valueStart = token + 1;
			}
			
			token++;
		}
		else if ( state == STATE_length )
		{
			if ( *token == ',' || *token == ')' )
			{
				//*token = 0;
				msg->set_branchlength(atof(valueStart));
				state = STATE_end;
			}
			else
			{
				token++;
			}
		}
	}
	
	return token;
}

void HarvestIO::writeNewickNode(std::ostream &out, const Harvest::Tree::Node & msg) const
{
	if ( msg.children_size() )
	{
		out << '(';
		
		for ( int i = 0; i < msg.children_size(); i++ )
		{
			writeNewickNode(out, msg.children(i));
		
			if ( i < msg.children_size() - 1 )
			{
				out << ',';
			}
		}
		
		out << ')';
		
		if ( msg.has_bootstrap() )
		{
			out << msg.bootstrap();
		}
	}
	else
	{
		out << '\'' << harvest.tracks().tracks(msg.track()).file() << '\'';
	}
	
	if ( msg.has_branchlength() )
	{
		out << ':' << msg.branchlength();
	}
}

char * removePrefix(char * string, const char * substring)
{
	size_t len = strlen(substring);
	
	if ( strncmp(string, substring, len) == 0 )
	{
		return string + len;
	}
	else
	{
		return 0;
	}
}

void ungap(string & gapped)
{
	int pos = 0;
	
	for ( int i = 0; i < gapped.length(); i++ )
	{
		if ( gapped[i] != '-' )
		{
			gapped[pos] = gapped[i];
			pos++;
		}
	}
}
