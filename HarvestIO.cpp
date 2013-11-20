
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
	Harvest::AnnotationList::Annotation * msgAnn;
	
	while ( in.getline(line, (1 << 20) - 1) )
	{
		if ( in.eof() )
		{
			printf("bad genbank\n");
			return;
		}
		
		if ( removePrefix(line, "FEATURES") )
		{
			break;
		}
	}
	
	while ( ! in.eof() )
	{
		in.getline(line, (1 << 20) - 1);
		
		char * token = line;
		char * suffix;
		
		while ( *token == ' ' )
		{
			token++;
		}
		
		if ( token == line + 5 && (suffix = removePrefix(token, "gene")) )
		{
			token = suffix;
			
			while ( *token == ' ' )
			{
				token++;
			}
			
			msgAnn = harvest.mutable_annotations()->add_annotations();
			msgAnn->set_sequence(0);
			suffix = removePrefix(token, "complement(");
			msgAnn->set_reverse(suffix);
			
			if ( suffix )
			{
				token = suffix;
			}
			
			suffix = removePrefix(token, "join(");
			
			if ( suffix )
			{
				token = suffix;
			}
			
			msgAnn->add_regions();
			msgAnn->mutable_regions(0)->set_start(atoi(strtok(token, ".")) - 1);
			msgAnn->mutable_regions(0)->set_end(atoi(strtok(0, ".,)")) - 1);
		}
		else
		{
			char * suffix;
			
			if ( (suffix = removePrefix(token, "/locus_tag=\"")) )
			{
				msgAnn->set_locus(strtok(suffix, "\""));
			}
			else if ( (suffix = removePrefix(token, "/gene=\"")) )
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
	
	findVariants(seqs);
}

void HarvestIO::loadNewick(const char * file)
{
	ifstream in(file);
	char line[1 << 20];
	
	bool useNames = ! harvest.has_tracks() || harvest.tracks().tracks_size() == 0;
	
	while ( in.getline(line, (1 << 20) - 1) )
	{
		loadNewickNode(line, harvest.mutable_tree()->mutable_root(), useNames);
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
			
			variant->set_sequence(atoi(strtok(line, "\t")));
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
	
	Harvest::Alignment * msgAlignment = harvest.mutable_alignment();
	Harvest::TrackList * msgTracks = harvest.mutable_tracks();
	Harvest::Alignment::Lcb * msgLcb;
	
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
			msgRegion->set_length(end - msgRegion->position());
			msgRegion->set_reverse(*strtok(0, " ") == '-');
			
			if ( track == 0 )
			{
				msgLcb->set_position(msgRegion->position());
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
				findVariants(seqs, msgLcb->position());
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
	}
	
	in.close();
}

void HarvestIO::writeHarvest(const char * file)
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

void HarvestIO::writeSnp(const char * file, bool indels)
{
	ofstream out(file);
	int wrap = 80;
	int col;
	
	for ( int i = 0; i < harvest.tracks().tracks_size(); i++ )
	{
		const Harvest::TrackList::Track & msgTrack = harvest.tracks().tracks(i);
		out << '>' << (msgTrack.has_name() ? msgTrack.name() : msgTrack.file()) << '\n';
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
	
	out.close();
}

void HarvestIO::findVariants(const vector<string> & seqs, int position)
{
	Harvest::Variation * msg = harvest.mutable_variation();
	char col[seqs.size() + 1];
	
	col[seqs.size()] = 0;
	
	for ( int i = 0; i < seqs[0].length(); i++ )
	{
		bool variant = false;
		
		col[0] = seqs[0][i];
		
		bool indel = col[0] == '-';
		
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
			Harvest::Variation::Variant * variant = msg->add_variants();
			
			variant->set_sequence(0);
			variant->set_position(position);
			variant->set_alleles(col);
			variant->set_filters(indel ? 1 : 0);
		}
		
		if ( col[0] != '-' )
		{
			position++;
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
					
					printf("leaf: %s\n", valueStart);
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
