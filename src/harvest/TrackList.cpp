// Copyright © 2014, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "harvest/TrackList.h"

using namespace::std;

TrackList::TrackList()
{
	trackReference = 0;
}

int TrackList::addTrack(const string & file, int size, const string & name, TrackType type)
{
	tracks.resize(tracks.size() + 1);
	
	tracks[tracks.size() - 1].file = file;
	tracks[tracks.size() - 1].name = name;
	tracks[tracks.size() - 1].size = size;
	tracks[tracks.size() - 1].type = type;
	
	tracksByFile[file] = tracks.size() - 1;
	
	return tracks.size() - 1;
}

void TrackList::clear()
{
	tracks.clear();
	tracksByFile.clear();
	trackReference = 0;
}

int TrackList::getTrackIndexByFile(const string & file) const
{
	try
	{
		return tracksByFile.at(file);
	}
	catch ( const out_of_range & e )
	{
		throw TrackNotFoundException(file);
		return -1;
	}
}

void TrackList::initFromCapnp(const capnp::Harvest::Reader & harvestReader)
{
	auto trackListReader = harvestReader.getTrackList();
	auto tracksReader = trackListReader.getTracks();
	
	tracks.resize(tracksReader.size());
	
	for ( int i = 0; i < tracks.size(); i++ )
	{
		auto trackReader = tracksReader[i];
		
		if ( trackReader.hasFile() )
		{
			tracks[i].file = trackReader.getFile();
		}
		
		if ( trackReader.hasName() )
		{
			tracks[i].name = trackReader.getName();
		}
		
		//if ( trackReader.hasType() )
		{
			tracks[i].type = (TrackType)trackReader.getType();
		}
		
		//if ( trackReader.hasSize() )
		{
			tracks[i].size = trackReader.getSize();
		}
	}
	
	setTracksByFile();
	
	//if ( trackListReader.hasVariantReference() )
	{
		trackReference = trackListReader.getVariantReference();
	}
}

void TrackList::initFromProtocolBuffer(const Harvest::TrackList & msg)
{
	tracks.resize(msg.tracks_size());
	
	for ( int i = 0; i < msg.tracks_size(); i++ )
	{
		if ( msg.tracks(i).has_file() )
		{
			tracks[i].file = msg.tracks(i).file();
		}
		
		if ( msg.tracks(i).has_name() )
		{
			tracks[i].name = msg.tracks(i).name();
		}
		
		if ( msg.tracks(i).has_type() )
		{
			tracks[i].type = (TrackType)msg.tracks(i).type();
		}
		
		if ( msg.tracks(i).has_size() )
		{
			tracks[i].size = msg.tracks(i).size();
		}
	}
	
	setTracksByFile();
	
	if ( msg.has_variantreference() )
	{
		trackReference = msg.variantreference();
	}
}

void TrackList::setTracksByFile()
{
	for ( int i = 0; i < tracks.size(); i++ )
	{
		tracksByFile[tracks[i].file] = i;
	}
}

void TrackList::writeToCapnp(capnp::Harvest::Builder & harvestBuilder) const
{
	auto trackListBuilder = harvestBuilder.initTrackList();
	auto tracksBuilder = trackListBuilder.initTracks(tracks.size());
	
	for ( int i = 0; i < tracks.size(); i++ )
	{
		const Track & track = tracks[i];
		auto trackBuilder = tracksBuilder[i];
		
		if ( track.file != "" )
		{
			trackBuilder.setFile(track.file);
		}
		
		if ( track.name != "" )
		{
			trackBuilder.setName(track.name);
		}
		
		if ( track.size != 0 )
		{
			trackBuilder.setSize(track.size);
		}
		
		if ( track.type != NONE )
		{
			trackBuilder.setType((capnp::Harvest::TrackList::Track::TrackType)track.type);
		}
	}
	
	trackListBuilder.setVariantReference(trackReference);
}

void TrackList::writeToProtocolBuffer(Harvest * msg) const
{
	for ( int i = 0; i < tracks.size(); i++ )
	{
		const Track & track = tracks[i];
		
		Harvest::TrackList::Track * msgTrack = msg->mutable_tracks()->add_tracks();
		
		if ( track.file != "" )
		{
			msgTrack->set_file(track.file);
		}
		
		if ( track.name != "" )
		{
			msgTrack->set_name(track.name);
		}
		
		if ( track.size != 0 )
		{
			msgTrack->set_size(track.size);
		}
		
		if ( track.type != NONE )
		{
			msgTrack->set_type((Harvest::TrackList::Track::Type)track.type);
		}
	}
	
	msg->mutable_tracks()->set_variantreference(trackReference);
}
