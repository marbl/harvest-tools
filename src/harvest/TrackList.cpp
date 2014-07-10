#include "harvest/TrackList.h"

using namespace::std;

int TrackList::addTrack(const char * file, int size, const char * name, TrackType type)
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
}

int TrackList::getTrackIndexByFile(const char * file) const
{
	return tracksByFile.at(file);
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
}

void TrackList::setTracksByFile()
{
	for ( int i = 0; i < tracks.size(); i++ )
	{
		tracksByFile[tracks[i].file] = i;
	}
}

void TrackList::writeToProtocolBuffer(Harvest * msg)
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
}
