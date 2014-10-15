// Copyright © 2014, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef TrackList_h
#define TrackList_h

#include "harvest/pb/harvest.pb.h"

#include <stdexcept>
#include <string>
#include <vector>
#include <map>

enum TrackType
{
	NONE = 0,
	GENOME = 1,
	CONTIGS = 2,
	SCAFFOLDS = 3,
	READS = 4,
};

class TrackList
{
public:

	struct Track
	{
		Track()
		{
			size = 0;
			type = NONE;
		}
	
		std::string file;
		std::string name;
		int size;
		TrackType type;
	};
	
	class TrackNotFoundException : public std::exception
	{
	public:
		
		TrackNotFoundException(const std::string & nameNew)
		{
			name = nameNew;
		}
		
		virtual ~TrackNotFoundException() throw() {}
		
		std::string name;
	};
	
	TrackList();
	
	int addTrack(const std::string & file, int size = 0, const std::string & name = "", TrackType type = NONE);
	void clear();
	const Track & getTrack(int index) const;
	int getTrackCount() const;
	int getTrackIndexByFile(const std::string & file) const;
	Track & getTrackMutable(int index);
	int getTrackReference() const;
	void initFromProtocolBuffer(const Harvest::TrackList & msg);
	void setTrackReference(int trackReferenceNew);
	void setTracksByFile();
	void writeToProtocolBuffer(Harvest * msg);
	
private:
	
	std::vector<Track> tracks;
	int trackReference;
	std::map<std::string, int> tracksByFile;
};

inline const TrackList::Track & TrackList::getTrack(int index) const { return tracks[index]; }
inline int TrackList::getTrackCount() const { return tracks.size(); }
inline TrackList::Track & TrackList::getTrackMutable(int index) { return tracks[index]; }
inline int TrackList::getTrackReference() const { return trackReference; }
inline void TrackList::setTrackReference(int trackReferenceNew) { trackReference = trackReferenceNew; }
#endif
