
#ifndef TrackList_h
#define TrackList_h

#include "harvest/pb/harvest.pb.h"

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
	
	int addTrack(const char * file, int size = 0, const char * name = "", TrackType type = NONE);
	const Track & getTrack(int index) const;
	Track & getTrackMutable(int index);
	int getTrackIndexByFile(const char * file) const;
	int getTrackCount() const;
	void initFromProtocolBuffer(const Harvest::TrackList & msg);
	void setTracksByFile();
	void writeToProtocolBuffer(Harvest * msg);
	
private:
	
	std::vector<Track> tracks;
	std::map<std::string, int> tracksByFile;
};

inline const TrackList::Track & TrackList::getTrack(int index) const { return tracks[index]; }
inline TrackList::Track & TrackList::getTrackMutable(int index) { return tracks[index]; }
inline int TrackList::getTrackCount() const { return tracks.size(); }
#endif
