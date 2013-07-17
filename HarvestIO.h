
#ifndef HarvestIO_h
#define HarvestIO_h

#include "harvest.pb.h"
#include <string>
#include <map>

class HarvestIO : public Harvest
{
public:

	HarvestIO();
	
	void loadFasta(const char * file);
	void loadGenbank(const char * file);
	void loadHarvest(const char * file);
	void loadNewick(const char * file);
	void loadVcf(const char * file);
	void loadXmfa(const char * file);
	
	void writeHarvest(const char * file);
	
private:
	
	enum ParseState
	{
		STATE_start,
		STATE_children,
		STATE_nameLeaf,
		STATE_nameInternal,
		STATE_length,
		STATE_end,
	};
	
	char * loadNewickNode(char * token, Harvest::Tree::Node * msg);
	
	std::map<std::string, google::protobuf::uint32> tracksByFile;
};

char * removePrefix(char * string, const char * substring);

#endif
