
#ifndef HarvestIO_h
#define HarvestIO_h

#include "harvest.pb.h"
#include <string>
#include <map>
#include <vector>

class HarvestIO
{
public:

	HarvestIO();
	
	void loadBed(const char * file);
	void loadFasta(const char * file);
	void loadGenbank(const char * file);
	bool loadHarvest(const char * file);
	void loadMFA(const char * file);
	void loadNewick(const char * file);
	void loadVcf(const char * file);
	void loadXmfa(const char * file, bool variants = false);
	
	void writeHarvest(const char * file) const;
	void writeNewick(std::ostream &out) const;
	void writeSnp(std::ostream &out, bool indels = false) const;
	void writeVcf(std::ostream &out, bool indels = false) const;
	void writeXmfa(std::ostream &out, bool split = false) const;
	
	Harvest harvest;
	
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
	
	void findVariants(const std::vector<std::string> & seqs, int position = 0);
	char * loadNewickNode(char * token, Harvest::Tree::Node * msg, bool useNames);
	void writeNewickNode(std::ostream &out, const Harvest::Tree::Node & msg) const;
	
	std::map<std::string, google::protobuf::uint32> tracksByFile;
};

char * removePrefix(char * string, const char * substring);
void ungap(std::string & gapped);

#endif
