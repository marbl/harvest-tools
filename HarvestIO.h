
#ifndef HarvestIO_h
#define HarvestIO_h

#include "harvest.pb.h"
#include <string>
#include <map>
#include <vector>

class HarvestIO
{
public:

	static const google::protobuf::uint32 undef = -1;
	
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
	void writeBackbone(std::ostream &out) const;
	
	Harvest harvest;
	
private:
	
	struct Variant
	{
		google::protobuf::uint32 sequence;
		google::protobuf::uint32 position;
		int offset;
		google::protobuf::string alleles;
		google::protobuf::uint64 filters;
		google::protobuf::uint32 quality;
	};
	
	enum ParseState
	{
		STATE_start,
		STATE_children,
		STATE_nameLeaf,
		STATE_nameInternal,
		STATE_length,
		STATE_end,
	};
	
	int getPositionFromConcatenated(int sequence, long int position) const;
	google::protobuf::uint32 getReferenceSequenceFromConcatenated(long int position) const;
	google::protobuf::uint32 getReferenceSequenceFromGi(long int gi) const;
	void findVariants(const std::vector<std::string> & seqs, std::vector<const Variant *> & vars, int sequence, int position = 0, bool lcbfilt = false);
	char * loadNewickNode(char * token, Harvest::Tree::Node * msg, bool useNames);
	void writeNewickNode(std::ostream &out, const Harvest::Tree::Node & msg) const;
	
	static bool variantLessThan(const Variant * a, const Variant * b)
	{
		if ( a->sequence == b->sequence )
		{
			if ( a->position == b->position )
			{
				return a->offset < b->offset;
			}
			else
			{
				return a->position < b->position;
			}
		}
		else
		{
			return a->sequence < b->sequence;
		}
	}
	
	std::map<std::string, google::protobuf::uint32> tracksByFile;
};

char * removePrefix(char * string, const char * substring);
void ungap(std::string & gapped);

#endif
