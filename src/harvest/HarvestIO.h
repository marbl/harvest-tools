// Copyright © 2014, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef HarvestIO_h
#define HarvestIO_h

#include "harvest/pb/harvest.pb.h"
#include <string>
#include <map>
#include <vector>

#include "harvest/ReferenceList.h"
#include "harvest/AnnotationList.h"
#include "harvest/PhylogenyTree.h"
#include "harvest/LcbList.h"
#include "harvest/VariantList.h"

static const char * capnpHeader = "Cap'n Proto";
static const int capnpHeaderLength = strlen(capnpHeader);

class HarvestIO
{
public:

	HarvestIO();
	
	void clear();
	
	void loadBed(const char * file, const char * name, const char * desc);
	void loadFasta(const char * file);
	void loadGenbank(const char * file, bool useSeq);
	bool loadHarvest(const char * file);
	bool loadHarvestCapnp(const char * file);
	bool loadHarvestProtocolBuffer(const char * file);
	void loadMaf(const char * file, bool findVariants, const char * referenceFileName);
	void loadMfa(const char * file, bool findVariants);
	void loadNewick(const char * file);
	void loadVcf(const char * file);
	void loadXmfa(const char * file, bool findVariants);
	
	void writeFasta(std::ostream &out) const;
	void writeHarvest(const char * file);
	void writeMfa(std::ostream &out) const;
	void writeFilteredMfa(std::ostream &out, std::ostream &out2) const;
	void writeNewick(std::ostream &out, bool useMult = false) const;
	void writeSnp(std::ostream &out, bool indels = false) const;
	void writeVcf(std::ostream &out, const std::vector<std::string> * trackNames = 0, const PhylogenyTreeNode * node = 0, bool indels = false, bool signature = false) const;
	void writeXmfa(std::ostream &out, bool split = false) const;
	void writeBackbone(std::ostream &out) const;
	
	ReferenceList referenceList;
	AnnotationList annotationList;
	PhylogenyTree phylogenyTree;
	TrackList trackList;
	LcbList lcbList;
	VariantList variantList;
	
private:
	
	void writeNewickNode(std::ostream &out, const Harvest::Tree::Node & msg) const;
};

int def(int fdSource, int fdDest, int level);
int inf(int fdSource, int fdDest);
void zerr(int ret);

#endif
