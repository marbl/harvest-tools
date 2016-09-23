// Copyright © 2014, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
// Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "HarvestIO.h"

#include <google/protobuf/io/gzip_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/coded_stream.h>
#include <capnp/message.h>
#include <capnp/serialize.h>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include "parse.h"
#include <sys/stat.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include <stdio.h>
#include <assert.h>

#define SET_BINARY_MODE(file)

#define CHUNK 16384

using namespace::std;
using namespace::google::protobuf::io;

HarvestIO::HarvestIO()
{
	GOOGLE_PROTOBUF_VERIFY_VERSION;
}

void HarvestIO::clear()
{
	referenceList.clear();
	trackList.clear();
	lcbList.clear();
	variantList.clear();
	annotationList.clear();
	phylogenyTree.clear();
}

void HarvestIO::loadBed(const char * file, const char * name, const char * desc)
{
	variantList.addFilterFromBed(file, name, desc);
}

void HarvestIO::loadFasta(const char * file)
{
	referenceList.initFromFasta(file);
}

void HarvestIO::loadGenbank(const char * file, bool useSeq)
{
	annotationList.initFromGenbank(file, referenceList, useSeq);
}

bool HarvestIO::loadHarvest(const char * file)
{
	ifstream in(file);
	
	char header[capnpHeaderLength];
	
	in.read(header, capnpHeaderLength);
	in.close();
	
	if ( false || strncmp(header, capnpHeader, capnpHeaderLength) == 0 )
	{
		return loadHarvestCapnp(file);
	}
	else
	{
		return loadHarvestProtocolBuffer(file);
	}
}

bool HarvestIO::loadHarvestCapnp(const char * file)
{
	// use a pipe to decompress input to Cap'n Proto
	
	int fds[2];
	int piped = pipe(fds);
	
	if ( piped < 0 )
	{
		cerr << "ERROR: could not open pipe for decompression\n";
		return 1;
	}
	
	int forked = fork();
	
	if ( forked < 0 )
	{
		cerr << "ERROR: could not fork for decompression\n";
		return 1;
	}
	
	if ( forked == 0 )
	{
		// read from zipped fd and write to pipe
		
		close(fds[0]); // other process's end of pipe
		
		int fd = open(file, O_RDONLY);
		
		if ( fd < 0 )
		{
			cerr << "ERROR: could not open " << file << " for reading.\n";
			exit(1);
		}
		
		char buffer[1024];
		
		read(fd, buffer, capnpHeaderLength);
		
		int ret = inf(fd, fds[1]);
		if (ret != Z_OK) zerr(ret);
		close(fd);
		exit(ret);
		
		gzFile fileIn = gzopen(file, "rb");
		
		int bytesRead;
		
		// eat header
		//
		//gzseek(fileIn, capnpHeaderLength, SEEK_SET);
		gzread(fileIn, buffer, capnpHeaderLength);
		
		printf("header: %s\n", buffer);
		while ( (bytesRead = gzread(fileIn, buffer, sizeof(buffer))) > 0)
		{
			printf("uncompressed: %s\n", buffer);
			write(fds[1], buffer, bytesRead);
		}
		
		gzclose(fileIn);
		close(fds[1]);
		exit(0);
	}
	
	// read from pipe
	
	close(fds[1]); // other process's end of pipe
	
	capnp::ReaderOptions readerOptions;
	
	readerOptions.traversalLimitInWords = 1000000000000;
	readerOptions.nestingLimit = 1000000;
	
	//char buffer[1024];
	//read(fds[0], buffer, 1024);
	//printf("data: %s\n", buffer);
	//return true;
	
	capnp::StreamFdMessageReader message(fds[0], readerOptions);
	
	capnp::Harvest::Reader harvestReader = message.getRoot<capnp::Harvest>();
	
	if ( harvestReader.hasReferenceList() )
	{
		referenceList.initFromCapnp(harvestReader);
	}
	
	if ( harvestReader.hasAnnotationList() )
	{
		annotationList.initFromCapnp(harvestReader, referenceList);
	}
	
	trackList.initFromCapnp(harvestReader);
	
	if ( harvestReader.hasTree() )
	{
		phylogenyTree.initFromCapnp(harvestReader);
	}
	
	if ( harvestReader.hasLcbList() )
	{
		lcbList.initFromCapnp(harvestReader);
	}
	
	if ( harvestReader.hasVariantList() )
	{
		variantList.initFromCapnp(harvestReader);
	}
	
	close(fds[0]);
	return true;
}

bool HarvestIO::loadHarvestProtocolBuffer(const char * file)
{
	Harvest harvest;
	
	int fd = open(file, O_RDONLY);
	
	if ( fd < 0 )
	{
		return false;
	}
	
	FileInputStream raw_input(fd);
	GzipInputStream gz(&raw_input);
	CodedInputStream coded_input(&gz);
	
	coded_input.SetTotalBytesLimit(INT_MAX, INT_MAX);
	
	if ( ! harvest.ParseFromCodedStream(&coded_input) )
	{
		printf("FAIL!\n");
		return false;
	}
	
	if ( harvest.has_reference() )
	{
		referenceList.initFromProtocolBuffer(harvest.reference());
	}
	
	if ( harvest.has_annotations() )
	{
		annotationList.initFromProtocolBuffer(harvest.annotations(), referenceList);
	}
	
	trackList.initFromProtocolBuffer(harvest.tracks());
	
	if ( harvest.has_tree() )
	{
		phylogenyTree.initFromProtocolBuffer(harvest.tree());
	}
	
	if ( harvest.has_alignment() )
	{
		lcbList.initFromProtocolBuffer(harvest.alignment());
	}
	
	if ( harvest.has_variation() )
	{
		variantList.initFromProtocolBuffer(harvest.variation());
	}
	
	harvest.Clear();
	
	close(fd);
	google::protobuf::ShutdownProtobufLibrary();
	return true;
}

void HarvestIO::loadMaf(const char * file, bool findVariants, const char * referenceFileName)
{
	lcbList.initFromMaf(file, &referenceList, &trackList, &phylogenyTree, findVariants ? &variantList : 0, referenceFileName);
}

void HarvestIO::loadMfa(const char * file, bool findVariants)
{
	lcbList.initFromMfa(file, &referenceList, &trackList, &phylogenyTree, findVariants ? &variantList : 0);
}

void HarvestIO::loadNewick(const char * file)
{
	if ( lcbList.getLcbCount() == 0 )
	{
		trackList.clear();
	}
	
	phylogenyTree.initFromNewick(file, &trackList);
}

void HarvestIO::loadVcf(const char * file)
{
	variantList.initFromVcf(file, referenceList, &trackList, &lcbList, &phylogenyTree);
}

void HarvestIO::loadXmfa(const char * file, bool findVariants)
{
	lcbList.initFromXmfa(file, &referenceList, &trackList, &phylogenyTree, findVariants ? &variantList : 0);
}

void HarvestIO::writeFasta(std::ostream &out) const
{
	referenceList.writeToFasta(out);
}

void HarvestIO::writeHarvest(const char * file)
{
	// use a pipe to compress Cap'n Proto output
	
	int fds[2];
	int piped = pipe(fds);
	
	if ( piped < 0 )
	{
		cerr << "ERROR: could not open pipe for compression\n";
		exit(1);
	}
	
	int forked = fork();
	
	if ( forked < 0 )
	{
		cerr << "ERROR: could not fork for compression\n";
		exit(1);
	}
	
	if ( forked == 0 )
	{
		// read from pipe and write to compressed file
		
		close(fds[1]); // other process's end of pipe
		
		int fd = open(file, O_CREAT | O_WRONLY | O_TRUNC, 0644);
		
		if ( fd < 0 )
		{
			cerr << "ERROR: could not open " << file << " for writing.\n";
			exit(1);
		}
		
		// write header
		//
		write(fd, capnpHeader, capnpHeaderLength);
		//close(fd);
		
		int ret = def(fds[0], fd, Z_DEFAULT_COMPRESSION);
		if (ret != Z_OK) zerr(ret);
        exit(ret);
        
		char buffer[1024];
		gzFile fileOut = gzopen(file, "ab");
		
		int bytesRead;
		
		while ( (bytesRead = read(fds[0], buffer, sizeof(buffer))) > 0)
		{
			printf("compressing: %s\n", buffer);
			gzwrite(fileOut, buffer, bytesRead);
		}
		
		gzclose(fileOut);
		close(fds[0]);
		exit(0);
	}
	
	// write to pipe
	
	close(fds[0]); // other process's end of pipe
	/*
	write(fds[1], "test\n", 5);
	close(fds[1]);
	return;
	*/
	
	capnp::MallocMessageBuilder message;
	capnp::Harvest::Builder harvestBuilder = message.initRoot<capnp::Harvest>();
	
	if ( referenceList.getReferenceCount() )
	{
		referenceList.writeToCapnp(harvestBuilder);
	}
	
	if ( annotationList.getAnnotationCount() )
	{
		annotationList.writeToCapnp(harvestBuilder, referenceList);
	}
	
	trackList.writeToCapnp(harvestBuilder);
	
	if ( phylogenyTree.getRoot() )
	{
		phylogenyTree.writeToCapnp(harvestBuilder);
	}
	
	if ( lcbList.getLcbCount() )
	{
		lcbList.writeToCapnp(harvestBuilder);
	}
	
	if ( variantList.getVariantCount() || variantList.getFilterCount() )
	{
		variantList.writeToCapnp(harvestBuilder);
	}
	
	//write(fds[1], capnpHeader, capnpHeaderLength);
	
	//write(fds[1], "hello", 5);
	writeMessageToFd(fds[1], message);
	
//	FileOutputStream stream(fd);
//	GzipOutputStream zip_stream(&stream);
	
//	if ( ! harvest.SerializeToZeroCopyStream(&zip_stream) )
//	{
//		printf("Failed to write.\n");
//	}
	
//	zip_stream.Close();
//	stream.Close();
	close(fds[1]);
	
//	harvest.Clear();
//	google::protobuf::ShutdownProtobufLibrary();
}

void HarvestIO::writeMfa(std::ostream &out) const
{
	lcbList.writeToMfa(out, referenceList, trackList, variantList);
}

void HarvestIO::writeFilteredMfa(std::ostream &out, std::ostream &out2) const
{
        lcbList.writeFilteredToMfa(out, out2, referenceList, trackList, variantList);
}

void HarvestIO::writeNewick(std::ostream &out, bool useMult) const
{
	if ( ! phylogenyTree.getRoot() )
	{
		printf("Cannot write Newick; no tree loaded.\n");
		exit(1);
	}
	
	phylogenyTree.writeToNewick(out, trackList, useMult);
}

void HarvestIO::writeXmfa(std::ostream &out, bool split) const
{
	lcbList.writeToXmfa(out, referenceList, trackList, variantList);
}

void HarvestIO::writeBackbone(std::ostream &out) const
{
	int i = 0;
	
	for ( i = 0; i < trackList.getTrackCount(); i++ )
	{
		const TrackList::Track & track = trackList.getTrack(i);
		out << (track.name.length() ? track.name : track.file) << "_start\t";
		out << (track.name.length() ? track.name : track.file) << "_end";
		
		if ( i == trackList.getTrackCount() - 1 )
		{
			out << endl;
		}
		else
		{
			out << '\t';
		}
	}
	
	// now iterate over alignments
	
	for ( int j = 0; j < lcbList.getLcbCount(); j++ )
	{
  	    const LcbList::Lcb & lcb = lcbList.getLcb(j);
            int r = 0;
            for (r= 0; r < lcb.regions.size(); r++)
	    {
	          const LcbList::Region & region = lcb.regions.at(r);
                  int start = region.position + 1;
                  int end = region.position + region.length;
                  if (!region.reverse)
		  {
		    out << start << "\t" << end;
		  }
                  else
		  {
		    out << -1*start << "\t" << -1*end;
		  }
		  
		  if ( r == lcb.regions.size() - 1 )
		  {
		  	out << endl;
		  }
		  else
		  {
		  	out << '\t';
		  }
	    }
	}

}

void HarvestIO::writeSnp(std::ostream &out, bool indels) const
{
	variantList.writeToMfa(out, indels, trackList);
}

void HarvestIO::writeVcf(std::ostream &out, bool indels) const
{
	variantList.writeToVcf(out, indels, referenceList, annotationList, trackList);
}


// The following functions are adapted from http://www.zlib.net/zpipe.c


/* Compress from file source to file dest until EOF on source.
   def() returns Z_OK on success, Z_MEM_ERROR if memory could not be
   allocated for processing, Z_STREAM_ERROR if an invalid compression
   level is supplied, Z_VERSION_ERROR if the version of zlib.h and the
   version of the library linked do not match, or Z_ERRNO if there is
   an error reading or writing the files. */
int def(int fdSource, int fdDest, int level)
{
    int ret, flush;
    unsigned have;
    z_stream strm;
    unsigned char in[CHUNK];
    unsigned char out[CHUNK];

    /* allocate deflate state */
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    ret = deflateInit(&strm, level);
    if (ret != Z_OK)
        return ret;

    /* compress until end of file */
    do {
        strm.avail_in = read(fdSource, in, CHUNK);
        if (strm.avail_in == -1) {
            (void)deflateEnd(&strm);
            return Z_ERRNO;
        }
        flush = strm.avail_in == 0 ? Z_FINISH : Z_NO_FLUSH;
        strm.next_in = in;

        /* run deflate() on input until output buffer not full, finish
           compression if all of source has been read in */
        do {
            strm.avail_out = CHUNK;
            strm.next_out = out;
            ret = deflate(&strm, flush);    /* no bad return value */
            assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
            have = CHUNK - strm.avail_out;
            if (write(fdDest, out, have) != have) {
                (void)deflateEnd(&strm);
                return Z_ERRNO;
            }
        } while (strm.avail_out == 0);
        assert(strm.avail_in == 0);     /* all input will be used */

        /* done when last data in file processed */
    } while (flush != Z_FINISH);
    assert(ret == Z_STREAM_END);        /* stream will be complete */

    /* clean up and return */
    (void)deflateEnd(&strm);
    return Z_OK;
}

/* Decompress from file source to file dest until stream ends or EOF.
   inf() returns Z_OK on success, Z_MEM_ERROR if memory could not be
   allocated for processing, Z_DATA_ERROR if the deflate data is
   invalid or incomplete, Z_VERSION_ERROR if the version of zlib.h and
   the version of the library linked do not match, or Z_ERRNO if there
   is an error reading or writing the files. */
int inf(int fdSource, int fdDest)
{
    int ret;
    unsigned have;
    z_stream strm;
    unsigned char in[CHUNK];
    unsigned char out[CHUNK];

    /* allocate inflate state */
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    ret = inflateInit(&strm);
    if (ret != Z_OK)
        return ret;

    /* decompress until deflate stream ends or end of file */
    do {
        strm.avail_in = read(fdSource, in, CHUNK);
        if (strm.avail_in == -1) {
            (void)inflateEnd(&strm);
            return Z_ERRNO;
        }
        if (strm.avail_in == 0)
            break;
        strm.next_in = in;

        /* run inflate() on input until output buffer not full */
        do {
            strm.avail_out = CHUNK;
            strm.next_out = out;
            ret = inflate(&strm, Z_NO_FLUSH);
            assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
            switch (ret) {
            case Z_NEED_DICT:
                ret = Z_DATA_ERROR;     /* and fall through */
            case Z_DATA_ERROR:
            case Z_MEM_ERROR:
                (void)inflateEnd(&strm);
                return ret;
            }
            have = CHUNK - strm.avail_out;
            if (write(fdDest, out, have) != have) {
                (void)inflateEnd(&strm);
                return Z_ERRNO;
            }
        } while (strm.avail_out == 0);

        /* done when inflate() says it's done */
    } while (ret != Z_STREAM_END);

    /* clean up and return */
    (void)inflateEnd(&strm);
    return ret == Z_STREAM_END ? Z_OK : Z_DATA_ERROR;
}

/* report a zlib or i/o error */
void zerr(int ret)
{
    fputs("zpipe: ", stderr);
    switch (ret) {
    case Z_ERRNO:
        if (ferror(stdin))
            fputs("error reading stdin\n", stderr);
        if (ferror(stdout))
            fputs("error writing stdout\n", stderr);
        break;
    case Z_STREAM_ERROR:
        fputs("invalid compression level\n", stderr);
        break;
    case Z_DATA_ERROR:
        fputs("invalid or incomplete deflate data\n", stderr);
        break;
    case Z_MEM_ERROR:
        fputs("out of memory\n", stderr);
        break;
    case Z_VERSION_ERROR:
        fputs("zlib version mismatch!\n", stderr);
    }
}
