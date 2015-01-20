# Copyright © 2014, Battelle National Biodefense Institute (BNBI);
# all rights reserved. Authored by: Brian Ondov, Todd Treangen, and
# Adam Phillippy
#
# See the LICENSE.txt file included with this software for license information.

@0x953fca0ab70f1586;

using Cxx = import "/capnp/c++.capnp";
$Cxx.namespace("capnp");

struct Harvest
{
	struct ReferenceList
	{
		struct Reference
		{
			tag @0 : Text;
			sequence @1 : Text;
		}
	
		references @0 : List(Reference);
	}

	struct TrackList
	{
		struct Track
		{
			enum TrackType
			{
				none @0;
				genome @1;
				contigs @2;
				scaffolds @3;
				reads @4;
			}

			file @0 : Text;
			name @1 : Text;
			type @2 : TrackType;
			size @3 : UInt32;
		}

		tracks @0 : List(Track);
		variantReference @1 : UInt32 = 0;
	}

	struct Tree
	{
		struct Node
		{
			children @0 : List(Node);
			branchLength @1 : Float32;
			track @2 : UInt32; # 0-index to 'TrackList.tracks' (leaves only)
			bootstrap @3 : Float32;
		}
	
		root @0 : Node;
		multiplier @1 : Float32;
	}

	struct LcbList
	{
		struct Lcb
		{
			struct Region
			{
				track @0 : UInt32; # 0-index to 'TrackList.tracks'
				position @1 : UInt32;
				length @2 : UInt32;
				reverse @3 : Bool;
			}
		
			regions @0 : List(Region);
			sequence @1 : UInt32;
			position @2 : UInt32;
			name @3 : Text;
			concordance @4 : Float32;
			tree @5 : Tree ;
			length @6 : UInt32;
		}
	
		lcbs @0 : List(Lcb);
	}

	struct VariantList
	{
		struct Filter
		{
			flag @0 : UInt64; # power of 2; used for 'Variant.filters'
			name @1 : Text;
			description @2 : Text;
		}
	
		struct Variant
		{
			sequence @0 : UInt32; # 0-index to 'ReferenceList.references'
			position @1 : UInt32; # ungapped ref pos; dupl. for ins.
			alleles @2 : Text; # column of multialignment
			filters @3 : UInt64; # bit field of 'Filter.flag'
			quality @4 : UInt32;
			reference @5 : UInt8; # char; if ref is not in alignment (eg VCF)
		}
	
		filters @0 : List(Filter);
		variants @1 : List(Variant);
		defaultFilters @2 : UInt64; # bit field of 'Filter.flag'
	}

	struct AnnotationList
	{
		struct Annotation
		{
			struct AnnotationRegion
			{
				start @0 : UInt32; # leftmost ref pos
				end @1 : UInt32; # rightmost ref pos
				reverse @2 : Bool;
			}
		
			sequence @0 : UInt32; # 0-index to 'ReferenceList.references'
			regions @1 : List(AnnotationRegion); 
			reverse @2 : Bool;
			locus @3 : Text;
			name @4 : Text;
			description @5 : Text;
			feature @6 : Text;
		}

		annotations @0 : List(Annotation);
	}

	referenceList @0 : ReferenceList;
	trackList @1 : TrackList;
	lcbList @2 : LcbList;
	tree @3 : Tree;
	variantList @4 : VariantList;
	annotationList @5 : AnnotationList;
}
