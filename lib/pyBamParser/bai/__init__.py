#Dan Blankenberg
import sys

from ..util.odict import odict
from ..util.packer import pack_int8, unpack_int8, pack_uint8, unpack_uint8, pack_int16, unpack_int16, pack_uint16, unpack_uint16, pack_int32, unpack_int32, pack_uint32, unpack_uint32, pack_int64, unpack_int64, pack_uint64, unpack_uint64

BAI_MAGIC = 'BAI\x01'
BAI_MAX_BINS = 37450 #(((1<<18)-1)/7) + 1 ; (8**6-1)/7+1
BAI_BINS = [None] * BAI_MAX_BINS
BAI_WINDOW_SIZE = 2 ** 14

class Reader( object ):
    
    def __init__( self, filename, bam_reader ):
        self._filename = filename
        self._bam_reader = bam_reader
        try:
            self._fh = open( filename, 'rb ')
            self._nonzero = True
        except IOError:
            self._nonzero = False
        if self:
            assert self._fh.read( 4 ) == BAI_MAGIC, "Not a BAM Index: %s" % ( filename )
            self._n_ref = unpack_int32( self._fh.read( 4 ) )[0]
            self._references = []
            for i in range( self._n_ref ):
                bins = odict()#do we care about order?
                raw_bins = {}
                n_bins = unpack_int32( self._fh.read( 4 ) )[0]
                for j in range( n_bins ):
                    bin = unpack_uint32( self._fh.read( 4 ) )[0]
                    n_chunk = unpack_int32( self._fh.read( 4 ) )[0]
                    bins[ bin ] = []
                    raw_bins[ bin ] = []
                    for k in range( n_chunk ):
                        chunk_beg = unpack_uint64( self._fh.read( 8 ) )[0]
                        raw_bins[ bin ].append( chunk_beg )
                        chunk_beg = ( chunk_beg >> 16, chunk_beg & 0xFFFF )
                        chunk_end = unpack_uint64( self._fh.read( 8 ) )[0]
                        raw_bins[ bin ].append( chunk_end )
                        chunk_end = ( chunk_end >> 16, chunk_end & 0xFFFF )
                        bins[ bin ].append( ( chunk_beg, chunk_end ) )
                n_intv = unpack_int32( self._fh.read( 4 ) )[0]
                intv = []
                for l in range( n_intv ):
                    offset = unpack_uint64( self._fh.read( 8 ) )[0]
                    if offset: #only append non-zero offsets
                        intv.append( ( offset >> 16, offset & 0xFFFF ) )
                self._references.append( { 'bins': bins, 'intv': intv } )
                # determine idxstats
                idxstats = raw_bins.get( BAI_MAX_BINS, [] )
                idxstats2 = bins.get( BAI_MAX_BINS, [] )
                if len( idxstats ) == 4:
                    idxstats = { 'unmapped_beg': idxstats2[0][0], 'unmapped_end': idxstats2[0][1], 'n_mapped': idxstats[2], 'n_unmapped': idxstats[3] }
                else:
                    idxstats = { 'unmapped_beg': None, 'unmapped_end': None, 'n_mapped': None, 'n_unmapped': None }
                self._references[-1].update( idxstats )
            #unaligned count is optional
            unaligned_count = self._fh.read( 8 )
            if unaligned_count:
                self._unaligned_count = unpack_uint64( unaligned_count )[ 0 ]
            else:
                self._unaligned_count = 0
            
            extra = self._fh.read()
            if extra:
                print("BAM Index appears to be malformed: %s." % ( repr( extra ) ), file=sys.stderr)
    def __bool__( self ):
        return self._nonzero 
    
    def _fix_region( self, seq_id, start, end ):
        ref = self._bam_reader.get_reference_by_id( seq_id )
        start = start or 0
        if end <= start or end > ref[1]:
            end = ref[1]
        return start, end
    
    def jump_to_region( self, seq_id, start, end ):
        if isinstance( seq_id, str ):
            seq_id = self._bam_reader.get_reference_id_by_name( seq_id )
        if seq_id is None or seq_id < 0 or seq_id >= len( self._references ):
            return False
        
        start, end = self._fix_region( seq_id, start, end )
        seq_bins = self._references[ seq_id ]

        try:
            linear = offset = seq_bins['intv'][start / BAI_WINDOW_SIZE ]
        except Exception:
            linear = offset = -1

        bin_list = self.reg2bins( start, end )
        for bin in seq_bins['bins']:
            if bin in bin_list:
                for left, right in seq_bins['bins'][bin]:
                    if right >= linear:
                        if left < offset:
                            offset = left

        if offset == -1:
            return False
        self._bam_reader.seek_virtual( offset )
        return True
    
    def reg2bins( self, beg, end ):
        #calculate the list of bins that may overlap with region [beg,end) (zero-based)
        bins = [0]
        if beg >= end:
            return bins
        for i in range( 1 + ( beg >> 26 ), 1 + ( end >> 26  ) ):
            bins.append( i )
        
        for i in range( 9 + ( beg >> 23 ), 9 + ( end >> 23 ) ):
            bins.append( i )
        
        for i in range( 73 + ( beg >> 20 ), 73 + ( end >> 20 ) ):
            bins.append( i )
        
        for i in range( 585 + ( beg >> 17 ), 585 + ( end >> 17 ) ):
            bins.append( i )
        
        for i in range( 4681 + ( beg >> 14 ), 4681 + ( end >> 14 ) ):
            bins.append( i )
        return bins
    
    def reg2bin( self, beg, end ):
        #calculate bin given an alignment covering [beg,end) (zero-based, half-close-half-open)
        end -= 1#end
        if (beg>>14 == end>>14):
            return ((1<<15)-1)/7 + (beg>>14)
        if (beg>>17 == end>>17):
            return ((1<<12)-1)/7 + (beg>>17)
        if (beg>>20 == end>>20):
            return ((1<<9)-1)/7 + (beg>>20)
        if (beg>>23 == end>>23):
            return ((1<<6)-1)/7 + (beg>>23)
        if (beg>>26 == end>>26):
            return ((1<<3)-1)/7 + (beg>>26)
        return 0
