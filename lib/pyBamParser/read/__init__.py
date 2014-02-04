#Dan Blankenberg
import struct

from ..util.packer import pack_int8, unpack_int8, pack_uint8, unpack_uint8, pack_int16, unpack_int16, pack_uint16, unpack_uint16, pack_int32, unpack_int32, pack_uint32, unpack_uint32, pack_int64, unpack_int64, pack_uint64, unpack_uint64
from ..util import NULL_CHAR

BAM_NO_QUAL = 0xFF #255
SEQ_4_BIT_TO_SEQ = list( '=ACMGRSVTWYHKDBN' )

CIGAR_OP = list( 'MIDNSHP=X' )

NULL_TERMINATED_TAGS = [ 'Z', 'H' ]
TAG_TYPE_TO_VALUE_LENGTH = {
                            'A': 1,
                            'c': 1,
                            'C': 1,
                            's': 2,
                            'S': 2,
                            'i': 4,
                            'I': 4,
                            'f': 4,
                            'Z': 1,
                            'H': 2,
                            'B': 1
                            }
TAG_TYPE_TO_STRUCT_TYPE = {
                            'A': 'c',
                            'c': 'b',
                            'C': 'B',
                            's': 'h',
                            'S': 'H',
                            'i': 'i',
                            'I': 'I',
                            'f': 'f',
                            'Z': 'c',
                            'H': 'c',
                            'B': 'c'
                            }

TAG_TYPE_TO_SAM_TYPE = {
                            'A': 'A',
                            'c': 'i',
                            'C': 'i',
                            's': 'i',
                            'S': 'i',
                            'i': 'i',
                            'I': 'i',
                            'f': 'f',
                            'Z': 'Z',
                            'H': 'H',
                            'B': 'B'
                            }

BAM_READ_BEGIN_UNPACKER = struct.Struct( "<iiIIiiii" ).unpack
READ_GROUP_RECORD_TAG = 'RG'

SEQ_UNPACKERS = {} #cache seq unpackers for reuse
QUAL_UNPACKERS = {}
CIGAR_UNPACKERS = {}

#cache seq bit to 4 bit lookup
SEQ_4_BIT_TO_SEQ_4 = [ SEQ_4_BIT_TO_SEQ[ i >> 4 ]  for i in range( 256 ) ]
SEQ_4_BIT_TO_SEQ_F = [ SEQ_4_BIT_TO_SEQ[ i & 0x0f ]  for i in range( 256 ) ]



#TODO: FIXME: make parsing occur on-demand, by attribute, not at initialization
class BAMRead( object ):
    
    def __init__( self, data, reader ): #TODO: instead of reader, can we take list of sequences? can we omit block size?
        self._block_size = len( data )
        self._reader = reader
        self._ref_id, self._pos, self._bin_mq_nl, self._flag_nc, self._l_seq, self._next_ref_id, self._next_pos, self._t_len = BAM_READ_BEGIN_UNPACKER( data[ :32 ] )
        self._bin = self._bin_mq_nl >> 16
        self._mapq = self._bin_mq_nl >> 8 & 0xff
        self._l_read_name = self._bin_mq_nl & 0xff
        self._flag = self._flag_nc >> 16
        self._n_cigar_op = self._flag_nc & 0xff
        
        self.__data = data
        self.__not_parsed_1 = True
        self.__not_parsed_2 = True
        self.__not_parsed_3 = True
        self.__not_parsed_4 = True
        self.__not_parsed_5 = True
        
        self.__zero_based_end_position = None
        self.__seq_string = None
        self.__qual_list = None
        self.__is_seq_reverse_complement = None
        self.__read_group = None
        self.__read_group_parsed = False
        self.__reference_name = None
        self.__reference = None

    def __parse_block_1( self ):
        if self.__not_parsed_1:
            self.__not_parsed_1 = False
            self._block_offset = 32 + self._l_read_name
            self._read_name = self.__data[ 32:self._block_offset].rstrip( NULL_CHAR )
            
    def __parse_block_2( self ):
        if self.__not_parsed_2:
            self.__parse_block_1()
            self.__not_parsed_2 = False
            n_cigar_op = self._n_cigar_op
            cigar_op_len = 4 * n_cigar_op
            self._new_block_offset = self._block_offset + cigar_op_len
            cigar_unpacker = CIGAR_UNPACKERS.get( n_cigar_op, None )
            if cigar_unpacker is None:
                cigar_unpacker = struct.Struct( "<" +"I" * n_cigar_op ).unpack
                CIGAR_UNPACKERS[ n_cigar_op ] = cigar_unpacker
            self._cigar = cigar_unpacker( self.__data[ self._block_offset:self._new_block_offset ] )
            self._cigar_list = []
            for cigar in self._cigar:
                op_len = cigar >> 4
                op = cigar & 0x07
                self._cigar_list.append( ( op_len, op ) )
            
    def __parse_block_3( self ):
        if self.__not_parsed_3:
            self.__parse_block_2()
            self.__not_parsed_3 = False
            seq_bin_len = ( self._l_seq + 1 ) / 2
            self._block_offset = self._new_block_offset
            self._new_block_offset += seq_bin_len
            seq_unpacker = SEQ_UNPACKERS.get( seq_bin_len, None )
            if seq_unpacker is None:
                seq_unpacker = struct.Struct( "<" + "B" * seq_bin_len ).unpack
                SEQ_UNPACKERS[ seq_bin_len ] = seq_unpacker #cache this unpacker for use later
            self._seq = seq_unpacker( self.__data[ self._block_offset: self._new_block_offset ] )
            half_seq_len = self._l_seq / 2
            seq = []
            for i, s in enumerate( self._seq, start=1 ):
                seq.append( SEQ_4_BIT_TO_SEQ_4[ s ] )
                if i <= half_seq_len:
                    seq.append( SEQ_4_BIT_TO_SEQ_F[ s ] )
                else:
                    break
            self._seq_list = seq
            self._block_offset = self._new_block_offset
            self._new_block_offset += self._l_seq
            
    def __parse_block_4( self ):
        if self.__not_parsed_4:
            self.__parse_block_3()
            self.__not_parsed_4 = False
            qual_unpacker = QUAL_UNPACKERS.get( self._l_seq, None )
            if qual_unpacker is None:
                qual_unpacker = struct.Struct( "<" + "c" * self._l_seq ).unpack
                QUAL_UNPACKERS[ self._l_seq ] = qual_unpacker
            self._qual = map( ord, qual_unpacker( self.__data[ self._block_offset: self._new_block_offset ] ) )
            if self._qual[0] == BAM_NO_QUAL:
                self._qual = None
            
    def __parse_block_5( self ):
        if self.__not_parsed_5:
            self.__parse_block_4()
            self.__not_parsed_5 = False
            self._aux_data = []
            while self._block_size > self._new_block_offset:
                self._block_offset = self._new_block_offset
                self._new_block_offset += 2
                tag = self.__data[ self._block_offset: self._new_block_offset ]
                self._block_offset = self._new_block_offset
                self._new_block_offset += 1
                val_type = self.__data[ self._block_offset: self._new_block_offset ]
                if val_type in NULL_TERMINATED_TAGS:
                    value = []
                    while True:
                        self._block_offset = self._new_block_offset
                        self._new_block_offset += 1
                        val = self.__data[ self._block_offset: self._new_block_offset ]
                        if val == NULL_CHAR:
                            break
                        else:
                            value.append( val )
                    value = ( "".join( value ), )
                else:
                    if val_type == 'B':
                        self._block_offset = self._new_block_offset
                        self._new_block_offset += 1
                        val_type = self.__data[ self._block_offset: self._new_block_offset ]
                        self._block_offset = self._new_block_offset
                        self._new_block_offset += 4
                        tag_length = unpack_int32( self.__data[ self._block_offset: self._new_block_offset ] )[ 0 ]
                    else:
                        tag_length = 1
                    val_size = TAG_TYPE_TO_VALUE_LENGTH[ val_type ] * tag_length
                    self._block_offset = self._new_block_offset
                    self._new_block_offset += val_size
                    value = self.__data[ self._block_offset: self._new_block_offset ]
                    value = struct.unpack( "<" + TAG_TYPE_TO_STRUCT_TYPE[ val_type ] * tag_length, value )
                self._aux_data.append( ( tag, val_type, value ) )
            self.__data = None
    
    def get_read_name( self ):
        self.__parse_block_1()
        return self._read_name
    def _get_bam_read_name( self ):
        return self.get_read_name() + NULL_CHAR
    
    def get_flag( self ):
        return self._flag
    
    def get_reference( self ):
        if self.__reference is None:
            self.__reference = self._reader.get_reference_by_id( self._ref_id )
        return self.__reference
        
    def get_reference_name( self  ):
        if self.__reference_name is None:
            self.__reference_name = self._reader.get_reference_name_by_id( self._ref_id )
        return self.__reference_name
    
    def get_reference_id( self ):
        return self._ref_id
    
    def _get_bam_ref_id( self ):
        return pack_int32( self._ref_id )
    
    def get_rnext( self ):
        return self._reader.get_reference_by_id( self._next_ref_id )
    
    def get_rnext_name( self ):
        return self._reader.get_reference_name_by_id( self._next_ref_id, self._ref_id  )
    
    def _get_bam_rnext_id( self ):
        return pack_int32( self._next_ref_id )
    
    def get_position( self, one_based=True ):
        return self._pos + one_based
    def get_position_zero_based( self ):
        return self._pos
    def _get_bam_pos( self ):
        return pack_int32( self.get_position( False ) )

    def indel_at( self, position, check_insertions=True, check_deletions=True, one_based=True ):
        """Does the read contain an indel at the given position?
        Return True if the read contains an insertion at the given position
        (position must be the base before the insertion event) or if the read
        contains a deletion where the base at position is deleted. Return False
        otherwise."""
        (insertions, deletions) = self.get_indels( one_based=one_based )
        if check_insertions:
            for insertion in insertions:
                if insertion[0] == position:
                    return True
        if check_deletions:
            for deletion in deletions:
                if deletion[0] < position < deletion[0] + deletion[1] + 1:
                    return True
        return False

    def get_indels( self, one_based=True ):
        """Return a data structure containing all indels in the read.
        Returns the tuple (insertions, deletions)
        insertions = [(pos1,ins1), (pos2,ins2)]
        posN = start position (preceding base, VCF-style)
        insN = length of inserted sequence (not including preceding base)
        deletions = [(pos1,del1), (pos2,del2)]
        posN = start position (preceding base, VCF-style)
        delN = length of deleted sequence (not including preceding base)
        """
        cigar = self.get_cigar() #CIGAR_OP = list( 'MIDNSHP=X' )
        insertions = []
        deletions = []
        position_offset = 0
        position_start = self.get_position( one_based=one_based )
        while cigar:
            cigar_size, cigar_op = cigar.pop( 0 )
            if cigar_op in [ 0, 7, 8 ]: #M alignment match (can be a sequence match or mismatch); = sequence match; x sequence mismatch
                position_offset += cigar_size
            elif cigar_op == 1: #I insertion
                insertions.append((position_start + position_offset - 1, cigar_size))
            elif cigar_op == 2: #D deletion from the reference
                deletions.append((position_start + position_offset - 1, cigar_size))
                position_offset += cigar_size
            elif cigar_op == 3: #N skipped region from the reference
                position_offset += cigar_size
            elif cigar_op == 4: #S soft clipping (clipped sequences present in SEQ)
                pass
            elif cigar_op == 5: #H hard clipping (clipped sequences NOT present in SEQ)
                position_offset += cigar_size
            elif cigar_op == 6: #P padding (silent deletion from padded reference)
                pass
            else: #unknown cigar_op
                print >>sys.stderr, 'unknown cigar_op', cigar_op, cigar_size
        return (insertions, deletions)
    
    def get_end_position( self, one_based=True ):
        if self.__zero_based_end_position is None:
            position_offset = 0
            self.__parse_block_2()
            for cigar_size, cigar_op in self._cigar_list:
                if cigar_op in [ 0, 7, 8 ]: #M alignment match (can be a sequence match or mismatch); = sequence match; x sequence mismatch
                    position_offset += cigar_size
                elif cigar_op == 1: #insertion
                    pass
                elif cigar_op == 2: #D deletion from the reference
                    position_offset += cigar_size
                elif cigar_op == 3: #N skipped region from the reference
                    #print >>sys.stderr, 'passing cigar_op N skipped', cigar_op, cigar_size
                    #print >>sys.stderr, self.to_sam()
                    position_offset += cigar_size
                elif cigar_op == 4: #S soft clipping (clipped sequences present in SEQ)
                    #print >>sys.stderr, 'passing cigar_op soft clipping', cigar_op, cigar_size
                    #print >>sys.stderr, self.to_sam()
                    pass
                elif cigar_op == 5: #H hard clipping (clipped sequences NOT present in SEQ)
                    #print >>sys.stderr, 'passing cigar_op hard clipping', cigar_op, cigar_size
                    #print >>sys.stderr, self.to_sam()
                    position_offset += cigar_size
                elif cigar_op == 6: #P padding (silent deletion from padded reference)
                    #print >>sys.stderr, 'passing cigar_op padding', cigar_op, cigar_size
                    #print >>sys.stderr, self.to_sam()
                    #position_offset += cigar_size
                    pass
                else: #unknown cigar_op
                    print >>sys.stderr, 'unknown cigar_op', cigar_op, cigar_size
                    #position_offset += cigar_size
            self.__zero_based_end_position = self._pos + position_offset
        return self.__zero_based_end_position + one_based
    
    def get_pnext( self, one_based=True ):
        return self._next_pos + one_based
    def _get_bam_next_pos( self ):
        return pack_int32( self.get_pnext( False ) )
    
    def get_mapq( self ):
        return self._mapq
    
    def get_cigar( self ):
        self.__parse_block_2()
        return self._cigar_list[:]
    def get_sam_cigar( self ):
        self.__parse_block_2()
        return "".join( "%s%s" % ( l, CIGAR_OP[o] ) for l, o in self._cigar_list )
    
    def get_t_len( self ):
        return self._t_len
    def _get_bam_t_len( self ):
        return pack_int32( self._t_len )
    
    
    def get_seq( self ):
        if self.__seq_string is None:
            self.__parse_block_3()
            self.__seq_string = "".join( self._seq_list )
        return self.__seq_string
    def _get_bam_seq( self ):
        self.__parse_block_3()
        return struct.pack( "<" + "B" * (  ( self._l_seq + 1 ) / 2 ), *self._seq )
    
    def get_l_seq( self ):
        return self._l_seq
    def _get_bam_seq_length( self ):
        self.__parse_block_3()
        return pack_int32( len( self._seq_list ) )
    
    def get_qual( self ):
        self.__parse_block_4()
        return self._qual
    def get_qual_list( self ):
        if self.__qual_list is None:
            self.__qual_list = self.get_qual()
            if not self.__qual_list:
                self.__qual_list = [ BAM_NO_QUAL for i in range ( self._l_seq ) ]
        return self.__qual_list[:]
    def get_sam_qual( self ):
        return "".join( chr( c + 33 ) for c in self.get_qual() )
    def get_qual_tuple( self ):
        return tuple( self.get_qual_list() )
    def _get_bam_qual( self ):
        return struct.pack( "<" + "c" * self._l_seq, *map( chr, self.get_qual_tuple() ) )
    
    def _get_bam_bin_mq_nl( self ):
        return pack_uint32( self._get_bin_mq_nl() )
    
    def _get_bin_mq_nl( self ):
        return self._bin_mq_nl
    
    def _get_bam_flag_nc( self ):
        return pack_uint32( self._get_flag_nc() )
    def _get_flag_nc( self ):
        #TODO: fix me to calculate unless set and parsed already
        return self._flag_nc
    
    def _get_bam_n_cigar_op( self ):
        return struct.pack( "<" +"I" * self._n_cigar_op, *self._get_cigar() )
    
    def _get_cigar( self ):
        self.__parse_block_2()
        return self._cigar
    
    def is_seq_reverse_complement( self ):
        if self.__is_seq_reverse_complement is None:
            self.__is_seq_reverse_complement = ( self.get_flag() & 0x0010 == 0x0010 )
        return self.__is_seq_reverse_complement
    
    def get_read_group( self ):
        if self.__read_group_parsed is False:
            self.__parse_block_5()
            for tag, val_type, value in self._aux_data:
                if tag == READ_GROUP_RECORD_TAG:
                    self.__read_group = value[0]
                    break
            self.__read_group_parsed = True
        return self.__read_group
    def get_sam_aux( self ):
        self.__parse_block_5()
        rval = ''
        for aux in self._aux_data:
            rval += "\t%s:%s:%s" % ( aux[0], TAG_TYPE_TO_SAM_TYPE[ aux[1] ], ",".join( map( str, aux[2] ) ) )
        return rval.strip( '\t' )
    def _get_bam_aux( self ):
        self.__parse_block_5()
        data = ''
        for tag, val_type, value in self._aux_data:
            data += tag
            tag_length = len( value )
            if tag_length > 1:
                data = '%sB%s' ( data, val_type )
                data += pack_int32( tag_length )
            else:
                data += val_type
                if val_type in NULL_TERMINATED_TAGS:
                    data += value[0]
                    data += NULL_CHAR
                    tag_length = None
            if tag_length:
                data += struct.pack( "<" + TAG_TYPE_TO_STRUCT_TYPE[ val_type ] * tag_length, *value )
        return data
    
    def get_bam_data( self ):
        #FIX ME: have these calculate from updatable properties
        rval = "%s%s%s%s%s%s%s%s%s%s%s%s%s" % ( self._get_bam_ref_id(), self._get_bam_pos(), self._get_bam_bin_mq_nl(), 
                                                self._get_bam_flag_nc(), self._get_bam_seq_length(), self._get_bam_rnext_id(),
                                                self._get_bam_next_pos(), self._get_bam_t_len(), self._get_bam_read_name(),
                                                self._get_bam_n_cigar_op(), self._get_bam_seq(), self._get_bam_qual(),
                                                self._get_bam_aux() )
        rval = "%s%s" % ( pack_int32( len( rval ) ), rval )
        return rval
    
    def to_sam( self ):
        rval = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % ( self.get_read_name(), self.get_flag(), self.get_reference_name(), self.get_position(), self.get_mapq(), self.get_sam_cigar(), self.get_rnext_name(), self.get_pnext(), self.get_t_len(), self.get_seq(), self.get_sam_qual() )
        aux = self.get_sam_aux()
        if aux:
            rval = "%s\t%s" % ( rval, aux )
        return rval
