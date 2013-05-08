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

#TODO: FIXME: make parsing occur on-demand, by attribute, not at initialization
class BAMRead( object ):
    def __init__( self, data, reader ): #TODO: instead of reader, can we take list of sequences? can we omit block size?
        block_size = len( data )
        self._reader = reader
        self._ref_id, self._pos, self._bin_mq_nl, self._flag_nc, self._l_seq, self._next_ref_id, self._next_pos, self._t_len = BAM_READ_BEGIN_UNPACKER( data[ :32 ] )
        self._bin = self._bin_mq_nl >> 16
        self._mapq = self._bin_mq_nl >> 8 & 0xff
        self._l_read_name = self._bin_mq_nl & 0xff
        self._flag = self._flag_nc >> 16
        self._n_cigar_op = self._flag_nc & 0xff
        
        block_offset = 32 + self._l_read_name
        self._read_name = data[ 32:block_offset].rstrip( NULL_CHAR )
        cigar_op_len = 4 * self._n_cigar_op
        new_block_offset = block_offset + cigar_op_len
        self._cigar = struct.unpack("<" +"I" * self._n_cigar_op, data[ block_offset:new_block_offset ] )
        self._cigar_list = []
        for cigar in self._cigar:
            op_len = cigar >> 4
            op = cigar & 0x07
            self._cigar_list.append( ( op_len, op ) )
        seq_bin_len = ( self._l_seq + 1 ) / 2
        block_offset = new_block_offset
        new_block_offset += seq_bin_len
        self._seq = struct.unpack( "<" + "B" * seq_bin_len, data[ block_offset: new_block_offset ] )
        seq = []
        for i, s in enumerate( self._seq ):
            seq.append( SEQ_4_BIT_TO_SEQ[ s >> 4 ] )
            if ( i + 1 ) * 2 <= self._l_seq:
                seq.append( SEQ_4_BIT_TO_SEQ[ s & 0x0f ] )
        self._seq_list = seq
        block_offset = new_block_offset
        new_block_offset += self._l_seq
        self._qual = map( ord, struct.unpack( "<" + "c" * self._l_seq, data[ block_offset: new_block_offset ] ) )
        if self._qual[0] == BAM_NO_QUAL:
            self._qual = None
        self._aux_data = []
        while block_size > new_block_offset:
            block_offset = new_block_offset
            new_block_offset += 2
            tag = data[ block_offset: new_block_offset ]
            block_offset = new_block_offset
            new_block_offset += 1
            val_type = data[ block_offset: new_block_offset ]
            if val_type in NULL_TERMINATED_TAGS:
                value = []
                while True:
                    block_offset = new_block_offset
                    new_block_offset += 1
                    val = data[ block_offset: new_block_offset ]
                    if val == NULL_CHAR:
                        break
                    else:
                        value.append( val )
                value = ( "".join( value ), )
            else:
                if val_type == 'B':
                    block_offset = new_block_offset
                    new_block_offset += 1
                    val_type = data[ block_offset: new_block_offset ]
                    block_offset = new_block_offset
                    new_block_offset += 4
                    tag_length = unpack_int32( data[ block_offset: new_block_offset ] )[ 0 ]
                else:
                    tag_length = 1
                val_size = TAG_TYPE_TO_VALUE_LENGTH[ val_type ] * tag_length
                block_offset = new_block_offset
                new_block_offset += val_size
                value = data[ block_offset: new_block_offset ]
                value = struct.unpack( "<" + TAG_TYPE_TO_STRUCT_TYPE[ val_type ] * tag_length, value )
            self._aux_data.append( ( tag, val_type, value ) )
    
    def get_read_name( self ):
        return self._read_name
    def _get_bam_read_name( self ):
        return self._read_name + NULL_CHAR
    
    def get_flag( self ):
        return self._flag
    
    def get_reference( self ):
        return self._reader.get_reference_by_id( self._ref_id )
        
    def get_reference_name( self  ):
        return self._reader.get_reference_name_by_id( self._ref_id )
    
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
    def _get_bam_pos( self ):
        return pack_int32( self.get_position( False ) )
    
    def get_end_position( self, one_based=True ):
        cigar = self.get_cigar() #CIGAR_OP = list( 'MIDNSHP=X' )
        position_offset = 0
        while cigar:
            cigar_size, cigar_op = cigar.pop( 0 )
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
                position_offset += cigar_size
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
        return position_offset + self.get_position( one_based=one_based )
    
    def get_pnext( self, one_based=True ):
        return self._next_pos + one_based
    def _get_bam_next_pos( self ):
        return pack_int32( self.get_pnext( False ) )
    
    def get_mapq( self ):
        return self._mapq
    
    def get_cigar( self ):
        return list( self._cigar_list )
    def get_sam_cigar( self ):
        return "".join( "%s%s" % ( l, CIGAR_OP[o] ) for l, o in self._cigar_list )
    
    def get_t_len( self ):
        return self._t_len
    def _get_bam_t_len( self ):
        return pack_int32( self._t_len )
    
    
    def get_seq( self ):
        return "".join( self._seq_list )
    def _get_bam_seq( self ):
        return struct.pack( "<" + "B" * (  ( self._l_seq + 1 ) / 2 ), *self._seq )
    
    def get_l_seq( self ):
        return self._l_seq
    def _get_bam_seq_length( self ):
        return pack_int32( len( self._seq_list ) )
    
    def get_qual( self ):
        return self._qual
    def get_sam_qual( self ):
        return "".join( chr( c + 33 ) for c in self.get_qual() )
    def get_qual_tuple( self ):
        if not self._qual:
            return tuple( [ BAM_NO_QUAL for i in range ( self._l_seq ) ] )
        return tuple( self._qual )
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
        #fix
        return self._cigar
    
    def get_read_group( self ):
        for tag, val_type, value in self._aux_data:
            if tag == 'RG':
                return value[0]
        return None
    def get_sam_aux( self ):
        rval = ''
        for aux in self._aux_data:
            rval += "\t%s:%s:%s" % ( aux[0], TAG_TYPE_TO_SAM_TYPE[ aux[1] ], ",".join( map( str, aux[2] ) ) )
        return rval.strip( '\t' )
    def _get_bam_aux( self ):
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
