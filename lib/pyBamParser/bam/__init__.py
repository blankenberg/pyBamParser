#Dan Blankenberg
from ..bgzf import Reader as BGZFReader
from ..bgzf import Reader as BGZFWriter
from ..bai import Reader as BAIReader
from ..read import BAMRead
from ..util import NULL_CHAR
from ..util.odict import odict
from ..util.packer import pack_int8, unpack_int8, pack_uint8, unpack_uint8, pack_int16, unpack_int16, pack_uint16, unpack_uint16, pack_int32, unpack_int32, pack_uint32, unpack_uint32, pack_int64, unpack_int64, pack_uint64, unpack_uint64

BAM_MAGIC = b'BAM\x01'

SAM_HEADER_NON_TAB_RECORDS = [ '@CO' ]
SAM_READ_GROUP_RECORD_CODE = '@RG'
SAM_READ_GROUP_ID_STR = 'ID'
SAM_READ_GROUP_SAMPLE_STR = 'SM'
SAM_NO_READ_GROUP_NAME = '__NONE__'

class Reader( object ):
    def __init__( self, filename, index_filename=None ):
        self._bgzf_reader = BGZFReader( filename )
        self._references_list = []
        self._buffer = next(self._bgzf_reader)
        self._read_groups = None
        magic = self.read( 4 )
        assert magic == BAM_MAGIC, "Bad BAM Magic (%s) in %s" % ( magic, self._filename )
        l_text = unpack_int32( self.read( 4 ) )[0]
        self._headers = self.read( l_text ).rstrip( NULL_CHAR ).decode()
        n_ref = unpack_int32( self.read( 4 ) )[0]
        for i in range( n_ref ):
            l_name = unpack_int32( self.read( 4 ) )[0]
            name = self.read( l_name ).rstrip( NULL_CHAR ).decode()
            l_ref = unpack_int32( self.read( 4 ) )[0]
            self._references_list.append( ( name, l_ref ) )
        self._bam_index = BAIReader( index_filename or "%s.bai" % self._filename, self )
    
    @property
    def _filename( self ):
        return self._bgzf_reader.filename
    
    def read( self, size ):
        while len( self._buffer ) < size:
            self._buffer += next(self._bgzf_reader)
        data = self._buffer[ 0:size ]
        self._buffer = self._buffer[ size: ]
        return data
    
    def seek_virtual( self, offset_tuple ):
        file_offset, block_offset = offset_tuple
        self._bgzf_reader.seek( file_offset )
        self._buffer = next(self._bgzf_reader)[block_offset:]
    
    def __next__( self ):
        block_size = unpack_int32( self.read( 4 ) )[ 0 ]
        return BAMRead( self.read( block_size ), self )
    
    def __iter__( self ):
        return self
    
    def get_references( self ):
        return self._references_list
    
    def get_reference_by_id( self, ref_id ):
        if ref_id < 0:
            return ( None, None )
        else:
            return self._references_list[ ref_id ]
    
    def get_reference_id_by_name( self, name ):
        for i, (ref_name, ref_len ) in enumerate( self._references_list ):
            if name == ref_name:
                return i
        return None
    
    def get_reference_name_by_id( self, ref_id, other_id=None ):
        if ref_id < 0:
            return '*'
        if ref_id == other_id:
            return '='
        return self._references_list[ ref_id ][0]
    
    def get_sam_header_dict( self ):
        rval = odict()
        for line in self._headers.split( '\n' ):
            line = line.rstrip( '\r' )
            if line:
                rec_code, value = line.split( '\t', 1 )
                if rec_code not in rval:
                    rval[rec_code] = []
                if rec_code in SAM_HEADER_NON_TAB_RECORDS:
                    rval[rec_code].append( value )
                else:
                    fields = value.split( '\t' )
                    attribs = odict()
                    for field in fields:
                        try:
                            tag, value = field.split( ':', 1 )
                        except ValueError:
                            continue
                        attribs[tag] = value
                    rval[rec_code].append( attribs )
        return rval
    
    def get_read_groups( self ):
        if self._read_groups is None:
            headers = self.get_sam_header_dict()
            if SAM_READ_GROUP_RECORD_CODE in headers:
                self._read_groups = [ rg[SAM_READ_GROUP_ID_STR] for rg in headers[SAM_READ_GROUP_RECORD_CODE] ]
            else:
                self._read_groups = []
        return self._read_groups
    
    def jump( self, seq_name, start, fetch_next=True ):
        assert self._bam_index, Exception( "You must provide a valid BAM index in order to use the jump.")
        if isinstance( seq_name, str ):
            seq_id = self.get_reference_id_by_name( seq_name )
        else:
            seq_id = seq_name  
        if self._bam_index.jump_to_region( seq_id, start, start+1 ):
            offset = self._bgzf_reader.tell()
            buffer = self._buffer
            read = next(self)
            read_ref_id = read.get_reference_id()
            while read_ref_id <= seq_id:
                if read.get_end_position( one_based=False ) > start and read_ref_id == seq_id:
                    break
                offset = self._bgzf_reader.tell()
                buffer = self._buffer
                read = next(self)
                read_ref_id = read.get_reference_id()
            if fetch_next:
                return read
            self._bgzf_reader.seek( offset )
            self._buffer = buffer
            return True
        else:
            if next:
                return None
            else:
                return False
    
    def get_sam_header_text( self ):
        header_dict = self.get_sam_header_dict()
        if '@SQ' not in header_dict:
            header_dict.insert( 0, '@SQ', [] )
        sq_in_dict = []
        for sq in header_dict.get( '@SQ' ):
            sq_in_dict.append( sq.get( 'SN' ) )
        rval = ''
        for ref_name, ref_len in self._references_list:
            if ref_name not in sq_in_dict:
                sn_dict = odict()
                sn_dict['SN'] = ref_name #SN comes before LN by convention
                sn_dict['LN'] = ref_len
                header_dict[ '@SQ' ].append( sn_dict )
        rval = ''
        for rec_code, values in header_dict.items():
            for value in values:
                if rec_code in SAM_HEADER_NON_TAB_RECORDS:
                    str_val = '\t' + value
                else:
                    str_val = ''
                    for tag, val in value.items():
                        str_val = "%s\t%s:%s" % ( str_val, tag, val )
                rval = "%s%s%s\n" % ( rval, rec_code, str_val )
        rval = rval.strip( '\n\r' )
        #TODO: add @PG header
        return rval
    
    
class Writer( object ):
    
    def __init__( self, filename, headers=None, references=None ):
        self._writer = BGZFWriter( filename )
        self._headers = headers or ''
        self._references = references
        bam_header = "%s%s%s%s" % (BAM_MAGIC, pack_int32( len( self._headers ) ), self._headers, pack_int32( "<i", len( self._references ) ) )
        for ref_name, ref_length in self._references:
            bam_header = "%s%s%s%s%s" % ( bam_header, pack_int32( "<i", len( ref_name ) + 1 ), ref_name, NULL_CHAR, pack_int32( "<i", ref_length ) )
        self._writer.write( bam_header )
        self._writer.flush() #bam header will have its own bgzf blocks, increases speed for replacing header later
        self._alignment_start_offset = self._writer._fh.tell() #alignment blocks start here, make it easy to replace header
    
    def write( self, read ):
        if isinstance( read, BAMRead ):
            self._writer.write( read.get_bam_data() )
        else:
            raise NotImplementedError( 'this write type is not implemented yet' )
        
    def flush( self ):
        self._writer.flush()
        
    def close( self ):
        self.flush()
        self._writer.close()
