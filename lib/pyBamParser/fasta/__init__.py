#Dan Blankenberg
import os
import sys

from ..util import StringIO
from ..util.odict import odict

class IndexedReferenceSequences( object ):
    
    UNKNOWN_SEQUENCE_CHAR = "N"
    FAI_FIELDS = [ 'name', 'len', 'offset', 'line_blen', 'line_len' ]
    FAI_TYPE_MAP = [ str, int, int, int, int ]
    FIELD_LEN = len( FAI_FIELDS )
    DEFAULT_SEQUENCE_LENGTH = 255000000
    
    def __init__( self, filename, index_filename=None, sequence_filter=None, default_sequence_length=None ):
        self._filename = filename
        self._sequence_filter = sequence_filter or str
        if default_sequence_length is None:
            default_sequence_length = self.DEFAULT_SEQUENCE_LENGTH
        self._default_sequence_length = default_sequence_length
        if filename:
            self._fh = open( filename )
        else:
            self._fh = StringIO( '' ) #can we handle this better?
        self._index = None
        if index_filename is None:
            self._index_filename = "%s.fai" % filename
        else:
            self._index_filename = index_filename
        if os.path.exists( self._index_filename ):
            self._parse_index()
        else:
            self._create_index_dictionary()
    
    def _parse_index( self ):
        self._index = odict()
        for line in open( self._index_filename ):
            line = line.rstrip( '\n\r' )
            fields = line.split( '\t' )
            assert len( fields ) == self.FIELD_LEN, "Bad number of FAI fields: %s" % fields
            info = dict( ( x[0], x[1]( x[2] ) ) for x in zip( self.FAI_FIELDS, self.FAI_TYPE_MAP, fields ) )
            info['name'] = info['name'].split()[0] #remove extended id, only use before first whitespace
            self._index[ info['name'] ] = info
    
    def _create_index_dictionary( self ):
        self._index = odict()
        if self._filename and os.path.exists( self._filename ):
            start_offset = self._fh.tell()
            self._fh.seek( 0 )
            seq_name = None
            last_line_len_mismatch = False
            checking_for_trailing_new_lines = False 
            while True:
                data = self._fh.readline()
                if not data:
                    break
                if data.startswith( '>' ):
                    seq_name = data.split()[0][1:].strip()
                    self._index[seq_name] = dict( name=seq_name, len=0, offset=self._fh.tell(), line_blen=None, line_len=None ) #TODO: store all the fai data here
                    last_line_len_mismatch = False
                elif seq_name:
                    line_len = len( data )
                    data = data.strip()
                    if not data:
                        #empty new lines between or after sequences are ignored
                        seq_name = None
                        continue
                    line_blen = len( data )
                    assert not last_line_len_mismatch, last_line_len_mismatch
                    if self._index[seq_name]['line_blen'] is None:
                        self._index[seq_name]['line_blen'] = line_blen
                    if self._index[seq_name]['line_len'] is None:
                        self._index[seq_name]['line_len'] = line_len
                    if self._index[seq_name]['line_len'] != line_len:
                        last_line_len_mismatch = "line_len mismatch: %s != %s" % ( line_len, self._index[seq_name]['line_len'] )
                    if self._index[seq_name]['line_blen'] != line_blen:
                        last_line_len_mismatch = "line_blen mismatch: %s != %s" % ( line_blen, self._index[seq_name]['line_blen'] )
                    self._index[seq_name]['len']  = self._index[seq_name]['len'] + line_blen
                else:
                    if data.strip():
                        raise Exception( "Unexpected characters found in FASTA at position %s: %s" % ( self._fh.tell() - len( data ), data ) )
            self._fh.seek( start_offset )
    
    def get_sequence_by_position( self, sequence_name, position, length=1, unknown_sequence_character=UNKNOWN_SEQUENCE_CHAR, die_on_error=True ):
        assert length > 0, "You must provide a positive length."
        if self._index and sequence_name in self._index:
            if position > self._index[ sequence_name ]['len']:
                if die_on_error:
                    raise ValueError( 'Requested position (%i) is greater than reference sequence length (%i) for "%s".' % ( position, self._index[ sequence_name ]['len'], sequence_name ) )
                else:
                    sys.stderr.write( 'Requested position (%i) is greater than reference sequence length (%i) for "%s".' % ( position, self._index[ sequence_name ]['len'], sequence_name ) )
                position = self._index[ sequence_name ]['len'] - 1
            total_length = length
            end_position = position + length
            if  end_position > self._index[ sequence_name ]['len']:
                if die_on_error:
                    raise ValueError( 'Requested position (%i) and length (%i) is greater than reference sequence length (%i) for "%s".' % ( position, length, self._index[ sequence_name ]['len'], sequence_name ) )
                else:
                    sys.stderr.write( 'Requested position (%i) and length (%i) is greater than reference sequence length (%i) for "%s".\n' % ( position, length, self._index[ sequence_name ]['len'], sequence_name ) )
                end_position = self._index[ sequence_name ]['len']
            length = end_position - position 
            self._fh.seek( self._index[ sequence_name ]['offset'] + ( ( position / self._index[ sequence_name ]['line_blen'] ) * self._index[ sequence_name ]['line_len'] ) + ( position % self._index[ sequence_name ]['line_blen'] ) )
            rval = ''
            while len( rval ) < length:
                data = self._fh.read( length - len( rval ) )
                if not data:
                    break
                rval += "".join( data.split() )#remove any white spaces
                #TODO: fix to remove line numbers etc
            if len( rval ) < total_length:
                rval += unknown_sequence_character * ( total_length - len( rval ) )
            return self._sequence_filter( rval )
        return unknown_sequence_character * length #TODO: make this configurable to raise an error
    
    def get_sequence_size_by_name( self, name ):
        if name in self._index:
            return self._index[ name ]['len']
        return self._default_sequence_length
    
    def get_sequence_names( self ):
        return self._index.keys()
    
    def sort_region_list( self, region_list ):
        region_list = sorted( region_list ) #make copy and sort positions within seq names
        seq_names = self.get_sequence_names()
        if not seq_names:
            return region_list
        rval = []
        for name in seq_names:
            if not region_list:
                break
            remove = []
            for i in xrange( len( region_list ) ):
                if region_list[i][0] == name:
                    rval.append( region_list[i] )
                    remove.append( i )
            for i in reversed( remove ):
                del region_list[i]
        return rval
