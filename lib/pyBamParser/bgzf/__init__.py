#Dan Blankenberg
import zlib
import struct
import sys
import os.path
import string

from ..util import get_filename_and_open
from ..util.packer import pack_int8, unpack_int8, pack_uint8, unpack_uint8, pack_int16, unpack_int16, pack_uint16, unpack_uint16, pack_int32, unpack_int32, pack_uint32, unpack_uint32, pack_int64, unpack_int64, pack_uint64, unpack_uint64

MAX_NEG_INT = -sys.maxsize

WBITS = -15 #Negative wbits removes zlib header usage. Larger numbers can always work on files that were compressed by smaller numbers, but not vis-versa, so always use the largest available number (15)
BGZF_MAGIC = b'\x1f\x8b\x08\x04'
BGZF_WRITE_HEADER = BGZF_MAGIC + b'\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00'
BGZF_EOF = BGZF_WRITE_HEADER + b'\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00'
BGZF_MAX_BLOCK_SIZE = 2**16

HEADER_BLOCK_SIZE = 25
CRC_PYTHON_COMPATIBILITY = 0xffffffff

MTIME_XFL_OS_XLEN_SI1_SI2_SLEN_BSIZE_UNPACKER = struct.Struct( "<IBBHBBHH" ).unpack

class Reader( object ):
    
    def __init__( self, filename ):
        self.filename, self.fh = get_filename_and_open( filename, mode='rb' )
    
    def seek( self, offset ):
        return self.fh.seek( offset )
    
    def tell( self ):
        return self.fh.tell()
    
    def __iter__(self):
        return self

    def __next__( self ):
        magic = self.fh.read( 4 )
        if magic:
            assert magic == BGZF_MAGIC, "Bad BGZF magic, are you sure this is a BGZF file (%s)?" % ( self.filename )
        else:
            self.fh.seek( -28, 2 )
            eof = self.fh.read( 28 )
            if eof != BGZF_EOF:
                print('BGZF EOF marker was not found. Confirm that the BGZF file is not truncated.', file=sys.stderr)
                raise StopIteration( 'End of file, but BGZF EOF marker was not found. Confirm that the BGZF file is not truncated.' )
            raise StopIteration( 'End of file.' )
        mtime, xfl, OS, xlen, si1, si2, slen, bsize = MTIME_XFL_OS_XLEN_SI1_SI2_SLEN_BSIZE_UNPACKER( self.fh.read( 14 ) )
        cdata = self.fh.read( bsize - xlen - 19 ) #Compressed DATA by zlib::deflate() uint8 t[BSIZE-XLEN-19]
        crc = unpack_uint32( self.fh.read( 4 ) )[0]
        isize = unpack_uint32( self.fh.read( 4 ) )[0]
        data = zlib.decompress( cdata, WBITS ) #decompress the data, no headers
        assert isize == len( data ), "Invalid decompressed data size"
        return data

class Writer( object ):
    
    def __init__( self, filename, compress_level=6 ):
        self._filename, self._fh = get_filename_and_open( filename, mode='wb' )
        self._compress_level = compress_level
        self._buffer = b''
        self._compressor = self._new_compressor()
        self._fh_is_open = True
        self.get_new_compressor = self._compressor.copy #supposedly faster than creating a new one? But need to test
    
    def _new_compressor( self ):
        return zlib.compressobj( self._compress_level, zlib.DEFLATED, WBITS, zlib.DEF_MEM_LEVEL, 0 )
    
    def _write_buffer( self, block_size=None ):
        if not block_size:
            block_size = BGZF_MAX_BLOCK_SIZE
        if len( self._buffer ) >= block_size:
            data = self._buffer[:block_size]
            self._buffer = self._buffer[block_size:]
            len_data = block_size
        else:
            data = self._buffer
            self._buffer = b''
            len_data = len( data )
        compressor = self.get_new_compressor()
        compressed_data = compressor.compress( data ) + compressor.flush()
        len_compressed = len( compressed_data )
        if len_compressed > BGZF_MAX_BLOCK_SIZE:
            #compressed too big, retry with half blocksize
            del compressed_data
            del len_compressed
            del len_data
            self._buffer = data + self._buffer
            del data
            return self._write_buffer( block_size = block_size/2 )
        
        crc = pack_uint32( zlib.crc32( data ) & CRC_PYTHON_COMPATIBILITY ) #Note: To generate the same numeric value across all Python versions and platforms use crc32(data) & 0xffffffff. If you are only using the checksum in packed binary format this is not necessary as the return value is the correct 32bit binary representation regardless of sign.
        bsize = pack_uint16( HEADER_BLOCK_SIZE + len_compressed )
        isize = pack_uint32( len_data )
        self._fh.write( BGZF_WRITE_HEADER + bsize + compressed_data + crc + isize ) #FIXME: adjust BGZF_WRITE_HEADER for compression level?
        
    def write( self, data ):
        self._buffer += data
        if len( self._buffer ) >= BGZF_MAX_BLOCK_SIZE:
            return self._write_buffer()
    
    def flush( self ):
        while self._buffer:
            self._write_buffer()
        self._fh.flush()
            
    def close( self ):
        if self._fh_is_open:
            self.flush()
            self._fh.write( BGZF_EOF )
            self._fh.close()
            self._fh_is_open = False
        
    def __del__( self ):
        self.close() #need to close file, so that EOF is written
