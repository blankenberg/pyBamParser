#Dan Blankenberg
import struct

ENDIAN = '<' #BAM is always little endian <

INT8 = '%sb' % ( ENDIAN )
UINT8 = '%sB' % ( ENDIAN )

INT16 = '%sh' % ( ENDIAN )
UINT16 = '%sH' % ( ENDIAN )

INT32 = '%si' % ( ENDIAN )
UINT32 = '%sI' % ( ENDIAN )

INT64 = '%sq' % ( ENDIAN )
UINT64 = '%sQ' % ( ENDIAN )


INT8_PACKER = struct.Struct( INT8 )
UINT8_PACKER = struct.Struct( UINT8 )
INT8_SIZE = INT8_PACKER.size
UINT8_SIZE = UINT8_PACKER.size

INT16_PACKER = struct.Struct( INT16 )
UINT16_PACKER = struct.Struct( UINT16 )
INT16_SIZE = INT16_PACKER.size
UINT16_SIZE = UINT16_PACKER.size

INT32_PACKER = struct.Struct( INT32 )
UINT32_PACKER = struct.Struct( UINT32 )
INT32_SIZE = INT32_PACKER.size
UINT32_SIZE = UINT32_PACKER.size

INT64_PACKER = struct.Struct( INT64 )
UINT64_PACKER = struct.Struct( UINT64 )
INT64_SIZE = INT64_PACKER.size
UINT64_SIZE = UINT64_PACKER.size

pack_int8 = INT8_PACKER.pack
pack_uint8 = UINT8_PACKER.pack
unpack_int8 = INT8_PACKER.unpack
unpack_uint8 = UINT8_PACKER.unpack
def unpack_int8_reader( reader ):
    return INT8_PACKER.unpack( reader.read( INT8_SIZE ) )
def unpack_uint8_reader( reader ):
    return UINT8_PACKER.unpack( reader.read( UINT8_SIZE ) )

pack_int16 = INT16_PACKER.pack
pack_uint16 = UINT16_PACKER.pack
unpack_int16 = INT16_PACKER.unpack
unpack_uint16 = UINT16_PACKER.unpack
def unpack_int16_reader( reader ):
    return INT16_PACKER.unpack( reader.read( INT16_SIZE ) )
def unpack_uint16_reader( reader ):
    return UINT16_PACKER.unpack( reader.read( UINT16_SIZE ) )

pack_int32 = INT32_PACKER.pack
pack_uint32 = UINT32_PACKER.pack
unpack_int32 = INT32_PACKER.unpack
unpack_uint32 = UINT32_PACKER.unpack
def unpack_int32_reader( reader ):
    return INT32_PACKER.unpack( reader.read( INT32_SIZE ) )
def unpack_uint32_reader( reader ):
    return UINT32_PACKER.unpack( reader.read( UINT32_SIZE ) )

pack_int64 = INT64_PACKER.pack
pack_uint64 = UINT64_PACKER.pack
unpack_int64 = INT64_PACKER.unpack
unpack_uint64 = UINT64_PACKER.unpack
def unpack_int64_reader( reader ):
    return INT64_PACKER.unpack( reader.read( INT64_SIZE ) )
def unpack_uint64_reader( reader ):
    return UINT64_PACKER.unpack( reader.read( UINT64_SIZE ) )
