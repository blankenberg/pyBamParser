#Dan Blankenberg
NULL_CHAR = '\x00'

try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO

def get_filename_and_open( filename, default=None, mode='rb'):
    if isinstance( filename, basestring ):
        fh = open( filename, mode )
    else:
        fh = filename
        try:
            filename = filename.name
        except AttributeError:
            filename = default
    return filename, fh
