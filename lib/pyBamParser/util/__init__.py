#Dan Blankenberg
NULL_CHAR = '\x00'

try:
    from io import StringIO
except ImportError:
    from io import StringIO

def get_filename_and_open( filename, default=None, mode='rb'):
    if isinstance( filename, str ):
        fh = open( filename, mode )
    else:
        fh = filename
        try:
            filename = filename.name
        except AttributeError:
            filename = default
    return filename, fh
