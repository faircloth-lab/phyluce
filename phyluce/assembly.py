import os
import re
import glob


class Read():
    """Fastq reads"""
    def __init__(self, dir, file):
        self.dir = dir
        self.file = file
        self.pth = os.path.join(dir, file)

    def __str__(self):
        return "{} fastq read".format(self.file)

    def __repr__(self):
        return "<{}.{} instance at {}>".format(self.file, self.__class__.__name__, hex(id(self)))


class Fastqs():
    """Container for fastq data"""
    def __init__(self):
        self.r1 = None
        self.r2 = None
        self.singleton = None
        self.type = None
        self.reads = ()

    def __str__(self):
        return "Fastq container of R1, R2, Singletons"

    def __repr__(self):
        return "<{}.{} instance at {}>".format("test", self.__class__.__name__, hex(id(self)))

    def set_read(self, read, dir, file):
        if read == 'r1':
            self.r1 = Read(dir, file)
            self.reads += ((self.r1),)
        elif read == 'r2':
            self.r2 = Read(dir, file)
            self.reads += ((self.r2),)
        elif read == 'singleton':
            self.singleton = Read(dir, file)
            self.reads += ((self.singleton),)


def get_fastq_input_files(dir, subfolder, log):
    log.info("Finding fastq files")
    types = ("*.fastq.gz", "*.fastq.gzip", "*.fq.gz", "*fq.gzip", "*.fq", "*.fastq")
    files = []
    for type in types:
        files.extend(glob.glob(os.path.join(dir, subfolder, type)))
    if not files:
        raise IOError("There are not appropriate files in {}".format(dir))
    fq = Fastqs()
    # get dirname of first file
    dir = os.path.dirname(files[0])
    ext = set()
    for f in files:
        # get file extension
        ext.add(os.path.splitext(f)[-1])
        # get file name
        fname = os.path.basename(f)
        # find which reach this is
        match = re.search("(?:.*)[_-](?:READ|Read|R)(\d)*[_-]*(singleton)*(?:.*)", fname)
        try:
            if match.groups()[0] == '1':
                assert fq.r1 is None
                fq.set_read('r1', dir, fname)
            elif match.groups()[0] == '2':
                assert fq.r2 is None
                fq.set_read('r2', dir, fname)
            elif match.groups()[1] == 'singleton':
                assert fq.singleton is None
                fq.set_read('singleton', dir, fname)
        except:
            raise IOError("The appear to be multiple files for R1/R2/Singleton reads")
    if len(ext) != 1:
        raise IOError("Files are of different types (e.g. gzip and fastq)")
    if '.gzip' in ext or '.gz' in ext:
        fq.type = 'gzip'
    else:
        fq.type = 'fastq'
    return fq
