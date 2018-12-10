import os
import re

# collection of useful snippets of code that's used frequently

nucleic = ['U', 'A', 'C', 'G', 'T']

class Singleton(type):
    """
        singleton metaclass that keeps track of instances
    """
    _instances = {}

    def __call__(self, *args, **kwargs):
        if self not in self._instances:
            #print "SINGLETON", self, type(self)
            self._instances[self] = {'obj': super(Singleton, self).__call__(*args, **kwargs),'count':0}
        self._instances[self]['count'] +=1
        return self._instances[self]['obj']

class Borg:
    __shared_state = {}
    def __init__(self):
        self.__dict__ = self.__shared_state
        # and whatever else you want in your class -- that's all!

def getNameExt(fname):
    """ extract name and extension from the input file, removing the dot
        filename.ext -> [filename, ext]
    """
    name, ext = os.path.splitext(fname)
    return name, ext[1:] #.lower()

def getLines(filename, doStrip = False, removeEmpty=False,
            removeCommentLines=False, removeAllComments=False):
    """
    doStrip             :   extra spaces
    removeEmpty         :   remove emtpy lines
    removeCommentLines  :   remove lines starting with "#"
    removeAllComments   :   truncate lines from the first occurrence of "#" on
    """
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()
    if doStrip:
        #lines = map(strip,lines)
        lines = [ x.strip() for x in lines ]
    if removeEmpty:
        #lines = removeEmptyLines(lines)
        lines = [ l for l in lines if l.strip() ]
    if removeCommentLines:
        lines = [ l for l in lines if not l.startswith("#") ]
    if removeAllComments:
        lines = [ l.split('#', 1)[0] for l in lines ]
    return lines

def writeList(filename, inlist, mode = 'w', addNewLine = False):
    if addNewLine: nl = "\n"
    else: nl = ""
    fp = open(filename, mode)
    for i in inlist:
        fp.write(str(i)+nl)
    fp.close()



def getResInfo(string):
    """ CHAIN:RESnum -> [ "CHAIN", "RES", num ]"""
    if ':' in string:
        chain, resraw = string.split(':')
    else:
        chain = ''
        resraw = string
    try:
        res = resraw[0:3]
        num = int(resraw[3:])
    except:
        # heuristic for nucleic acids
        regex = r'[UACGT]+\d'
        match = re.search(regex, resraw)
        if match == None:
            print "WARNING! Unknown residue naming scheme"
            return chain, "X", "X"
        res = resraw[0]
        num = int(resraw[1:])
        #print "NUCLEIC:",  chain, res, num
    return chain, res, num

def getRaccoonPath():
    import Raccoon
    return os.path.dirname(Raccoon.__file__)


def get_data_file(file_handle, dir_name, data_file):
    module_dir, module_fname = os.path.split(file_handle)
    DATAPATH = os.path.join(module_dir, dir_name, data_file)
    return DATAPATH

