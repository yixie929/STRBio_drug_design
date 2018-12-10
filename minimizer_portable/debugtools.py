#       
#           AutoDock | Raccoon2
#
#       Copyright 2013, Stefano Forli
#          Molecular Graphics Lab
#  
#     The Scripps Research Institute 
#           _  
#          (,)  T  h e
#         _/
#        (.)    S  c r i p p s
#          \_
#          (,)  R  e s e a r c h
#         ./  
#        ( )    I  n s t i t u t e
#         '
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import inspect
from datetime import datetime
import sys
import logging
from generic import Singleton

# XXX TODO make it a singleton?
# http://stackoverflow.com/questions/8237390/can-a-python-singleton-class-be-inherited

# Use the standard Python logging:
# https://fangpenlin.com/posts/2012/08/26/good-logging-practice-in-python/

class DebugMessageDispatcher:
    """ Gods' messenger"""
    __metaclass__ = Singleton
    def __init__(self, caller=None, outfile=None):
        """ """
        self.caller = caller
        self.outfile = outfile
        if not outfile == None:
            try:
                self.fileStream = open(outfile, 'a', 0)
                self.dispatch("\n=[ %s %s log file opened ]===========================================" % (datetime.now(), self.outfile) )
            except:
                print "FATAL ERROR! Impossible to open the filestream [%s]: %s" % (self.outfile, sys.exc_info()[1])
                sys.exit(1)
        else:
            self.fileStream = None
        self.dispatch("DEBUG[%s] %s() initialized in DEBUG MODE" % (datetime.now(), self.caller))

    def __del__(self):
        """ destructor, called when garbage collecting
            to close the filestream
        """
        #print "DEBUG MERCURY CALLED"
        #print "TESTING THIS", self._instances_count
        #print ">>>", self.__metaclass__._instances.keys()
        #print "TEST", self.__class__.__name__ in self.__metaclass__._instances.keys()
        #print "TEST", type(self).__name__ in self.__metaclass__._instances.keys()
        #print "TEST", type(self) in self.__metaclass__._instances.keys()
        if self.__metaclass__._instances[type(self)]['count'] > 1:
        #if self.__metaclass__._instances[type(self)]['count'] > 1:
            # if we're not the last one, we're good...
            #print "THIS IS NOT THE LAST",
            #print self.__metaclass__._instances[type(self)]['count']
            self.__metaclass__._instances[type(self)]['count'] -= 1
            #self._instances[type(self)]['count'] > 1:
            #self._instances_count -= 1
            return
        if self.fileStream == None:
            return
        self.dispatch("=[ %s %s log file closed ]===========================================" % (datetime.now(), self.outfile) )
        self.fileStream.flush()
        self.fileStream.close()
        #self.__metaclass__.__del__()

    def dispatch(self, string):
        """ performs the actual dispatch of the message """
        print string
        if not self.fileStream == None:
            self.fileStream.write( "%s\n" % string)

class DebugObj:
    """
    General class that can be subclassed by every object
    to provide logging/debugging functionalities

    The class provides superpowers to inspect the call
    stack and identify who called a method, etc...

    OPTIONS
        debug = False       : debug messages disabled
              = True        : debug messages enabled
              = string      : debug log file

        logfile = None      : debug message printed on screen
                = string    : interpreted as filename, opened, messages written to file
                = file object: : messages written to file

    dprint()  : debug-level dependent print function
              return the string "function [CALLER:function_caller]> message"

    The message delivery (printing, file logging) is delegated to
    DebugDispatcher, which is a singleton which opens a file stream
    and prints messages
    
    """
    def __init__(self, debug=False, logfile=None):
        """ """
        # XXX when initialize the singleton,
        # the calling function should be used in the callerinitialization
        self.debug = debug
        self._dispatcher = None
        self._dbug_spacer = "%s\n" % ("="*70)
        if self.debug == False:
            return
        caller = self.caller()
        logfile = None
        #if isinstance(debug, str):
        #    logfile = debug
        self._dispatcher = DebugMessageDispatcher(caller, logfile)
        # XXX SOURCE: http://stackoverflow.com/questions/582056/getting-list-of-parameter-names-inside-python-function
        # XXX SOURCE: http://stackoverflow.com/questions/8315389/how-do-i-print-functions-as-they-are-called
        # XXX SOURCE: http://stackoverflow.com/questions/2654113/python-how-to-get-the-callers-method-name-in-the-called-method
        # XXX PERFORMANCE: http://stackoverflow.com/questions/5067604/determine-function-name-from-within-that-function-without-using-traceback
        """
        for s in inspect.stack(4):
            a,v,k,l =  inspect.getargvalues(s[0])
            print "AVKL", a,k,v
            #,v,k,l
        """
        """
        frame = inspect.currentframe()
        print "ANALYSIS", type(frame)
        # XXX useful: print "ARGV", inspect.getargvalues(frame)
        a,v,k,l =  inspect.getargvalues(frame)
        print "ARG", a
        print "VAL", v
        print "KWD", k
        print "LOC", l
        #print "\n>>>>>> ARGV", inspect.getargvalues(frame)
        stack = inspect.stack(10)
        print "________________"
        for s in stack:
            print "====", inspect.getargvalues(s[0])
        print "\n\n"
        for a in range(10):
            prev = inspect.getouterframes(frame, a)
            print "\t\tPREV", inspect.getargvalues(prev)
        print "-----------------"
        print "GOT THIS", inspect.getframeinfo(frame)
        args, _, _, values = inspect.getargvalues(frame)
        print "ARGS VALS", args, values
        # XXX arguments passed to the caller should go here
        """

        """
   V:  F   T  's'
 D:  ----------------
    |
  F |  0 start start
    |
  T | stop 0  start
    |
 's'|
    |
    # SOURCE: http://stackoverflow.com/questions/2825452/correct-approach-to-validate-attributes-of-an-instance-of-class
    @property
    def debug(self):
        return self._debug
    @debug.setter
    def debug(self, value):
        if not self._debug:
            if value:
                self.startDebugger(value)
            self._debug = value
            return
        if isinstance(value
        self._debug = value
    """

    def __del__(self):
        """ """
        #print "DEBUG OBJ destroier"
        if not self._dispatcher == None:
            self._dispatcher.__del__()

    def caller(self, level=1):
        try:
            """
            currFrame = inspect.currentframe()
            for a in range(0,5):
                calframe = inspect.getouterframes(currFrame, a)
                print "MASTER", calframe
                for x in calframe:
                    print "[%d] =] " % a, x[0:4]
                print "********************"
            print "------------------------------"
            print "XXXX", inspect.stack()[0][3]
            print "YYYY", inspect.stack()[1][3]
            print "ZZZZ", inspect.stack()[2][3]
            print "WWWW", inspect.stack()[3][3]
            print "EEEE", inspect.stack()[4][3]
            print "RRRR", inspect.stack()[5][3]
            #print "CACACA", inspect.stack()[1][3]
            """
            stack = inspect.stack()
            #print "STACK", stack
            _class = stack[1][0].f_locals["self"].__class__
            _method = stack[3][0].f_code.co_name
            #print "\n\t\t", the_class, the_method
            #return inspect.stack()[2+level][3]
            return "%s.%s" % (_class, _method)
        except:
            return "(unknown)"

    def warning(self, msg):
        """ using the logger"""
        self.logger.warn(msg)

    def info(self, msg):
        """ """
        self.logger.info(msg)

    def dprint(self, msg, new=False, args=False, section=False):
        """
            this is where the magic happens
            TODO: implement indentation?
            TODO: indentation should be made by looking at the length of the stack

            it should be written as a recursive function using getouterframes
            until we hit __main__ or something similar...
        """
        # DOC: https://pymotw.com/2/inspect/
        if not self.debug: return

        # if level == 'info':
        #   self.logger.info(msg)
        # elif level == 'debug':
        #   self.logger.debug(msg)
        # elif level == 'warning':
        #   self.logger.warn(msg)
        # elif level == 'critical'
        #   self.logger.critical(msg)

        

        # inspect the stack frame
        curframe = inspect.currentframe()
        # frame object, the filename where the code exists, the line number in 
        # that file for the current line being run, the function name being called, 
        # a list of lines of context from the source file, and the index into that 
        # list of the current line

        # measure the depth of the function
        depth = 0
        cf = curframe
        while not cf == None:
            cf = cf.f_back
            depth+=1
        depth = max(0, depth-2)
        #    print "DONE recursing", c, sys.exc_info()[1]
        #    print cf, type(cf)
            
        #print "CURRENTFRAME", curframe.f_lineno, curframe.f_back,curframe.f_locals, dir(curframe)
        calframe = inspect.getouterframes(curframe, 2)
        fname = calframe[1][3]
        if args == True:
            _arg, _, _, _val = inspect.getargvalues(curframe)
            a_v = ",".join(['%s = %s' % (x, _val[x]) for x in _arg])
        else:
            a_v = ""
        # add newline if necessary
        if new == True: 
            nl = '\n\n'
        else:
            nl = ''
        # composing the message
        bar = ''
        if section:
            bar = self._dbug_spacer
        msg = "%s%sDEBUG[%s] %s(%s) [CALLER:%s()] >>%s%s" % (nl,bar, datetime.now(), fname, a_v, self.caller(), "--"*depth, msg)
        self._dispatcher.dispatch(msg)

    def dwritepdbfromcoord(self, coord, atype=[], title=None, fname=None):
        """ write coordinates as pdb entries"""
        buff = []
        if len(atype):
            if not len(coord) == len(atype):
                print "ERROR! coord count [%d] mismatch with atype count [%d]" % (len(coord), len(atype))
        else:
            atype == None
        for i in range(len(coord)):
            c = coord[i]
            if not atype == None:
                a = atype[i]
            pdb = self.makePdb(c, keyw='HETATM', at_index =i+1, atype=a)
            buff.append(pdb)
        if not fname == None:
            fp = open(fname, 'w')
            for b in buff:
                fp.write(b + '\n')
            fp.close()
        else:
            if not title:
                title = 'PDB coordinates'
            self.dprint("==============[ DEBUG: %s ]===========" % title)
            for b in buff:
                self.dprint(b)

    def debug_atomToPdb(self, atom):
        """ used by objects that manage OBMol and OBAtoms"""
        coord = atom.GetX(), atom.GetY(), atom.GetZ()
        return self.makePdb(coord)


    def makePdb(self, coord, keyw = "ATOM  ", at_index = 1, res_index = 1, atype = 'X', elem = None,
                res = "CNT", chain  ="Z", bfactor = 10,pcharge = 0.0):
        """ function that can be used for debugging code generating coordinates,
            which can be written as PDB structures and visualized
        """
        if not elem: elem = atype
        # padding bfactor
        bfactor = "%2.2f" % bfactor
        if len(bfactor.split(".")[0]) == 1:
            bfactor = " "+bfactor
        if len(atype) == 1:
            atype = atype + " "
        atom = "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 %s  %8.3f %1s" % (keyw,
                at_index, elem, res, chain, res_index, 
                coord[0], coord[1], coord[2], bfactor, pcharge, atype)
        return atom
