#!/usr/bin/env python

import numpy as np
import math
import string
import sys
import argparse 
import os
import itertools
import multiprocessing as mp

# import psutil

import pybel
ob = pybel.ob
from scrub import Scrub
from chiralenumerator import ChiralEnumerator
from debugtools import DebugObj
OPTIMALAMIDEDIHE = math.pi
ZEROISH = 0.005
# format used to generate string of molecule in the processing queue 
"""
ScrubMT and MolWriter work in separate threads:
 - ScrubMT receives data from queueIn, generates optimal coordinates and puts 
    (molString, fname) in queueOut. Multiple instances of ScrubMT are created
    in multiprocessing. Each time ScrubMT receives a poison pill, it will pass 
    it to the MolWriter, then will exit.

 - MolWriter reads 'molString' and writes it in fname; if MolWriter is initialized 
    with a fname in the constructor, then fhane specified in the queue is ignored 
    and all outputs are written in the constructor single fname (self.fp).
    When initialized, the number of ScrubMT instances is defined, so that when
    the same number of poison pills is received, MolWriter will flush the data
    in the output file, then will exit.
"""

class ScrubMT(mp.Process, Scrub):
    def __init__(self, queueIn, queueOut, nice=None, inType='mol2', outType='mol2', scrubOpts={}):
        mp.Process.__init__(self)
        self.inType = inType
        self.outType = outType
        if 'auto' in scrubOpts:
            scrubOpts['auto'] = True
        Scrub.__init__(self, **scrubOpts)
        self.queueIn = queueIn
        self.queueOut = queueOut
        self.nice = nice
        self._c = 0

    def run(self):
        # XXX psutil missing in MGLTools
        # nice from here: https://andrewbolster.info/2014/05/multiprocessing-niceness-in-python
        if not self.nice == None:
            p = psutil.Process(os.getpid())
            niceness = os.nice(0)
            os.nice(self.nice - niceness)
        while True:
            # TODO
            # this should become: name, string, fname
            string, fname = self.queueIn.get()
            if string == None:
                # poison pill!
                # print "Poison pill is gone"
                self.queueIn.task_done()
                self.queueOut.put((None, None))
                break
            # this should become (name, string, self.inType, self.outType)
            # stringOut = self.process(string, self.inType, self.outType)
            stringOut = self.process(string, self.inType, self.outType, uniqueAtomNames=True)
            self.queueOut.put((stringOut, fname))
            self.queueIn.task_done()
            self._c +=1
        # DEBUG: 
        # print "Structures processed [%d]" % self._c
        return

# TODO make sure that debug is in place

class MolWriter(mp.Process):
    def __init__(self, queue, threads, fname=None, ext=None):
        """ molecular writer process, used to write 
            molecules from the queue.

            by default, the queue provides tuples as
            (mol_string, fname), which will be written in
            the provided filename.

            If fname is provided at the constructor, 
            a file pointer is created when the class is 
            instanciated (self.fp). The 'fname' from the queue will
            be ignored and all results will be appended in
            self.fp. 

            if 'ext' is provided, it will be used as extension
            for the output files from the queue
        """
        mp.Process.__init__(self)
        self.queue = queue
        self.threads = threads
        self.fp = None
        if not fname == None:
            self.fp = open(fname, 'wb', 0)
        self.ext = ext
        self._c=0

    def run(self):
        """ """
        while True:
            # TODO XXX handle error messages here 
            string, fname = self.queue.get()
            if string == None:
                # poison pills are coming!
                self.threads -= 1
                if self.threads == 0:
                    break
                else:
                    #print "Pill received...", self.threads
                    pass
            else:
                self._c+=1
                if not (self.fp == None):
                    self.fp.write(string)
                else:
                    fp = open("%s.%s" % (fname, self.ext), 'wb', 0)
                    fp.write(string)
                    fp.close()
        if not self.fp == None:
            self.fp.close()
        # DEBUG
        #print "Structures written [%d]" % self._c
        return 

class MolecularHound:
    """ a(n allegedly) efficient multi-SMARTS matcher"""
    def __init__(self, patterns):
        """ """
        self.patterns = patterns

    def sniffMolecule(self, mol):
        """ """
        matchers = self.createMatchers()
        c = 1
        for m in matchers:
            if m.Match(mol):
                return c
            c += 1
        return False

    def createMatchers(self):
        """ generator of OB matchers"""
        for p in self.patterns:
            m = ob.OBSmartsPattern()
            m.Init(p)
            yield m

class MiniMee(DebugObj):

    # XXX update documentation: usemolname + outfile creates prefix_name
    # XXX add a property filter?

    def __init__(self, debug=False):
        """ """
        DebugObj.__init__(self, debug)
        self.progname = 'Raccoon Scrub'
        self.setDefaults()
        # init stuff
        self.initOBChoices()
        self.initDocs()
        self.initOptions()
        self.initOptParser()
        self.parseOptions()
        self.start()

###    def vprint(self, string, newline=True):
###        """ verbose level printing"""
###        if self.verbose:
###            msg = "VERBOSE: %s" % string
###            if newline:
###                msg += "\n"
###            print msg,


    def setDefaults(self):
        """ set all predefined options"""
        # in/out
        self.default_out = 'mol2'
        self.default_informat = None
        self.default_outformat = None
        self.default_usemolname = None
        self.default_usemolnamesafe = None

        # slicing
        self.default_split = False
        self.default_only = False
        self.default_begin = 1
        self.default_end = None
            
        # mol processing
        self.default_pH = 7.4
        self.default_nopH = False
        self.default_charges = 'gasteiger'
        self.default_noflipamide = False
        self.default_nostripsalts = False
        self.default_nocheckhydro = False
        self.default_noprocess = False
        self.default_enumchiral = 'off' # all, undefined, protomer, off
        self.default_maxenumchiral = 10
        self.default_nomini = False
        self.default_noextra = False
        self.default_norotamer = False
        self.default_exclude = None

        # minimizer
        self.default_sdconv = 1e-5
        self.default_sdsteps = 300
        self.default_cgconv = 1e-6
        self.default_cgsteps = 300
        self.default_forcefield = 'mmff94s'

        # extra steps
        #self.default_sdconv = 1e-9
        #self.default_sdsteps = 5000000
        #self.default_cgconv = 1e-10
        #self.default_cgsteps = 10000000

        # rotamer search
        self.default_rotamer_conf = 10
        self.default_rotamer_steps = 5

        # minimizer (post-processing)
        self.default_sdconv_extra = 1e-5
        self.default_sdsteps_extra = 1000
        self.default_cgconv_extra = 1e-6
        self.default_cgsteps_extra = 1000
        self.default_forcefield_extra = 'uff'

        self.default_strict = False # abort minimization if no parameters?
        # logging and info
        self.default_log = None
        self.default_verbose = False

        # multiprocessor
        self.default_multiproc_mode = 'all'
        print "XXX URGENT: add psutil to the packages"
        #self.default_multiproc_max = psutil.cpu_count()
        self.default_multiproc_max = 1
        self.multiproc_max = self.default_multiproc_max
        # NICE XXX TODO 
        # psutil
        # source: https://stackoverflow.com/questions/1023038/change-process-priority-in-python-cross-platform
        #self.default_multiproc_nice = psutil.

        self._allowedOutFormats = ['mol2', 'pdb', 'pdbqt', 'sdf', 'ent', 
            'gpr', 'mdl', 'mol', 'mop', 'sd']
        self._validFnameChars = "-_.()%s%s" % (string.ascii_letters, string.digits)


    def initOptParser(self):
        """ initialize the parser"""
        # http://stackoverflow.com/questions/3853722/python-argparse-
        # how-to-insert-newline-the-help-text
        class SmartFormatter(argparse.HelpFormatter):
            def _split_lines(self, text, width):
                # this is the RawTextHelpFormatter._split_lines
                if text.startswith('R|'):
                    return text[2:].splitlines()  
                return argparse.HelpFormatter._split_lines(self, text, width)
        self.parser = argparse.ArgumentParser(description=self.desc,
            usage=self.usage, epilog=self.epilog,
            formatter_class = argparse.RawDescriptionHelpFormatter)
        for gName, gDesc in self.groupsOrder:
            group = self.parser.add_argument_group(gName, gDesc)
            for name in self._argsorder[gName]:
                opts = self._argsDict[name]
                group.add_argument(name, **opts)
        # help advanced checking
        if '--help_advanced' in sys.argv:
            self.showAdvancedHelp()
            sys.exit(0)
            
        self.args = self.parser.parse_args()
        # activate verbose as soon as possible, if requested
        self.verbose = self.args.verbose
        self._data = []


    def printMsg(self, msg, mtype = 'e', q=None):
        """ utility to report messages """
        if mtype == 'e':
            buff = '*** ERROR *** '
        elif mtype == 'w':
            buff = '*** WARNING *** '
        qs = self._getString(q)
        print buff + msg + qs 
        if mtype == 'e':
            sys.exit(1)
        
            
    def _getString(self,q):
        """ Base/64"""
        return ""
        return self._data[q]



    def parseOptions(self):
        """ perform the parsing of options"""
        self.parseLogging()
        self.parseMolNames()
        self.parseSlicing()
        self.parseInputOutput()
        #self.updateOutputFname()
        self.parseMinimizer()
        self.parseMolProcessing()
        self.parseExclude()
        self.parseMultiProc()
        
    def parseLogging(self):
        """ set verbose mode and logging"""
        #self.verbose = self.args.verbose
        if self.verbose:
            print "==========VERBOSE MODE================\nOPTS"
            print "Called with following settings:"
            for gName, gDesc in self.groupsOrder:
                for i in self._argsorder[gName]:
                    print i, "=>", getattr(self.args,i[2:])
            print "-" * 30
            print "Supported forcefields:",
            print self.fflist_names
            print "Supported charge models:",
            print self.chargelist_names
            print "======================"
        # XXX INIT LOG FILE POINTER HERE
        # logging
        log = self.args.log
        #if log == None:
        #    self.log = self.outname + '.log'

    def parseMolNames(self):
        """ check what's the policy for molecule names
        
            if using molecule names is requested, split is 
            switched on, and the dictionary to keep track
            of molecules is created
        """
        self.usemolname =  self.args.usemolname
        self.safenames =  self.args.usemolnamesafe
        if self.safenames:
            self.usemolname = True
        self.molNames = {}
        if self.usemolname:
            self.args.split = True
            self.dprint("[usemolname] output file names will use molecule name if found")
            if self.safenames:
                self.dprint("[safenames] molecule names will be sanitized")
            self.dprint("[split] split is forced to be ON by usemolname")

    def parseInputOutput(self):
        """ parse input and output file management"""
        self.parseInput()
        self.parseMulti()
        self.parseOutput()

    def parseInput(self):
        """ manage the input file name options"""
        # infile
        self.infile = self.args.infile
        self.inname, self.inext = self.getNameExt(self.infile)
        self.informat = self.args.informat
        self.intype = None

        if (self.inext == '') and (self.informat == None):
            msg = ('The input file has no extension and no input '
                   'format has been specified.\nUse \'--informat\' to '
                   'specify a supported format type\.\n'
                   '"One goes to the right, the other to the left; ' 
                   'both are wrong, but in different directions." '
                   '-- Quintus Horatius Flaccus, Satires (II,3,5)'
                   )
            print msg
            sys.exit(1)

        if self.informat:
            self.intype = self.informat
            self.dprint("[parseInput] input file type defined as: %s" % self.informat)
        else:
            self.intype = self.inext
            msg = ("[parseInput] input file name \"%s\" (guessed type: %s)") 
            self.dprint( msg %  (self.infile, self.inext))
            self.informat = self.intype


    def parseMulti(self):
        """ perform the multi-structure check"""
        # check for multi-structure 
        try:
            self.multi = self.isMulti(self.infile)
            if self.multi == True:
                msg = ("[parseMulti] Input file is a multi-structure file")
                self.dprint(msg)
            elif self.multi == False:
                self.args.split = False
                msg = ('(parseMulti) Input file is not a multi-structure '
                       'file (forcing --split off)')
                self.dprint(msg)
            elif self.multi == None:
                print "*** WARNING ***\nand now?"
                    
        except IOError as e:
            msg =  'An error occurred while reading file [%s]:' % self.infile
            msg += 'I/O error (%d) : %s' % (e.errno, e.strerror)
            self.printMsg(msg)

    
    def parseOutput(self):
        """ parse output information to determine output filename
            and format.

        Priorities in the assignments are:

        Filename:
            outname and usemolname specified ->  %outname_%usemolname
            outname only specified           -> %outname_(auto_progressive naming)
            no outname specified             -> %inname 

        File type:
            outformat specified -> %outformat
            outfile specified   -> %outfile_extension
            no outfile specified-> %informat (if valid); default otherwise

        File extension:
            (non-empty) -> use specified file extension
            (empty) + outformat -> keep empty
            (empty)             -> use informat if valid; default otherwise
        """
        # check if the output file is actually a directory
        self.outdir = None
        self.outfile = self.args.outfile
        if os.path.isdir(self.outfile):
            self.outdir = self.outfile
            self.outfile = ''
        self.outformat = self.args.outformat.lower()

        self.outtype = ''
        self.outname, self.outext = self.getNameExt(self.outfile)

        #### file name
        if not len(self.outname) and not self.usemolname:
                # "%inname_processed"
                self.outname = '%s_processed' % self.inname 
        ### file type
        if self.outformat:
            # requested format
            if self.outformat in self._allowedOutFormats:
                self.outtype = self.outformat
                #self.outext = self.outformat
            else:
                msg = ('Specified format [%s] is not an acceptable output format')
                self.printMsg(msg)
        else:
            if len(self.outext):
                # format from output file extension
                if self.outext.lower() in self._allowedOutFormats:
                    self.outtype = self.outext.lower()
                # format from input type (if acceptable)
            elif self.informat in self._allowedOutFormats:
                    self.outtype = self.informat
            else:
                # default format type
                self.dprint('[parseOutput] output set to default (%s)' % self.default_out)
                self.outtype = self.default_out

        # file extension
        if len(self.outext) == 0:
            if self.outformat:
                self.outext = self.outformat
                pass
            elif self.informat in self._allowedOutFormats:
                self.outext = self.informat
            else:
                self.outext = self.default_out
        self.updateOutputFname()


    def parseSlicing(self):
        """ check that all splicing options are well formed"""
        self.single = single = self.args.single
        begin = self.args.begin
        end = self.args.end
        self.split = self.args.split
        # single
        if self.single:
            if begin or end:
                msg =  "Use either --single or --begin/--end"
                self.printMsg(msg, q=1)
            self.begin = single
            self.end = single
            msg = "[parseSlicing] processing only molecule #%d (%s,%s)"
            self.dprint(msg % (single,self.begin,self.end))
        # begin/end
        else:
            if begin or end:
                self.dprint("[parseSlicing] begin/end detected (%s,%s)" % (begin, end))
            self.begin = ( begin or self.default_begin)
            self.end = end
            if self.end and (self.begin > self.end):
                msg = ('The value specified in --begin should be '
                         'smaller than the one used in --end')
                self.printMsg(msg, q=2)
            self.dprint("[parseSlicing] processing molecule range:(%s,%s)" % (begin,end))

    
    def updateOutputFname(self):
        """ check various option combinations to provide
            the appropriate file name
        """
        self.outfile = self.outname
        # updating naming scheme (redundant 'usemolname', but good reminder)
        if self.split or self.single or self.usemolname: 
            if len(self.outfile):
                self.outfile += "_%s"
            else:
                self.outfile = "%s"
        # create output filename
        if len(self.outext):
            self.outfile = "%s.%s" % (self.outfile, self.outext)
        if not self.outdir == None:
            self.outfile = self.outdir + os.path.sep + self.outfile

        self.dprint("[updateOutputFname] the output file is \"%s\"" % self.outfile)


    def parseMinimizer(self):
        """ minimizer parameters (forcefield, charges, minimizer steps)"""
        # get options 
        self.nomini = self.args.nomini
        self.noextra = self.args.noextra
        self.norotamer = self.args.norotamer
        self.rotamer_conf = self.args.rotamer_conf
        self.rotamer_steps = self.args.rotamer_steps

        # minimizer
        self.sdsteps = self.args.sdsteps
        self.sdconv = self.args.sdconv
        self.cgsteps = self.args.cgsteps
        self.cgconv = self.args.cgconv
        self.forcefield = self.args.forcefield.lower()
        # minimizer extra
        self.sdsteps_extra = self.args.sdsteps_extra
        self.sdconv_extra = self.args.sdconv_extra
        self.cgsteps_extra = self.args.cgsteps_extra
        self.cgconv_extra = self.args.cgconv_extra
        self.forcefield_extra = self.args.forcefield_extra.lower()
        self.chargemodel = self.args.chargemodel.lower()

        # charges (check them here, they're going to be used always)
        if not self.chargemodel in self.chargelist:
            msg = ('Specified charge model [%s] is not recognized. '
                   'Supported charge models are: %s.') % (self.chargemodel,
                   self.chargelist_str)
            self.printMsg(msg)

        self.checkMiniParms()
        self.checkMiniParmsExtra()
        self.checkRotamerParms()

    def checkRotamerParms(self):
        """ check rotamer search parameters"""
        if self.norotamer:
            anyparm = (self.rotamer_conf > 0 or self.rotamer_steps > 0)
            if anyparm:
                msg = ('Rotameric search disabled (--norotamer) but rotamer '
                      'parameters defined (--rotamer_conf/--rotamer_steps) requested')
                self.printMsg(msg, q=5)
            self.rotamer_conf = 0
            self.rotamer_steps = 0
        else:
            if self.rotamer_conf == None:
                self.rotamer_conf = self.default_rotamer_conf
            if self.rotamer_steps == None:
                self.rotamer_steps = self.default_rotamer_steps
            if (self.rotamer_conf == 0) and (self.rotamer_steps == 0):
                self.dprint("[parseRotamer] rotameric search disabled")
            if (self.rotamer_conf < 0) or (self.rotamer_steps < 0):
                msg = ('Number of rotamer conformers or steps negative!')
                self.printMsg(msg, q=5)

    def checkMiniParms(self):
        # disable minimization
        if self.nomini:
            anysteps = (self.sdsteps > 0 or self.cgsteps > 0)
            if anysteps:
                msg = ('Minimization disabled (--nomini) but minimizer '
                      'steps defined (--sdsteps/--cgsteps) requested')
                self.printMsg(msg, q=5)
                # printmsg will exit by now
            self.sdsteps = 0
            self.cgsteps = 0
            self.dprint('[parseMinimizer] minimization disabled (--nomini)')
            return
        else:
            if self.sdsteps == None:
                self.sdsteps = self.default_sdsteps
            if self.sdconv == None:
                self.sdconv = self.default_sdconv
            if self.cgsteps == None:
                self.cgsteps = self.default_cgsteps
            if self.cgconv == None:
                self.cgconv = self.default_cgconv           
        # steepest descent
        if self.sdsteps == 0:
            self.dprint("[parseMinimizer] SD disabled")
        if self.sdsteps < 0:
            msg = ('SD steps are negative! (seriously? A negative number of steps '
                    'for the minimizer? Try again)')
            self.printmsg(msg)
        # conjugate gradient
        if self.cgsteps == 0:
            self.dprint("[parseMinimizer] CG disabled (0 steps)")
        if self.cgsteps < 0:
            msg = ('CG steps are negative! (seriously? A negative number of steps '
                    'for the minimizer? Try again)')
            self.printmsg(msg)
        if self.sdsteps == self.cgsteps == 0:
            msg = ('(SD and CG steps have been set to 0: '
                    'use --nomini to disable minimizer)')
            print msg
        # forcefield
        self.dprint("[parseMinimizer] forcefield set to [%s]"% self.forcefield.upper())
        if not self.forcefield in [x.lower() for x in self.fflist.keys()]:
            msg = ('Specified forcefield [%s] is not recognized. '
                   'Supported forcefields are: %s.') % (self.forcefield,
                   self.fflist_str)
            self.printMsg(msg)

    def checkMiniParmsExtra(self):
        # disable minimization
        if self.noextra:
            anysteps = (self.sdsteps_extra > 0 or self.cgsteps_extra > 0)
            if anysteps:
                msg = ('Extra minimization disabled (--noextra) but minimizer '
                      'steps defined (--sdsteps_extra/--cgsteps_extra) requested')
                self.printMsg(msg, q=5)
                # printmsg will exit by now
            self.sdsteps_extra = 0
            self.cgsteps_extra = 0
            self.dprint('[parseMinimizer] extra minimization disabled (--noextra)')
            return
        else:
            if self.sdsteps_extra == None:
                self.sdsteps_extra = self.default_sdsteps_extra
            if self.sdconv_extra == None:
                self.sdconv_extra = self.default_sdconv_extra
            if self.cgsteps_extra == None:
                self.cgsteps_extra = self.default_cgsteps_extra
            if self.cgconv_extra == None:
                self.cgconv_extra = self.default_cgconv_extra
        # steepest descent
        if self.sdsteps_extra == 0:
            self.dprint("[parseMinimizer] extra minimization SD disabled (0 steps)")
        if self.sdsteps_extra < 0:
            msg = ('SD steps are negative! (seriously? A negative number of steps '
                    'for the minimizer? Try again)')
            self.printmsg(msg)
        # conjugate gradient
        if self.cgsteps_extra == 0:
            self.dprint("[parseMinimizer] extra minimization CG disabled (0 steps)")
        if self.cgsteps_extra < 0:
            msg = ('CG steps are negative! (seriously? A negative number of steps '
                    'for the minimizer? Try again)')
            self.printmsg(msg)
        if self.sdsteps_extra == self.cgsteps_extra == 0:
            msg = ('(SD and CG steps have been set to 0: '
                    'use --nomini to disable minimizer)')
            print msg
        # forcefield
        self.dprint("[parseMinimizer] forcefield set to [%s]"% self.forcefield_extra.upper())
        if not self.forcefield_extra in [x.lower() for x in self.fflist.keys()]:
            msg = ('Specified forcefield [%s] is not recognized. '
                   'Supported forcefields are: %s.') % (self.forcefield_extra,
                   self.fflist_str)
            self.printMsg(msg)
    
    def parseMolProcessing(self):
        """ parse molecule processing options"""
        # mol processing 

        self.pH = self.args.pH
        self.nopH = self.args.nopH
        self.flipamide = not self.args.noflipamide
        self.stripsalts = not self.args.nostripsalts
        self.checkhydro = not self.args.nocheckhydro

        if self.args.noprocess:
            self.nopH = self.noflipamide = True
            self.stripsalts = self.checkhydro = self.flipamide = False

        # check that --nopH and --pH are not used at the same time
        if (not self.nopH == None) and (not self.pH == None):
            msg = ('--pH value for protonation specified together with '
                   '--nopH option: only one can be used at the time')
            self.printMsg(msg, q=6)
        if (self.pH == None) and (self.nopH == None):
            self.pH = self.default_pH
            self.nopH = self.default_nopH
            self.dprint ("[parseMolProcessing) pH and protonation settings set to default")
        if self.nopH:
            self.pH = None
            self.dprint ("[parseMolProcessing) protonation disabled")

        if not self.args.enumchiral.lower() in ['all', 'undefined', 'protomer', 'off']:
            msg = ('incorrect value for --enumchiral option '
                    '[%s], allowed:' % self.args.enumchiral.lower())
                      #  "|".join(['all', 'undefined', 'None']))
            self.printMsg(msg)
        #else:
        #    self.chirality = self.args.enumchiral
        #print "SELFCHIRALITY", self.chirality

    def parseExclude(self):
        """ parse and check that SMARTS patterns are valid"""
        # parse options
        self.exclude = self.args.exclude
        self.excludefromfile = self.args.excludefromfile
        self.blacklist = []
        # exclude
        if len(self.exclude):
            msg = "[parseExclude] %d SMARTS pattern(s) specified with --exclude:\n" % len(self.exclude)
            patterns = ["   |%s|"%x for x in self.exclude]
            msg += "\n".join(patterns)
            self.dprint(msg)
            result = self.validateSMARTS(self.exclude)
            if not result == True:
                msg = ('Pattern defined in --exclude (%d) is not a valid '
                        'SMARTS pattern |%s|') % (result[0], result[1])
                self.printMsg(msg)
            self.blacklist += self.exclude
        # exclude from file
        if self.excludefromfile:
            msg = None
            patterns = self.getLinesFromFile(self.excludefromfile)
            if len(patterns) == 0:
                msg = '[%s] is an empty file'
            result = self.validateSMARTS(patterns)
            if not result == True:
                msg = ('Pattern at line %d in %s is not a valid SMARTS pattern '
                       '|%s|') % (result[0], self.excludefromfile, result[1])
            if not msg == None:
                self.printMsg(msg)
            self.blacklist += patterns
        self.initMolHound(self.blacklist)

    def parseMultiProc(self):
        self.multiproc_max = self.args.multicore
        # if multiproc_max > num_cpu:
        # raise warning
        self.nice = self.args.nice
        print "NICE NOT USED"


    def initOBChoices(self):
        """ initialize choices depending on OB features"""
        # forcefields list
        fflist =  ob.vectorString()
        ob.OBPlugin.ListAsVector('forcefields', None, fflist)
        self.fflist = {}
        self.fflist_str = []
        self.fflist_names = []
        for k, v in [ x.split(None, 1) for x in fflist ]:
            self.fflist[k] = v[:-1]
            self.fflist_str.append( "\"%s\" (%s)" % (k, v[:-1]))
            self.fflist_names.append('"%s"' % k.lower())
        self.fflist_str = ", ".join(self.fflist_str)
        self.fflist_names = ', '.join(self.fflist_names)

        # charges list
        chargelist = ob.vectorString()
        ob.OBPlugin.ListAsVector('charges', None, chargelist)
        self.chargelist = {}
        self.chargelist_str = []
        self.chargelist_names = []
        for k, v in [x.split(None, 1) for x in chargelist]:
            self.chargelist[k] = v #[:-1]
            self.chargelist_str.append( "\"%s\" (%s)" % (k, v)) # [:-1]))
            self.chargelist_names.append('"%s"' % k)
        self.chargelist_str = ", ".join(self.chargelist_str)
        self.chargelist_names = ', '.join(self.chargelist_names)

    def isMulti(self, fname, firstonly=True):
        """ check if there is more than one structure
            NOTE (2014.9.15): this function is never used 
            for actual counting, that could turn out to be
            very expensing when processing large data files

            The function should be re-written to check 
            that empty files are recognized!
        """
        string_patterns = { 'smi': '\n', 'can': '\n', 
                'sdf' : '$$$$',
              'mol2': '@<TRIPOS>MOLECULE',
              'pdb': 'MODEL',
              }

        binary_patterns = {'cdx' : None,
                        }
        #name, ext = self.getNameExt(fname)
        ext = self.informat
        if ext in string_patterns:
            patt = string_patterns.get(ext, None)
            if patt == None:
                print ('Warning! The [%s] format is a not '
                       'supported multi-structure file '
                       '(assuming multi-structure)' % ext)
                return True
            c = 0
            with open(fname, 'r') as f:
                for l in f:
                    if l.find(patt) >= 0:
                        c +=1
                        if c==2 and firstonly:
                            f.close()
                            return True
                # full count
                f.close()
                return c
            f.close()
        elif ext in binary_patterns:
            self.dprint("[isMulti] found CDX/binary")
            patt = binary_patterns.get(ext, None)
            # XXX this is an approximation
            # by definition, a CDX is considered a
            # multi-structure file
            return True
            
        return False

    def getNameExt(self, fname):
        """ extract name and extension from the input file"""
        #name, ext = fname.rsplit('.', 1)
        name, ext = os.path.splitext(fname)
        return name, ext[1:].lower()

    def getLinesFromFile(self, fname):
        """ parse stripped lines from file"""
        try:
            fp = open(fname, 'r')
            lines = fp.readlines()
            fp.close()
            lines = [ x.strip() for x in lines if x.strip() ]
            return lines
        except IOError as e:
            msg = ('Error while reading [%s]: I/O error (%d) : %s') % (fname, e.errno, e.strerror)
            self.printMsg(msg)


    def validateSMARTS(self, plist):
        """ check that every smarts pattern is """
        sys.stdout.flush()
        for i, p in enumerate(plist):
            self.dprint('[validateSMARTS]: checking pattern[%d]:"%s"' % (i,p))
            #print "[validateSMARTS]: checking %s\n" %p * 100
            validator = ob.OBSmartsPattern()
            validator.Init(p)
            if not validator.IsValid():
                self.dprint('[validateSMARTS]: invalid pattern found')
                
                return i,p
        self.dprint('[validateSMARTS]: all patterns are valid')
        return True

    def initMolHound(self, blacklist):
        """ initialize molHound SMARTS matcher"""
        self.molHound = MolecularHound(blacklist)

    def filterSMARTS(self, mol):
        """ basic SMARTS pattern matcher """
        if len(self.blacklist) == 0:
            return True
        result = self.molHound.sniffMolecule(mol)
        if result == False:
            self.dprint("[filterSMARTS] all filters passed")
            return True
        self.dprint("[filterSMARTS] [%s] rejected by SMARTS pattern #%d" % (self.currentname, result))
        self._rejected+=1
        return False

    def initDocs(self):
        """ initialize documentation"""
        # args/opt initialization
        self.desc = ('Minimize molecular structure(s) in input file.\n\n'
                'By default, steepest descent minimization is performed, \n'
                'followed by conjugate gradients minimization.\nWhen necessary, 3D structure'
                ' is generated automatically using MMFF94.\n'
                'All (reasonable) OpenBabel input/output file formats are supported.\n\n'
                '      _(\-/)_\n'
                '     {(#b^d#)} \n'
                '     `-.(Y).-` \n') # JGS

        #usage = 'This. Is. Sparta!'
        #self.usage = 'This. Is. Sparta!'
        self.usage = None
        self.epilog = None

    def showAdvancedHelp(self):
        """ do what it says..."""
        
        cc = ''
        for k in sorted(self.chargelist.keys()):
            v = self.chargelist[k]
            cc += '        {0:9s} : {1:30s}\n'.format(k,v)
            #cc += '        %s\t: %s\n' % (k, v)

        ff = ''
        for k in sorted(self.fflist.keys()):
            v = self.fflist[k]
            ff += '        {0:9s} : {1:30s}\n'.format(k,v)

        self.parser.print_help()
        self._advancedHelp = """
    
An example usage is the generation of ligand structures to be used for docking, 
with AutoDock, starting from SMILES files or ChemDraw schemes:

  %(progname)s --infile ligand_library.smi --outfile ligand_library.mol2
   
This will generate automatically 3D coordinates and save them in a multi-structure
Mol2 file.

Another usage could be to split a Mol2 file that already has 3D coordinates and 
generate single ligands suitable for docking with AutoDock (i.e., minimized 
conformation, protonation state, etc...):

  %(progname)s --infile nci_p0.0.mol2 --outfile dockings.pdbqt \
       --usemolname --end 10
   
This will convert the first 10 molecules found in the file Mol2, 
saving ZINC79106062.pdbqt, ZINC33512584.pdbqt, etc...

   

   %(progname)s is a general tool...
   PROTOCOL   When a molecule is read from a file, the following protocol is
   applied by default:

                                    [ START ]
                                        |
                                        |
                                  SMARTS filter
                                        |
                                        |
                               Strip salts/fragments
                                        |
                ________________________|_______________________
               |                        |                       |
               |                        |                       |
          1D formats               2D formats                3D formats
        (SMI, CDX,...)           (MOL2, SDF...)           (MOL2, SDF, PDB...)
               |                        |                       |
               |                        |                       |
         canonical SMI                  |                       |
               |                        |                       |
               |                        |                       |
           3D coordinates generation (mmff94)             check hydrogens
                         |                                      |
                         |______________________________________|
                                        |
                                        |
                      Check primary amide trans-conformation
                                        |
                                        |
                                 protonate for pH
                               (default: %(pH)1.1f)
                                        |
                                        |
                            Calculate partial charges
                         (default charge model: %(charge)s)
                                        |
                                        |
                           Steepest descent minimization
                          (default forcefield: %(forcefield)s)
                                        |
                                        |
                          Conjugate gradient minimization
                          (default forcefield: %(forcefield)s)
                                        |
                       Quick weighted conformational search 
                                        |
                          Post-optimization using UFF ???
                                        |
                                        |
                                 Save output file
                              (default format: %(outformat)s)
                                        |
                                        |
                                     [ END ]            


   In the default non-verbose mode, processing output is shown is presented as following:

                      .- Steepest-descent minimization
                      |
           molname:SD[x],CG[x]
                            |
                            '- Conjugate gradient minimization

   where 'x' can be 0 ( = converged) or  1 ( = max.iterations reached)
   Ideally, all minimizations should converge to guarantee that a true minimum
   conformation is reached. This is particularly important for docking.

   3D coordinates are generated automatically when 1D or 2D structures are found, then
   a series of structural clean-ups are performed.   

   SUPPORTED FORCEFIELDS
   The following forcefields are supported:
%(ff)s   
    SUPPORTED PARTIAL CHARGE MODELS
    The following partial charge models are available:
%(cc)s

EXAMPLES

   %(scriptname)s --infile ligand_library.smi --outfile ligand_library.mol2
       #  Generate 3D structures for all molecules found in the input file, adding hydrogens
       #  protonating groups at pH 7.4, minimizing coordinates and performing all default 
       #  structural cleanups, and saving the result in a single multi-structure file
        
   %(scriptname)s --infile nci_p0.0.mol2 --outfile docking.pdbqt --noprocess --nomini --end 100 --usemolname
       #  Convert first 100 molecules in the input file to separate PDBQT files ready 
       #  to be docked and using molecule names contained in the input file; no minimization is
       #  performed; all processing steps (salt stripping, hydrogens/protonation, trans-amide check)
       #  are disabled; ligands will be saved as: 
       #       docking_ZINC04783481.pdbqt
       #       docking_ZINC04786808.pdbqt
       #       docking_ZINC04786811.pdbqt
       #       ...
              
====================================================	 
%(progname)s (C)2014 Stefano Forli, MGL, TSRI 
        """ % { 'progname' :self.progname,
        'pH': self.default_pH, 'charge':self.default_charges, 'forcefield':self.default_forcefield,
        'outformat': self.default_out, 'cc': cc, 'ff':ff, 'scriptname': os.path.basename(sys.argv[0]),
        }  
        print self._advancedHelp


    def initOptions(self):
        """ initialize command line options and information to properly
            pack the parser respecting proper order
        """
        print "\n\n ADD SHORT AND LONG VERSION! for INITIOPTIOHNS!\n\n"

        self._argsDict = { '--infile' : { 'help':'input file; '
                                  'format is guessed from file extension; ',
                                 #'this is the only required argument', 
                                 'action': 'store',
                                 'metavar' : 'INPUT_FILE[.EXT]', 'required' : True, 'type': str},

        '--single': {  'help':'process only Nth structure in the file' ,
                        'action' :  'store', 'metavar': 'STRUCTURE_NUMBER', 'type' : int},

        '--begin' : { 'help':'start processing from Nth molecule', # (default: first molecule)', 
                       'metavar' : 'FIRST_STRUCTURE_NUMBER', 'type': int, 'action':'store', 
                       'default': None},

        '--end' : { 'help': 'stop processing at Nth molecule ', # (default: last molecule)', 
                    'action' : 'store', 'metavar' : 'LAST_STRUCTURE_NUMBER', 'type': int, 
                    'default': None},

        '--byname' : { 'help' : 'process only molecules with specified name (note: cAse sEnSiTivE,'
                    ' spaces are not allowed; ...; )',
                    'metavar' : 'MOL_NAME', 'action' : 'store', 'type' : str },

        '--outfile' :  { 'help':('set output file; format is guessed from file extension; '
                                   'by default, output filename is \'INFILE_mini\' and same '
                                   'type as input; if input file format is not a 3D format '
                                   '(e.g. SMI) default %s is used; if an existing  directory '
                                   'is specified, output will be written in there') % self.default_out.upper(),
                        'metavar': 'OUTPUT_FILE[.EXT] ', 'action': 'store', 'type' :str, 'default':''},

        '--informat' : {'help':('force input file to be parsed as specified'),
                             'action':'store', 'metavar':'[mol2|sdf|smi|pdb|...]', 'type':str},

        '--outformat' : {'help':('force output file format to be written as specified'),
                             'action':'store', 'metavar' : '[mol2|sdf|pdb|...]', 'type':str, 'default':''},

        # Molecule naming schemes
        '--usemolname':{'help':('enable \'--split\' and use molecule names for saving splitted output files; '
                           'when names are not found or duplicated, progressive numbering is used; if --outfile\n'
                           'is used, names will have the outname string as prefix'
                           ),
                           'action':'store_true', 'default': self.default_usemolname},

        '--usemolnamesafe':{'help':('same as \'--usemolname\', but sanitize names to be valid/clean '
                           'filenames'),
                           'action':'store_true', 'default': self.default_usemolnamesafe},

        '--usefieldname' : {'help' : ('use specified field as molecule name (only for SDF, Mol2 '
                             'and PDB input files)')
                             #; NOTE: spaces must specified with HTML notation, '
                             #'i.e. \'Cpd ID\' -> \'Cpd%%20ID\' (thank you sys.argv!)')
                             , 'metavar': 'FIELD', 'action':'store', 
                             'default' : None},

        # Splitting
        '--split' : { 'help':('split each molecule as separate file, progressively numbered '
                                '(see \'--usemolname\' for using molecule names); '
                                'disabled by default if output format supports multi-structure'),
                                'default': self.default_split, 'action': 'store_true'},
        
        # minimizing 
        # XXX TODO XXX TODO XXX TODO
        # provide a minimum set of values good for everything
        '--mini_quick' : { 'help':('steepest descent max. iterations; set to 0 to disable '
                        '(default %d).') % self.default_sdsteps, 'action': 'store', 'metavar':'[ %4g ]' % self.default_sdsteps,
                        'type' : int, 'default' : None},
        '--mini_std' : {},
        '--mini_aggressive' : {},
        # XXX TODO XXX TODO XXX TODO




        '--sdsteps' : { 'help':('steepest descent max. iterations; set to 0 to disable '
                        '(default %d).') % self.default_sdsteps, 'action': 'store', 'metavar':'[ %4g ]' % self.default_sdsteps,
                        'type' : int, 'default' : None},

        '--sdconv' : {  'help':('set steepest descent convergence criterion. (default: %s)') % self.default_sdconv, 
                        'action': 'store', 'metavar':  '[ %4g ]' % self.default_sdconv, 'type': float,
                        'default': self.default_sdconv},

        '--cgsteps': {  'help': ('conjugate gradients max. iterations; set to 0 to disable '
                         '(default %d).' % self.default_cgsteps), 'action': 'store', 'metavar':'[ %3g ]' % self.default_cgsteps,
                         'type' :  int, 'default': None},

        '--cgconv' : { 'help':('set conjugate gradients descent convergece criterion. '
                             '(default: %s)' % self.default_cgconv), 'action': 'store',
                             'metavar': '[ %4g ]' % self.default_cgconv, 'type':float, 'default': None},

        '--sdsteps_extra' : { 'help':('steepest descent max. iterations (post-processing minimization with UFF?); set to 0 to disable '
                        '(default %d).') % self.default_sdsteps_extra, 'action': 'store', 'metavar':'[ %4g ]' % self.default_sdsteps_extra,
                        'type' : int, 'default' : None},

        '--sdconv_extra' : {  'help':('set steepest descent convergence criterion (post-processing minimization with UFF?). (default: %s)') % self.default_sdconv_extra, 
                        'action': 'store', 'metavar':  '[ %4g ]' % self.default_sdconv_extra, 'type': float,
                        'default': self.default_sdconv_extra},

        '--cgsteps_extra': {  'help': ('conjugate gradients max. iterations (post-processing minimization); set to 0 to disable '
                         '(default %d).' % self.default_cgsteps_extra), 'action': 'store', 'metavar':'[ %4g ]' % self.default_cgsteps_extra,
                         'type' :  int, 'default': None},

        '--cgconv_extra' : { 'help':('set conjugate gradients descent convergece criterion (post-processing minimization). '
                             '(default: %s)' % self.default_cgconv_extra), 'action': 'store',
                             'metavar': '[ %4g ]' % self.default_cgconv_extra, 'type':float, 'default': None},

        '--rotamer_conf' : { 'help':('number of random conformers to use in weighted conformational search '
                             '(default: %s)' % self.default_rotamer_conf), 'action': 'store',
                             'metavar': '[ %d ]' % self.default_rotamer_conf, 'type':int, 'default': None},
        
        '--rotamer_steps' : { 'help':('number of geometry minimization steps to per form during weighted '
                              'conformational search; this value should be small, otherwise calculation times '
                              'will increase considerably) (default: %s)' % self.default_rotamer_steps), 'action': 'store',
                             'metavar': '[ %s ]' % self.default_rotamer_steps, 'type':int, 'default': None},

        '--forcefield' : {'help': ('set forcefield for minimization (default: \'%s\'; see --help_advanced for '
                            'choices)') % (self.default_forcefield) , 
                         'default': self.default_forcefield, 'metavar':'[ %s ]' % self.default_forcefield,
                         'type': str}, 

        '--forcefield_extra' : {'help': ('set forcefield for post-minimization (default: \'%s\'; see --help_advanced for '
                            'choices)') % (self.default_forcefield_extra) , 
                         'default': self.default_forcefield_extra, 'metavar':'[ %s ]' % self.default_forcefield_extra,
                         'type': str}, 

        '--strict'     : {'help':'do not process molecules with atoms missing forcefield/charges '
                          'parameters (default: %s)' % self.default_strict, 'action':'store_true' },

        '--nomini' :{ 'help': 'disable minimization',
                    'default':self.default_nomini, 'action': 'store_true'},

        '--noextra' :{ 'help': 'disable extra minimization',
                    'default':self.default_noextra, 'action': 'store_true'},

        '--norotamer' :{ 'help': 'disable rotameric minimum search',
                    'default': self.default_norotamer, 'action': 'store_true'},

        '--log' :  { 'help': 'XXASADAD Asave minimization logs in specified file; by default the output '
                     'filename is used  (OUTFILE.log)', 'action': 'store', 'metavar':'LOGFILENAME'},

        '--verbose' : {'help': ('enable verbose mode'),
                                'default' : self.default_verbose, 'action': 'store_true'},

        '--noflipamide' :  { 'help': 'disable forcing cis conformation for non-cyclic primary/secondary amides', 
                             'action': 'store_true', 'default': self.default_noflipamide},

        '--nostripsalts' :  { 'help': 'disable removing salts/ions and small fragments (i.e. '
                            'keep largest fragment only)', 'action': 'store_true', 
                             'default': self.default_nostripsalts},

        '--nocheckhydro' :  { 'help': 'disable check for missing hydrogens', 'action': 'store_true', 
                            'default': self.default_nocheckhydro },

        '--noprocess' :  { 'help': 'disable all structure processing (equivalent to (--noflipamide) + '
                            '(--nocheckhydro) + (--nostripsalt) + (--nopH)', 'action': 'store_true', 
                            'default': self.default_noprocess },

        '--pH' :  { 'help': 'generate protonation state for requested pH (default: %1.1f)' % self.default_pH, 
                    'action': 'store', 'type': float, 'default': None, 'metavar':'[ 7.4 ]'}, # default set to None to perform check later

        '--nopH' :  { 'help': 'disable protonation state modifications', 
                    'action': 'store_true', 'default': None}, # default set to None to perform check later

        '--chargemodel' : { 'help': ('set charge model to compute atom partial '
                            'charges (default \'%s\'; see --help_advanced for options)')
                            % (self.default_charges), 'metavar': '[ %s ]' % self.default_charges,
                            'action': 'store', 'type': str, 'default': self.default_charges}, 
        '--enumchiral' : { 'help' : ('generate multiple stereoisomers for tetrahedral chiral'
                                'centers (note: input geometry will be ignored); \'undefined\': '
                                'enumerate only undefined chiral centers; '
                                '\'all\': enumerate all chiral centers; default: \'off\' (no enumeration)'),
                                'default': self.default_enumchiral, 'action' : 'store',
                                'metavar' : '[undefined|protomer|all|off]', 'type' : str, },

        '--maxenumchiral' : { 'help': ('maximum number of enantiomers to generate if --enumchiral'
                        ' is active; default: %d' % self.default_maxenumchiral),
                'default': self.default_maxenumchiral, 'action' : 'store', 
                'metavar': 'MAX_ENANTIOMERS', 'type':int},

        '--exclude' :  { 'help': ('skip molecules matching specified SMARTS pattern; '
                                 'multiple patterns can be repeated by using multiple \'--exclude\' '
                                'or using \'--excludefromfile\''), 
                    'action': 'append', 'metavar': 'SMARTS', 'type': str, 'default': []},

        '--excludefromfile' :  { 'help': ('skip molecules matching any of the SMARTS patterns'
                                        ' in each line of FILENAME'), 
                    'action': 'store', 'metavar': 'FILENAME', 'type': str, 'default': None},

        '--multicore' :  { 'help':('specify how many cores/cpu to use '
                             '(default: %d)' % self.default_multiproc_max), 'action': 'store',
                             'metavar': '[ %d ]' % self.default_multiproc_max, 'type':int, 
                            'default': self.default_multiproc_max},


        '--nice' : {'help': 'set nice level (not working: requires psutil)'},

        '--help_advanced' : {'help': 'show the advanced help with methods, etc...', 'action': 'store_true',
                    'default': False},

                    }

        # enforce the proper arguments order for the help message
        self._argsorder = {  # in/out
                'INPUT/OUTPUT' : [ '--infile', '--outfile',  '--informat', 
                        '--outformat', '--usemolname', '--usemolnamesafe', '--usefieldname' ],
                # slicing
                'SLICING' : ['--single', '--begin', '--end', '--split', '--byname'],
                # mol processing 
                'CLEANUP': [ '--pH', '--nopH', '--noflipamide', '--nostripsalts', '--nocheckhydro', 
                            '--noprocess', '--enumchiral', '--maxenumchiral', '--exclude', '--excludefromfile',],
                # minimizer 
                'ENERGY MINIMIZATION AND FORCEFIELD PARAMETERS' : [ '--sdsteps', '--sdconv', '--cgsteps', '--cgconv', '--forcefield', 
                             '--sdsteps_extra', '--sdconv_extra', '--cgsteps_extra', '--cgconv_extra', '--forcefield_extra', 
                             '--rotamer_conf', '--rotamer_steps',

                '--nomini', '--noextra', '--norotamer', '--chargemodel', '--strict'],

                # multiprocessing
                'MULTIPROCESSING' : ['--multicore', '--nice'],

                # logging and info
                'LOGGING': ['--log', '--verbose', '--help_advanced'],
                     }

        self.groupsOrder = [
            ('INPUT/OUTPUT', ('Control input/output data files and types.')),
                             #'%s supports (in theory)  most  of the OpenBabel  file formats.\n'
                             #'3D structures are generated automatically for molecules and files\n'
                             #'that do not have 3D information (i.e., SMI, CDX).')),

            ('SLICING', ('Multi-structure files management')),
                        #'For files containing multiple molecules, it is possible to specify\n'
                        #'the range of molecules to be processed, split a multi-structure file in\n'
                        #'separate files.' )),
            ('CLEANUP', 'Pre- and post-processing cleanups performed on molecules'),

            ('ENERGY MINIMIZATION AND FORCEFIELD PARAMETERS',
                        'Options for the structure optimization'),

            ('MULTIPROCESSING', 'Manage multi-processing and parallel computation, nice, etc...'),

            ('LOGGING', 'Control logging, message reports and extra info'),
                        ]


    def initMolParser(self):
        """ initialize the molecule parser"""
        self.molParser = ob.OBConversion()
        """
        self.molParser.SetInAndOutFormats(self.intype,self.outtype)
        self.dprint("[initMolParser] initialized molParser")
        if self.split or self.single:
            # files will be read and written separately
            return
        self.dprint("[initMolParser] Open combined in/output parser")
        self.molParser.OpenInAndOutFiles(self.infile, self.outfile)
        """
        self.molParser.SetInAndOutFormats(self.intype,self.intype) #outtype)
        self.dprint("[initMolParser] initialized molParser: %s,%s"% (self.intype, self.intype))
        #if self.split or self.single:
        #    # files will be read and written separately
        #    return
        #self.dprint("[initMolParser] Open combined in/output parser")
        #self.molParser.OpenInAndOutFiles(self.infile, self.outfile)

    def readMolecule(self, first=False):
        """ initialze the molecule parser and get the first molecule"""
        """ read the next molecule from the molecular parser
            or report when the end is reached
        """
        mol = ob.OBMol()
        # XXX This part is not very clear...
        if first:
            #if (not self.split) and (not self.single):
            #    more = self.molParser.Read(mol)
            #else:
            #    more = self.molParser.ReadFile(mol, self.infile)
            if self.split or (not self.single == None):
                more = self.molParser.ReadFile(mol, self.infile)
            else:
                more = self.molParser.ReadFile(mol, self.infile)
        else:
            more = self.molParser.Read(mol)
        if not more:
            return False
        self._counter += 1
        return mol

    def getMolString(self, mol):
        """ generate string to submit to the queue"""
        return self.molParser.WriteString(mol) #, self.intype)
    
    def getMolTitle(self, mol):
        """ get the molecule name stored in the OBMol obj
            and sanitize it if required: spaces will be converted to
            underscores, then unacceptable characters will be
            discarded
        """
        title = mol.GetTitle()
        if not self.safenames:
            return title
        title = title.replace(' ', '_')
        validname = ''.join(c for c in title if c in self._validFnameChars)
        return validname

    def getMolName(self, mol, data={}): #, field=None):
        """ find the appropriate unique name for the molecule
            accordingly to the naming policy selected;
            if the molecule lacks a proper name or the
            progressive numbering is asked, a new name
            is genearated and accounting performed

            perform also the filter by name, returning None
            if no match with byname value is found
        """
        # get the molecule name
        mname = ''
        #print "FFFFIEL", self.args.usefieldname
        if not self.args.usefieldname == None:
            self.args.usefieldname = self.args.usefieldname.replace('%20', ' ')
            self.dprint('[getMolName] Search field name :"%s"' % self.args.usefieldname)
            # by field 
            try:
                mname = data[self.args.usefieldname].strip()
                self.dprint('[getMolName] Found name :"%s"' % mname)
                mol.SetTitle(mname)
            except KeyError:
                mname = ''
                self.dprint('[getMolName] WARNING: field name not found')
        elif self.usemolname:
            # by mol name
            mname = self.getMolTitle(mol)
            if mname == '':
                mname = 'MOL'
        else:
            # by default
            mname = "MOL_%d" % self._counter
        if mname.strip() == '':
            mname = 'MOL'
        # filter and accept only molecule by name
        if self.args.byname:
            if not self.getMolTitle(mol).strip() == self.args.byname:
                msg = ("[getMolName] Skipping molecule by name filter (%s) [%s]")
                self.dprint(msg % (self.getMolTitle(mol), self.args.byname) )
                return None
        if not mname in self.molNames.keys():
            self.molNames[mname] = -1
        self.molNames[mname] += 1
        if self.molNames[mname] == 0:
            name = mname
        else:
            name = "%s_%d" % (mname, self.molNames[mname])
        self.dprint("[getMolName] Name is now [%s]" % name)
        return name

    def getMolData(self, mol):
        """ parse generic data (e.g., SDF, Mol2 fields...)
            [ this code is shamelessly borrowed from pybel ]
        """
        data = mol.GetData()
        answer = [ x for x in data if
            x.GetDataType() == ob.PairData or
            x.GetDataType() == ob.CommentData]    
        raw = [ ( x.GetAttribute(), x.GetValue() ) for x in answer ]
        return dict(raw)        
        
    def writeOutputFile(self, mol):
        """ write processed molecule accordingly to the file
            output policy: update the multi-structure file
            or create a new file 
        """
        if self.split or self.single:
            fname = self.outfile % (self.currentname)
            self.dprint( "[writeOutputFile] Saving new file: %s"% fname)
            self.molParser.WriteFile(mol, fname)
        else:
            self.dprint( "[writeOutputFile] Updating file: %s"% self.outfile)
            self.molParser.Write(mol)

    def closeMolParser(self):
        """ close output file used for multi-structure savings"""
        if not self.split:
            self.molParser.CloseOutFile()

    def guessParameters(self, mol):
        """ try to guess """
        rot = mol.NumRotors()
        atoms = mol.NumHvyAtoms()
        bonds = mol.NumBonds()
        print "[[[[[[[[[[ %s R:%d  A:%d  B:%d"% (mol.GetTitle(), rot, atoms, bonds),
        sdsteps = max( ( (bonds/10.) * (atoms/10.) * rot), 10)
        cgsteps = sdsteps / 2
        sdsteps_extra = sdsteps * 0.75
        cgsteps_extra = cgsteps * 0.75
        print "==>", sdsteps, cgsteps, sdsteps_extra, cgsteps_extra

        self.scrubOptions['sdsteps'] = int(sdsteps)
        self.scrubOptions['cgsteps'] = int(cgsteps)
        self.scrubOptions['sdsteps_extra'] = int(sdsteps_extra)
        self.scrubOptions['cgsteps_extra'] = int(cgsteps_extra)
        self.scrubOptions['rotamerConf'] = int(rot)


    def initThreads(self):
        """ """
        self.scrubOptions = { 'ff':self.forcefield, 'sdsteps':self.sdsteps,
            'sdconv':self.sdconv, 'cgsteps':self.cgsteps, 
            'cgconv':self.cgconv, 'ff_extra':self.forcefield_extra, 
                'sdsteps_extra':self.sdsteps_extra, 
                'sdconv_extra':self.sdconv_extra, 
                'cgsteps_extra':self.cgsteps_extra,
                'cgconv_extra':self.cgconv_extra, 
                'rotamerConf' : self.rotamer_conf,
                'rotamerGeomSteps' : self.rotamer_steps,
                'pH':self.pH, 'chargemodel':self.chargemodel,
                'flipamide': self.flipamide,
                'stripsalts':self.stripsalts, 
                'checkhydro':self.checkhydro,   'name': None, #self.currentname, 
                'verbose':self.verbose, 'auto':True
            }
        # feeding queue
        self.queueIn = mp.JoinableQueue(maxsize=self.multiproc_max)
        # writer queue
        self.queueOut = mp.Queue(maxsize=-1)
        # create processes and start them
        for i in xrange(self.multiproc_max):
           s = ScrubMT(self.queueIn, self.queueOut, nice=None, 
                inType=self.intype, outType=self.outtype,  scrubOpts = self.scrubOptions)
           s.start()
        # XXX this has to be fixed for split, etc...
        if self.split or self.single:
            fname = None
        else:
            fname = self.outfile
        writer = MolWriter(queue=self.queueOut, threads=self.multiproc_max,
                    fname=fname, ext=self.outtype)
        writer.start()
        print "%d threads initialized" % self.multiproc_max


    def closeThreads(self):
        """ """
        for i in xrange(self.multiproc_max):
            self.queueIn.put((None,None))
        #self.queueOut.put((None, None))
        
    def initializeLoop(self):
        """ initialize variables, molecule input file parser and threads"""
        self._counter = 1 # general counter
        self._rejected = 0
        self._processed = 0
        self.initMolParser()
        self.initThreads()
        
    
    def start(self):
        """ main loop does the actual job"""
        self.dprint('====================[ START ]====================')
        # start!
        self.initializeLoop()
        rawMol = self.readMolecule(first=True)
        if not self.begin == 1:
            print "Seeking begin molecule %d... " % (self.begin)
        while rawMol:
            # getting the data is more efficient (before enumerating chirals)
            # check that there are atoms (useful for CDX files)
            if not rawMol.NumAtoms():
                self.dprint("No atoms, skipping molecule")
                rawMol = self.readMolecule()
                continue
            # skip by counter
            if self._counter < self.begin:
                rawMol = self.readMolecule()
                self.dprint("[start] molecule[%s] skipped by counting" % self._counter)
                continue
            molData = self.getMolData(rawMol) 
            #print "========================== MOLE READ", rawMol.NumAtoms()
            # enumerate chiral structures
            for mol in ChiralEnumerator(rawMol, self.args.enumchiral,
                                self.args.maxenumchiral, debug=self.verbose):
                # TODO consider if grouping these functions in an "EvaluateMol" wrapper?
                # naming
                self.currentname = self.getMolName(mol, molData)
                if self.currentname == None:
                    # skip by name
                    continue
                msg = "-------------[processing mol.%d: %s]---------" 
                if 0:
                    print "\n\nGUESSING PARMS!!!!\n\n"
                    self.guessParameters(mol)
                    
                args = (self._counter, self.currentname)
                self.dprint(msg % args)
                # skip by SMARTS
                if not self.filterSMARTS(mol):
                    molRaw = self.readMolecule()
                    self.dprint("[start] molecule[%s] skipped by SMARTS" % self.currentname)
                    continue
                self.addMolToQueue(mol, self.currentname)
                self._processed += 1
            # terminate by counter
            if self.end and (self._counter == self.end):
                break
            rawMol = self.readMolecule()
        print "\n\n\n\n\n[ DONE ]Total structures processed: %d" % (self._processed)
        # poison pill
        self.closeThreads()

    def addMolToQueue(self, mol, name):
        """ """
        string = self.getMolString(mol)
        self.queueIn.put((string,name), block=True)

    def setStopCriterion(self):
        for i in xrange(self.multiproc_max): 
            self.queueIn.put((None,None))
        self.queueOut.put(None)
        #self.queueIn.close()
        #self.queueIn.join_thread()

    def debugPrint(self, mol):
        """ """
        c = 0
        for a in ob.OBMolAtomIter(mol):
            c += a.GetPartialCharge()
        print "DEBUGGO: charge    ", mol.GetTotalCharge()
        print "DEBUGGO: charge (P)", c
            
         
###### XXX XXX HIC SUNT LEONES
if __name__ == '__main__':
    MiniMee()
    """
    print "- COUNT THE MOLECULES IN THE FILE? TO HAVE A PROPER ZERO-PADDING?"
    print "- ADD EXAMPLES"
    print "- ADD INCLUDE (opposite than EXCLUDE"
    print "PH RANGES..."
    print "TAUTOMERS"
    print "SULFONAMIDES FLIP?"
    print "NONPRIMARY AMIDE FLIPS?"
    print "- ADD AUTOFILTER FOR REACTIVE GROUPS"
    print "preserve properties  (GETDATA)"
    print "RENAME atoms to have unique names"
    print "\n\n\n"


    print "GET IDEAS:http://openbabel.org/docs/dev/Command-line_tools/babel.html#append-option " 
    print "\n Implement UNIQUE ATOM NAMES"
    #print "-- enumerate chirality"
    """
