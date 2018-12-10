#!/usr/bin/env python

import pybel
ob = pybel.ob
# import shlex
import numpy as np
import math, string
import sys
import os

#from CopyMol import MoleculeDuplicator
import itertools

OPTIMALAMIDEDIHE = math.pi
ZEROISH = 0.005

###
### IMPORTANT: get example to suppress log from Dalke
#              self.lvl = ob.obErrorLog.GetOutputLevel()
#              ob.obErrorLog.SetOutputLevel(-1)

###

from debugtools import DebugObj
# XXX refactor this class to have multiple set methods
# XXX so objects are initialized once, performance will be better
class Scrub(DebugObj):
    """ This class is devoted to generate a minimum conformation
            for a ligand.

        [ http://www.merriam-webster.com/dictionary/scrub ]
        2scrub verb
                : to rub (something) hard with a rough object or 
                  substance and often with soap in order to clean it
    """
    # TODO
    #   - add tautomeric search
    #   - add quick conformational search
    #   - include (as opposed to exclude)
    #   - add PAINS support
    #   - read from compressed formats
    def __init__(self, ff='MMFF94s', ff_extra='UFF',
        sdsteps=250, sdconv= '1e-6', 
                cgsteps=100, cgconv='1e-8', 
                sdsteps_extra=100, sdconv_extra= '1e-6', 
                cgsteps_extra=100, cgconv_extra='1e-8', 
                rotamerConf = 10, rotamerGeomSteps=5,
                extra3Dmini = False,
                flipamide = True,
                checkhydro = True,
                pH = 7.4, stripsalts=True,
                chargemodel='gasteiger',
                name = None,
                verbose=False,
                strict=False,
                inFormat=None,
                outFormat='mol2',
                debug=False,
                auto=False):
        DebugObj.__init__(self, debug)
        self.verbose = verbose
        self.stripsalts = stripsalts
        self.ready = True
        self.log = []
        self.mol = None
        #self.name = mol.GetTitle()
        #self.vprint('(__init__): name found [%s]' % self.name)
        #if name == None:
        #    name = "MultiProcessMolecule"
        self.molName = name
        #if not name == None:
        #    self.name = name
        #    if self.verbose:
        #        self.vprint( "(__init__): molecule name set to user defined [%s]" % self.name)
        self.ff = ff.lower()
        self.ff_extra = ff_extra.lower()
        if not pH == None:
            self.pH = self._double(pH)
        else:
            self.pH = None
        self.sdsteps = sdsteps
        self.sdconv = sdconv
        self.cgsteps = cgsteps
        self.cgconv= cgconv
        self.sdsteps_extra = sdsteps_extra
        self.sdconv_extra = sdconv_extra
        self.cgsteps_extra = cgsteps_extra
        self.cgconv_extra = cgconv_extra
        self.rotamerConf = rotamerConf
        self.rotamerGeomSteps = rotamerGeomSteps
        self.chargemodel = chargemodel
        self.strict = strict
        self.inFormat = inFormat
        self.outFormat = outFormat
        self.flipAmide = flipamide
        self.checkHydro = checkhydro
        self.initEssentials() # XXX <- this should help when refactoring into reusable class
        #if auto:
        #    self.process()

    def initEssentials(self):
        """ initialize OB object used during the processing"""
        # 3D structure builder
        self.builder = ob.OBBuilder()
        # molecule writer to create string
        self.parser = ob.OBConversion()
        self.parser.SetOutFormat(self.outFormat)
        if not self.inFormat == None:
            self.parser.SetInFormat(self.inFormat)

    def vprint(self, string, addnewline=True, flush=True, rewind=False):
        """ debugging printing"""
        if not self.verbose:
            return
        msg = ''
        if rewind:
            msg += '\r'
        msg += "VERBOSE: [SCRUB] %s" % string
        if not addnewline:
            print msg,
        else:    
            print msg
        if flush:
            sys.stdout.flush()

    def _double(self, value):
        """ helper function to generate double floats"""
        return ob.double_array([value])[0]

    def process(self, molString, inType, outType, uniqueAtomNames=True):
        """ perform the operations requested"""
        if not inType == None:
            self.parser.SetInAndOutFormats(inType,outType)
        self.mol = ob.OBMol()
        self.parser.ReadString(self.mol, molString)
        if self.molName == None:
            self.molName = self.mol.GetTitle()

        self.initForcefield()
        self.prepareMolecule()
        #if not self.ready:
        #    print "*** Processing aborted, mol not ready ***"
        #    return
        # perform minimization
        self.optimize()
        # check if amides are in cis-form
        if self.flipAmide:
            self.checkAmide()
        # perform UFF post-processing
        self.postOptimize()
        self.ready = True
        if uniqueAtomNames:
            #print "X"
            self.setUniqueAtomNames()
        out = self.parser.WriteString(self.mol)
        return ("%s\n" % out)
        

    def prepareMolecule(self):
        """ initialize presettings"""
        #if not self.ready:
        #    return
        # strip fragments, salts...
        if self.stripsalts:
            self.mol.StripSalts()
        # check if structure needs 3D coordinates
        self.checkStructure()
        if not self.ready:
            return
        # check hydrogens
        self.fixPH()
        # calculate charges
        self.calculateCharges()
        

    def checkStructure(self):
        """ check if the molecule has a 3D structure
            otherwise generate it
            if the structure is already 3D, 
            check that hydrogens are properly added
        """
        # 3D + hydrogens
        dimension = self.mol.GetDimension()
        if dimension < 3:
            canonicize = False
            self.vprint("[checkStructure] Generating 3D coords...")
            if dimension < 2:
                canonicize=True
            self.make3D(canonicize)
            self.vprint("[checkStructure] [ DONE ]")
        # hydrogens only
        else:
            if self.checkHydro:
                self.checkMissingHydrogens()

    def checkMissingHydrogens(self):
        """ check that implicit and explicit valences
            match, otherwise add appropriate hydrogens
        """
        #self.mol.PerceiveBondOrders()
        #print "SKIPP CHECKISS"
        return 
        for a in ob.OBMolAtomIter(self.mol):
            valence = a.GetValence()
            implicit = a.GetImplicitValence()
            if not (valence == implicit):
                self.mol.AddHydrogens(a)

    def fixPH(self, pH=None):
        """ performs pH magics"""
        if pH == None:
            pH = self.pH
        if self.pH == None:
            return
        #self.verbose=True
        self.vprintCan('fixPH: START')
        self.mol.DeleteHydrogens()
        self.mol.UnsetFlag(ob.OB_PH_CORRECTED_MOL)
        self.vprintCan('fixPH: UNSETFLAG   ')
        for a in ob.OBMolAtomIter(self.mol):
            a.SetFormalCharge(0)
        self.vprintCan('fixPH: FORMALCHARGE')
        self.mol.SetAutomaticFormalCharge(True)
        self.vprintCan('fixPH: AUTOFORMAL  ')
        self.mol.AddHydrogens(False, True, pH)
        self.vprintCan('fixPH: ADDYDRO  ')
        #self.verbose=False


    def calculateCharges(self):
        """ calculate partial charges using selected charge model"""
        # XXX partial charges are queried to trigger something in the 
        #     openbabel lib that fixes charge bug 
        # TODO this function should be disabled if an option to use original
        #      charges is used
        print "WARNING: destroying any input charges!"
        self.vprint("Charge assignment; model [ %s ]" % self.chargemodel )
        charger = ob.OBChargeModel.FindType(self.chargemodel)
        self.mol.UnsetPartialChargesPerceived()
        self.mol.SetAutomaticPartialCharge(False)
        report = charger.ComputeCharges(self.mol)
        c = sum([a.GetPartialCharge() for a in ob.OBMolAtomIter(self.mol)])
        self.vprint("Total partial charges sum: %2.3f"  % c)
        if not report:
            msg = ('[calculateCharges] WARNING: charge model is '
                     'missing parameters for some atoms')
            if self.verbose:
                print msg
            if self.strict:
                msg = ('MISSING CHARGES (--strict mode on) charge model '
                       'not available for molecule -> rejecting')
                print msg
                self.ready = False
                return
        partialCharges = charger.GetPartialCharges()
        if self.verbose:
            #self.vprint("[calculateCharges] Charges assigned:")
            totalCharge = sum(partialCharges)
            buff = '[calculateCharges] Charges assigned, total charge: %2.3f\n' % totalCharge
            indices = range(len(partialCharges))
            elements = [ self.mol.GetAtom(i+1).GetAtomicNum() for i in indices]
            charges = ['%3.2f' % x for x in partialCharges]

            idxString = ",".join([str(x) for x in indices])
            elementString = ",".join([str(x) for x in elements])
            chargeString = ",".join(charges)

            buff += "  AtomIdx\t:" + idxString + "\n"
            buff += "  Element\t:" + elementString + "\n"
            buff += "  Charge\t:" + chargeString + "\n"
            self.vprint(buff)

    def checkAmide(self):
        """ check if amide is in cis conformation, and flips it if necessary"""
        for b in ob.OBMolBondIter(self.mol):
            if b.IsPrimaryAmide() or b.IsSecondaryAmide():
                if not b.IsInRing():
                    self.forceTransAmide(b)

    def forceTransAmide(self, bond):
        """ measure amide torsion angle and force it to be 180"""
        C, N, O, H = 6,7,8,1

        optimal = OPTIMALAMIDEDIHE
        tol = math.radians(36)
        r180 = math.radians(180)
        begin = bond.GetBeginAtom()
        end = bond.GetEndAtom()
        if begin.GetAtomicNum() == C:
            carbon = begin
            nitro = end
        else:
            carbon = end
            nitro = begin
        oxy = self._findAttached(carbon, O)
        if not oxy:
            print "WARNING! MISSING O FOR AMIDE!"
            return
        hydro = self._findAttached(nitro, H)
        if not hydro:
            #print "WARNING! MISSING H FOR AMIDE!"
            return
        dihe = self.calcDihedral(oxy, carbon, nitro, hydro)
        if abs(dihe) <= r180:
            deviation = abs(dihe) 
        else:
            deviation = abs(dihe) - r180
        if deviation <= tol:
            if self.verbose: 
                print ('[VERBOSE] Amide TRANS conformation detected'
                        ' ( %3.2f deg, deviation: %3.3f): OK' % (dihe, deviation  ))
            return
        if self.verbose: 
                print ('[VERBOSE] Amide CIS conformation '
                       'detected ( %2.2f deg): FIXING!' % math.degrees(dihe) )
        angleFix = ob.double_array([-dihe])[0]
        self.mol.SetTorsion(oxy,carbon,nitro,hydro, angleFix)
        dihe2 = self.calcDihedral(oxy, carbon, nitro, hydro)


    def _flipDihedral(self, atoms=[], bond=None, angle=None):
        """ rotate a bond by the required angle"""
        idx = [ x.GetIdx() for x in atoms]
        rotor = ob.OBRotor()
        rotor.SetBond(bond)
        rotor.SetDihedralAtoms(idx)
        return
            

    def calcDihedral(self, a1, a2, a3, a4): 
        """ given 4 OBAtom return the dihedral 
            angle between them
        """
        v1 = self._atToVec(a1, a2)
        v2 = self._atToVec(a3, a2)
        v3 = self._atToVec(a3, a4)

        v4 = np.cross(v1, v2)
        v5 = np.cross(v2, v4)
        try:
            dihe = math.atan2(np.dot(v3,v4), np.dot(v3,v5) * math.sqrt(np.dot(v2,v2)))
        except ZeroDivisionError:
            dihe = 0.
        return dihe

    def _atToVec(self, a1, a2):
        """ return the vector between two atoms"""
        c1 = self._atcoord(a1)
        c2 = self._atcoord(a2)
        return c2-c1
        
    def _atcoord(self, a):
        """ given an atom return its coords as float"""
        coord = [ a.GetX(), a.GetY(), a.GetZ() ]
        return np.array(coord, 'f')

    def _findAttached(self, atom, anum):
        """ return the first neighbor of atom that is of type anum"""
        for n in ob.OBAtomAtomIter(atom):
            if n.GetAtomicNum() == anum:
                return n

    def setFFoutlev(self):
        """ set reduced output for the forcefield"""
        pass


    def initForcefield(self):
        """ retrieve known forcefields"""
        # XXX change to INIT FORCEFIELD
        forcefields =  ob.vectorString()
        ob.OBPlugin.ListAsVector('forcefields', None, forcefields)
        self.knownForcefields = [x.split()[0].lower() for x in forcefields]
        if not self.ff in self.knownForcefields:
            print "ERROR: Unknown forcefield! [%s]" % self.ff
            self.ready = False
            return False
        self.forcefield = ob.OBForceField.FindType(self.ff)
        # setup the outlevel
        self.setFFoutlev()
        
        outlev = ob.OBFF_LOGLVL_NONE
        if self.verbose:
            outlev = ob.OBFF_LOGLVL_LOW
        self.forcefield.SetLogLevel(outlev)

    def setForcefield(self, ff):
        """ """
        ff= ob.OBForceField.FindType(ff)
        outlev = ob.OBFF_LOGLVL_NONE
        if self.verbose:
            outlev = ob.OBFF_LOGLVL_LOW
        ff.SetLogLevel(outlev)
        # XXX ff.SetLogFile() here!
        # setup molecule
        setup = ff.Setup(self.mol)

        
        if self.verbose:
            buff = '[setForceField] Atom types assigned:\n'
            idx, atype = [], []
            for a in ob.OBMolAtomIter(self.mol):
                idx.append(a.GetIdx())
                atype.append(a.GetType())
            idx = [str(x) for x in idx]
            buff += '  AtomIdx \t' + ','.join(idx) + '\n'
            buff += '  AtomType\t' +  ','.join(atype) + '\n'
            self.vprint(buff)

        if not setup:
            if self.strict:
                msg =  ('MISSING FF PARAMETERS (--strict mode on) forcefield parameters'
                       'not available for molecule -> rejecting')
                self.ready = False
                return False
            print "WARNING: Missing parameters for molecule [%s] !" % self.molName
        return ff

    def vprintCan(self, msg='VPRINTCAN'):
        """ print canonical SMI string"""
        if not self.verbose:
            return
        canonicForm = ob.OBMol()
        obc = ob.OBConversion()
        obc.SetInAndOutFormats('smi', 'smi')
        can = obc.WriteString(self.mol)
        self.vprint('[VPRINTCAN] Generated canonical form: %s' %can.strip())
        print '[%s] %s' % (msg, can.strip())
        

    def make3D(self, canonicize=False, extra3Dmini=False):
        """ make 3D on request 
            by default, a canonical SMILES form is generated
            to avoid twisted/interconnected groups of 
            atoms ? USEFUL?
        """
        #self.verbose = True
        #canonicize = False

        #obc = ob.OBConversion()
        #obc.SetInAndOutFormats('mol2', 'can')
        #can = obc.WriteFile(self.mol, 'pre_all.can')

        if canonicize or False: # possibly useless...
            canonicForm = ob.OBMol()
            obc = ob.OBConversion()
            obc.SetInAndOutFormats('smi', 'smi')
            can = obc.WriteString(self.mol)
            obc.ReadString(self.mol,can)
            self.vprint('[make3D] Generated canonical form: %s' %can.strip())
        self.vprint('[make3D] Generating 3D structure')
        #obc = ob.OBConversion()
        #obc.SetInAndOutFormats('mol2', 'mol2')
        #can = obc.WriteFile(self.mol, 'pre_builder.mol2')
        #print "BUILDER"
        outcome = self.builder.Build(self.mol)
        #bbb = ob.OBBuilder()
        #outcome = bbb.Build(self.mol)
        # check if there are hydrogens
        self.mol.SetDimension(3)
        self.mol.AddHydrogens(False, False)
        # as from Pybel, generation of initial 3D coords is made with blessed mmff94s
        #print "\n\n\nNO 3D GEN MINI!!!!"
        #return 
        ff = self.setForcefield('mmff94s')
        ff = self.setForcefield('mmff94s')
        setup_out = ff.Setup(self.mol)
        if not setup_out:
            print "ERROR SETTING UP THE MOLECULE!"
            self.ready=False
            return
        #if not self.ready:
        #    return
        #print "SETUP", setup_out
        ff.EnableCutOff(True);
        ff.SetVDWCutOff(10.0);
        ff.SetElectrostaticCutOff(20.0);
        ff.SetUpdateFrequency(10)
        steps = 500
        if extra3Dmini:
            steps = 1000
        self.minimize(ff=ff, gradient='sd', steps=steps)
        #self.verbose = False
        if self.debug:
            obc = ob.OBConversion()
            obc.SetInAndOutFormats('mol2', 'mol2')
            can = obc.WriteFile(self.mol, 'builder.mol2')


    def optimize(self):
        """ generate optmized structure"""
        ff = self.setForcefield(self.ff)
        if not self.ready:
            return
        # XXX Add logging here
        #print "%s:" % self.name, 
        if self.cgsteps == self.sdsteps == 0:
            print " (minimization skipped)"
        # steepest descent
        if self.sdsteps:
            self.runGradient(ff, 'sd', steps=self.sdsteps, convergence=self.sdconv)

        # rotameric search
        if self.rotamerConf: # and self.rotamerGeomSteps:
            self.rotamerSearch(ff, numConf=self.rotamerConf, geomSteps=self.rotamerGeomSteps)

        # conjugated gradients
        if self.cgsteps:
            self.runGradient(ff, 'cg', steps=self.cgsteps, convergence=self.cgconv)

    def rotamerSearch(self, ff, numConf=10, geomSteps=1):
        """ perform rotameric search"""
        self.vprint("[rotamerSearch] numConf=%d, geomSteps=%d" % (numConf, geomSteps))
        coord = self.mol.GetAtom(1).GetX()
        ff.WeightedRotorSearch(numConf, geomSteps)
        ff.UpdateCoordinates(self.mol)
        #print "COORDA PRE POST", coord, self.mol.GetAtom(1).GetX()

    def postOptimize(self):
        """ perform minimzation post-processing (by default using UFF)
            to overcome limitations of MMMFF94  (aromatic amines, etc...)
        """
        #print "POSTOPT"
        ff = self.setForcefield(self.ff_extra)
        if not self.ready:
            return
        # add logging here
        #print "%s:" % self.name, 
        self.vprint("Extra minimization called:")
        if self.sdsteps_extra:
            self.vprint("Steepest descent (extra)")
            self.runGradient(ff, 'sd', steps=self.sdsteps_extra, convergence=self.sdconv_extra)
        if self.cgsteps_extra:
            self.vprint("Conjugated gradients (extra)")
            self.runGradient(ff, 'cg', steps=self.cgsteps_extra, convergence=self.cgconv_extra)
        if self.cgsteps == self.sdsteps == 0:
            print " (minimization skipped)"


    def runGradient(self, ff, gradient='sd', steps=0, convergence=0):
        """ generalized SD/CG gradient runner"""
        messages = { 'sd' : { False: 'SD', True: "Steepest descent..."},
                     'cg' : { False: 'CG', True: "Conjugated gradients..."}
                    }
        #print messages[gradient][self.verbose]
        stepsDone,convergenceDone = self.minimize(ff, gradient, steps=steps, convergence=convergence)
        if convergenceDone:
            if self.verbose: 
                print "converged", stepsDone
            if not self.verbose:
                #print "[0],",
                pass
        else:
            if not self.verbose:
                #print "[1],",
                pass
        sys.stdout.flush()

    def minimize(self, ff=None, gradient='sd', steps=None, convergence=1e-5, pace=10):
        """ minimize energy by combining multiple minimizers """
        convergence = ob.double_array([float(convergence)])[0]
        if self.verbose:
            if gradient == 'sd':
                init = ff.SteepestDescentInitialize
                _minimizer = ff.SteepestDescentTakeNSteps
            else:
                init = ff.ConjugateGradientsInitialize
                _minimizer = ff.ConjugateGradientsTakeNSteps
            init(steps, convergence)
            keep = False
            pool = steps
            while pool > 0:
                pool -= pace
                pool = max(pool, 0)
                if not _minimizer(pace):
                    if pool:
                        break
                msg = "[minimize] %s: steps %d/%d     "
                args = (gradient.upper(), steps-pool, steps)
                self.vprint(msg % args, addnewline=False, flush=True, rewind=True)
        else:
            if gradient == 'sd':
                gradientObj = ff.SteepestDescent
            else:
                gradientObj = ff.ConjugateGradients
            gradientObj(steps, convergence)
            
        self.vprint(" ")
        #ff.GetCoordinates(self.mol)
        ff.UpdateCoordinates(self.mol)
        if self.verbose:
            converged = bool(pool)
            if converged:
                msg = "[minimize] %s: convergence reached (%d steps)" 
                args = (gradient.upper(),steps-pool)
            else: 
                msg = '[minimize] %s: max. steps reached (%d steps)'
                args = (gradient.upper(), steps-pool)
            self.vprint(msg % args)
            return (steps-pool, converged )
        else:
            return (0, True)


    def setUniqueAtomNames(self):
        """ create unique names for each atom in the molecule"""
        etable = ob.OBElementTable()
        self.mol._atomNames = {}
        resList = [x for x in ob.OBResidueIter(self.mol)]
        if len(resList) == 0:
            res = self.mol.CreateResidue()
            for a in ob.OBMolAtomIter(self.mol):
                res.AddAtom(a)
            resList = [ res ]

        for res in resList:
            for atom in ob.OBResidueAtomIter(res):
                name = res.GetAtomID(atom)
                if name == "":
                    name = etable.GetSymbol( atom.GetAtomicNum() )
                # print "NAME", name,
                if not name in self.mol._atomNames:
                    self.mol._atomNames[name] = 1
                    #newName = name
                else:
                    self.mol._atomNames[name] += 1
                newName = "%s%d" % (name,  self.mol._atomNames[name])
                #print newName
                res.SetAtomID(atom, newName)

