# PhosphorylationFit_netModel.py
#
# Bryan Daniels
# 06.15.2009 - 06.16.2009
#
# 

import os
import scipy.io
import resource,time # for measuring cpu usage and wall clock time
import PhosphorylationFit_MakeBNGL as writeBNGL
reload(writeBNGL)

class netModel():
    """
    Class that implements a BioNetGen model of a single protein
    species that can be activated on n different sites.
    
    To do: clean up BioNetGen files automatically.
    """
    def __init__(self, n, rulesList, endTime, nSteps, filename=None,            \
        verbose=True, BNGpath = "~/Downloads/BioNetGen_2.0.46/",                \
        MichaelisMenten=True):
        
        if os.getcwd()[:3] == '/nh': # for cnls cluster
            BNGpath = "~/BioNetGen_2.0.46/"
        
        self.n = n
        self.rulesList = rulesList
        self.endTime = endTime
        self.nSteps = nSteps
        self.filename = filename
        self.verbose = verbose
        self.BNGpath = BNGpath
        self.MichaelisMenten = MichaelisMenten
        
        if filename is None:
            self.filename = self._setupFilename()
        else:
            self.filename = filename
        
        self.namesList = self._createNetwork(verbose)
        self.numParams = len(self.namesList)
        
        self.recentOutput = None
    
    def _setupFilename(self):
        """
        Finds a filename that won't be overwriting anything (hopefully).
        """
        try:
            os.mkdir('./.netModel')
        except:
            pass # hope it's already there...
        filenames = os.listdir('./.netModel')
        configNum = 1
        i = 0
        configNumString = '%(c)04d' % {'c':configNum}
        while i < len(filenames):
            configNumString = '%(c)04d' % {'c':configNum}
            if (filenames[i][:4]==configNumString):
                configNum += 1
                i = 0
            else:
                i += 1
        return os.path.realpath('.')+'/.netModel/'+configNumString
    
    def _createNetwork(self,verbose):
        """
        Creates BioNetGen network 'filename.net'.
        
        Returns list of names of network parameters.
        """
        filename,n,rulesList = self.filename,self.n,self.rulesList
        if self.verbose:
            mult = 2
            if self.MichaelisMenten:
              mult = 4
            start,startWall = cpuTime(),wallTime()
            print ""
            print "Creating network with "+str(n)+" activation sites"
            print "  and "+str(len(rulesList))+" additional rules ("                \
                  +str(mult*(n+len(rulesList)))+" parameters)."
        
        namesList = writeBNGL.writeBNGLnetwork(n,rulesList,filename,                \
            MichaelisMenten=self.MichaelisMenten)
        self._runBNGLfile(filename)
        
        if self.verbose:
            print "Network creation took "+bothTimeStr(start,startWall)
        
        return namesList

    def _readModelOutput(self,outputFilename):
        """
        Returns model output as array: 
        output[0] = time series for species 1,
        output[1] = time series for species 2, 
          ...
        """
        modelData = scipy.io.read_array(outputFilename)
        modelOutput = scipy.transpose(modelData[:,1:])
        return modelOutput
        
    def _runBNGLfile(self,filename):
        error = os.system(self.BNGpath+"Perl2/BNG2.pl "+filename+".bngl &> "        \
            +filename+"_messages.txt")
        if error:
            stdoutFile = open(filename+"_messages.txt")
            stdout = stdoutFile.read()
            print "PhosphorylationFit_netModel._runBNGLfile: Error calling BioNetGen."
            print stdout
            raise Exception, "Error calling BioNetGen."
        
    def modelOutput(self,params):
        """
        Returns model output as array: 
        output[0] = time series for species 1,
        output[1] = time series for species 2, 
          ...
        """
        netFile,namesList = self.filename,self.namesList
        endTime,nSteps = self.endTime,self.nSteps
        
        # use BNG to integrate the model
        writeBNGL.writeModifiedNet(netFile,namesList,params,"_modified")
        writeBNGL.writeBNGLsimulate(netFile+"_modified",endTime,nSteps)
        self._runBNGLfile(netFile+"_modified_simulate")
        
        # or use (slow)
        #writeBNGL.writeBNGLsimulateSlow(filename,namesList,params,endTime,nSteps)
        #os.system(BNGpath+"Perl2/BNG2.pl "+filename+"_simulate.bngl > "            \
        #       +filename+"_simulate_messages.txt")
        
        # read in output of model
        output = self._readModelOutput(netFile+"_modified_simulate.cdat")
        self.recentOutput = output
        return output

    def modelOutputSlow(self,params):
        filename,namesList = self.filename,self.namesList
        endTime,nSteps = self.endTime,self.nSteps
        
        # use BNG to integrate the model
        writeBNGL.writeBNGLsimulate(filename,namesList,params,endTime,nSteps)
        os.system(BNGpath+"Perl2/BNG2.pl "+filename+"_simulate.bngl > "             \
                +filename+"_simulate_messages.txt")
        
        # read in output of model
        output = self._readModelOutput(filename+"_simulate.cdat")
        self.recentOutput = output
        return output
        
    def writeSBML(self,params):
        writeBNGL.writeBNGL_SBML(self.filename,self.namesList,params)
        self._runBNGLfile(self.filename[:-4]+"sbml")
        return self.filename[:-4]+"sbml.xml"

# time stuff
def cpuTime():
    return resource.getrusage(resource.RUSAGE_SELF)[0]
    
def wallTime():
    return time.time()

def cpuTimeStr(startTime):
    return str(round(cpuTime()-startTime,1))+" sec cpu time"

def wallTimeStr(startTime):
    return str(round(wallTime()-startTime,1))+" sec wall time"

def bothTimeStr(cpuStart,wallStart):
    return cpuTimeStr(cpuStart)+" and "+wallTimeStr(wallStart)+"."


