# matlab.py
# 
# Bryan Daniels
# 7.18.2011, 8.19.2011
#
# For calling MATLAB code from Python.
#

from subprocess import call, STDOUT # for callMATLAB
import scipy, os

# 7.18.2011 generalizing from what was in inverseIsing
def callMATLAB(callString,outputFilename=None,codeDirList=[],               \
    tmpStdoutFilename="MATLAB_tmpStdout.txt"):
    """
    Calls MATLAB code given in callString. 
    
    outputFilename (None)   : Looks for output and returns it if 
                              possible.  If no output is there, raises 
                              error with text from MATLAB that hopefully 
                              expresses the problem.
                              Can also be a list of outputFilenames.
    codeDir ([])            : An optional list of code directories
                              for MATLAB to look in for code.
    
    
    """
    codeDirStr = ""
    for codeDir in codeDirList:
        codeDirStr += "path(path,'"+codeDir+"');"
    
    stdoutFile = open(tmpStdoutFilename,'w')
    
    # 2.29.2012 remember stty settings to fix weird Terminal problem
    sttyFilename = "temp_stty_args_for_MATLABpy.txt"
    sttyFile = open(sttyFilename,'w')
    call(["stty","-g"],stdout=sttyFile)
    sttyFile.close()
    sttyFile = open(sttyFilename,'r')
    sttyArgs = sttyFile.readline()
    sttyFile.close()
    os.remove(sttyFilename)
    
    # set MATLAB path to matlabCodeDirectories and call the MATLAB code
    call(["matlab","-nodesktop","-nosplash","-nojvm","-nodisplay",          \
          "-r","try,"+codeDirStr+callString+                                \
          ";catch exception,getReport(exception),end,exit;"],               \
          stderr=stdoutFile,stdout=stdoutFile)
    stdoutFile.close()
    
    # 2.29.2012 fix weird Terminal problem 
    # see also http://www.mathworks.com/matlabcentral/newsreader/view_thread/288861
    call(["stty",sttyArgs])

    if outputFilename is not None:
        if type(outputFilename) == str: outputFilenameList = [outputFilename]
        else: outputFilenameList = outputFilename
        try:
            outputList = []
            for outputFilename in outputFilenameList:
                outputList.append( scipy.loadtxt(outputFilename) )
                os.remove(outputFilename)
        except IOError:
            # there was a MATLAB error
            print("callMATLAB MATLAB error:")
            stdoutFile = open(tmpStdoutFilename)
            stdout = stdoutFile.read()
            print(stdout)
            print("callMATLAB: Current directory is",os.getcwd())
            raise Exception("callMATLAB MATLAB error.")
        if len(outputList)==1: return outputList[0]
        return outputList
