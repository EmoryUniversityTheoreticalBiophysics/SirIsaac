# outputTag.py
#
# Bryan Daniels
# 1.29.2010
#
# Use: python outputTag.py codeToRun.py
#
# Runs python script in background, passing a unique identifier tag
# to the script and piping stdout and stderr to a file
# named with the same tag.
#
# Finds next available tag and runs, eg,
#     python codeToRun.py 0001 &> 0001_stdout_stderr.txt &
#
# Advanced use: python outputTag.py [-o outputDirectory] [-n niceLevel]
#               [-m multipleRunNumber] [-d runDirectory] codeToRun.py"
#

import sys, os, getopt

def nextFileNumString(directory='.',returnConfigNum=False,configNum=None):
    # (make sure you don't start runs in rapid succession, or you might
    #  start before the old one has written a file)
    if configNum is None:
        filenames = os.listdir(directory)
        configNum = 1
        i = 0
        while i < len(filenames):
            configNumString = '%(c)04d' % {'c':configNum}
            if (filenames[i][:4]==configNumString                               \
             or filenames[i][1:5]==configNumString):
                configNum += 1
                i = 0
            else:
                i += 1
    configNumString = '%(c)04d' % {'c':configNum}
    if os.uname()[1][:4] == 'node': # emory machines
        if os.getcwd().startswith('/star'):
            configNumString = 's' + configNumString
        if os.getcwd().startswith('/spark'):
            configNumString = 'k' + configNumString
    if returnConfigNum: return configNum
    else: return configNumString

def die():
    print "Use: python outputTag.py [-o outputDirectory] [-n niceLevel] "       \
          "[-m multipleRunNumber] [-d runDirectory] codeToRun.py"
    sys.exit()

if __name__ == '__main__':
    
    # see example http://www.faqs.org/docs/diveintopython/kgp_commandline.html
    
    outputDirectory,niceLevel = None,None
    multipleRunNumber = 1
    runDirectory = None
    try:
        optsList, otherArgs = getopt.getopt(sys.argv[1:], "o:n:t:m:d:")
    except:
        die()
    for opt, arg in optsList:
        if opt == '-o':
            outputDirectory = arg
            if outputDirectory[-1] != "/":
                outputDirectory = outputDirectory + "/"
        if opt == '-n':
            niceLevel = arg
        if opt == '-m':
            multipleRunNumber = int(arg)
        if opt == '-d':
            runDirectory = arg
    if outputDirectory is None:
        outputDirectory = './'
    
    if (len(sys.argv) < 2):
      die()
    if otherArgs[0][-3:] != '.py':
      die()
    
    if os.uname()[1][:8] == 'vader': pythonName = "python2.7 "
    else: pythonName = "python "
    
    configNumInit = nextFileNumString(outputDirectory,returnConfigNum=True)
    
    for i in range(multipleRunNumber):

        tag = nextFileNumString(outputDirectory,configNum=configNumInit+i)

        #command = pythonName+otherArgs[0]+" '"+outputDirectory+tag+"' &> '"     \
        #    +outputDirectory+"stdout_stderr_"+tag+".txt' &"
        command = pythonName+otherArgs[0]+" '"+outputDirectory+tag+"' "+str(i)+ \
            " &> '"+outputDirectory+"stdout_stderr_"+tag+".txt' &"
        if niceLevel is not None:
            command = "nice -n "+str(niceLevel)+" "+command
        if runDirectory is not None:
            os.chdir(runDirectory)
        os.system(command)
        print command
