# powerLawYeastOscillator.py
# 
# Bryan Daniels
# 2.6.2013
#
# Sets up yeast oscillator as a 19-dimensional power-law network.

from FittingProblem import *

class PowerLawFittingModel_yeastOscillator(PowerLawFittingModel_FullyConnected):
    """
    Models yeast oscillator as a 19-dimensional power-law network.
    
    See notes 2.6.2013 - 2.13.2013 and Ruoff_model_original.m.
    """
    
    # (the first in each list is the one that gets the initial condition)
    definitionDict = { 'S2': [('S2A',1),('S2B',1)],
                       'S3': [('S3A','xi_S3A'),('S3B','xi_S3B')],
                       'S4': [('S4A','xi_S4A'),('S4B','xi_S4B'),('S4C','xi_S4C')],
                       'N2': [('N2A',1),('N2B',1)],
                       'A3': [('A3A','xi_A3A'),('A3B','xi_A3B'),('A3C','xi_A3C')],
                       'v1': [('v1A','xi_v1A'),('v1B',1),('v1C',1),('v1D',1),('v1E',1)] }
            

    def __init__(self,temperature=286.5,prune=True,**kwargs):
        """
        Models yeast oscillator as a 19-dimensional power-law network.
        
        See notes 2.6.2013 - 2.13.2013 and Ruoff_model_original.m.
        
        temperature (286.5)     : In Kelvin.  Rates are dependent on 
                                  temperature according to an Arrhenius
                                  law.
                                  
        10.16.2013 Default temperature changed to 286.5 K.
        """
        
        self.speciesNames = scipy.array(['S1','S2A','S2B','S3A','S3B',          \
            'S4A','S4B','S4C','N2A','N2B','A3A','A3B','A3C','S4ex',             \
            'v1A','v1B','v1C','v1D','v1E'])
        self.ICnames = [ name+"_init" for name in self.speciesNames ]
        
        # () store typical ranges for initial conditions
        # taken from SchValJen11 Table 2
        # (ones are for extra recasting variables)
        #ICranges = scipy.array(                                                 \
        #    [[0.15,1.60],                                                       \
        #     [0.19,2.16],[1.,1.],                                               \
        #     [0.04,0.20],[1.,1.],                                               \
        #     [0.10,0.35],[1.,1.],[1.,1.],                                       \
        #     [0.08,0.30],[1.,1.],                                               \
        #     [0.14,2.67],[1.,1.],[1.,1.],                                       \
        #     [0.05,0.10],                                                      \
        #     [1.,1.],[1.,1.],[1.,1.],[1.,1.],[1.,1.]] ) # mM
        #self.indepParamRanges = ICranges
        # taken from SchValJen11 Table 2
        ICranges = scipy.array(                                                 \
                   [[0.15,1.60],[0.19,2.16],                                    \
                    [0.04,0.20],[0.10,0.35],                                    \
                    [0.08,0.30],[0.14,2.67],[0.05,0.10]] ) # mM
        self.indepParamRanges = ICranges
        self.typicalTimeRange = [0.,3] # minutes
        
        # () set up a fully-connected 19-dimensional model
        PowerLawFittingModel_FullyConnected.__init__(self,19,                   \
            indepParamNames=self.ICnames,outputNames=self.speciesNames,         \
            includeRegularizer=False,logParams=False,useDeltaGamma=False,       \
            **kwargs)
        
        # () set the constant original model parameters
        # parameter names are J0,k1,k2,k3,k4,k5,k6,k,kappa,q,K1,psi,N,A
        net = self.net
        
        # temperature-independent parameters q,psi,N,A
        q, psi = 4., 0.1
        N, A = 1., 4.
        net.addParameter('q',q,isOptimizable=False)
        net.addParameter('psi',psi,isOptimizable=False)
        net.addParameter('N',N,isOptimizable=False)
        net.addParameter('A',A,isOptimizable=False)
        
        # temperature-dependent parameters J0,k1,k2,k3,k4,k5,k6,k,kappa,K1
        # (same as in Ruoff_model_original.m)
        # energetic/enthalpic barriers in kJ/mol
        #   from Table 4 from the Ruoff paper RuoChrWol03
        # rates in mM/min (except for K1, in mM)
        # (I think k1,k2,k3,k4,k6 really have units 1/mM/min)
        #   k1,k2 from "Fit of model to Hemker et al. data"
        #   at end of Fig.3 caption; rest from Table 1
        barriersAndRates = {'J0': (16.2, 2.5 ),     'k1': (13.8, 25. ),         \
                            'k2': (60.7, 2.0 ),     'k3': (44.9, 16.0),         \
                            'k4': (58.7, 100.),     'k5': (41.2, 1.28),         \
                            'k6': (31.4, 12.0),     'k':  (24.3, 1.8 ),         \
                            'kappa': (15.9, 13.0),  'K1': (47.0, 0.52)}
        T,Tref = temperature, 286.5 # K
        R = 0.0083144 # kJ/K/mol
        tempDepRates = dict( [ (k,br[1]*scipy.exp(br[0]/(R*Tref)-br[0]/(R*T)))  \
                               for k,br in barriersAndRates.items() ] )
        for name,rate in tempDepRates.items():
            net.addParameter(name,rate)
        
        # () add assignment rules and composite species
        comp = net.compartments[0].id
        net.addParameter('y',1.,isOptimizable=False)
        net.addAssignmentRule('y','q/(k1*K1**q)')
        for name,defList in self.definitionDict.items():
            # also set initial condition
            net.addParameter(name+"_init",1.0,isOptimizable=False)
            net.addSpecies(name,comp,name+"_init")
            if False: # old as of 2.26.2013
                productStr = str(defList).replace("', '","*")[2:-2]
                #print "productStr =",productStr
                net.addAssignmentRule(name,productStr)
            else: # define using odes
                rhs = ''
                # derivative is sum of derivative of one factor
                # times the other factors
                for index in range(len(defList)):
                  diName,diExp = defList[index]
                  rhs += ' + (' + net.rateRules.get(diName) + ')'
                  rhs += '*' + str(diExp) + '*' + diName + '**-1'
                  for otherIndex in range(len(defList)):
                    djName,djExp = defList[otherIndex]
                    rhs += '*' + djName + '**' + str(djExp)
                #print "Rate for",name,"=",rhs
                net.addRateRule(name,rhs)
        
        # 3.4.2013 add variable definition parameters
        # (tuned to make species stay closer to O(1))
        net.addParameter('xi_S3A',13.,isOptimizable=True)
        net.addParameter('xi_S3B',13.,isOptimizable=True)
        net.addParameter('xi_v1A',10.,isOptimizable=True)
        net.addParameter('xi_S4C',10.,isOptimizable=True)
        net.addParameter('xi_S4B',10.,isOptimizable=True)
        net.addParameter('xi_A3C',10.,isOptimizable=True)
        net.addParameter('xi_S4A',10.,isOptimizable=True)
        net.addParameter('xi_A3A',10.,isOptimizable=True)
        net.addParameter('xi_A3B',10.,isOptimizable=True)
                    
        
        # set initial conditions
        for name in self.speciesNames:
            if name.endswith("A"):
                pooledVarName = name[:-1]
                pooledVarInitName = pooledVarName+"_init"
                #net.addParameter(pooledVarInitName,1.0,isOptimizable=False)
                defExp = self.definitionDict[pooledVarName][0][1]
                initVal = pooledVarInitName+'**(1./'+str(defExp)+')'
                net.setInitialVariableValue(name,initVal)
            else:
                net.setInitialVariableValue(name,name+"_init")
        v1initStr = "(k1*S1_init*A3_init / (1. + (A3_init/K1)**q ))"
        defExp = self.definitionDict['v1'][0][1]
        net.setInitialVariableValue("v1A",v1initStr+'**(1./'+str(defExp)+')')
        self.indepParamNames = [ 'S1_init','S2_init','S3_init','S4_init',       \
                                 'N2_init','A3_init','S4ex_init' ]
        
        # () set the (sparse) nonzero structural parameters
        # see notes 2.6.2013
        self._setTerm('S1', +1,'J0',[])
        self._setTerm('S1', -1,'1',[('v1',1)])
        
        self._setTerm('S2A',+1,'k2',[('N2',1),('S2A',1)])
        self._setTerm('S2A',-1,'k6',[('N2',1),('S2A',1)])
        self._setTerm('S2B',+1,'2.',[('v1',1),('S2A',-1)])
        self._setTerm('S2B',-1,'k2*N',[('S2B',1)])
        
        self._setTerm('S3A',+1,'k3/xi_S3A',[('A3',1),('S3A',1)])
        self._setTerm('S3A',-1,'k2/xi_S3A',[('N2',1),('S2',1),('S3A','1-xi_S3A'),('S3B','-xi_S3B')])
        self._setTerm('S3B',+1,'k2*N/xi_S3B',[('S2',1),('S3A','-xi_S3A'),('S3B','1-xi_S3B')])
        self._setTerm('S3B',-1,'k3*A/xi_S3B',[('S3B',1)])
        
        #self._setTerm('N2A',+1,'0.',[])
        #self._setTerm('N2A',-1,'k2+k6',[('S2',1),('N2A',1)])
        #self._setTerm('N2B',+1,'k2*N',[('S2',1),('N2A',-1)])
        #self._setTerm('N2B',-1,'k4',[('S4',1),('N2B',1)])
        self._setTerm('N2A',+1,'k2*N/2.',[('S2',1),('N2B',-1)])
        self._setTerm('N2A',-1,'k2+k6',[('S2',1),('N2A',1)])
        self._setTerm('N2B',+1,'k2*N/2.',[('S2',1),('N2A',-1)])
        self._setTerm('N2B',-1,'k4',[('S4',1),('N2B',1)])

        #self._setTerm('A3A',+1,'0.',[])
        #self._setTerm('A3A',-1,'k5',[('A3A',1)])
        #self._setTerm('A3B',+1,'2.*k3*A',[('S3',1),('A3A',-1),('A3C',-1)])
        #self._setTerm('A3B',-1,'2.*k3',[('A3B',1),('S3',1)])
        #self._setTerm('A3C',+1,'0.',[])
        #self._setTerm('A3C',-1,'2.',[('v1',1),('A3A',-1),('A3B',-1)])
        self._setTerm('A3A',+1,'1.*k3*A/3./xi_A3A',[('S3',1),('A3A','1-xi_A3A'),('A3B','-xi_A3B'),('A3C','-xi_A3C')])
        self._setTerm('A3A',-1,'k5/xi_A3A',[('A3A',1)])
        self._setTerm('A3B',+1,'1.*k3*A/3./xi_A3B',[('S3',1),('A3A','-xi_A3A'),('A3B','1-xi_A3B'),('A3C','-xi_A3C')])
        self._setTerm('A3B',-1,'2.*k3/xi_A3B',[('A3B',1),('S3',1)])
        self._setTerm('A3C',+1,'4.*k3*A/3./xi_A3C',[('S3',1),('A3A','-xi_A3A'),('A3B','-xi_A3B'),('A3C','1-xi_A3C')])
        self._setTerm('A3C',-1,'2./xi_A3C',[('v1',1),('A3A','-xi_A3A'),('A3B','-xi_A3B'),('A3C','1-xi_A3C')])

        net.addParameter('theta_S4',0.5,isOptimizable=True)
        #self._setTerm('S4A',+1,'kappa',[('S4ex',1),('S4B',-1),('S4C',-1)])
        #self._setTerm('S4A',-1,'kappa',[('S4A',1)])
        #self._setTerm('S4B',+1,'0.',[])
        #self._setTerm('S4B',-1,'k4',[('S4B',1),('N2',1)])
        self._setTerm('S4A',+1,'(theta_S4)*kappa/xi_S4A',[('S4ex',1),('S4A','1-xi_S4A'),('S4B','-xi_S4B'),('S4C','-xi_S4C')])
        self._setTerm('S4A',-1,'k3/xi_S4A',[('A3',1),('S3',1),('S4A','1-xi_S4A'),('S4B','-xi_S4B'),('S4C','-xi_S4C')])
        self._setTerm('S4B',+1,'(1.-theta_S4)*kappa/xi_S4B',[('S4ex',1),('S4A','-xi_S4A'),('S4B','1-xi_S4B'),('S4C','-xi_S4C')])
        self._setTerm('S4B',-1,  'kappa/xi_S4B',[('S4B',1)])
        self._setTerm('S4C',+1,'k3*A/xi_S4C',[('S3',1),('S4A','-xi_S4A'),('S4B','-xi_S4B'),('S4C','1-xi_S4C')])
        self._setTerm('S4C',-1,  'k4/xi_S4C',[('S4C',1),('N2',1)])
        
        self._setTerm('S4ex',+1,'psi*kappa',[('S4',1)])
        self._setTerm('S4ex',-1,'psi*kappa+k',[('S4ex',1)])
        
        #self._setTerm('v1A',+1,'2.*y',[('v1',2),('v1A',1),('S1',-1),('A3',q-2)])
        #self._setTerm('v1A',-1,'2.',[('v1',1),('v1A',1),('A3',-1)])
        #self._setTerm('v1B',+1,'2.*k3*A',[('v1B',1),('S3',1),('A3',-1)])
        #self._setTerm('v1B',-1,'2.*y*k3*A',[('v1',1),('v1B',1),('S3',1),('S1',-1),('A3',q-2)])
        self._setTerm('v1A',+1,'2.*y/xi_v1A',[('v1',2),('v1A',1),('S1',-1),('A3',q-2)])
        self._setTerm('v1A',-1,    '2.*y*k3*A/xi_v1A',[('v1',1),('v1A',1),('S3',1),('S1',-1),('A3',q-2)])
        self._setTerm('v1B',+1,'2.*k3*A',[('v1B',1),('S3',1),('A3',-1)])
        self._setTerm('v1B',-1,    '2.',[('v1',1),('v1B',1),('A3',-1)])
        self._setTerm('v1C',+1,'2.*y*k3',[('v1',1),('v1C',1),('S3',1),('S1',-1),('A3',q-1)])
        self._setTerm('v1C',-1,'2.*k3',[('v1C',1),('S3',1)])
        self._setTerm('v1D',+1,'y*k5',[('v1',1),('v1D',1),('S1',-1),('A3',q-1)])
        self._setTerm('v1D',-1,'k5',[('v1D',1)])
        self._setTerm('v1E',+1,'J0',[('v1E',1),('S1',-1)])
        self._setTerm('v1E',-1,'1.',[('v1',1),('v1E',1),('S1',-1)])

        if prune: self.prune()
                


# 2.14.2013
class PowerLawFittingModel_stirredTank(PowerLawFittingModel_FullyConnected):
    """
    Example s-system power law network from SavVoi87 (p. 98).
    """
    
    def __init__(self,prune=True,**kwargs):
        """
        Example s-system power law network from SavVoi87 (p. 98).
        """
        
        self.speciesNames = scipy.array(['X1','X2','X3','X4','X5','X6','X7','X8'])
        self.ICnames = [ name+"_init" for name in self.speciesNames ]
        
        # () set up a fully-connected 8-dimensional model
        PowerLawFittingModel_FullyConnected.__init__(self,8,                    \
            indepParamNames=self.ICnames,outputNames=self.speciesNames,         \
            includeRegularizer=False,logParams=False,useDeltaGamma=False,       \
                                                     **kwargs)
        
        # () set the constant original model parameters
        # parameter names are A,B,C,b,c
        net = self.net
        
        # temperature-independent parameters q,psi,N,A
        A = 0.103
        B = 19.*A
        c = 50.
        b = 3.
        a = 1.
        C = a*b + 2.*(1.+b)
        net.addParameter('A',A,isOptimizable=True)
        net.addParameter('B',B,isOptimizable=True)
        net.addParameter('C',C,isOptimizable=True)
        net.addParameter('b',b,isOptimizable=True)
        net.addParameter('c',c,isOptimizable=True)

        # () add assignment rules
        comp = net.compartments[0].id
        net.addSpecies('Z1',comp)
        net.addAssignmentRule('Z1','X1*X2 - 1.')
        net.addSpecies('Z2',comp)
        net.addAssignmentRule('Z2','X3*X4 - 2.')

        # () set initial conditions
        Z1_init,Z2_init = 1.,1.
        net.setInitialVariableValue('X1',1)
        net.setInitialVariableValue('X2',Z1_init+1)
        net.setInitialVariableValue('X3',1)
        net.setInitialVariableValue('X4',Z2_init+2)
        net.setInitialVariableValue('X5',1)
        net.setInitialVariableValue('X6',scipy.exp(c*Z2_init/(c+Z2_init)))
        net.setInitialVariableValue('X7',1)
        net.setInitialVariableValue('X8',(c+Z2_init)/c)

        # () set the (sparse) nonzero structural parameters
        self._setTerm('X1', +1,'1.',[('X2',-1)])
        self._setTerm('X1', -1,'1.',[('X1',1)])

        self._setTerm('X2', +1,'2.*A',[('X1',-1),('X5',1),('X6',1)])
        self._setTerm('X2', -1,'A',[('X2',1),('X5',1),('X6',1)])
        
        self._setTerm('X3', +1,'C',[('X4',-1)])
        self._setTerm('X3', -1,'(1.+b)',[('X3',1)])
        
        self._setTerm('X4', +1,'2.*B',[('X3',-1),('X5',1),('X6',1)])
        self._setTerm('X4', -1,'B',[('X1',1),('X2',1),('X3',-1),('X5',1),('X6',1)])
        
        self._setTerm('X5', +1,'C',[('X5',1),('X7',-2),('X8',-2)])
        self._setTerm('X5', -1,'(1.+b)',[('X3',1),('X4',1),('X5',1),('X7',-2),('X8',-2)])
        
        self._setTerm('X6', +1,'2.*B',[('X5',1),('X6',2),('X7',-2),('X8',-2)])
        self._setTerm('X6', -1,'B',[('X1',1),('X2',1),('X5',1),('X6',2),('X7',-2),('X8',-2)])
        
        self._setTerm('X7', +1,'C/c',[('X8',-1)])
        self._setTerm('X7', -1,'(1.+b)/c',[('X3',1),('X4',1),('X8',-1)])
        
        self._setTerm('X8', +1,'2.*B/c',[('X5',1),('X6',1),('X7',-1)])
        self._setTerm('X8', -1,'B/c',[('X1',1),('X2',1),('X5',1),('X6',1),('X7',-1)])

        
        if prune: self.prune()



