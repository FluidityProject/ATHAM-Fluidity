try:
    import _FortranMicrophysics
    FM=_FortranMicrophysics.fw_auto
except ImportError:
    print "No _FortranMicrophysics module"
    print "Module must be built"
    print "See python/Fortran section of Microphysics.pdf"
import numpy as np
from FortranMicrophysicsWrapper import MakeWrapperFiles

"""This module provides convenience wrappers for the python side of the python/Fortan interface"""


### Objects defined for testing purposes only

fields={'gas_tracers'   :[('WaterVapour','MassFraction')],
        'incomp_tracers':[],
        'other_fields'  :[('Bulk','IteratedEOSDensity')]
 }
call_string='microphysics_main_pointwise'

class M(object):
    def __init__(self,s):
        self.scalar_fields={}
        for S,k in s:
            self.scalar_fields[S]=k
class N(object):
    def __init__(self,R):
        self.val=R
    

States={'Bulk':M((('IteratedEOSDensity',N(np.ones(5,order='F'))),)),
       'WaterVapour':M((('MassFraction',N(2.0*np.ones(5,order='F'))),))}
    
### end of test objects

def UnwrapFields(states=States,fields_dict=fields):
    Setup(fields_dict)
    for key,fields in fields_dict.items():
        for n,(state,field) in enumerate(fields):
            Fnew=states[state].scalar_fields[field].val
            try:
                FOld=states[state].scalar_fields['Old'+field].val
            except KeyError:
                FOld=np.zeros(0,order='F')
            try:
                FSource=states[state].scalar_fields[field+'Source'].val
            except KeyError:
                FSource=np.zeros(0,order='F')
            getattr(FM,'set_'+key)(n+1,F,FOld,FSource)

def Run(current_time,dt):
    FM.run_microphysics(current_time,dt)

def Setup(fields_dict=fields):
    nm=[]
    for key,fields in fields_dict.items():
        nm.append(len(fields))
    FM.allocate_storage(nm)

def Finish():
    FM.finalize()

if __name__=="__main__":
    MakeWrapperFiles(fields,call_string,pointwise=True)
    
