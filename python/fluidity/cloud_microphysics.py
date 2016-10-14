import numpy as np
import fluidity_tools
import __builtin__

import warnings,sys

def myWarning(*args):
    ewrite(3,warnings.formatwarning(*args))

def myNumpyWarning(errstr,errflag):
    ewrite(3,errstr)

warnings.showwarning=myWarning
np.seterrcall(myNumpyWarning)
np.seterr(all='call')

if 'current_debug_level' not in __builtin__.__dict__:
    __builtin__.current_debug_level=0

def ewrite(priority_level,output):
    if priority_level<=__builtin__.current_debug_level:
        print output

class microphysics_limit(dict):
    
    def __init__(self,value,arg,logic,fargs,name=''):
        dict.__init__(self)
        self.value=value
        self.name=name
        self.arg=arg
        self.logic=logic
        for larg in fargs:
            try:
                self[larg[0]]=larg[1]
            except TypeError:
                self[larg]=1.0

class microphysics_forcing(dict):
    
    def __init__(self,forcing,args,bounds=(),alpha=0.0,name=''):
        dict.__init__(self)
        self.forcing=forcing
        self.bounds=bounds
        self.alpha=alpha
        self.name=name
        for arg in args:
            try:
                self[arg[0]]=arg[1]
            except TypeError:
                self[arg]=1.0

    def bound(self,f,alpha):
        if not self.bounds==():
            lbound=self.bounds[0]
            ubound=self.bounds[1]
            return np.where(f>0.0,f*alpha[ubound],f*alpha[lbound])
        else:
            return f

class slip_velocity(object):
    
    def _scale(self,**kwargs):
        return 1.0

    def __init__(self,field,function,scale=None):
        self.field=field
        self.function=function
        if scale==None:
            self.scale=self._scale
        else:
            self.scale=scale

    def __call__(self,kwargs):
        return self.function(kwargs,**kwargs)*self.scale(**kwargs)



class microphysics_model(object):

    """ An example factory class for python microphysics models.
    The state passed as input has the source terms of the forced 
    fields set to match the forcings provided, subject to imposed 
    constraints on the water phases, namely that they lie between 0 and 1."""

    def __init__(self,state,
                 prescribed_fields={},
                 forced_fields={},
                 python_fields=tuple(),
                 limits=tuple(),
                 forcings=tuple(),
                 slip_vels=tuple(),
                 dt=None):
        self.state=state
        self.prescribed_fields=prescribed_fields.keys()
        self.forced_fields=forced_fields.keys()
        self.slip_vels=slip_vels
        self.dt=dt
        self.forcings=forcings
        self.limits=limits
        self.fields={'old':False}
        self.oldfields={'old':True}
        self._get_fields(prescribed_fields,forced_fields)
        self._add_python_fields(python_fields)
        self._calculateLimits()
        self._calculateForcings()
        self._calculateSlipVels()
        self._addSponges()
    
    def _get_fields(self,prescribed_fields,forced_fields):
        for field,name in prescribed_fields.items():
            try:
                self.fields[field]=self.state.scalar_fields['Iterated'+name].val[:]
                self.fields['old_'+field]=self.state.scalar_fields['Old'+name].val[:]
                self.oldfields[field]=self.state.scalar_fields['Old'+name].val[:]
            except KeyError:
                ewrite(3,'KeyError! %s'%field)
                self.prescribed_fields.remove(field)
        for field,name in forced_fields.items():
            try:
                self.fields[field]=self.state.scalar_fields['Iterated'+name].val[:]
                self.oldfields[field]=self.state.scalar_fields['Old'+name].val[:]
                self.fields['old_'+field]=self.state.scalar_fields['Old'+name].val[:]
            except KeyError:
                ewrite(3,'KeyError! %s'%field)
                self.forced_fields.remove(field)
        for field,name in forced_fields.items():
            try:
                self.fields['delta'+field]=(
                    self.state.scalar_fields[name+'Source'].val)
                self.fields['delta'+field][:]=0.0
            except KeyError:
                try:
                    self.fields['delta'+field]=(
                    self.state.scalar_fields[name+'MicrophysicsSource'].val)
                    self.fields['delta'+field][:]=0.0
                except KeyError:
                    ewrite(3,'KeyError! %s'%field)
                    pass
            try:
                self.fields[field+'slip']=(
                    self.state.scalar_fields[name+'SinkingVelocity'].val)
                self.fields[field+'slip'][:]=0.0
            except KeyError:
                pass
            try:
                self.fields[field+'sponge']=(
                    self.state.scalar_fields[name+'Sponge'].val)
            except KeyError:
                pass
        self.fields['dt']=self.dt
        self.oldfields['dt']=self.dt


    def _add_python_fields(self,python_fields):
        for field_name,func in python_fields:
            self.fields[field_name]=func(**self.fields)
            self.fields['old_'+field_name]=func(**self.oldfields)
            self.oldfields[field_name]=func(**self.oldfields)

    def _calculateLimits(self):
        for F in self.limits:
            lim=F.value(**self.fields)
            dlim=lim-self.fields[F.arg][:]
            mask=F.logic(lim,self.fields[F.arg][:],**self.fields)
            lim=np.where(mask,lim,self.fields[F.arg][:])
            dlim=lim-self.fields[F.arg][:]
            self.fields[F.arg][:]=lim

#            print 'limiter:', lim.min(), lim.max()
#            print mask

            for fname,scale in F.items():
                try:
                    if fname[0]=='q':
                        self.fields[fname][:]+=scale*dlim
                    else:
                        self.fields['delta'+fname][:]+=scale*dlim/self.dt
                except KeyError:
                    pass

            lim=F.value(**self.oldfields)
            dlim=lim-self.oldfields[F.arg][:]
            mask=F.logic(lim,self.oldfields[F.arg][:],**self.oldfields)
            lim=np.where(mask,lim,self.oldfields[F.arg][:])
            dlim=lim-self.oldfields[F.arg][:]

 #           print 'limiter:', lim.min(), lim.max()
            
            self.oldfields[F.arg][:]=lim
            for fname,scale in F.items():
                try:
                    self.oldfields[fname][:]+=scale*dlim
                except KeyError:
                    pass
            

    def _calculateForcings(self):

        a=0.0

        alpha={}
        dF={}

        for F in self.forcings:
#            f=F.forcing(**self.fields)
            oldf=(F.alpha*F.forcing(**self.fields)
                  +(1.0-F.alpha)*F.forcing(**self.oldfields))
            ewrite(3,(F.name, np.min(oldf), np.max(oldf),'First'))
            for field_name,scale in F.items():
                if field_name in dF:
                    dF[field_name][:]=dF[field_name][:]+scale*oldf
                else:
                    dF[field_name]=scale*(oldf)

#        print 'DF', dF['q_v'].max(),dF['q_v'].min()

        for field in dF:
            try:
                fv=self.fields[field]
                dF[field]=np.where(abs(dF[field])<1.0e-8,fv/self.dt,dF[field]) 
                A=-fv/(self.dt*dF[field])
                A=np.where(np.isnan(A),1.0,A)
                alpha[field]=np.where((A>1.0)*(fv>0.0),1.0,A)
                alpha[field]=np.where((A<0.0)*(fv>0.0),1.0,alpha[field])
                alpha[field]=np.where((A>0.0)*(fv<0.0),1.0,alpha[field])
                alpha[field]=np.where((A<0.0)*(fv<0.0),
                                      np.minimum(-A,1.0e6),alpha[field])
            except KeyError:
                pass

        for F in self.forcings:
#            f=F.bound(F.forcing(**self.fields),alpha)
#            oldf=F.bound(F.forcing(**self.oldfields),alpha)
            oldf=(F.bound(F.alpha*F.forcing(**self.fields)
                  +(1.0-F.alpha)*F.forcing(**self.oldfields),alpha))
            ewrite(3,(F.name, np.min(oldf), np.max(oldf),'second'))
            for field_name,scale in F.items():
#                try:
                    if field_name[0]=='q':
                        q_w=self.oldfields['q_w'][:]
                        q_lw=self.oldfields['q_lw'][:]
                        q_v=abs(self.oldfields['q_v'][:])
#                        print F.name, field_name, oldf.max(), oldf.min(), 
                        if field_name=='q_v':
                            self.fields['delta'+field_name][:]=(self.fields['delta'+field_name][:]+fix(scale*(oldf),-q_v,q_lw))
                        else:
                            self.fields['delta'+field_name][:]=(self.fields['delta'+field_name][:]+fix(scale*(oldf),-q_lw,q_v))
                    elif field_name in('theta','e_i'):
                        try:
                            self.fields['delta'+field_name][:]=(self.fields['delta'+field_name][:]+scale*(oldf))
                        except:
                            pass
                    else:
                        self.fields['delta'+field_name][:]=(self.fields['delta'+field_name][:]+scale*(oldf))
#                except KeyError:
#                    pass

        for field_name in self.forced_fields:
                ewrite(3,(field_name,
                          np.min(self.fields['delta'+field_name][:]),
                   np.max(self.fields['delta'+field_name][:])))


    def _addSponges(self):
        for field_name in self.forced_fields:
            if (field_name+'sponge') in self.fields:
                self.fields['delta'+field_name][:]=(self.fields['delta'+field_name][:]+self.fields[field_name+'sponge'][:])


    def _calculateSlipVels(self):
        for S in self.slip_vels:
            try:
                self.fields[S.field+'slip'][:]=(0.5*S(self.fields)+
                                                0.5*S(self.oldfields))
            except KeyError:
                ewrite(3,'No slip velocity for %s!'%S.field)

L_v=2500000.0
#m=0.018
m=3.0e-3
rho_w=1020.0
R_v=461.5
D_v=2.0e-5
K_a=0.03
N_0r=10.0e7
a_r=201
nu=20.0
rho_0=1.0
k_c=1e-3
a_c=5e-4
cp=1000.0
g=9.81
cv=714.285714
R=cp-cv
cp_v=1840.0
Rv=461.5
cv_v= cp_v-Rv
cw=4100.0

def fix(x,a,b):
    return np.maximum(np.minimum(x,b),a)

def testing(states,dt,parameters):


    print "python Microphysics"

    if len(states)>1:
        state=states['Bulk']
        
        for a in ('WaterVapour','CloudWater','RainWater'):
            state.scalar_fields['Iterated'+a+'Fraction']=states[a].scalar_fields['IteratedMassFraction']
            if 'MassFractionMicrophysicsSource' in states[a].scalar_fields:
                state.scalar_fields[a+'FractionSource']=states[a].scalar_fields['MassFractionMicrophysicsSource']
            else:
                state.scalar_fields[a+'FractionSource']=states[a].scalar_fields['MassFractionSource']
            state.scalar_fields[a+'FractionSinkingVelocity']=states[a].scalar_fields['MassFractionSinkingVelocity']
            state.scalar_fields['Old'+a+'Fraction']=states[a].scalar_fields['OldMassFraction']
    else:
        state=states.items()[0]

    M=hot_moist_microphysics_model(state,dt,parameters)

    


        


        
def get_gas_density(rho=None,q_g=None,q_v=None,q_c=None,q_r=None,**kwargs):
    ewrite(3,'Gas Density %f %f %f %f %f'%(rho.min(),rho.max(),
                                           q_g.min(),q_c.max(),q_r.max()) )
    rho_g=rho*q_g/(1.0-rho*(q_c+q_r)/rho_w)
    return np.where(rho_g>0,rho_g,np.inf)

def get_T(q_v=None,q_c=None,q_r=None,e_i=None,theta=None,T=None,**kwargs):
    c_v=cv+(cv_v-cv)*q_v+(cw-cv)*(q_c+q_r)
    if not e_i==None:
        return e_i/c_v
    elif not T==None:
        return T
    else:
        return 280.0

def get_p(rho_g=None,T=None,q_v=None,q_g=None,**kwargs):
    R_g=(R*q_g+(Rv-R)*q_v)/q_g
    return rho_g*R_g*T

def get_q_g(q_v=None,q_c=None,q_r=None,**kwargs):
    return np.minimum(1.0-q_c-q_r,1.0)

def get_q_w(q_v=None,q_c=None,q_r=None,**kwargs):
    return q_v+q_c+q_r

def get_q_lw(q_v=None,q_c=None,q_r=None,**kwargs):
    return q_c+q_r

def p_sat(T):
    return 611.2*np.exp(17.62*(T-273.0)/(T-30.0))

def get_q_sat(rho=None,p=None,T=None,rho_g=None,q_v=None,
              q_c=None,q_r=None,**kwargs):
    q_g=get_q_g(q_v,q_c,q_r)
    ewrite(3,('T : ', np.min(T), np.max(T)))
    T=np.where(T<273.0,273.0,T)
    q_sat =p_sat(T)/R_v*(q_g/(rho_g*T))

    return q_sat
    
def av(v1,v2,alpha=1.0):
    return (alpha*v1+(1.0-alpha)*v2)


def negative_fix_logic(f,v,**kwargs):
    return v<0.0

def negative_fix(q_v=None,**kwargs):
    return 0.0*q_v

def supersaturation_fix_logic(r,v,**kwargs):
    return r<v

def supersaturation_fix(rho=None,p=None,T=None,rho_g=None,q_v=None,
              q_c=None,q_r=None,old_q_v=None,
                        old_rho=None,old_p=None,old_T=None,old_rho_g=None,
              old_q_c=None,old_q_r=None,dt=None,old=None,**kwargs):
    q_sat=get_q_sat(rho,p,T,rho_g,q_v,
              q_c,q_r,**kwargs)
    return q_sat
    if old==False:
#        q_sat=get_q_sat(rho=av(rho,old_rho),p=None,T=av(T,old_T),
#                        rho_g=av(rho_g,old_rho_g),q_v=av(q_v,old_q_v),
#              q_c=av(q_c,old_q_c),q_r=av(q_r,old_q_r),**kwargs)
        q_v=av(q_v,old_q_v)
        ewrite(3,('RH : ', np.min(q_v/q_sat), np.max(q_v/q_sat)))
        ConN=np.maximum((q_v-q_sat)/dt,0.0)
        ewrite(3,('ConN :', np.min(ConN), np.max(ConN)))
        return ConN
    else:
        return 0


def cloudwater_condensation(rho=None,p=None,T=None,
                            q_v=None,q_c=None,q_r=0.0,
                            q_sat=None,rho_g=None,**kwargs):

    from numpy import pi

    S_w=np.minimum(q_v/q_sat,1.0)
    q_g=1.0-q_c-q_r
    r_c=10.0e-6
    N_c=3.0*q_c/(4.0*rho_w*pi*r_c**3)
#    N_c=10000

    F_k=(L_v/(R_v*T)-1.0)*(L_v/(K_a*T))

    F_d=(R_v*T)/(p_sat(T)*D_v)

    ConC=4.0*np.pi*(q_g/rho_g)*N_c*r_c*(S_w-1.0)/(F_k+F_d)

    return np.minimum(ConC,0)

    
def rainwater_condensation(rho=None,p=None,T=None,q_r=None,q_c=None,
                           q_v=None,q_sat=None,rho_g=None,q_g=None,
                           dt=None,**kwargs):
    
    from scipy.special import gamma

    P_sat=611.2*np.exp(17.62*(T-273.0)/(T-30.0))
    S_w=np.minimum(q_v/q_sat,1.0)

    ilambda_r=(np.maximum((rho_g*q_r)/(np.pi*q_g*rho_w*N_0r),0))**0.25
    F_k=(L_v/(R_v*T)-1.0)*(L_v/(K_a*T))
    F_d=(R_v*T)/(P_sat*D_v)


    ewrite(3,'rain condensation in: %f %f'%(q_g.min(),rho_g.min()))
    ConR=(2.0*np.pi*q_g/rho_g*N_0r*(S_w-1.0)/(F_k+F_d)*
            (ilambda_r**2.0
            +0.22*gamma(2.75)*np.sqrt(a_r/nu)
                  *(rho_0/rho_g)**0.25*ilambda_r**(11.0/12.0)))
    ewrite(3,'rain condensation out %f %f'%(ilambda_r.max(),ConR.max()))

    return ConR

def water_droplet_autoconversion(q_c=None,q_r=None,rho_g=None,dt=None,**kwargs):

    q_g=1.0-q_c-q_r

    AutC=k_c*(q_c-q_g/rho_g*a_c)
    return np.maximum(AutC,0.0)


def warm_accreation(q_g=None,q_c=None,q_r=None,rho_g=None,rho_0=1.0,**kwargs):

    from scipy.special import gamma

    ewrite(3,'accreation in')
    ilambda_r=(np.maximum((rho_g*q_r)/(np.pi*q_g*rho_w*N_0r),0))**0.25
    Acc=np.pi/4.0*N_0r*a_r*np.sqrt(rho_0/rho_g)*gamma(3.5)*ilambda_r**3.5*q_c
    ewrite(3,'accreation out')

    return Acc

def w_r(field,rho_g=None,q_r=None,q_g=None,w_r=None,**kwargs):

    if w_r==None:
        ilambda_r=(np.maximum((rho_g*q_r)/(np.pi*q_g*rho_w*N_0r),0))**0.25
        r_r=ilambda_r/2.0
        ewrite(3,'r_max=%f, q_g=(%f,%f)'%(r_r.max(),q_g.max(),q_g.min()))
        field['w_r']=a_r*np.sqrt(r_r*rho_0/rho_g)
        return field['w_r']
    else:
        return w_r

def CloudScale(q_r=None,**kwargs):
    return -q_r/(1.0-q_r)

def VapourScale(q_r=None,**kwargs):
    return -q_r/(1.0-q_r)

class hot_moist_microphysics_model(microphysics_model):

    prescribed_fields={'T':"InsituTemperature",
                    'rho':"EOSDensity",
                    'p':"EOSPressure"}
    forced_fields={'q_v':"WaterVapourFraction",
                   'q_c':"CloudWaterFraction",
                   'q_r':"RainWaterFraction",
                   'e_i':"InternalEnergy",
                   'theta':"PotentialTemperature"}
    python_fields=[
        ('T',get_T),
        ('q_g',get_q_g),
        ('rho_g',get_gas_density),
        ('p',get_p),
        ('q_sat',get_q_sat),
        ('q_w',get_q_w),
        ('q_lw',get_q_lw)
        ]


    # First build your microphysical source terms:


    ConN=microphysics_limit(supersaturation_fix,'q_v',
                            supersaturation_fix_logic,
                            (('q_c',-1.0),
                             ('e_i',-L_v),
                             ('theta',-L_v/cp)),
                            name='ConN')
    ConC=microphysics_forcing(cloudwater_condensation,(('q_c',1.0),
                                                       ('q_v',-1.0),
                                                       ('e_i',L_v),
                                                       ('theta',L_v/cp)),
                               ('q_c','q_v'),alpha=0.0,name='ConC')
    ConR=microphysics_forcing(rainwater_condensation,(('q_r',1.0),
                                                      ('q_v',-1.0),
                                                      ('e_i',L_v),
                                                      ('theta',L_v/cp)),
                              ('q_r','q_v'),alpha=0.0,name='ConR')
    Aut=microphysics_forcing(water_droplet_autoconversion,(('q_r',1.0),
                                                           ('q_c',-1.0)),
                             ('q_r','q_c'),alpha=0.0,name='Aut')
    Acc=microphysics_forcing(warm_accreation,(('q_r',1.0),
                                              ('q_c',-1.0)),
                             ('q_r','q_c'),alpha=0.0,name='Acc')


    

    RainSlip=slip_velocity('q_r',w_r)
    CloudSlip=slip_velocity('q_c',w_r,CloudScale)
    VapourSlip=slip_velocity('q_v',w_r,VapourScale)

    qvlim=microphysics_limit(negative_fix,'q_v',
                            negative_fix_logic,
                             tuple(),
                            name='qvlim')

    qclim=microphysics_limit(negative_fix,'q_c',
                            negative_fix_logic,
                            tuple(),
                            name='qvlim')


    limits=(ConN,)
    sources=(ConC,ConR,Aut,Acc)



    slip_vels=(RainSlip,CloudSlip,VapourSlip)

    def __init__(self,state,dt,parameters):
        microphysics_model.__init__(self,state, 
                                  self.prescribed_fields,
                                  self.forced_fields,
                                  self.python_fields,
                                  self.limits,
                                  self.sources,
                                  self.slip_vels,
                                  dt
                                  )
