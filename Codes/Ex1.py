'''
Created on Oct 22, 2013

@author: Chandrasekar
'''
import os
import sys
import numpy as NP
from numpy import *
import matplotlib.pyplot as plt
import zipfile
import time
import shutil
import scipy.optimize

try:
  # JModelica
  from pyfmi.common.core import unzip_unit, get_files_in_archive
  from pymodelica import compile_fmux
  #from pyjmi import CasadiModel
  #from pyjmi.common import xmlparser
  #from pyjmi.common.xmlparser import XMLException  
  from pymodelica import compile_jmu
  from pyjmi import JMUModel
  import pymodelica
  use_precompiled = False
except:
  print "No jmodelica installation, falling back to precompiled XML-files"
  use_precompiled = True

# CasADi
from casadi import *
import os
   
def run_Ex1(with_plots=True):

#    os.chdir("C:\JModelica.org-SDK-1.8.1\install\Python\pyjmi\examples")
    curr_dir = os.path.dirname(os.path.abspath(__file__));
    casadi_name = compile_fmux("OPT.Ex1_Opt",(curr_dir+"/OPT.mop",curr_dir+"/OPT_Ctl_DT.mo"))
    tmpdir = unzip_unit(casadi_name)
    fmux_files = get_files_in_archive(tmpdir)
    xmlFile = fmux_files['model_desc']
#    #xmldoc = xmlparser.ModelDescription(xmlFile)
    
    ocp = casadi.SymbolicOCP()
    options = {}
    options["verbose"] = True
    options["sort_equations"] = False
    options["eliminate_dependent"] = False
    ocp.parseFMI(xmlFile,options)
#    
#    
    print ocp.ode
    
    # Get RHS of ode
    ode_rhs = SXMatrix(casadi.der(ocp.x)) - ocp.ode
    simplify(ode_rhs)
    print ode_rhs
    
    dLdu = casadi.jacobian(ode_rhs[0],SXMatrix(casadi.var(ocp.u)))
    simplify(dLdu)
    print dLdu
    
    dfdu = casadi.jacobian(ode_rhs[1:ocp.ode.size()],SXMatrix(casadi.var(ocp.u)))
    simplify(dfdu)
    print dfdu
    
    Cu = ssym('Cu',ocp.x.size()-1,degree+1) # single input
    x = SXMatrix.ones(ocp.x.size()-1,degree+1)
    x[:,0] = SXMatrix.ones(ocp.x.size()-1)
    x[:,1] = SXMatrix(casadi.var(ocp.x[1:ocp.x.size()]))
    for i in range(2,degree+1):
        x[:,i] = pow(x[:,1],i)
    print x

    u_x = mul(Cu,x.T)
    print u_x
    dLdu = casadi.substitute(dLdu,SXMatrix(casadi.var(ocp.u)),u_x)
    print dLdu
    
    dU = DT_coeffs(dLdu,SXMatrix(casadi.var(ocp.x[1:ocp.x.size()])))
    print dU
    
    V = ssym('V',ocp.x.size()-1,degree+1)
    V_x = mul(V,x.T)
    dV_dx = jacobian(V_x,x[:,1])
    print dV_dx
    dVdx = DT_coeffs(dV_dx,x[:,1])
    print dVdx*dfdu
#    dVdx = ssym('dVdx',ocp.x.size()-1,degree+1)
#    dVdx[ocp.x.size()] = SX(0)

    
    expr1 = dU + dVdx.T*dfdu #'''dHdu = 0 is the first algebraic difference equation'''
#    print expr1
    
    Lag = substitute(ode_rhs[0],SXMatrix(casadi.var(ocp.u)),u_x)
    L = DT_coeffs(Lag,SXMatrix(casadi.var(ocp.x[1:ocp.x.size()])))
#    print L
    
    fun = substitute(ode_rhs[1:ocp.ode.size()],SXMatrix(casadi.var(ocp.u)),u_x)
    f = DT_coeffs(fun,SXMatrix(casadi.var(ocp.x[1:ocp.x.size()])))
#    print f

    expr2 = L + dVdx.T*f
#    print expr2
#    
    global dVdx, Cu, ocp
    return expr1, expr2
    
    
def F(x):
    global expr11,expr22,dVdx,Cu
    x1 = x[0]
    z1 = casadi.substitute(expr11,Cu,SXMatrix([x1[0:3]]))
    y1 = casadi.substitute(z1,dVdx,SXMatrix([x1[3:6]]))
    
    z2 = casadi.substitute(expr22,Cu,SXMatrix([x1[0:3]]))
    y2 = casadi.substitute(z2,dVdx,SXMatrix([x1[3:6]]))
    
    return [double(y1),double(y2)]
    
    
def DT_coeffs(f,x):
# 1D DT of a scalar function f about the scalar variable x
    global degree, H
    X = ssym('X',degree+1);
    f1 = f
    X[0] = f1#casadi.jacobian(f1,pow(x,0))
    for i in range(1, degree+1):
        f1 = casadi.jacobian(f1,x)
        X[i] = pow(H,i)/math.factorial(i)*f1
        
    return X 
        
if __name__ == '__main__':
    global degree, H, z, dVdx, expr11, expr22, ocp
    degree = 2
    H = 1.0
    expr1, expr2 = run_Ex1()
#    expr11 = casadi.substitute(expr1,SXMatrix(casadi.var(ocp.x[1:ocp.x.size()])),SXMatrix.zeros(ocp.x.size()-1,1))
#    expr22 = casadi.substitute(expr2,SXMatrix(casadi.var(ocp.x[1:ocp.x.size()])),SXMatrix.zeros(ocp.x.size()-1,1))
##    print F([[1,1,1,1,1,1]])
#    
#    print expr11
#    print expr22
#    
#    
##    z = run_Ex1()
#    x = scipy.optimize.root(F,[[1,1,1,1,1,1]], method='krylov', options={'disp': True})
##    x = scipy.optimize.broyden1(F, [[1,1,1,1,1,1]], f_tol=1e-14)
    print x
