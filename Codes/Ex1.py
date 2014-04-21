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
import math

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
    par_dir = os.path.split(curr_dir)[0]
    
    # Define model directory
    mo_dir = par_dir+"/Models/Ex1"
    
    # Define model
    mo_name = mo_dir+"/OPT_Ctl_DT.mo"
    
    # Define optimization file
    mop_file =  mo_dir+"/OPT.mop"
    
    # Define optimization problem
    mop_name = "OPT.Ex1_Opt"
    
    # Algorithm parameters
    degree = 2
    xe = [0,1,2,3]

    # Call optimization algorithm & solve
    eq = runDT(mo_dir, mo_name, mop_file, mop_name, degree)
    print eq
#     print casadi.getFree(eq[0])
    solveEqns(eq,xe)

def runDT(mo_dir, mo_name, mop_file, mop_name, degree):    
#     casadi_name = compile_fmux("OPT.Ex1_Opt",(mo_dir+"/OPT.mop",mo_dir+"/OPT_Ctl_DT.mo"))
    casadi_name = compile_fmux(mop_name,(mop_file,mo_name))  
    tmpdir = unzip_unit(casadi_name)
    fmux_files = get_files_in_archive(tmpdir)
    xmlFile = fmux_files['model_desc']
#    #xmldoc = xmlparser.ModelDescription(xmlFile)
#     
    ocp = casadi.SymbolicOCP()
    options = {}
    options["verbose"] = True
    options["sort_equations"] = False
    options["eliminate_dependent"] = False
    ocp.parseFMI(xmlFile,options)
# #    
# #    
#     print ocp.ode
#     print ocp.u
#     print ocp.x
#     print ocp.y
    
    # Get states
    x = SXMatrix(casadi.var(ocp.x[1:ocp.x.size()]))
#     print x
    
    u = ocp.u
    
    # Set expansion points as a varible
    xe = ssym("xe",x.size())
    
    # Get size of input and state
    num_inp = ocp.u.size()
    num_state = ocp.x.size() - 1 # Exclude cost from the state
         
    # Get RHS of ode
    ode_rhs = SXMatrix(casadi.der(ocp.x)) - ocp.ode
    simplify(ode_rhs)
    cost = ode_rhs[0]
    f = ode_rhs[1:num_state + 1]
    
    F = DT(f,x,degree,xe)
    L = DT(cost,x,degree,xe)
    
#   Derivative of cost w.r.t inputs  
#     derCost_derU = casadi.jacobian(cost,SXMatrix(casadi.var(ocp.u)))   
    derCost_derU = casadi.jacobian(L,SXMatrix(casadi.var(ocp.u)))
    simplify(derCost_derU)
#     print derCost_derU

#   Derivative of system dynamics w.r.t inputs     
#     dfdu = casadi.jacobian(f,SXMatrix(casadi.var(ocp.u)))   
    dfdu = casadi.jacobian(F,SXMatrix(casadi.var(ocp.u)))
    simplify(dfdu)
#     print dfdu

#   Define coefficients of value function as: V(x) = Cv where Cv is an array.    
    Cv0 = ssym('Cv0') # constant term
    Cv = ssym('Cv',num_state,degree) # single input, multiple states
    
    # Define value function expression
    V = Cv0
    for i in range(0,degree):
        V = V + Cv[:,i]*pow(x-xe,i+1)
    
    simplify(V)
#     print V
    
    dVdx = casadi.jacobian(V,x)
    simplify(dVdx)
#     print dVdx
    
    # Define the equations for solving
    eq1 = derCost_derU + dVdx*dfdu
    eq2 = L + dVdx*F
    
    simplify(eq1)
    simplify(eq2)
    
#     print eq1
#     print eq2
    
    return [eq1,eq2,V]
#     print casadi.getFree(eq1)
#     print casadi.substitute(f,x,xe)
    
#     x = SXMatrix.ones(ocp.x.size()-1,degree+1)
#     x[:,0] = SXMatrix.ones(ocp.x.size()-1)
#     x[:,1] = SXMatrix(casadi.var(ocp.x[1:ocp.x.size()]))
#     for i in range(2,degree+1):
#         x[:,i] = pow(x[:,1],i)
#     print x
# 
#     u_x = mul(Cu,x.T)
#     print u_x
#     dLdu = casadi.substitute(dLdu,SXMatrix(casadi.var(ocp.u)),u_x)
#     print dLdu
#     
#     dU = DT_coeffs(dLdu,SXMatrix(casadi.var(ocp.x[1:ocp.x.size()])))
#     print dU
#     
#     V = ssym('V',ocp.x.size()-1,degree+1)
#     V_x = mul(V,x.T)
#     dV_dx = jacobian(V_x,x[:,1])
#     print dV_dx
#     dVdx = DT_coeffs(dV_dx,x[:,1])
#     print dVdx*dfdu
# #    dVdx = ssym('dVdx',ocp.x.size()-1,degree+1)
# #    dVdx[ocp.x.size()] = SX(0)
# 
#     
#     expr1 = dU + dVdx.T*dfdu #'''dHdu = 0 is the first algebraic difference equation'''
# #    print expr1
#     
#     Lag = substitute(ode_rhs[0],SXMatrix(casadi.var(ocp.u)),u_x)
#     L = DT_coeffs(Lag,SXMatrix(casadi.var(ocp.x[1:ocp.x.size()])))
# #    print L
#     
#     fun = substitute(ode_rhs[1:ocp.ode.size()],SXMatrix(casadi.var(ocp.u)),u_x)
#     f = DT_coeffs(fun,SXMatrix(casadi.var(ocp.x[1:ocp.x.size()])))
# #    print f
# 
#     expr2 = L + dVdx.T*f
# #    print expr2
# #    
#     global dVdx, Cu, ocp
#     return expr1, expr2
#     
#     
# def F(x):
#     global expr11,expr22,dVdx,Cu
#     x1 = x[0]
#     z1 = casadi.substitute(expr11,Cu,SXMatrix([x1[0:3]]))
#     y1 = casadi.substitute(z1,dVdx,SXMatrix([x1[3:6]]))
#     
#     z2 = casadi.substitute(expr22,Cu,SXMatrix([x1[0:3]]))
#     y2 = casadi.substitute(z2,dVdx,SXMatrix([x1[3:6]]))
#     
#     return [double(y1),double(y2)]
#     
#     
# def DT_coeffs(f,x):
# # 1D DT of a scalar function f about the scalar variable x
#     global degree, H
#     X = ssym('X',degree+1);
#     f1 = f
#     X[0] = f1#casadi.jacobian(f1,pow(x,0))
#     for i in range(1, degree+1):
#         f1 = casadi.jacobian(f1,x)
#         X[i] = pow(H,i)/math.factorial(i)*f1
#         
#     return X 
#      

def DT(f,x,degree,xe):
        
    vars = casadi.getFree(f)
    out = casadi.substitute(f,x,xe)
    coeffs = casadi.jacobian(f,x);

    for i in range(0,degree):
        out = out + casadi.substitute(coeffs,x,xe)*pow(x-xe,i+1)/math.factorial(i+1)
        coeffs = casadi.jacobian(coeffs,x);
    
    simplify(out)
    return out

def solveEqns(f,xe1):
    BC = 0
    var = casadi.getFree(f[0])
    for i in range(1,2): #len(xe):
        expr = casadi.substitute(f[0],var[3],xe1[i-1])
        simplify(expr)
#         x = scipy.optimize.root([expr[0],expr[1],expr[2]-BC],[[1,1,1,1,1,1]], method='krylov', options={'disp': True})
        print var
        print f[0]
        print expr
    return 0
#     for i = 
   
if __name__ == '__main__':
#     global degree, H, z, dVdx, expr11, expr22, ocp
#     degree = 2
#     H = 1.0
#     expr1, expr2 = run_Ex1()
# #    expr11 = casadi.substitute(expr1,SXMatrix(casadi.var(ocp.x[1:ocp.x.size()])),SXMatrix.zeros(ocp.x.size()-1,1))
# #    expr22 = casadi.substitute(expr2,SXMatrix(casadi.var(ocp.x[1:ocp.x.size()])),SXMatrix.zeros(ocp.x.size()-1,1))
# ##    print F([[1,1,1,1,1,1]])
# #    
# #    print expr11
# #    print expr22
# #    
# #    
# ##    z = run_Ex1()
# #    x = scipy.optimize.root(F,[[1,1,1,1,1,1]], method='krylov', options={'disp': True})
# ##    x = scipy.optimize.broyden1(F, [[1,1,1,1,1,1]], f_tol=1e-14)
    run_Ex1()
    print 0
