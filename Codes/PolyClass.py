'''
Created on Apr 22, 2014

@author: Chandrasekar
'''

from casadi import *
class Atom:
    # Constructor: allows Atom to be created from the coefficient, coefficient & name or coefficient, name & index.
    def __init__(self, c = 0, v = ssym("@"), i=0):
        v1 = ssym("@")
        self.c = c # coefficient
        self.v = v  # variable name
        if casadi.isEqual(self.v,v1):
            self.i = 0  # exponent
        else:
            self.i = 1 # exponent
            
    def setCoeff(self,c=0):
        self.c = c
#         self = Atom(c)
        return self
        
    def setVar(self,v = ssym("@")):
        self.v = v
#         self = Atom(,v,)
        return self
            
    def setExp(self,i=0):
        self.i = i
#         self = Atom(,,i)
        return self
        
    def isEqual(self,a):
        if (self.i == a.i and self.c == a.c and casadi.isEqual(self.v,a.v)):
            return 1
        else:
            return 0
        
    def multInv(self):
        try:
            self.setCoeff(1/self.c)
            if self.i == 0:
                self.setExp(0)
            else:
                self.setExp(1/self.i)
        except:
            print "Multiplicative Inverse cannot be computed"
        
        return self
    
    def addInv(self):
        if self.c != 0:
            self.setCoeff(-self.c)
        return self
    
    def printAtom(self):
#         print str(self.c) + " " + self.v + "^" + str(self.i)
        print self.c*pow(self.v,self.i)
        
    def evalAtom(self,val):
        var = casadi.getFree(self.v)
        return self.c*pow(casadi.substitute(self.v,var,val),self.i)
            
                                
        
if __name__ == '__main__':
    x = Atom()
    y0 = Atom()
    y1 = Atom()
    y2 = Atom()
    y3 = Atom()
    
    y = ssym("y")
    y0 = Atom(1,y,)
    y0.printAtom()
    
    x = ssym("x")
    y1.setVar(x)
    y1.setCoeff(3.0)
    y1.printAtom()

    y2.setVar(x)
    y2.setCoeff(4.0)
    y2.printAtom()
    
    y3.setExp(3.0)
    y3.printAtom()

    print y0.isEqual(y0), y0.isEqual(y1), y0.isEqual(y2), y0.isEqual(y3) 
    
    y2.setCoeff(4.0)
    y2.printAtom()
    
    y2.multInv()
    y2.printAtom()
    
    y2.addInv()
    y2.printAtom()
    print y2.evalAtom(4)