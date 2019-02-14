import numpy as np

class Equation:
    """
    Python class to conveniently represent Reluplex Equation
    """
    def __init__(self):
        """
        Construct empty equation
        """
        self.addendList = []
        self.participatingVariables = set()
        self.scalar = None
        self.auxVar = None

    def setScalar(self, x):
        """
        Set scalar of equation
        Arguments:
            x: (float) scalar RHS of equation
        """
        self.scalar = x
        
    

    def addAddend(self, c, x):
        """
        Add addend to equation
        Arguments:
            c: (float) coefficient of addend
            x: (int) variable number of variable in addend
        """
        self.addendList += [(c, x)]
        self.participatingVariables.update([x])

    def markAuxiliaryVariable(self, aux):
        """
        Mark variable as auxiliary for this equation
        Arguments:
            aux: (int) variable number of variable to mark
        """
        self.auxVar = aux

    def getParticipatingVariables(self):
        """
        Returns set of variables participating in this equation
        """
        return self.participatingVariables

    def participatingVariable(self, var):
        """
        Check if the variable participates in this equation
        Arguments:
            var: (int) variable number to check
        """
        return var in self.getParticipatingVariables()

    def replaceVariable(self, x, xprime, c):
        """
        Replace x with xprime + c
        Arguments:
            x: (int) old variable to be replaced in this equation
            xprime: (int) new variable to be added, does not participate in this equation
            c: (float) difference between old and new variable
        """
        assert self.participatingVariable(x)
        assert not self.participatingVariable(xprime)
        assert self.auxVar != x and self.auxVar != xprime
        for i in range(len(self.addendList)):
            if self.addendList[i][1] == x:
                coeff = self.addendList[i][0]
                self.addendList[i] = (coeff, xprime)
                self.setScalar(self.scalar + coeff*c)
                self.participatingVariables.remove(x)
                self.participatingVariables.update([xprime])

def addEquality(network, vars, coeffs, scalar):
    """
    Function to conveniently add equality constraint to network
    \sum_i vars_i*coeffs_i - scalar = 0
    Arguments:
        network: (MarabouNetwork) to which to add constraint
        vars: (list) of variable numbers
        coeffs: (list) of coefficients
        scalar: (float) representing RHS of equation
    """
    assert len(vars)==len(coeffs)
    e = Equation()
    aux = network.getNewVariable()
    network.setLowerBound(aux, 0.0)
    network.setUpperBound(aux, 0.0)
    e.markAuxiliaryVariable(aux)
    for i in range(len(vars)):
        e.addAddend(coeffs[i], vars[i])
    e.setScalar(-scalar)
    network.addEquation(e)

def addInequality(network, vars, coeffs, scalar):
    """
    Function to conveniently add inequality constraint to network
    \sum_i vars_i*coeffs_i -scalar <= 0
    Arguments:
        network: (MarabouNetwork) to which to add constraint
        vars: (list) of variable numbers
        coeffs: (list) of coefficients
        scalar: (float) representing RHS of equation
    """
    assert len(vars)==len(coeffs)
    e = Equation()
    aux = network.getNewVariable()
    network.setUpperBound(aux, 0.0)
    e.markAuxiliaryVariable(aux)
    for i in range(len(vars)):
        e.addAddend(coeffs[i], vars[i])
    e.setScalar(-scalar)
    network.addEquation(e)

def addComplementOutputSet(network, LB, UB, x):
    """
    Function to convert an output specification of staying within a set defined 
    by a lower bound and upper bound to its complement appropriate for Marabou.
    Arguments:
        network: (MarabouNetwork) to which to add constraint
        LB: (float) specifying the lower bound
        UB: (float) specifying the upper bound
        x: (int) specifying the variable
    """
    def add_aux_var(network, equation):
        aux = network.getNewVariable()
        network.setLowerBound(aux, 0.0)
        network.setUpperBound(aux, 0.0)
        equation.markAuxiliaryVariable(aux)
    # define x_l = l - x
    x_l = network.getNewVariable()
    eq = Equation()
    add_aux_var(network, eq)
    eq.addAddend(1.0, x_l)
    eq.addAddend(1.0, x)
    eq.setScalar(LB)
    network.addEquation(eq)
    # define x_u = x - u
    x_u = network.getNewVariable()
    eq1 = Equation()
    add_aux_var(network, eq1)
    eq1.addAddend(1.0, x_u)
    eq1.addAddend(-1.0, x)
    eq1.setScalar(-UB)
    network.addEquation(eq1)
    #
    # For a validity interface we would want both x_l and x_u to be negative, but
    # for Marabou's satisfiability interface, we assert that at least on of the
    # constraints must be in violation, meaning one on them is greater than zero
    # max(x_l, x_u) > 0
    Y = network.getNewVariable()
    network.addMaxConstraint({x_l,x_u}, Y)
    network.setLowerBound(Y, 0.0)
    if network.outputVars is None:
        network.outputVars = np.array([[Y]])
    else: 
        network.outputVars = np.vstack([network.outputVars, np.array([[Y]]) ])
