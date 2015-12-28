#!/usr/bin/env python
# 
# Copyright (c) 2001 Vivake Gupta (vivakeATomniscia.org).  All rights reserved.
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
#
# This software is maintained by Vivake (vivakeATomniscia.org) and is available at:
#     http://www.omniscia.org/~vivake/python/Simplex.py
 
 
""" Simplex - a regression method for arbitrary nonlinear function minimization
 
Simplex minimizes an arbitrary nonlinear function of N variables by the
Nedler-Mead Simplex method as described in:
 
Nedler, J.A. and Mead, R. "A Simplex Method for Function Minimization." 
    Computer Journal 7 (1965): 308-313.
 
It makes no assumptions about the smoothness of the function being minimized.
It converges to a local minimum which may or may not be the global minimum
depending on the initial guess used as a starting point.
"""
 
import math
import copy
 
class Simplex:
    def __init__(self, testfunc, guess, increments, kR = -1, kE = 2, kC = 0.5):
        """Initializes the simplex.
        INPUTS
        ------
        testfunc      the function to minimize
        guess[]       an list containing initial guesses
        increments[]  an list containing increments, perturbation size
        kR            reflection constant
        kE            expansion constant
        kC            contraction constant
        """
        self.testfunc = testfunc
        self.guess = guess
        self.increments = increments
        self.kR = kR
        self.kE = kE
        self.kC = kC
        self.numvars = len(self.guess)
        self.simplex = []
 
        self.lowest = -1
        self.highest = -1
        self.secondhighest = -1
 
        self.errors = []
        self.currenterror = 0
 
        # Initialize vertices
        for vertex in range(0, self.numvars + 3): # Two extras to store centroid and reflected point
            self.simplex.append(copy.copy(self.guess))
        # Use initial increments
        for vertex in range(0, self.numvars + 1):
            for x in range(0, self.numvars):
                if x == (vertex - 1):
                    self.simplex[vertex][x] = self.guess[x] + self.increments[x]
            self.errors.append(0)
        self.calculate_errors_at_vertices()
 
    def minimize(self, maxiters = 250):
        """Walks to the simplex down to a local minima.
        INPUTS
        ------
        epsilon       convergence requirement
        maxiters      maximum number of iterations
        monitor       if non-zero, progress info is output to stdout  
 
        OUTPUTS
        -------
        an array containing the final values
        lowest value of the error function
        number of iterations taken to get here
        """
        iter = 0
        for iter in range(0, maxiters):
            self.simplexStep()
        # Either converged or reached the maximum number of iterations.
        # Return the lowest vertex and the currenterror.
        for x in range(0, self.numvars):
            self.guess[x] = self.simplex[self.lowest][x]
        self.currenterror = self.errors[self.lowest]
        return self.guess, self.currenterror, iter
 
    def simplexStep(self, epsilon = 0.0001):
        # Identify highest, second highest, and lowest vertices
        self.highest = 0
        self.lowest = 0
        for vertex in range(0, self.numvars + 1):
            if self.errors[vertex] > self.errors[self.highest]:
                self.highest = vertex
            if self.errors[vertex] < self.errors[self.lowest]:
                self.lowest = vertex
        self.secondhighest = 0
        for vertex in range(0, self.numvars + 1):
            if vertex == self.highest:
                continue
            if self.errors[vertex] > self.errors[self.secondhighest]:
                self.secondhighest = vertex
        # Test for convergence
        S = 0.0
        S1 = 0.0
        for vertex in range(0, self.numvars + 1):
            S = S + self.errors[vertex]
        F2 = S / (self.numvars + 1)
        for vertex in range(0, self.numvars + 1):
            S1 = S1 + (self.errors[vertex] - F2)**2
        T = math.sqrt(S1 / self.numvars)
        # Optionally, print progress information
        #if monitor:
        #    print 'Iteration = %d   Best = %f   Worst = %f' % (iter,self.errors[self.lowest],self.errors[self.highest])
            
        if T <= epsilon:   # We converged!  Break out of loop!
            return True, T, self.errors[self.lowest],self.errors[self.highest]
        else:                   # Didn't converge.  Keep crunching.
            # Calculate centroid of simplex, excluding highest vertex
            for x in range(0, self.numvars):
                S = 0.0
                for vertex in range(0, self.numvars + 1):
                    if vertex == self.highest:
                        continue
                    S = S + self.simplex[vertex][x]
                self.simplex[self.numvars + 1][x] = S / self.numvars
 
            self.reflect_simplex()
 
            self.currenterror = self.testfunc(self.guess)
 
            if self.currenterror < self.errors[self.lowest]:
                tmp = self.currenterror
                self.expand_simplex()
                self.currenterror = self.testfunc(self.guess)
                if self.currenterror < tmp:
                    self.accept_expanded_point()
                else:
                    self.accept_reflected_point()
 
            elif self.currenterror <= self.errors[self.secondhighest]:
                self.accept_reflected_point()
 
            elif self.currenterror <= self.errors[self.highest]:
                self.accept_reflected_point()
 
                self.contract_simplex()
                self.currenterror = self.testfunc(self.guess)
                if self.currenterror < self.errors[self.highest]:
                    self.accept_contracted_point()
                else:
                    self.multiple_contract_simplex()
 
            elif self.currenterror >= self.errors[self.highest]:
                self.contract_simplex()
                self.currenterror = self.testfunc(self.guess)
                if self.currenterror < self.errors[self.highest]:
                    self.accept_contracted_point()
                else:
                    self.multiple_contract_simplex()
        return False, T, self.errors[self.lowest],self.errors[self.highest]
 
 
    def contract_simplex(self):
        for x in range(0, self.numvars):
            self.guess[x] = self.kC * self.simplex[self.highest][x] + (1 - self.kC) * self.simplex[self.numvars + 1][x]
        return
 
    def expand_simplex(self):
        for x in range(0, self.numvars):
            self.guess[x] = self.kE * self.guess[x]                 + (1 - self.kE) * self.simplex[self.numvars + 1][x]
        return
 
    def reflect_simplex(self):
        for x in range(0, self.numvars):
            self.guess[x] = self.kR * self.simplex[self.highest][x] + (1 - self.kR) * self.simplex[self.numvars + 1][x]
            self.simplex[self.numvars + 2][x] = self.guess[x] # REMEMBER THE REFLECTED POINT
        return
 
    def multiple_contract_simplex(self):
        for vertex in range(0, self.numvars + 1):
            if vertex == self.lowest:
                continue
            for x in range(0, self.numvars):
                self.simplex[vertex][x] = 0.5 * (self.simplex[vertex][x] + self.simplex[self.lowest][x])
        self.calculate_errors_at_vertices()
        return
 
    def accept_contracted_point(self):
        self.errors[self.highest] = self.currenterror
        for x in range(0, self.numvars):
            self.simplex[self.highest][x] = self.guess[x]
        return
 
    def accept_expanded_point(self):
        self.errors[self.highest] = self.currenterror
        for x in range(0, self.numvars):
            self.simplex[self.highest][x] = self.guess[x]
        return
 
    def accept_reflected_point(self):
        self.errors[self.highest] = self.currenterror
        for x in range(0, self.numvars):
            self.simplex[self.highest][x] = self.simplex[self.numvars + 2][x]
        return
 
    def calculate_errors_at_vertices(self):
        for vertex in range(0, self.numvars + 1):
            if vertex == self.lowest:
                continue
            for x in range(0, self.numvars):
                self.guess[x] = self.simplex[vertex][x]
            self.currenterror = self.testfunc(self.guess)
            self.errors[vertex] = self.currenterror
        return
 
def objective_function(args):
    #return abs(args[0] * args[0] * args[0] * 5 - args[1] * args[1] * 7 + math.sqrt(abs(args[0])) - 118)
    return abs(args[0] **3 * 5 - args[1]**2 * 7 + math.sqrt(abs(args[2])) - 118)
 
def main():
    s = Simplex(objective_function, [1, 1, 1], [2, 4, 6])
    values, err, iter = s.minimize()
    print 'args = ', values
    print 'error = ', err
    print 'iterations = ', iter
 
if __name__ == '__main__':
    main()
