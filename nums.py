#Numerical methods: Find Numerical solution to given function e^x
#to n signifigant digits

import os, sys
import math
import time

def exponential():
    #Use taylor series to converge on a value within an acceptle estimated error for an expoential function 
    sigFigs = float(input("Enter number of signifigant figures: "))
    xVal = float(input("Enter x-value: "))
    e_s = 0.5 * 10**(2 - sigFigs)
    e_a = 100
    denom = 1
    ans = 1
    while e_a > e_s:
        coeff = (xVal**(denom))/(math.factorial(denom)) 
        newAns = ans + coeff
        e_a = ((newAns - ans)/newAns)*100
        denom += 1
        
        ans = newAns
        print("The ans: " + str(ans))
        print(e_a)
    print("The final value is " + str(ans))
    print("The final e_a is " + str(e_a))
    
#exponential()


#-----------------------------------------#
# Section 5: Finding roots of polynomials #
#-----------------------------------------#

# Functions use "bracketing" methods, which chosses a bound above and below the root and systematically
# "narrows" the gap between the bounds until the root is found. Method is said to be convergent since
# as the method progresses, it will converge on a root.

def bolzano():
    #Numerically find roots for polynomal equation, "Brute Force" method 
    x_u = float(input("Enter upper bound x val: "))
    x_l = float(input("enter lower bound x val: "))
    sigFigs = float(input("Enter number of signifigant figures: "))
    t0 = time.clock()
    func_u = 2*x_u**2 + 7*x_u + 2
    func_l = 2*x_l**2 + 7*x_l + 2
    check = func_u * func_l
    e_s = 0.5 * 10**(2 - sigFigs)
    e_a = 100
    if check > 0:
        while check > 0:
            x_l -= 1
            func_l = 2*x_l**2 + 7*x_l + 2
            check = func_u * func_l
    while e_a > e_s:
        newX_r = (x_l + x_u)/2
        func_r = 2*newX_r**2 + 7*newX_r + 2
        check_l = func_l * func_r
        e_a = abs((x_u - x_l)/(x_u + x_l))*100
        if check_l < 0:
            x_u = newX_r
        elif check_l > 0:
            x_l = newX_r
    deltat = time.clock() - t0
    print("The final value is " + str(newX_r))
    print("The final e_a is " + str(e_a))
    print("The bolzano method took " + str(deltat) + " seconds") 
        
            

#bolzano()

def falsePosition():
    #Numerically find roots using the false-postition formula, or linear interpolation
    x_u = float(input("Enter upper bound x val: "))
    x_l = float(input("enter lower bound x val: "))
    sigFigs = float(input("Enter number of signifigant figures: "))
    t0 = time.clock()
    e_s = 0.5 * 10**(2 - sigFigs)
    func_u = 2*x_u**2 + 7*x_u + 2
    func_l = 2*x_l**2 + 7*x_l + 2
    check = func_u * func_l
    e_a = 100
    if check > 0:
        while check > 0:
            x_l -= 1
            func_l = 2*x_l**2 + 7*x_l + 2
            check = func_u * func_l
    while e_a > e_s:
        x_r = x_u - (func_u * (x_l - x_u))/(func_l - func_u)
        func_r = 2*x_r**2 + 7*x_r + 2
        check_l = func_l * func_r
        e_a = abs((x_u - x_l)/(x_u + x_l))*100
        if check_l < 0:
            x_u = x_r
        elif check_l > 0:
            x_l = x_r
        
    deltat = time.clock() - t0 
    print("The final value is " + str(x_r))
    print("The final e_a is " + str(e_a))
    print("The false position method took " + str(deltat) + " seconds") 
    
#falsePosition()

#-------------------------#
# Chapter 6: Open Methods #
#-------------------------#

# open methods begin with a single number or a numbers that dont neccesarily bracket the root, and can
# diverge away from the root.


    
def fixedPoint():
    # Transform equation so that f(x) = 0 is in the form x = g(x) (i.e. algebraically re-arrange eqn. to get x = somthing with x^n)
    
    x_o = float(input("Enter inital x-value guess (suggest x = 0): "))
    sigFigs = float(input("Enter number of signifigant figures: "))
    t0 = time.clock()
    e_s = 0.5 * 10**(2 - sigFigs)
    e_a = 100
    while e_a > e_s:
        xPlus = math.exp(-x_o)
        e_a = abs((xPlus - x_o)/xPlus)*100
        x_o = xPlus
    deltat = time.clock() - t0 
    print("The final x value is " + str(xPlus))
    print("The final e_a value is " + str(e_a))
    print("The fixed point method found a root in " + str(deltat) + " seconds.") 
        
#fixedPoint()

def newtonRaphson():
    #Using known derivatives, the Newton-Raphson method will converge on a root much faster on the fixed point method. Timing indicates this converges
    #about twice as fast as the fixed point method

    # Problems with this method is that it can (as in the case of 10^10 - 1) converge to a root very very slowly. Also cases where multiple roots exist
    # it will also preform slowly. 
    
    sigFigs = float(input("Enter number of signifigant figures: "))
    t0 = time.clock()
    x_o = 0
    e_s = 0.5 * 10**(2 - sigFigs)
    e_a = 100
    while e_a > e_s:
        xPlus = x_o - (math.exp(-x_o) - x_o)/(-math.exp(x_o) - 1)
        e_a = abs((xPlus - x_o)/xPlus)*100
        x_o = xPlus
    deltat = time.clock() - t0
    check = math.exp(-xPlus) - xPlus
    print(str(check))
    print("The final x value is " + str(xPlus))
    print("The final e_a value is " + str(e_a))
    print("The fixed point method found a root in " + str(deltat) + " seconds.") 
    
#newtonRaphson()

def secant():
    #Sometimes the derivative is difficult to compute. The secan tmethod eliminates the need to find the derivaitve and uses the finite difference approach
    #to calculate the derivative

    # Important to note that since this is an incremental method of finding the root. This means that the root can lie below the first two guesses, and will
    # Will diverge.

    # When the secant method converges, it will do so at a faster rate than the false position method 
    sigFigs = float(input("Enter number of signifigant figures: "))
    t0 = time.clock()
    x_1 = 0
    x_i = 1
    e_s = 0.5 * 10**(2 - sigFigs)
    e_a = 100
    while e_a > e_s:
        func_i = math.exp(-x_i) - x_i
        func_i1 =(x_1 - x_i) 
        func_1 = math.exp(-x_1) - x_1
        xPlus = x_i - (func_i * func_i1)/(func_1 - func_i)
        e_a = abs((xPlus - x_i)/xPlus)*100
        x_1 = x_i
        x_i = xPlus
    deltat = time.clock() - t0
    check = math.exp(-xPlus) - xPlus
    print(str(check))
    print("The final x value is " + str(xPlus))
    print("The final e_a value is " + str(e_a))
    print("The fixed point method found a root in " + str(deltat) + " seconds.")    

#secant()

def modNewtonRaphson():
    # The Newton Raphson method is modified to include the possibilty of there being mulitple roots for a particular function.

    # Funtion evaluated is f(x) = x^3 - 5x^2 + 7x - 3 (this will have two roots)

    # The modification is made by using a ratio of the function and its derivative and taking the derivative of this ratio.

    # x_i+1 = x_i - f(x_i)f'(x_i)/[f'(x_i)]^2 - f(x_i)f"(x_i) 
    x_i = float(input("Enter inital x-value guess (suggest x = 0): "))
    sigFigs = float(input("Enter number of signifigant figures: "))
    t0 = time.clock()
    e_s = 0.5 * 10**(2 - sigFigs)
    e_a = 100
    while e_a > e_s:
        xPlus = x_i - ((x_i**3 - 5*x_i**2 + 7*x_i - 3)*(3*x_i**2 -10*x_i + 7))/((3*x_i**2 - 10*x_i + 7)**2 - (x_i**3 - 5*x_i**2 + 7*x_i - 3)*(6*x_i - 10))
        e_a = abs((xPlus - x_i)/xPlus)*100
        x_i = xPlus
    deltat = time.clock() - t0
    check = x_i**3 - 5*x_i**2 + 7*x_i - 3
    print(str(check))
    print("The final x value is " + str(xPlus))
    print("The final e_a value is " + str(e_a))
    print("The fixed point method found a root in " + str(deltat) + " seconds.") 
    

    

#modNewtonRaphson()

def systemsNR():
    #Systems of equations for the Newton Rahpson method similar method, but uses partial deriviatives with respect to x and y
    #instead of using the derivative of the eqautoin with respect to x.

    #Example equations will be u(x, y) = x^2 + xy - 10 = 0
    # and                      v(x, y) = y + 3xy^2 - 57 = 0

    #The iterative equation becomes x_i+1 = x_i - (u_i(dvi/dy) - vi(dui/dy))/((dui/dx)(dvi/dy) - (dui/dy)(dvi/dx))
    # and                           y_i+1 = y_i - (v_i(dui/dx) - ui(dvi/dx))/((dui/dx)(dvi/dy) - (dui/dy)(dvi/dx)
    
    x_i = float(input("Enter initial x-value guess: "))
    y_i = float(input("Enter initial y-value guess: ")) 
    sigFigs = float(input("Enter number of signifigant figures: "))
    t0 = time.clock()
    e_s = 0.5 * 10**(2 - sigFigs)
    e_ax = 100
    e_ay = 100
    checku = abs((2*x_i + y_i)) + abs(x_i)
    checkv = abs((3*y_i**2)) + abs(1 + 6*y_i)
    if (checku and checkv) > 1:
        i = True
    while i is True:
        pp = False
        mm = False
        mp = False
        pm = False 
        inchecku = abs(checku) + 1 
        incheckv = abs(checkv) + 1 
        #check for convergence of the u function
        while inchecku >= checku:
            x_i -= 1
            y_i -= 1
            newChecku = abs(2*x_i + y_i) + abs(x_i) 
            if newChecku > inchecku:
                break
            if inchecku < checku:
                mm = True 
            inchecku = newChecku
            print("C1")
            print(inchecku) 
            if inchecku < checku:
                mm = True
        while inchecku >= checku:
            x_i += 1
            y_i += 1
            newChecku = abs(2*x_i + y_i) + abs(x_i)
            if newChecku > inchecku:
                break
            inchecku = newChecku
            print("C2")
            print(inchecku)
            if inchecku < checku:
                pp = True
        while inchecku >= checku:
            x_i += 1
            y_i -= 1
            newChecku = abs(2*x_i + y_i) + abs(x_i)
            if newChecku > inchecku:
                break
            inchecku = newChecku
            print("C3")
            print(inchecku)
            if inchecku < checku:
                pm = True
        while inchecku >= checku:
            x_i -= 1
            y_i += 1
            newChecku = abs(2*x_i + y_i) + abs(x_i)
            if newChecku > inchecku:
                break
            inchecku = newChecku
            print("C4")
            print(inchecku)
            if inchecku < checku:
                mp = True
        #Check for convergence of the v function
        while incheckv >= checkv:
            x_i += 1
            y_i += 1
            newCheckv = (3*y_i**2) + (1 + 6*y_i) 
            if newCheckv > incheckv:
                break
            incheckv = newCheckv
            print("C5")
            print(incheckv)
            if incheckv < checkv:
                pp = True
        while incheckv >= checkv:
            x_i -= 1
            y_i -= 1
            newCheckv = abs(3*y_i**2) + abs(1 + 6*y_i)
            if newCheckv > incheckv:
                break
            incheckv = newCheckv
            print("C6")
            print(incheckv)
            if incheckv < checkv:
                mm = True
        while incheckv >= checkv:
            x_i += 1
            y_i -= 1
            newCheckv = abs(3*y_i**2) + abs(1 + 6*y_i)
            if newCheckv > incheckv:
                break
            incheckv = newCheckv
            print("C7")
            print(incheckv)
            if incheckv < checkv:
                pm = True
        while incheckv >= checkv:
            x_i -= 1
            y_i += 1
            newCheckv = abs(3*y_i**2) + abs(1 + 6*y_i)
            if newCheckv > incheckv:
                break
            incheckv = newCheckv
            print("C8")
            print(incheckv)
            if incheckv < checkv:
                mp = True
        if pp is True:
            checku = abs(2*x_i + y_i) + abs(x_i)
            checkv = abs(3*y_i**2) + abs(1 + 6*y_i)
            while (checku and checkv) > 1: 
                x_i += 1
                y_i += 1
                checku = abs(2*x_i + y_i) + abs(x_i)
                checkv = abs(3*y_i**2) + abs(1 + 6*y_i)
            i = True
        elif mm is True:
            checku = abs(2*x_i + y_i) + abs(x_i)
            checkv = abs(3*y_i**2) + abs(1 + 6*y_i)
            while (checku and checkv) > 1: 
                x_i -= 1
                y_i -= 1
                checku = abs(2*x_i + y_i) + abs(x_i)
                checkv = abs(3*y_i**2) + abs(1 + 6*y_i)
            i = True
        elif mp is True:
            checku = abs(2*x_i + y_i) + abs(x_i)
            checkv = abs(3*y_i**2) + abs(1 + 6*y_i)
            while (checku and checkv) > 1: 
                x_i -= 1
                y_i += 1
                checku = abs(2*x_i + y_i) + abs(x_i)
                checkv = abs(3*y_i**2) + abs(1 + 6*y_i)
            i = True
        elif pm is True:
            checku = abs(2*x_i + y_i) + abs(x_i)
            checkv = abs(3*y_i**2) + abs(1 + 6*y_i)
            while (checku and checkv) > 1: 
                x_i += 1
                y_i -= 1
                checku = abs(2*x_i + y_i) + abs(x_i)
                checkv = abs(3*y_i**2) + abs(1 + 6*y_i)
            i = True 
            
    while e_ax > e_s and e_ay > e_s: # Execute root finding iteration 
        xPlus = x_i - ((x_i**2 + x_i*y_i -10)*(1 + 6*x_i*y_i) - (y_i + 3*x_i*y_i**2 - 57)*(x_i))/((2*x_i + y_i)*(1 + 6*x_i*y_i) - (x_i)*(3*y_i**2))  
        yPlus = y_i - ((y_i + 3*x_i*y_i**2 - 57)*(2*x_i + y_i) - (x_i**2 + x_i*y_i -10)*(3*y_i**2))/((2*x_i + y_i)*(1 + 6*x_i*y_i) - (x_i)*(3*y_i**2))
        e_ax = abs((xPlus - x_i)/xPlus)*100
        e_ay = abs((yPlus - y_i)/yPlus)*100
        x_i = xPlus
        y_i = yPlus 
    deltat = time.clock() - t0
    print("The final x value is " + str(xPlus))
    print("The final y value is " + str(yPlus))
    print("The final e_a value is " + str(e_ax))
    print("The final e_a value is " + str(e_ay)) 
    print("The fixed point method found a root in " + str(deltat) + " seconds.") 
    
systemsNR()

        
        
