#!/usr/bin/env python
import numpy as np
import numpy.linalg as npl

def bancroft(B_pass):
    """

    Calculate a preliminary position
 
    """

    v_light = 299792458.0
    pos = np.matrix( np.zeros(4) )

    for iter in range(1,3):
        B = B_pass
        [m,n] = size(B)

        for i in ranage(1,m+1) :
            x = B(i,0)
            y = B(i,1)
            if iter == 1:
                traveltime = 0.072
            else:
                z = B(i,3);
                rho = (x-pos(0))**2 + (y-pos(1))**2 + (z-pos(2))**2
                traveltime = np.sqrt(rho)/v_light

            angle = traveltime * 7.292115147e-5
            cosa = np.cos(angle)
            sina = np.sin(angle)
            B(i,0) =  cosa*x + sina*y
            B(i,1) = -sina*x + cosa*y

        if m > 4:
            BBB = npl.pinv(B.T*B)*B.T
        else:
            BBB = npl.pinv(B)

        e = np.ones(m,1)
        alpha = np.zeros((m,1))
        #for i in range(1,m+1): 
        #    alpha(i) = lorentz(B(i,:).T,B(i,:).T)/2.
        
        BBBe = BBB*e
        BBBalpha = BBB*alpha
        a = lorentz(BBBe,BBBe)
        b = lorentz(BBBe,BBBalpha) - 1
        c = lorentz(BBBalpha,BBBalpha)
        root = np.sqrt(b*b-a*c)
        r(1) = (-b-root)/a
        r(2) = (-b+root)/a
        possible_pos = np.zeros((4,2))

        for i in range(1,3):
            possible_pos(:,i) = r(i) * BBBe + BBBalpha
            possible_pos(3,i) = -possible_pos(3,i)
        
        for j in range(1,m+1): 
            for i in range(1,3):
                c_dt = possible_pos(3,i)
                calc = norm(B(j,0:2).T - possible_pos(0:2,i))+c_dt
                omc = B(j,3)-calc
                abs_omc(i) = np.abs(omc)

        # discrimination between roots
        if abs_omc(1) > abs_omc(2):
            pos = possible_pos(:,2)
        else
            pos = possible_pos(:,1)


if __name__ == "__main__":

B_pass =[ -11716227.778 -10118754.628 21741083.973 22163882.029,
          -12082643.974 -20428242.179  11741374.154 21492579.823,
           14373286.650 -10448439.349    19596404.858 21492492.771,
           10278432.244 -21116508.618  -12689101.970 25284588.982]]

# Solution:    595025.053  -4856501.221  4078329.981
    bancroft(B_pass)

