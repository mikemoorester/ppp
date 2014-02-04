#!/usr/bin/env python
import numpy as np
import numpy.linalg as npl

def lorentz(x,y):

    """ LORENTZ  Calculates the Lorentz inner product
%          of the two 4 by 1 vectors x and y

% Written by Kai Borre
% April 22, 1995

%  M = diag([1 1 1 -1]);
%  p = x'*M*y;
    """
    p = x[0]*y[0] + x[1]*y[1] + x[2]*y[2] - x[3]*y[3]
    return p

def bancroft(B_pass):
    """

    Calculate a preliminary GPS user's position by using a direct 
      (non-iterative) GPS solution method proposed by Bancroft [1] [2], 
       when at least four satellite positions and the corresponding 
       pseudoranges are known. 

%  Usage:   [nsol,upos,ucbias] = uposdg(svpos,rho) 

   Description of parameters: 
                svpos  -  input, at least four ECEF satellite positions, each row  
                          contains the three component for a satellite,  
                          all components are in meters 
                rho    -  input, columnwise, at least four pseudoranges,  
                          the pseudoranges are in meters 
                nsol   -  output, number of solutions (0, 1 or 2) 
                upos   -  output, ECEF user's position, the components are  
                          in meters 
                ucbias -  output, user's clock bias (the difference between user's 
                          time and GPS time) measured in units of distance (meters) 
   References: 
       [1] Bancroft, S., An algebraic solution of the GPS equations. 
           IEEE Transactions on Aerospace and Electronic Systems,  
           vol. AES-21, No. 7, 1985, pp. 56-59. 
       [2] Chaffee, J. W., Abel, J. S., Bifurcation of pseudorange  
           equations. Institute of Navigation, Proceedings of the  
           1993 National Technical Meeting, San Francisco, CA, Jan. 
           20-22, 1993, pp. 203-211. 

    """

    v_light = 299792458.0
    pos = np.matrix( np.zeros(4) )

    for iter in range(1,3):
        B = B_pass
        [m,n] = np.shape(B)

        for i in range(0,m) :
            print(i)
            x = B[i][0]
            y = B[i][1]
            if iter == 1:
                traveltime = 0.072
            else:
                z = B[i][3]
                rho = (x-pos[0])**2 + (y-pos[1])**2 + (z-pos[2])**2
                traveltime = np.sqrt(rho)/v_light

            angle = traveltime * 7.292115147e-5
            cosa = np.cos(angle)
            sina = np.sin(angle)
            B[i][0] =  cosa*x + sina*y
            B[i][1] = -sina*x + cosa*y

        if m > 4:
            BBB = npl.pinv(B.T*B)*B.T
        else:
            BBB = npl.pinv(B)

        e = np.ones(m)
        #e = np.ones(m,1)
        alpha = np.zeros((m,1))
        #for i in range(1,m+1): 
        #    alpha(i) = lorentz(B(i,:).T,B(i,:).T)/2.
        
        BBBe = BBB*e
        BBBalpha = BBB*alpha
        a = lorentz(BBBe,BBBe)
        b = lorentz(BBBe,BBBalpha) - 1
        c = lorentz(BBBalpha,BBBalpha)
        root = np.sqrt(b*b-a*c)
        print("root",root)

        # Look at the two possible solutions
        r = np.zeros(2)
        print("r",r,np.shape(r))
        print("b",b)
        print("a",a)
        r[0] = (-b-root)/a
        r[1] = (-b+root)/a
        possible_pos = np.zeros((4,2))

        for i in range(0,2):
            possible_pos[:,i] = r[i] * BBBe + BBBalpha
            possible_pos[3,i] = -possible_pos[3,i]
        
        for j in range(1,m+1): 
            for i in range(1,3):
                c_dt = possible_pos(3,i)
                calc = norm(B[j,0:2].T - possible_pos[0:2,i])+c_dt
                omc = B[j,3]-calc
                abs_omc[i] = np.abs(omc)

        # discrimination between roots
        if abs_omc(1) > abs_omc(2):
            pos = possible_pos[:,2]
        else:
            pos = possible_pos[:,1]


if __name__ == "__main__":

    B_pass =[ [-11716227.778, -10118754.628, 21741083.973, 22163882.029],
          [-12082643.974, -20428242.179, 11741374.154, 21492579.823],
          [ 14373286.650, -10448439.349, 19596404.858, 21492492.771],
          [ 10278432.244, -21116508.618, -12689101.970, 25284588.982]]

# Solution:    595025.053  -4856501.221  4078329.981
    bancroft(B_pass)

