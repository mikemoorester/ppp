#!/usr/bin/env python
import numpy as np
import numpy.linalg as npl

def noApriori_position():
    """
    Calculate a position, with no apriori coordinates from pseudorange observation
 
    position()


    """

    # Receiver position
    #Aguess = np.array([ 0, 0, 0])
    # should I have an initial guess for the receiver clock offset
    Aguess = np.array([ -730000., -5440000., 3230000.])

    # satellite coordinates in ECEF XYZ from ephemeris parameters and SV time
    SatPos = np.array([
                      [15524471.175, -16649826.222, 13512272.387],
                      [-2304058.534, -23287906.465, 11917038.105],
                      [16680243.357, -3069625.561, 20378551.047],
                      [-14799931.395, -21425358.240, 6069947.224]
                     ])

    # satellite pseudorange in meters, from C/A code epochs in milliseconds
    # Range + Receiver clock bias
    pseudorange = np.array([89491.971, 133930.500, 283098.754, 205961.742])

    predictedRange = np.array(np.zeros(SatPos.shape[0]))

    # predicted ranges from receiver position to SV nav position
    ctr = 0
    for row in SatPos:
        predictedRange[ctr] = np.sqrt( (row[0] - Aguess[0])**2 + 
                                       (row[1] - Aguess[1])**2 + 
                                       (row[2] - Aguess[2])**2  )

        #predictedRange[ctr] = np.sqrt(np.sum(np.square(row - pseudorange[ctr])))
        ctr += 1    
   
#    print("Predicted Range",predictedRange,predictedRange.shape)
#    print("Aguess",Aguess,Aguess.shape) 
 
    # Observed minus Computed
    #OMC = np.matrix( (pseudorange / 299792.458) - predictedRange )
    #OMC = np.matrix( np.mod(pseudorange, 299792.458) - predictedRange )
    #OMC = np.matrix( (predictedRange / 299792.458) - pseudorange )
    OMC = np.matrix( np.mod(predictedRange, 299792.458) - pseudorange )

    # Solve for Correction to Receiver position estimates
    A = np.zeros((pseudorange.size,4))
    ctr = 0
    for row in A:
        #A[ctr,2] = (SatPos[ctr,2] - Aguess[2])/ pseudorange[ctr]
        A[ctr,0] = (SatPos[ctr,0] - Aguess[0])/ predictedRange[ctr]
        A[ctr,1] = (SatPos[ctr,1] - Aguess[1])/ predictedRange[ctr]
        A[ctr,2] = (SatPos[ctr,2] - Aguess[2])/ predictedRange[ctr]
        A[ctr,3] = -1.
        ctr += 1

    dR = ((npl.pinv((A.T) * A) ) * A.T) * OMC.T

    Rx = Aguess[0] + dR[0]
    Ry = Aguess[1] + dR[1]
    Rz = Aguess[2] + dR[2]
    Rt = dR[3]

    #print("A:",A)
    print("Rx",Rx)
    print("Ry",Ry)
    print("Rz",Rz)
    print("Rt",Rt)

    print(dR)

if __name__ == "__main__":
    
    noApriori_position()

