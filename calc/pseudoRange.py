#!/usr/bin/env python
import numpy as np
import numpy.linalg as npl

def noApriori_position():
    """
    Calculate a position, with no apriori coordinates from pseudorange observation
 
    position()


    """

    # Receiver position
    Aguess = np.array([ 0, 0, 0])

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

    # predicted ranges from receiver position to SVs
    ctr = 0
    for row in SatPos:
        predictedRange[ctr] = np.sqrt(np.sum(np.square(row - Aguess)))
        ctr += 1    
   
    print("Predicted Range",predictedRange,predictedRange.shape)
    print("Aguess",Aguess,Aguess.shape) 
 
    # Observed minus Computed
    OMC = np.abs(pseudorange * 299792.458) - predictedRange
    #OMC = np.abs(predictedRange * 299792.458) - Aguess.T

    print("OMC",OMC,"shape:",OMC.shape)

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

    print("A",A,A.shape)
    print("OMC",OMC,OMC.shape)
    #dR = (npl.pinv(A.T * A) * A.T) * OMC.T
    dR = npl.pinv(A.transpose() * A) * A.transpose() * OMC
    print("dR",dR,dR.shape)
    

if __name__ == "__main__":
    
    noApriori_position()

