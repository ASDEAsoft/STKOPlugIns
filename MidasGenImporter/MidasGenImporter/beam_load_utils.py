from typing import List

def compute_equivalent_constant_load(values : List[float]) -> float:
    '''
    Temporary function, used until we have polygonal load support.
    This is the command in Midas Gen:
    *BEAMLOAD    ; Element Beam Loads
    ; ELEM_LIST, CMD, TYPE, DIR, bPROJ, [ECCEN], [VALUE], GROUP
    ; ELEM_LIST, CMD, TYPE, TYPE, DIR, VX, VY, VZ, bPROJ, [ECCEN], [VALUE], GROUP
    ; [VALUE]       : D1, P1, D2, P2, D3, P3, D4, P4
    ; [ECCEN]       : bECCEN, ECCDIR, I-END, J-END, bJ-END
    ; [ADDITIONAL]  : bADDITIONAL, ADDITIONAL_I-END, ADDITIONAL_J-END, bADDITIONAL_J-END

    In this example:
    521, FLOOR  , UNILOAD, GZ, NO , NO, LY, , , , 0, -0, 0.583894, -16.3125, 1, -16.3125, 0, 0, , NO, 0, 0, NO,

    [VALUE] = 0, -0, 0.583894, -16.3125, 1, -16.3125, 0, 0
    where D1 to D4 are relative distances along the beam length (from 0 to 1)
    and P1 to P4 are the load values at those distances.

    This function computes the equivalent constant load value as:
    N = integral from 0 to 1 of the load distribution
    Q = N / L
    but L = 1, so Q = N
    '''

    if len(values) != 8:
        raise Exception('Invalid number of load values: {}, expecting 8 values'.format(len(values)))
    D1, P1, D2, P2, D3, P3, D4, P4 = values
    # Compute the integral using the trapezoidal rule
    # NOTE: in Midas, if D = 1.0 is reached before D4, the remaining D and P are 0.
    # case: D3 = 1.0, D4 = 0.0 -> we have to make sure D4-D3 is set to 0, so we use max(D2-D1, 0.0)
    N = 0.0
    N += max(D2 - D1, 0.0) * (P1 + P2) / 2.0
    N += max(D3 - D2, 0.0) * (P2 + P3) / 2.0
    N += max(D4 - D3, 0.0) * (P3 + P4) / 2.0
    return N