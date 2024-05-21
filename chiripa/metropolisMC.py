import random
import math

def metropolis_MC(E_current, E_next, T, iseed=None):

    """

    Args:
        E_current:
        E_next:
        T:
        iseed:

    Returns:

    """

    if iseed is None:
        random.seed()
    else:
        random.seed(iseed)

    R = 1.987 * 1e-03 # kcal/molK
    RT = T * R # kcal/mol

    if E_next < E_current:
        return True, RT
    else:
        u = random.uniform(0.0, 1.0)
        DE = E_next - E_current
        f = math.exp(-DE/RT)
        if u <= f:
            return True, RT
        else:
            return False, RT
