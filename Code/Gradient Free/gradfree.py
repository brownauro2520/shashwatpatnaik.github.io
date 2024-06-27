import numpy as np

def simplex(n, l, x0, func):

    X_vec = [x0]
    f0 = func(x0)[0]
    F_vec = [f0]
    
    for j in range(n):
        s_vec = np.zeros(n)

        for i in range(n):
            if i == j:
                s_vec[i] = (l / (n*np.sqrt(2))) * (np.sqrt(n + 1) - 1) + (l / np.sqrt(2))
            else:
                s_vec[i] = (l / (n*np.sqrt(2))) * (np.sqrt(n + 1) - 1)

        x_new = x0 + s_vec
        X_vec = np.vstack([X_vec, x_new])

        f_new = func(x_new)[0]
        F_vec = np.append(F_vec, f_new)
    
    return X_vec, F_vec

def sort(X_vec, F_vec):
    # Create a list of tuples (function value, corresponding point)
    simplex_tuples = list(zip(F_vec, X_vec))
    # Sort the list of tuples based on function values
    simplex_tuples.sort(key=lambda x: x[0])
    # Recreate X_vec and F_vec using the sorted list of tuples
    X_vec = np.array([point for _, point in simplex_tuples])
    F_vec = np.array([value for value, _ in simplex_tuples])
    return simplex_tuples, X_vec, F_vec


def delx(X_vec, n):

    delx = 0
    for i in range(n):
        norm = np.linalg.norm(np.subtract(X_vec[i], X_vec[-1]))
        delx += norm

    return delx

def delf(F_vec, n):

    sum_term = 0
    f_bar = np.sum(F_vec) / (n + 1)
    for i in range(n+1):
        square_term = (F_vec[i] - f_bar)**2
        sum_term += square_term
    delf = np.sqrt(sum_term / (n + 1))

    return delf



def NDM(x0, l, func, tau_x, tau_f):
    
    n = len(x0)

    X_vec, F_vec = simplex(n, l, x0, func)
    
    del_x = delx(X_vec, n)
    del_f = delf(F_vec, n)

    # history vectors
    normx = [del_x]
    normf = [del_f]
    simplex_hist = [X_vec.copy()]  

    k = 1
    while del_x > tau_x or del_f > tau_f:

        simplex_dict, X_vec, F_vec = sort(X_vec, F_vec)

        x_c = np.sum(X_vec[:-1], axis=0) / n
        x_r = x_c + (x_c - X_vec[-1])

        f_r = func(x_r)[0]

        if f_r < F_vec[0]:
            x_e = x_c + 2 * (x_c - X_vec[-1])
            f_e = func(x_e)[0]

            if f_e < F_vec[0]:
                X_vec[-1] = x_e
                F_vec[-1] = func(X_vec[-1])[0]
            else:
                X_vec[-1] = x_r
                F_vec[-1] = func(X_vec[-1])[0]
        
        elif f_r <= F_vec[-2]:
            X_vec[-1] = x_r
            F_vec[-1] = func(X_vec[-1])[0]
        
        else:
            if f_r > F_vec[-1]:
                x_ie = x_c - 0.5 * (x_c - X_vec[-1])
                f_ie = func(x_ie)[0]

                if f_ie < F_vec[-1]:
                    X_vec[-1] = x_ie
                    F_vec[-1] = func(X_vec[-1])[0]
                else:
                    for j in range(n+1):
                        X_vec[j] = X_vec[0] + 0.5 * (X_vec[j] - X_vec[0])
                        F_vec[j] = func(X_vec[j])[0]
            else:
                x_oc = x_c + 0.5 * (x_c - X_vec[-1])
                f_oc = func(x_oc)[0]
                if f_oc < f_r:
                    X_vec[-1] = x_oc
                    F_vec[-1] = func(X_vec[-1])[0]
                else:
                    for j in range(n+1):
                        X_vec[j] = X_vec[0] + 0.5 * (X_vec[j] - X_vec[0])
                        F_vec[j] = func(X_vec[j])[0]

        #updating error
        del_x = delx(X_vec, n)
        del_f = delf(F_vec, n)
        
        #storing optmization path
        normx.append(del_x)
        normf.append(del_f)
     
        simplex_hist.append(X_vec.copy())

        k += 1

        #print('Iteration ', k, ' Error delx = ', del_x, ' Error delf = ', del_f)

    simplex_dict, X_vec, F_vec = sort(X_vec, F_vec)
    xopt = X_vec[0]
    fopt = F_vec[0]
    return xopt, fopt, normx, normf, simplex_hist