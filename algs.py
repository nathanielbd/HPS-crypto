import math
import pandas as pd
import random

# TODO: Pollard's rho algorithm, Pohlig-Hellman algorithm

def gcd(a, b, verbose=False):
    """
    Compute the GCD using the Euclidean algorithm described in
    theorem 1.7
    """
    if verbose:
        print(f'Computing GCD({a}, {b})...')
        a, b = min(a, b), max(a, b)
        tex = r'\begin{align*}' + '\n'
        while a != 0:
            tex += f'\t{b} &= {b//a} \cdot {a} + {b%a} \\\\ \n'
            a, b = b%a, a
        tex += r'\end{align*}'
        print(tex)
        return b
    else:
        return a if b == 0 else gcd(b, a%b)

def bezout(a, b, verbose=False, lenstra=False):
    """
    Compute the Bezout coefficients using the extended Euclidean algorithm
    described in theorem 1.11
    i.e. compute u,v where au + bv = gcd(a, b)
    """
    a, b = min(a, b), max(a, b)
    # d holds the quotients
    d = [b//a]
    # e holds the remainders
    e = [a, b%a]
    while e[-1] != 0:
        d.append(e[-2]//e[-1])
        e.append(e[-2]%e[-1])
    if e[-2] > 1 and lenstra:
        if verbose:
            gcd(a, b, verbose=True)
        raise NonUnityGCDEvent(e[-2])
    # z holds the intermediates to the Bezout coefficients
    z = [(0, 1), (1, -d[0])]
    for i in range(1, len(d)-1):
        x = z[i-1][0] - d[i]*z[i][0]
        y = z[i-1][1] - d[i]*z[i][1]
        z.append((x,y))
    if verbose:
        c = gcd(a, b, verbose=True)
        # express the remainders in terms of a, b, and the intermediate coefficients
        print(f'Computing Bezout coefficients for {a} and {b}...')
        tex = r'\begin{align*}' + '\n'
        for i, (x,y) in enumerate(z[1:]):
            if i == 0:
                tex += f'\t{e[i+1]} &= {x} \cdot {b} + {y} \cdot {a} \\\\ \n'
            else:
                tex += f'\t{e[i+1]} &= {z[i-1][0]} \cdot {b} + {z[i-1][1]} \cdot {a} - {d[i]} ({z[i][0]} \cdot {b} + {z[i][1]} \cdot {a}) = {x} \cdot {b} + {y} \cdot {a} \\\\ \n'
        tex += r'\end{align*}'
        print(tex)
    return z[-1][1], z[-1][0]

def fast_pow(g, A, N, verbose=False):
    """
    Exponentiating by squaring aka square-and-multiply
    aka binary exponentiation. Computes
        g^A (mod N)
    as described in Exercise 1.25
    """
    # hold intermediate computations in the lists
    a = [g]
    b = [1]
    As = [A]
    while As[-1] > 0:
        if As[-1]%2==1:
            b.append((b[-1]*a[-1])%N)
        else:
            b.append(b[-1])
        a.append((a[-1]**2)%N)
        As.append(As[-1]//2)
    if verbose:
        print(f'Computing modular exponentiation {g}^{a} (mod {N})...')
        print(pd.DataFrame({'$a$': a, '$b$': b, '$A$': As}).to_latex(index=False, escape=False))
    return b[-1]

def mr_test(n, a, verbose=False):
    """
    Miller-Rabin test for composite numbers described in
    table 3.2. Failure indicates that a is not a Miller-Rabin
    witness for the compositeness of n
    """
    if n%2 == 0 or 1 < gcd(a, n, verbose=verbose) < n:
        return 'composite'
    k = 0
    while ((n-1)/(2**k))%2==0:
            k += 1
    q = (n-1)/(2**k)
    # write n-1 = (2^k)*q
    a = (a**q)%n
    tex = r'\begin{align*}' + '\n'
    tex += f'\tn - 1 &= {n-1} = 2^{k} \cdot {q} \\\\ \n'
    if a%n == 1:
            return 'fail'
    for i in range(0,k):
            if a%n == n-1:
                    return 'fail'
            tex += f'\t2^{{{2**i}\cdot{q}}} &\equiv {a} \pmod{{{n}}} \\\\ \n'
            a = (a**2)%n
    tex += r'\end{align*}'
    if verbose:
        print(f'Performing Miller-Rabin test for {n} with possible witness {a}...')
        print(tex)
    return 'composite'

def pollard_pminus(N, k=100, verbose=False):
    """
    Pollard's p-1 factorization algorithm described in
    table 3.3. Finds a nontrivial factor of N by searching
    up to j=k
    """
    a = 2
    tex = r'\begin{align*}' + '\n'
    for j in range(2,k):
        a = (a**j)%N
        d = gcd(a-1,N)
        tex += f'\t2^{{{j}!}} - 1 &\equiv {a-1} \pmod{{{N}}}, \quad \mathrm{{gcd}}({a-1}, {N}) &= {d} \\\\ \n'
        if 1 < d < N:
            if verbose:
                tex += r'\end{align*}'
                print(f'Performing Pollard\'s p-1 algorithm to factorize {N}...')
                print(tex)
            return d
    if verbose:
        tex += r'\end{align*}'
        print(f'Performing Pollard\'s p-1 algorithm to factorize {N}...')
        print(tex)
    return 'fail'

def shanks_bsgs(g, h, p, verbose=False):
    """
    Shanks's Babystep-Giantstep algorithm to solve the 
    discrete logarithm problem
        g^x = h (mod p)
    described in proposition 2.21
    """
    n = math.floor(math.sqrt(p - 1)) + 1
    tbl = {pow(g, i, p): i for i in range(n)}
    inv = pow(g, n*(p-2), p)
    h_uk = []
    for j in range(n):
        y = (h*pow(inv, j, p))%p
        h_uk.append(y)
        if y in tbl:
            # found collision
            if verbose:
                print(f'Performing Shanks\'s algorithm to solve DLP {g}^x = {h} (mod {p})...')
                print(pd.DataFrame({'$k$': range(1, n), '$g^k$': list(tbl.keys())[1:], '$h \cdot u^k$': h_uk + ['-']*(n-len(h_uk)-1)}).to_latex(index=False, escape=False))
            return j*n + tbl[y]
    if verbose:
        print(pd.DataFrame({'$k$': range(1, n), '$g^k$': list(tbl.keys())[1:], '$h \cdot u^k$': h_uk}).to_latex(index=False, escape=False))
    return None

# Elliptic Curve Cryptography

def pt_add(P, Q, p, verbose=False, lenstra=False):
    """
    Elliptic curve addition algorithm described in theorem 6.6:
    Given
        P = (x_1, y_1),
        Q = (x_2, y_2),
        base p of the finite field,
    return
        R = P + Q = (x,y)
    """
    # 0 represents the identity element of E(F_p)
    if P==0:
        return Q
    if Q==0:
        return P
    tex = r'\begin{align*}' + '\n'
    try:
        # compute modular multiplicative inverse with extended Euclidean algorithm
        inv = bezout((Q[0]-P[0])%p, p, verbose=verbose, lenstra=lenstra)[0]
        l = ((Q[1]-P[1])%p)*inv
        tex += f'\t\lambda &= \\frac{{{Q[1]}-{P[1]}}}{{{Q[0]}-{P[0]}}} = \\frac{{{(Q[1]-P[1])%p}}}{{{(Q[0]-P[0])%p}}} \pmod{{{p}}} \\\\ \n'
        tex += f'\t        &= {(Q[1]-P[1])%p} \cdot {inv} \equiv {l%p} \pmod{{{p}}} \\\\ \n'
    except NonUnityGCDEvent as e:
        # will be raised by `bezout` if this computation is part of a call to
        # `lenstras` and a GCD>1 is found -- will be caught later up in the
        # stack
        raise e
    except:
        # could alternatively been another case, but this can be caught here
        if verbose:
            print('$x_1 = x_2$ and $y_1 = -y_2$, so $P + Q = \mathcal{O}$')
        return 0
    x = (l**2 - P[0] - Q[0])%p
    y = (l*(P[0]-x)-P[1])%p
    tex += f'\tx       &= \lambda^2 - x(P) - x(Q) \equiv {x} \pmod{{{p}}} \\\\ \n'
    tex += f'\ty       &= \lambda(x(P)-x)-y(P) \equiv {y} \pmod{{{p}}} \\\\ \n'
    tex += r'\end{align*}'
    if verbose:
        print(f'Adding {P} + {Q}...')
        print(tex)
    return (x,y)

def pt_double(P, p, a, verbose=False, lenstra=False):
    """
    Optional method for the special case of `pt_add` when
    Q=P. Alternatively, call `pt_add(P, P, p)`
    """
    if P==0:
        if verbose:
            print('$P = 2P = \mathcal{O}$')
        return P
    tex = r'\begin{align*}' + '\n'
    try:
        # compute modular multiplicative inverse with extended Euclidean algorithm
        inv = bezout((2*P[1])%p, p, verbose=verbose, lenstra=lenstra)[0]
        l = ((3*(P[0]**2)+a)%p)*inv
        tex += f'\t\lambda &= \\frac{{3({P[0]})^2+{a}}}{{2({P[1]})}} = \\frac{{{(3*(P[0]**2)+a)%p}}}{{{(2*P[1])%p}}} \pmod{{{p}}} \\\\ \n'
        tex += f'\t        &= {(3*(P[0]**2)+a)%p} \cdot {inv} \equiv {l%p} \pmod{{{p}}} \\\\ \n'
    except NonUnityGCDEvent as e:
        # will be raised by `bezout` if this computation is part of a call to
        # `lenstras` and a GCD>1 is found -- will be caught later up in the
        # stack
        raise e
    except:
        # could alternatively been another case, but this can be caught here
        return 0
    x = (l**2 - P[0] - P[0])%p
    y = (l*(P[0] - x) - P[1])%p
    tex += f'\tx       &= \lambda^2 - x(P) - x(Q) \equiv {x} \pmod{{{p}}} \\\\ \n'
    tex += f'\ty       &= \lambda(x(P)-x)-y(P) \equiv {y} \pmod{{{p}}} \\\\ \n'
    tex += r'\end{align*}'
    if verbose:
        print(f'Doubling {P}...')
        print(tex)
    return (x,y)

def double_add(P, p, a, b, n, verbose=False, lenstra=False):
    """
    Elliptic curve cryptography analogue to `fast_pow` described in
    table 6.3. Compute nP in field with base p and elliptic curve
        E: Y^2 = X^3 + aX + b
    where
        P = (x,y)
    """
    # hold intermediate computations in the lists
    Q = [P]
    R = [0]
    i = [0]
    ns = [n]
    while ns[-1]>0:
        if ns[-1]%2==1:
            R.append(pt_add(R[-1], Q[-1], p, verbose=verbose, lenstra=lenstra))
        else:
            R.append(R[-1])
        Q.append(pt_double(Q[-1], p, a, verbose=verbose, lenstra=lenstra))
        ns.append(ns[-1]//2)
        i.append(i[-1]+1)
    if verbose:
        print(f'Computing {n}{P}...')
        print(pd.DataFrame({'Step $i$': i, '$n$': ns, '$Q = 2^i P$': Q, '$R$': R}).to_latex(index=False, escape=False))
    return R[-1]

class NonUnityGCDEvent(Exception):
    """
    Exception to be raised during computations in Lenstra's algorithm.

    Attributes:
        d -- integer GCD 1 < d <= N resulting from computation of a modular multiplicative inverse
    """
    def __init__(self, d, message="Found non-unity GCD d"):
        self.d = d
        self.message = message
        super().__init__(self.message)

class RepeatLenstrasEvent(Exception):
    """
    Exception to be raised if d=N in a NonUnityGCDEvent
    """
    def __init__(self, message="GCD = N. Repeat with new curve and point"):
        self.message = message
        super().__init__(self.message)

def lenstras_determ(N, A, a, b, k=100, verbose=False):
    """
    Lenstra's elliptic curve factorization algorithm as
    described in table 6.8:
        Factor integer N
        using point P = (a, b)
        in elliptic curve E: Y^2 = X^3 + AX + B
        where B = b^2 - a^3 - A*a (mod N)
    """
    P = (a, b)
    B = (b**2 - a**3 - A*a)%N
    pts = [P]
    for j in range(2, k):
        try:
            Q = double_add(P, N, A, B, j, verbose=verbose, lenstra=True)
            pts.append(Q)
        except NonUnityGCDEvent as e:
            # possible nontrivial factor
            if e.d < N:
                if verbose:
                    print(f'Found nontrivial factor {e.d} for {N} with curve E: Y^2 = X^3 + {A}X + {B}')
                    print(pd.DataFrame({'$i$': range(1, len(pts)+1), '$i! \cdot P \pmod{N}$': pts}).to_latex(index=False, escape=False))
                return e.d
            else:
                print(f'Failed to find nontrivial factor for {N} with curve E: Y^2 = X^3 + {A}X + {B}')
                # should be caught in `lenstras` or user should change their values of A, a, b
                raise RepeatLenstrasEvent
        P = Q
    if verbose:
        print(f'Failed to find nontrivial factor for {N} with curve E: Y^2 = X^3 + {A}X + {B}')
        print(pd.DataFrame({'$i$': range(1, len(pts)+1), '$i! \cdot P \pmod{N}$': pts}).to_latex(index=False, escape=False))
    return None

def lenstras(N, k=100, verbose=False):
    """
    Randomized wrapper for Lenstra's algorithm. Search
    up to j=k
    """
    for j in range(k):
        A = random.randint(1, N)
        a = random.randint(1, N)
        b = random.randint(1, N)
        try:
            return lenstras_determ(N, A, a, b, verbose=verbose)
        except RepeatLenstrasEvent:
            return lenstras(N, k=k, verbose=verbose)

# print help, function signatures, and docstrings upon import or running 
help(__name__)