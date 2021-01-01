# HPS-crypto

Python implementations of selected algorithms presented in [*An Introduction to Mathematical Cryptography* by Hoffstein, Pipher and Silverman](https://github.com/isislovecruft/library--/blob/master/cryptography%20%26%20mathematics/An%20Introduction%20to%20Mathematical%20Cryptography%20(2014)%20-%20Hoffstein%2C%20Pipher%2C%20Silverman.pdf).

Not intended to be efficient. Merely for learning purposes and for self-learners doing the exercises in LaTeX. Core feature is LaTeX output which shows intermediate steps in computations.

## Examples

Example 1.10

```python
>>> help(bezout)
Help on function bezout in module algs:

bezout(a, b, verbose=False, lenstra=False)
    Compute the Bezout coefficients using the extended Euclidean algorithm
    described in theorem 1.11
    i.e. compute u,v where au + bv = gcd(a, b)
>>> bezout(2024, 748)
(19, -7)
>>> bezout(2024, 748, verbose=True)
Computing GCD(748, 2024)...
\begin{align*}
        2024 &= 2 \cdot 748 + 528 \\
        748 &= 1 \cdot 528 + 220 \\
        528 &= 2 \cdot 220 + 88 \\
        220 &= 2 \cdot 88 + 44 \\
        88 &= 2 \cdot 44 + 0 \\
\end{align*}
Computing Bezout coefficients for 748 and 2024...
\begin{align*}
        528 &= 1 \cdot 2024 + -2 \cdot 748 \\
        220 &= 0 \cdot 2024 + 1 \cdot 748 - 1 (1 \cdot 2024 + -2 \cdot 748) = -1 \cdot 2024 + 3 \cdot 748 \\
        88 &= 1 \cdot 2024 + -2 \cdot 748 - 2 (-1 \cdot 2024 + 3 \cdot 748) = 3 \cdot 2024 + -8 \cdot 748 \\
        44 &= -1 \cdot 2024 + 3 \cdot 748 - 2 (3 \cdot 2024 + -8 \cdot 748) = -7 \cdot 2024 + 19 \cdot 748 \\
\end{align*}
(19, -7)
```

Example 6.16

```python
>>> help(double_add)
Help on function double_add in module algs:

double_add(P, p, a, b, n, verbose=False, lenstra=False)
    Elliptic curve cryptography analogue to `fast_pow` described in
    table 6.3. Compute nP in field with base p and elliptic curve
        E: Y^2 = X^3 + aX + b
    where
        P = (x,y)
>>> double_add((6, 730), 3623, 14, 19, 947)
(3492, 60)
>>> double_add((6, 730), 3623, 14, 19, 947, verbose=True)
... depth-first stacktrace of intermediate computations will print
Computing 947(6, 730)...
\begin{tabular}{rrll}
\toprule
 Step $i$ &  $n$ &   $Q = 2^i P$ &           $R$ \\
\midrule
        0 &  947 &      (6, 730) &             0 \\
        1 &  473 &  (2521, 3601) &      (6, 730) \\
        2 &  236 &   (2277, 502) &   (2149, 196) \\
        3 &  118 &   (3375, 535) &   (2149, 196) \\
        4 &   59 &  (1610, 1851) &   (2149, 196) \\
        5 &   29 &  (1753, 2436) &  (2838, 2175) \\
        6 &   14 &  (2005, 1764) &   (600, 2449) \\
        7 &    7 &  (2425, 1791) &   (600, 2449) \\
        8 &    3 &  (3529, 2158) &  (3247, 2849) \\
        9 &    1 &  (2742, 3254) &   (932, 1204) \\
       10 &    0 &  (1814, 3480) &    (3492, 60) \\
\bottomrule
\end{tabular}
```

## Requirements

- Python 3
- Pandas

## Usage

- Integrate into your LaTeX workflow
  - Make sure to include `\usepackage{booktabs}` and `\usepackage{amsmath}`
  - Either copy-paste or pipe `stdout` to `file.tex` and read with `\input{file.tex}`

OR

- Run in debugger with breakpoints while reading the textbook's pseudocode

OR

- Run in python REPL:

```bash
$ git clone https://github.com/nathaniel/HPS-crypto.git
$ cd HPS-crypto/
$ python
>>> from algs import *
... help, function signatures, docstrings will print
```

## Future

- Implement Pohlig-Hellman algorithm
- Implement Pollard's rho algorithm
- Create Flask/KaTeX web application