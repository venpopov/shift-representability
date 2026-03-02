## Definitions

### Family of random variables

Let $X=\{X_1,X_2,\dots,X_k\}$ be a finite set of random variables defined on a common interval DF, which is well-defined for continuous strictly increasing CDFs.

------------------------------------------------------------------------

### Shift-representability

The set $X$ is said to be **shift-representable** if there exists an increasing homeomorphism $g:I\to\mathbb R$, a base CDF $H$ on $\mathbb R$, and constants $\mu_1, \dots, \mu_k \in \mathbb R$ such that

$$
F_n(x) = H(g(x) - \mu_n).
$$

Equivalently, after transforming $Y_n=g(X_n)$, the set $Y = \{Y_1, \dots, Y_k\}$ belongs to a common shift family.

------------------------------------------------------------------------

### Transport maps

A quantile-matching transport map is a function $\psi: I \to I$ that transports the quantiles of one random variable onto the scale of another. For the set $X$ define the pairwise transport maps

$$
\psi_{n \to m} := F_m^{-1} \circ F_n, \qquad n, m \in \{1, \cdots, k\},
$$

$I \subset \mathbb R$. Assume that each $X_n$ has a continuous, strictly increasing cumulative distribution function $F_{X_n}: I \to [0,1]$. For ease of notation, write $F_n\equiv F_{X_n}$. For the quantile function we use the notation $F_n^{-1}$, i.e. the inverse of the Cwhere $f \circ g$ denotes the composition of the functions $f$ and $g$: i.e. $(f \circ g)(x) := f(g(x))$. For ease of notation we sometimes omit $\circ$ and simply write $fg := f \circ g$.

We can choose one variable, such as $X_1$, to be the reference scale, in which case we refer to the reference transport maps

$$
\psi_n := \psi_{n \to 1} = F_1^{-1}F_n.
$$

------------------------------------------------------------------------

### Commutativity

Two functions $f$ and $g$ are said to be commutative if the order in which they are applied doesn't matter:

$$
fg = gf
$$