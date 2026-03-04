# Shift representability: The UVSD case study

## I. The Two-Distribution Case

Consider two continuous, strictly increasing distribution functions $F_1$ and $F_2$ on a common interval $I$, and define the transport map

$$
\psi(x) = F_1^{-1}(F_2(x)).
$$

This map is itself continuous and strictly increasing. If the two distributions do not cross, then $\psi(x)\neq x$ for all $x\in I$. Because $\psi$ is continuous and increasing, this implies that either $\psi(x)>x$ for all $x\in I$ or $\psi(x)<x$ for all $x\in I$. In other words, $\psi$ is fixed-point-free.

For such maps, Abel's functional equation provides the relevant representation theorem: under mild regularity conditions, there exists a strictly monotone function $g$ and a constant $c$ such that

$$
g(\psi(x)) = g(x) + c.
$$

Equivalently, after the monotone reparameterization $y=g(x)$, the transport map becomes a pure translation. This implies that any pair of continuous strictly non-crossing distributions is shift-representable. In the two-distribution case, there is no commutativity constraint, because there is only one transport map to linearize.

The substantive point is that with only two distributions, any fixed-point-free geometric distortion between them can always be absorbed into a suitable monotone warping of the measurement axis.

---

## II. Unequal-Variance Signal Detection: Global Failure and Local Recovery

Consider the unequal-variance signal detection model with

$$

N \sim \operatorname{Normal}(0,1), \qquad S \sim \operatorname{Normal}(\mu,\sigma^2),
$$

where $\sigma\neq 1$.

### Global failure due to crossing on the full real line

The transport map from signal to noise is

$$
\psi(x)=F_N^{-1}(F_S(x)) = (x-\mu) / \sigma
$$

A fixed point $\psi(x) = x$ thus satisfies

$$
x = \mu + \sigma x,
$$

so the unique fixed point is

$$

x^* = \frac{\mu}{1-\sigma}.
$$

This is exactly the point at which the two normal cdfs intersect. Because $\psi$ has a fixed point, it cannot be conjugate on all of $\mathbb R$ to a nontrivial translation. Indeed, if one had

$$
g(\psi(x))=g(x)+c,
$$

then evaluating at $x=x^*$ would give

$$

g(x^*) = g(x^*) + c,
$$

forcing $c=0$. Thus, on the full real line, the unequal-variance Gaussian pair is not shift-representable.

### Recovery on a one-sided domain

If one restricts attention to one side of the crossing point, the fixed point is removed from the domain. Write

$$

y = x - x^*.
$$

Then the transport map becomes

$$
\psi(y) = y/\sigma.
$$

On either half-line $y>0$ or $y<0$, this map is fixed-point-free. Abel's equation now takes the form

$$
g(y / \sigma)=g(y)+c,
$$

whose standard monotone solution is logarithmic. On $y>0$, one may take

$$

g(y)=\log y,
$$

and on $y<0$,

$$

g(y)=\log(-y).
$$

Translating back to the original coordinate, the corresponding transformation is

$$

g(x)=\log(x-x^*) \qquad \text{for } x>x^*,
$$

or

$$

g(x)=\log(x^*-x) \qquad \text{for } x<x^*.
$$

Thus the unequal-variance Gaussian pair is not globally shift-representable, but it is shift-representable on either side of the crossing point separately. The logarithmic transformation introduces a vertical asymptote at the crossing boundary and converts the affine dilation about that point into an additive shift.

### Definition: **weak shift-representability**

A set of distributions will be called weakly shift-representable, if they are not shift-representable on their entire domain, but they become shift-representable when the values they could take are restricted to fall only one side of the fixed crossing point. A set of distributions that are shift-representable on their entire support satisfy strong shift-representability, or for short, are shift-representable.

---

## IV. Multiple Signal Conditions and the Common-Crossing Constraint

Now suppose we have one noise distribution and several signal conditions, all normally distributed:

$$

N \sim \operatorname{Normal}(0,1), \qquad S_i \sim \operatorname{Normal}(\mu_i,\sigma_i^2), \quad i=1,\dots,k.
$$

From the two-variable case we know that this set of distributions is not strongly shift-representable. Are there conditions under which they are weakly shift-representable? 

The transport maps are

$$
\psi_i(x)=(x-\mu_i)/\sigma_i.
$$

To obtain a single transformation $g$ that shift-represents the entire family, these maps must be simultaneously conjugate to translations. A necessary condition is therefore that they commute pairwise, just as for full shift-representability:

$$
\psi_i\circ \psi_j = \psi_j\circ \psi_i \qquad \text{for all } i,j.
$$

Computing both compositions gives

$$
\psi_i(\psi_j(x)) = ((x - \mu_j)/\sigma_j - \mu_i)/\sigma_i = \frac{x}{\sigma_j\sigma_i} - \frac{\mu_j}{\sigma_j\sigma_i} - \frac{\mu_i}{\sigma_i}
$$

and

$$
\psi_j(\psi_i(x)) = \frac{x}{\sigma_i\sigma_j} - \frac{\mu_i}{\sigma_i\sigma_j} - \frac{\mu_j}{\sigma_j}
$$

Thus commutativity holds if and only if

$$
\frac{\mu_i}{(1-\sigma_i)}=\frac{\mu_j}{(1-\sigma_j)}.
$$

in other words - this ratio is a constant that doesn’t depend on $i$ or $j$. Hence, there exists a constant $a$ such that

$$

\sigma_i = 1 + a\mu_i

\qquad\text{for all } i.
$$

This is the linear relation between standard deviation and mean.

Under this condition, every transport map has the same fixed point:

$$

x_i^*=\frac{\mu_i}{1-\sigma_i}

=\frac{\mu_i}{-a\mu_i}

=-\frac{1}{a},
$$

independent of $i$. Thus all distributions cross the noise distribution at the same coordinate, and the same logarithmic transformation can be used for all conditions:

$$

g(x)=\log\!\left(x+\frac{1}{a}\right)
$$

on the appropriate half-line.

Hence, under the weaker one-sided notion of shift-representability, a multi-condition unequal-variance Gaussian family is admissible only when all signal distributions share a common crossing point with the noise distribution. For affine Gaussian transports, this is equivalent to requiring that the standard deviations vary linearly with the means, with intercept $1$.

As an example, consider the common finding that the standard deviation of the signal in recognition memory is $\sigma_s = 1.25$. For the UVSD model to be one-sidedly shift-representable, whenever we measure multiple conditions $i, j$ such that $\mu_j > \mu_i > 0$, then the standard deviation cannot be shared in the signal conditions. If we choose $\mu_i = 1, \mu_j = 2$, and that $\sigma_i = 1.25$, then we must have $\sigma_j = 1 + 0.25\mu_j = 1.5$. In such case the distributions will cross at $x^*=-4$. Under this model if the standard full support $x \in (-\infty, \infty)$ is considered, then $1 - \Phi(-4) = 0.99997$ proportion of the probability mass falls entirely to the right of $x^*$ and the model could be reasonably restricted to the interval $(-4, \infty)$ without substantially changing its predictions, but making it shift-representable.

However this weak shift-representability comes at a stringent cost - if the measured standard deviation in recognition memory tasks doesn’t scale *exactly* as the linear function $1 + a \text{d'}$ of sensitivity $d'$, then the model is not shift-representable, neither strongly nor weakly.

---

## Summary

The two-distribution case is structurally simple: any fixed-point-free transport map can be converted into a translation by a suitable monotone reparameterization, so any pair of continuous non-crossing distributions is shift-representable.

The unequal-variance Gaussian model fails globally because its transport map has a fixed point, namely the crossing point of the two cdfs. However, once one restricts attention to one side of that point, the map becomes conjugate to a translation via a logarithmic transformation.

With multiple signal conditions, the problem changes qualitatively. One must find a single transformation that linearizes several transport maps simultaneously. This introduces a commutativity constraint, and for Gaussian affine transports that constraint is satisfied only when all distributions share the same crossing point. Algebraically, this occurs exactly when the signal standard deviations are linear functions of their means with intercept $1$.