### Transport-invariance identification from exemplar CDFs

#### 1) Objects and assumptions

- Let $\mathcal{F}=\{F_\theta:\theta\in\Theta\subset\mathbb{R}^d\}$ be a parametric family of strictly increasing CDFs on a common support $I\subseteq\mathbb{R}$.

- You are given $k\ge 2$ exemplar members from this family:

   $$
   F_1,\dots,F_k,
   $$

   where each $F_i:I\to(0,1)$ is strictly increasing and invertible with quantile function $Q_i(u):=F_i^{-1}(u)$ for $u\in(0,1)$.

- Goal: determine whether there exists a parameter $\theta_*\in\Theta$ such that the corresponding CDF $F_{\theta_*}$ is **invariant under the exemplar-induced transport dynamics**, and if so, whether $\theta_*$ is (numerically) unique.

---

#### 2) Transport maps and the invariance condition

**Probability transport (on $(0,1)$).**\
For each exemplar $F_i$, define the probability-scale transport from a chosen reference $F_1$ by

$$
\phi_i(u) := F_i(Q_1(u)),\qquad u\in(0,1).
$$

This is a strictly increasing bijection $\phi_i:(0,1)\to(0,1)$.

**Candidate CDF in $u$\-coordinates.**\
For any $\theta\in\Theta$, represent $F_\theta$ in the reference probability coordinate $u=F_1(x)$ by

$$
H_\theta(u) := F_\theta(Q_1(u)),\qquad u\in(0,1).
$$

So $H_\theta:(0,1)\to(0,1)$ is the “CDF of $F_\theta$ expressed on the $F_1$ probability scale.”

**Transport invariance / commutation criterion.**\
The candidate $F_\theta$ is transport-invariant (compatible with the exemplar planar transport dynamics) if and only if $H_\theta$ commutes with every $\phi_i$, i.e.

$$
H_\theta(\phi_i(u))=\phi_i(H_\theta(u))\quad \text{for all }u\in(0,1)\text{ and for }i=2,\dots,k.
$$

Define the commutator residual

$$
c_{i,\theta}(u) := H_\theta(\phi_i(u))-\phi_i(H_\theta(u)).
$$

Then invariance is equivalent to $c_{i,\theta}(u)=0$ for all $u$ and all $i\ge 2$.

---

#### 3) Computational problem (finite discretization)

Choose a mesh $U=\{u_j\}_{j=1}^m\subset(0,1)$ (avoid endpoints; e.g. $u_j\in[\varepsilon,1-\varepsilon]$).

Define the discrete objective

$$
J(\theta):=\sum_{i=2}^k \sum_{j=1}^m w_j\,\rho\!\left(c_{i,\theta}(u_j)\right),
$$

where:

- $w_j>0$ are quadrature weights (often $w_j\equiv 1$ is fine),

- $\rho(z)=z^2$ (least squares) or $\rho(z)=|z|$ (robust absolute loss).

**Existence test (numerical).**\
A transport-invariant member exists within the parametric family if

$$
\min_{\theta\in\Theta} J(\theta)\le \tau,
$$

for a chosen tolerance $\tau>0$.

**Uniqueness test (numerical).**\
The solution is treated as unique if repeated optimizations from diverse initializations return the same $\theta$ (within tolerance), and local curvature at the minimizer indicates an isolated minimum (e.g. a well-conditioned approximate Hessian for squared loss).

---

#### 4) How to evaluate the pieces numerically

You need to compute, for many $u$ values:

- $Q_1(u)$ (quantiles of the reference),

- $F_i(x)$ for exemplars,

- $F_\theta(x)$ for candidate parameters.

Then:

- $\phi_i(u)=F_i(Q_1(u))$,

- $H_\theta(u)=F_\theta(Q_1(u))$,

- $c_{i,\theta}(u)=H_\theta(\phi_i(u))-\phi_i(H_\theta(u))$ where

   - $H_\theta(\phi_i(u)) = F_\theta(Q_1(\phi_i(u)))$,

   - $\phi_i(H_\theta(u)) = F_i(Q_1(H_\theta(u)))$.

**Implementation note.** If $Q_1$ is not available in closed form, approximate it by monotone interpolation over a dense precomputed table of $(u,Q_1(u))$ values.

---

#### 5) Minimal algorithm (planning-level pseudocode)

**Inputs**

- Exemplar CDF routines $\{F_i\}_{i=1}^k$ and reference quantile routine $Q_1$.

- Parametric candidate routine $F_\theta$ (and optionally its gradients in $\theta$).

- Parameter domain $\Theta$, mesh $U=\{u_j\}_{j=1}^m$, loss $\rho$, weights $\{w_j\}$.

- Optimizer settings: multi-start count $M$, tolerance $\tau$.

**Precompute**

- For each $j$, compute $x_j := Q_1(u_j)$.

- For each $i=2,\dots,k$ and $j$, compute $\phi_i(u_j)=F_i(x_j)$.

**Objective evaluation at $\theta$**

- For each $j$:

   - compute $H_\theta(u_j)=F_\theta(x_j)$.

- For each $i=2,\dots,k$ and $j$:

   - compute $H_\theta(\phi_i(u_j)) = F_\theta(Q_1(\phi_i(u_j)))$,

   - compute $\phi_i(H_\theta(u_j)) = F_i(Q_1(H_\theta(u_j)))$,

   - set $c_{i,\theta}(u_j)$ as the difference,

- Return $J(\theta)=\sum_{i,j} w_j\rho(c_{i,\theta}(u_j))$.

**Optimization**

- For $m=1,\dots,M$:

   - initialize $\theta^{(m)}$ (random or grid-based in $\Theta$),

   - run constrained optimization to obtain $\hat\theta^{(m)}=\arg\min_\theta J(\theta)$.

- Choose $\theta_*=\hat\theta^{(m)}$ with smallest achieved objective.

**Decision**

- Existence: declare “invariant member found” if $J(\theta_*)\le\tau$.

- Uniquenessdeclare “unique” if all near-optimal runs return parameters within a small radius of $\theta_*$ and there are no distinct separated minima with objective $\le\tau$.

**Outputs**

- Estimated invariant parameter $\theta_*$ (if found),

- Objective value $J(\theta_*)$,

- Diagnostics: distribution of best objectives across starts; stability of $\theta$ across starts.

---

#### 6) Practical default choices (minimal)

- Mesh: $u_j=\varepsilon + (1-2\varepsilon)\frac{j-1}{m-1}$ with $\varepsilon\in[10^{-4},10^{-2}]$ and $m$ in the low thousands.

- Loss: squared loss $\rho(z)=z^2$.

- Multi-start: $M\ge 10$ for low-dimensional $\theta$, larger if $d$ is bigger or the objective is nonconvex.

- Tolerance: set $\tau$ based on expected numerical error from inverse-CDF approximation and floating-point noise (empirically by checking the objective for known commuting cases, if available).

This document defines the terminology, invariance criterion, and a complete computational workflow for searching within a parametric family for a member whose induced $H_\theta$ commutes with the exemplar probability transports $\{\phi_i\}$.
