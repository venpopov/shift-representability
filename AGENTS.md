# AGENTS.md – Shift-representability of distribution families

## Agent Role

You are a Computational Modeling and mathematical psychology expert working on simulation research project. Your priorities:
1. Efficient code that is easy to understand and maintain
2. System-levels thinking over quick solutions
3. Help develop deep insight into models' behaviors

## Project Overview

This project develops theory and numerical tools for **shift-representability**: given a finite family of continuous random variables $X_1, \dots, X_k$, does there exist a single strictly monotone transformation $g$ such that the transformed variables $g(X_n)$ form a pure location (shift) family? The question arises in signal detection theory and measurement theory, where shift-representability implies an additive conjoint measurement structure on the latent scale.

**Core mathematical framework.** Fix $X_1$ as reference and define quantile-matching transport maps $\psi_n = F_1^{-1} \circ F_n$. Shift-representability is equivalent to solving the simultaneous Abel equation $g(\psi_n(x)) = g(x) + c_n$ for a common increasing $g$. Theorem 2 gives necessary conditions (pairwise commutativity and non-crossing of transport maps) and a sufficient condition (disjointness of the generated group with $L(G) = \mathbb{R}$).

**Computational goal.** Build an extensible R toolkit that, given a family of distributions (specified analytically via CDF/quantile functions, or via simulation), numerically evaluates these conditions—commutativity, non-crossing, and disjointness—either confirming shift-representability or quantifying the degree of violation. The system should support both interactive single-family checks and batch sweeps over parameter grids.

**Status.** Early stage — foundational code in `R/transport_maps.R`; theory documented in `meta/`.

## Key Commands

```bash
# Environment setup
Rscript -e 'renv::restore()'   # Restore all package dependencies

# Quarto website
quarto preview                   # Live preview of the website
quarto render                    # Render full website to docs/
```

## Code Style

### ✅ Always
- Write self-documenting code with clear, descriptive names
- Refactor unclear code rather than adding explanatory comments
- Extract complex logic into well-named functions
- Use early returns to simplify conditional flow
- Use implicit returns at the end of functions (no explicit `return()`)
- Prefer vector operations over loops
- Use functional programming where possible
- Use `with_cache()` wrapper for long computations in R
- In notebooks, label code chunks via `#| label:` (do not put labels in the chunk header)

Example:

```r
recovery_est_bayes <- with_cache(
  "output/algo5_recovery_bayes.rds",
  purrr::pmap_dfr(recovery_true, fit_one_bayes, .progress = TRUE)
)
```


### 🚫 Never
- Add comments explaining what code does (refactor instead)

## Boundaries

### ✅ Always
- Write computation artifacts to `output`
- Source shared functions from `R/` rather than redefining them inline in notebooks
- Use project-root-relative paths (Quarto will run with `execute-dir: project`)
- Run `renv::restore()` before attempting to install packages

### 🚫 Never
- Commit secrets or `.env` files
- Modify content within `<!-- BEGIN USER-SPECIFIED -->` blocks
- Modify files in `data-raw`, `archive`, `_freeze`, `docs`, `renv`
- Force push to main
- Use `setwd()` or `rm(list=ls())` in function files

- Quarto notebooks are the primary analysis documents; they `source()` from `R/`

## User-Specified Content

Use `<!-- BEGIN USER-SPECIFIED -->` and `<!-- END USER-SPECIFIED -->` to mark
sections that AI agents must not modify or contradict. Useful for protecting
project-specific decisions that might conflict with general best practices

## Directory Structure

```
R/                   # shared functions (source from notebooks); no computations
notebooks/           # Quarto notebooks reporting computational results
output/              # computed artifacts (.csv, .rds, .h5) — write here
meta/                # project logs, planning, other admin
img/                 # static images
docs/                # AUTO-GENERATED quarto site — do not edit
_freeze/             # AUTO-GENERATED quarto cache — do not edit
_quarto.yml          # Quarto website project configuration
```
