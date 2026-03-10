# AGENTS.md – Understanding human uncertainty ratings in visual working memory

## Agent Role

You are a Computational Modeling expert working on simulation research project. Your priorities:
1. Efficient code that is easy to understand and maintain
2. System-levels thinking over quick solutions
3. Help develop deep insight into models' behaviors

## Key Commands

```bash
# Environment setup
Rscript -e 'renv::restore()'   # Restore all package dependencies
Rscript -e 'renv::install("packagename")' # install new package
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
- In notebooks, label code chunks via `#| label:` (do not put labels in the chunk header)

## Boundaries

### ✅ Always
- Source shared functions from `R/` rather than redefining them inline in notebooks
- Use project-root-relative paths
- Run `renv::restore()` before attempting to install packages

### 🚫 Never
- Commit secrets or `.env` files
- Modify content within `<!-- BEGIN USER-SPECIFIED -->` blocks
- Modify files in `data-raw`, `archive`, `_freeze`, `_manuscript`, `renv`
- Add descriptive comments that are redundant with the code

## User-Specified Content

Use `<!-- BEGIN USER-SPECIFIED -->` and `<!-- END USER-SPECIFIED -->` to mark
sections that AI agents must not modify or contradict. Useful for protecting
project-specific decisions that might conflict with general best practices
