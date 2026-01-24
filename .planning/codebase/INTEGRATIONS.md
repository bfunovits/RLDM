# External Integrations

**Analysis Date:** 2026-01-24

## APIs & External Services

**Package Management:**
- GitHub - Source code hosting and package distribution via `remotes::install_github()`
- CRAN - Standard R package repository for dependencies

**Documentation Hosting:**
- GitHub Pages - Package documentation website at https://bfunovits.github.io/RLDM/

## Data Storage

**Databases:**
- None - This is a statistical modeling package, not a database application
- Data is processed in memory using R data structures

**File Storage:**
- Local filesystem only - Package data stored in `data/` directory
- Example datasets: `BQdata` (Blanchard-Quah), `RSdata` (Ramey/Shapiro)

**Caching:**
- None - No external caching services detected

## Authentication & Identity

**Auth Provider:**
- None - No authentication required for package functionality
- GitHub authentication only for package installation via `remotes`

## Monitoring & Observability

**Error Tracking:**
- None - Standard R error handling via `tryCatch()` and `stop()`
- GitHub Issues for bug reporting (configured in `DESCRIPTION`)

**Logs:**
- Console output only - No external logging services
- Debugging via R's built-in debugging tools

## CI/CD & Deployment

**Hosting:**
- GitHub - Source code repository and issue tracking
- GitHub Pages - Documentation hosting

**CI Pipeline:**
- GitHub Actions - Single workflow at `.github/workflows/pkgdown.yaml`
  - Builds pkgdown documentation site on push to main/master branches
  - Deploys to GitHub Pages automatically
  - Uses r-lib/actions for R environment setup

**Package Distribution:**
- GitHub releases - Source package distribution
- No CRAN submission detected (version 0.0.0.9006 indicates development version)

## Environment Configuration

**Required env vars:**
- `GITHUB_PAT` - GitHub Personal Access Token for CI/CD (configured as GitHub secret)
- No other environment variables required for package functionality

**Secrets location:**
- GitHub repository secrets for CI/CD
- No local `.env` files detected

## Webhooks & Callbacks

**Incoming:**
- None - No webhook endpoints or API servers

**Outgoing:**
- None - No external API calls or webhook notifications

## Package Dependencies Integration

**Critical External Package:**
- `rationalmatrices` - Sister package from GitHub
  - Required for rational matrix computations
  - Linked via `LinkingTo:` in `DESCRIPTION`
  - Header-only C++ implementation included via `#include <rationalmatrices_lyapunov.h>`
  - Installed via `remotes::install_github("bfunovits/rationalmatrices")` in CI pipeline

**Build System Integration:**
- LAPACK/BLAS - System linear algebra libraries linked via `Makevars`
- OpenMP - Parallel computing framework for C++ code

## Data Source Integration

**Built-in Datasets:**
- `BQdata` - Blanchard-Quah economic dataset (stored in package)
- `RSdata` - Ramey/Shapiro government spending dataset (stored in package)
- No external data fetching or API calls for data retrieval

---

*Integration audit: 2026-01-24*