test_that("riccati works", {
  for (isim in (1:10)) {
    m = sample(1:3, 1) # number of outputs
    s = sample(1:5, 1) # number of states

    # generate random model ("full" parametrization, but not necessarily miniphase)
    model = r_model(tmpl_stsp_full(m, m, s), bpoles = 1.1, sd = 0.25)
    # scale sigma
    model$sigma_L = model$sigma_L / sqrt(sum(diag(model$sigma_L %*% t(model$sigma_L))))

    A = model$sys$A
    B = model$sys$B
    C = model$sys$C
    sigma = model$sigma_L %*% t(model$sigma_L)

    # solve Lyapunov equation: P = A P A' + B sigma B'
    P = riccati(A, B %*% sigma, NULL, sigma, only.X = TRUE)
    # check that P solves the Lyapunov equation
    expect_equivalent(P - A %*% P %*% t(A), B %*% sigma %*% t(B))

    M = A %*% P %*% t(C) + B %*% sigma
    G = C %*% P %*% t(C) + sigma

    # check that P solves the Riccati equation
    expect_equivalent(P - A %*% P %*% t(A),
                      (M - A %*% P %*% t(C)) %*% solve(G - C %*% P %*% t(C)) %*% t(M - A %*% P %*% t(C)))

    # solve Riccati equation: P = APA' + (M -APC')(G-CPC')^{-1}(M-APC')'
    out = riccati(A, M, C, G, only.X = FALSE)

    P2 = out$X
    # check that P2 solves the Riccati equation
    expect_equivalent(P2 - A %*% P2 %*% t(A),
                      (M - A %*% P2 %*% t(C)) %*% solve(G - C %*% P2 %*% t(C)) %*% t(M - A %*% P2 %*% t(C)))

    # eigenvalues of (A-BC)
    lambda = eigen(A - B %*% C, only.values=TRUE)$values

    if (max(abs(lambda))<1) {
      # if the system is miniphase, then P = P2 must hold
      expect_equivalent(P,P2)
    }

    # construct a miniphase model
    sigma2 = out$sigma # G - C %*% P2 %*% t(C)
    B2 = out$B         # (M - A %*% P2 %*% t(C)) %*% solve(sigma2)

    # check that P2 solves the Lyapunov equation P2 = A P2 A' + B2 sigma2 B2'
    expect_equivalent(P2 - A %*% P2 %*% t(A), B2 %*% sigma2 %*% t(B2))

    # model (A, B2, C, sigma2) is miniphase!
    lambda2 = eigen(A - B2 %*% C, only.values = TRUE)$values
    expect_lt(max(abs(lambda2)),1)

    }

  })
