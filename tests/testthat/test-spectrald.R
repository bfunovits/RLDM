# test-spectrald #########################################
# test freqresp and spectral density for MA(q) processes

freqresp_R = function(b, n.f = 128) {
  m = dim(b)[1]
  n = dim(b)[2]
  q = dim(b)[3] -1

  f = (0:(n.f-1))/n.f
  z = exp(complex(imaginary = -2*pi*f))
  frr = array(complex(real = numeric(m*n*n.f)), dim = c(m,n,n.f))

  for (i in (1:n.f)) {
    for (j in (0:q)) {
      frr[,,i] = frr[,,i] + b[,,j+1] * (z[i]^j)
    }
  }


  attr(frr, 'z') = z
  return(frr)
}

autocov_R = function(b, sigma_L, lag.max = (dim(b)[3] -1) ) {
  m = dim(b)[1]
  n = dim(b)[2]
  q = dim(b)[3] -1
  for ( i in (0:q) ) b[,,i + 1] = matrix(b[,,i+1], nrow = m, ncol = n) %*% sigma_L

  acf = array(0, dim = c(m,m,lag.max +1))
  for (i in (0:lag.max)) {
    for (j in iseq(0, q-i)) {
      acf[,,i+1] = acf[,,i+1] + matrix(b[,,i+j+1], nrow = m, ncol = n) %*% t(matrix(b[,,j+1], nrow = m, ncol = n))
    }
  }

  return(acf)
}

spectrald_R = function(acf, n.f) {
  m = dim(acf)[1]
  lag.max = dim(acf)[3] - 1
  acf = acf / (2*pi)

  f = (0:(n.f-1))/n.f
  z = exp(complex(imaginary = -2*pi*f))
  spd = array(complex(real = numeric(m*m*n.f)), dim = c(m,m,n.f))

  for (i in (1:n.f)) {
    spd[,,i] = acf[,,1]
    for (j in iseq(1, lag.max)) {
      spd[,,i] = spd[,,i] + acf[,,j+1] * (z[i]^j) + t(acf[,,j+1]) * Conj((z[i]^j))
    }
  }


  attr(spd, 'z') = z
  return(spd)

}

set.seed(12345)

test_that("freqresp methods work", {
  for (irun in (1:20)) {
    m = sample(1:5, 1)
    q = sample(1:5, 1)
    n.f = max(q, sample(1:256, 1))
    n.obs = sample((1:10)*5, 1)

    model = test_armamod(dim = c(m,m), degrees = c(0,q))
    y = sim(model, n.obs = n.obs)$y

    b = unclass(model$sys$b)
    sigma_L = model$sigma_L

    frr = unclass(freqresp(model, n.f = n.f)$frr)
    frr_R = freqresp_R(b, n.f = n.f)
    expect_equal(frr, frr_R)

    acf = unclass(autocov(model, lag.max = q)$acf)
    acf_R = autocov_R(b, sigma_L, lag.max = q)
    expect_equal(acf, acf_R)

    spd = unclass(spectrald(model, n.f = n.f)$spd)
    spd_R = spectrald_R(acf_R, n.f = n.f)
    expect_equal(spd, spd_R)

    sacf = autocov(y, lag.max = n.obs - 1)

    per = spectrald(sacf, n.f = n.obs)
    per2 = spectrald(y, n.f = n.obs)
    expect_equal(per, per2)

    gamma0 = apply(unclass(per$spd), MARGIN = c(1,2), FUN = mean)*(2*pi)
    expect_equal(max(abs(Im(gamma0))), 0)
    gamma0 = Re(gamma0)
    expect_equal(matrix(gamma0, nrow = m, ncol = m),
                 matrix(sacf$gamma[,,1], nrow = m, ncol = m))
  }
})

set.seed(NULL)
