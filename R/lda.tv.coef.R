#' @title Longitudinal data analysis
#'
#' @description This function provide regression analysis of mixed sparse synchronous and asynchronous longitudinal covariates with time-varying coefficients.
#'
#' @param  data_res  An object of class tibble. The structure of the tibble must be: tibble(id_y=ID, ty=measurement time for response, y=observation for response, x=matrix(observation for synchronous covariates), x_add=matrix(observation for uninterested synchronous covariates)).
#'
#' @param  data_cov  An object of class tibble. The structure of the tibble must be: tibble(id_z=ID, tz=measurement time for response, z=matrix(observation for asynchronous covariates)).
#'
#' @param  time An object of class vector. The interest time.
#'
#' @param  N An object of class integer. The sample size.
#'
#' @param  bd An object of class vector. If use auto bandwidth selection, the structure of the vector must be: bd=c(the maximum bandwidth for h1, the minimum bandwidth for h1, the maximum bandwidth for h2, the minimum bandwidth for h2, the fold of cross-validation, the number of bandwidth divided). If use fixed bandwidth, bd=c(the chosen bandwidth).
#'
#' @param  method An object of class integer indicating the method used to do estimation for asynchronous covariates. If use one-stage method, method=1; if use two-stage method with centering method for the first stage, method=1; if use two-stage method with time-varying method for the first stage, method=2.
#'
#' @param  scb An object of class vector. If need to construct the simultaneous confidence band, the structure of the vector must be: c(alpha=desirable confidence level, B=repeat times). Otherwise, scb=0.
#'
#' @return a list with the following elements:
#' \item{est.b}{The estimation for the parameter of synchronous covariates.}
#' \item{est.g}{The estimation for the parameter of asynchronous covariates.}
#' \item{se.b}{The estimation of standard error for the parameter of synchronous covariates.}
#' \item{se.g}{The estimation of standard error for the parameter of asynchronous covariates.}
#' \item{c_alpha_x}{The empirical percentile used to construct the simultaneous confidence band for the parameter of synchronous covariates.}
#' \item{c_alpha_z}{The empirical percentile used to construct the simultaneous confidence band for the parameter of asynchronous covariates.}
#'
#' @import  dplyr  tibble
#' @importFrom MASS mvrnorm
#' @importFrom  stats quantile runif
#' @examples
#'
#'library(dplyr)
#'library(MASS)
#'library(tibble)

#'N=400
#'ty=tz=y=x=x1=z=id_y=id_z=list()
#'beta<-function(t){
#'  0.3*(t-0.4)^2
#'}
#'gamma<-function(t){
#'  sin(2*pi*t)
#'}
#'ny=rpois(N,5)+1
#'nz=rpois(N,5)+1
#'for(i in 1:N){
#'  ty[[i]]=as.matrix(runif(ny[i]))
#'  tz[[i]]=as.matrix(runif(nz[i]))
#'  t.temp=rbind(tz[[i]],ty[[i]])
#'  n.temp=nz[i]+ny[i]
#'  corr=exp(-abs(rep(1,n.temp)%*%t(t.temp)-t.temp%*%t(rep(1,n.temp))))
#'  corr.e=2^(-abs(rep(1,n.temp)%*%t(t.temp)-t.temp%*%t(rep(1,n.temp))))
#'  MX=rep(0, n.temp)
#'  MZ= 2*(t.temp-0.5)^2
#'  x.temp=mvrnorm(1,MX,corr)
#'  z.temp=mvrnorm(1,MZ, corr)
#'  z[[i]]=as.matrix(z.temp[1:nz[i]])
#'  x[[i]]=as.matrix(x.temp[-(1:nz[i])])
#'  id_z[[i]]=rep(i,nz[i])
#'  id_y[[i]]=rep(i,ny[i])
#'  y.temp=gamma(t.temp)*z.temp+beta(t.temp)*x.temp+as.matrix(mvrnorm(1,rep(0,n.temp),corr.e))
#'  y[[i]]=as.matrix(y.temp[-(1:nz[i])])
#'}
#' data_cov=tibble(id_z=unlist(id_z),tz=unlist(tz),z=matrix(unlist(z),length(unlist(z))))
#' data_res=tibble(id_y=unlist(id_y),ty=unlist(ty),x=matrix( unlist(x),length(unlist(x))), y=unlist(y))
#' ldatv(data_res,data_cov,time=0.3,N,bd=c(N^(-0.5),N^(-0.5)),method=1,scb=0)
#' @export




ldatv <- function(data_res, data_cov, time, N, bd, method, scb) {
  K1 = function(x, h) {
    0.75 * (1 - (x / h) ^ 2) / h * (abs(x) < h)
  }

  dim_x = dim(data_res$x)[2]
  dim_z = dim(data_cov$z)[2]
  id_y=data_res$id_y
  id_z=data_cov$id_z
  nt = length(time)

  Qn.x = Qn.z = 0
  c_alpha.x =  matrix(0, dim_z )
  c_alpha.z =  matrix(0, dim_z )
  if (length(scb) == 2) {
    B = scb[2]
    alpha = scb[1]
    c_alpha.x = matrix(0, dim_x, B)
    c_alpha.z = matrix(0, dim_z, B)

  }

  if (length(bd) == 6) {
    cv_num = bd[5]
    bd_num = bd[6]
    h1.can = seq(bd[1], bd[2], length = bd_num)
    h2.can = seq(bd[3], bd[4], length = bd_num)
    gr = sample(1:cv_num, N, replace = T)


    if (method == 1) {
      mse = matrix(nrow = bd_num, ncol = bd_num)
      Mse = 0
      for (i in 1:nt) {
        s = time[i]
        for (hh1 in 1:bd_num) {
          for (hh2 in 1:bd_num) {
            h1.temp = h1.can[hh1]
            h2.temp = h2.can[hh2]
            aspe = aspe1 = 0
            for (cv in 1:cv_num) {
              cv_class = which(gr != cv)
              n_cv = length(cv_class)
              cv_res = data_res %>% group_by(id_y) %>% filter(id_y %in% cv_class)
              cv_cov = data_cov %>% group_by(id_z) %>% filter(id_z %in% cv_class)

              x = y = z = ty = tz = list()
              ny = nz = rep(0, n_cv)
              for (i in 1:n_cv) {
                x[[i]] = matrix(cv_res$x[which(cv_res$id_y == cv_class[i]), ], ncol = dim_x)
                y[[i]] = cv_res$y[which(cv_res$id_y == cv_class[i])]
                z[[i]] = matrix(cv_cov$z[which(cv_cov$id_z == cv_class[i]), ], ncol =
                                  dim_z)
                ty[[i]] = cv_res$ty[which(cv_res$id_y == cv_class[i])]
                tz[[i]] = cv_cov$tz[which(cv_cov$id_z == cv_class[i])]
                ny[i] = length(which(cv_res$id_y == cv_class[i]))
                nz[i] = length(which(cv_cov$id_z == cv_class[i]))
              }



              U.dot = 0
              U = 0
              for (i in 1:n_cv) {
                for (j in 1:length(ty[[i]])) {
                  for (k in 1:length(tz[[i]])) {
                    temp = c(x[[i]][j, ], x[[i]][j, ] * (ty[[i]][j] - s), z[[i]][k, ], z[[i]][k, ] *
                               (tz[[i]][k] - s))
                    k1 = K1(tz[[i]][k] - s, h2.temp) * K1(ty[[i]][j] - s, h1.temp)
                    U = U + k1 * y[[i]][j] * temp
                    U.dot = U.dot + k1 * temp %*% t(temp)
                  }
                }
              }
              est = solve(U.dot) %*% U

              rsd = fn = 0
              for (i in 1:n_cv) {
                for (j in 1:length(ty[[i]])) {
                  for (k in 1:length(tz[[i]])) {
                    temp = c(x[[i]][j, ], x[[i]][j, ] * (ty[[i]][j] - s), z[[i]][k, ], z[[i]][k, ] *
                               (tz[[i]][k] - s))
                    k1 = K1(tz[[i]][k] - s, h2.temp) * K1(ty[[i]][j] - s, h1.temp)
                    rsd = rsd + k1 * (y[[i]][j] - temp %*% est) ^ 2
                    fn = fn + k1
                  }
                }
              }
              aspe = aspe + rsd / fn


            }
            mse[hh1, hh2] = aspe
          }
        }
        Mse = Mse + mse
      }

      ind = as.numeric(which(Mse == min(Mse), arr.ind = TRUE))
      h1 = h1.can[ind[1]]
      h2 = h2.can[ind[2]]


    }

    if (method > 1) {
      Mse = mse = rep(0, bd_num)
      for (m in 1:nt) {
        s = time[m]
        for (hh1 in 1:bd_num) {
          h1.temp = h1.can[hh1]
          aspe = 0
          for (cv in 1:cv_num) {
            cv_class = which(gr != cv)
            n_cv = length(cv_class)
            cv_res = data_res %>% group_by(id_y) %>% filter(id_y %in% cv_class)
            cv_cov = data_cov %>% group_by(id_z) %>% filter(id_z %in% cv_class)

            x = y = z = ty = list()
            ny = rep(0, n_cv)
            for (i in 1:n_cv) {
              x[[i]] = matrix(cv_res$x[which(cv_res$id_y == cv_class[i]), ], ncol = dim_x)
              y[[i]] = cv_res$y[which(cv_res$id_y == cv_class[i])]
              ty[[i]] = cv_res$ty[which(cv_res$id_y == cv_class[i])]
              ny[i] = length(y[[i]])
            }

            x.center = x
            y.center = y

            x.uls = cv_res$x
            t.uls = cv_res$ty
            y.uls = cv_res$y
            Ex = list()
            Ey = list()
            for (i in 1:n_cv) {
              Ex[[i]] = x[[i]]
              Ey[[i]] = y[[i]]

              for (k in 1:ny[i]) {
                temp = K1(ty[[i]][k] - t.uls, h1.temp) / sum(K1(ty[[i]][k] - t.uls, h1.temp))
                Ex[[i]][k, ] = t(x.uls) %*% temp
                Ey[[i]][k] = t(y.uls) %*% temp
              }
              x.center[[i]] = x[[i]] - Ex[[i]]
              y.center[[i]] = y[[i]] - Ey[[i]]
            }


            X = NULL
            U.dot = U = 0
            for (i in 1:n_cv) {
              for (j in 1:length(ty[[i]])) {
                temp = c(x.center[[i]][j, ], x.center[[i]][j, ] * (ty[[i]][j] - s))
                k1 = K1(ty[[i]][j] - s, h1.temp)
                U = U + k1 * temp * y.center[[i]][j]
                U.dot = U.dot + k1 * temp %*% t(temp)

              }
            }
            est.ts = solve(U.dot) %*% U

            rsd = fn = 0
            for (i in 1:n_cv) {
              for (j in 1:length(ty[[i]])) {
                k1 = K1(ty[[i]][j] - s, h1.temp)
                rsd = rsd + k1 * (y.center[[i]][j] - x.center[[i]][j, ] %*%
                                    est.ts[c(1:dim_x)]) ^ 2
                fn = fn + k1
              }
            }
            aspe = aspe + rsd / fn

          }
          mse[hh1] = aspe
        }
        Mse = Mse + mse
      }
      h1 = h1.can[which(Mse == min(Mse))]

      Mse = mse = rep(0, bd_num)
      for (m in 1:nt) {
        s = time[m]

        x = y = z = ty = tz = list()
        ny = rep(0, N)
        for (i in 1:N) {
          x[[i]] = matrix(data_res$x[which(data_res$id_y == i), ], ncol = dim_x)
          y[[i]] = data_res$y[which(data_res$id_y == i)]
          ty[[i]] = data_res$ty[which(data_res$id_y == i)]
          ny[i] = length(y[[i]])
        }

        x.center = x
        y.center = y

        x.uls = data_res$x
        t.uls = data_res$ty
        y.uls = data_res$y
        Ex = list()
        Ey = list()
        for (i in 1:N) {
          Ex[[i]] = x[[i]]
          Ey[[i]] = y[[i]]

          for (k in 1:ny[i]) {
            temp = K1(ty[[i]][k] - t.uls, h1) / sum(K1(ty[[i]][k] - t.uls, h1))
            Ex[[i]][k, ] = t(x.uls) %*% temp
            Ey[[i]][k] = t(y.uls) %*% temp
          }
          x.center[[i]] = x[[i]] - Ex[[i]]
          y.center[[i]] = y[[i]] - Ey[[i]]
        }

        U.dot = U = 0
        for (i in 1:N) {
          for (j in 1:length(ty[[i]])) {
            temp = c(x.center[[i]][j, ], x.center[[i]][j, ] * (ty[[i]][j] - s))
            k1 = K1(ty[[i]][j] - s, h1)
            U = U + k1 * temp * y.center[[i]][j]
            U.dot = U.dot + k1 * temp %*% t(temp)

          }
        }
        est.full = solve(U.dot) %*% U
        est.full.b = est.full[c(1:dim_x)]



        for (hh2 in 1:bd_num) {
          h2.temp = h2.can[hh2]
          aspe = 0
          for (cv in 1:cv_num) {
            cv_class = which(gr != cv)
            n_cv = length(cv_class)
            cv_res = data_res %>% group_by(id_y) %>% filter(id_y %in% cv_class)
            cv_cov = data_cov %>% group_by(id_z) %>% filter(id_z %in% cv_class)

            x = y = z = ty = tz = list()
            ny = nz = rep(0, n_cv)
            for (i in 1:n_cv) {
              x[[i]] = matrix(cv_res$x[which(cv_res$id_y == cv_class[i]), ], ncol = dim_x)
              y[[i]] = cv_res$y[which(cv_res$id_y == cv_class[i])]
              z[[i]] = matrix(cv_cov$z[which(cv_cov$id_z == cv_class[i]), ], ncol =
                                dim_z)
              ty[[i]] = cv_res$ty[which(cv_res$id_y == cv_class[i])]
              tz[[i]] = cv_cov$tz[which(cv_cov$id_z == cv_class[i])]
            }

            x.center = x
            y.center = y

            x.uls = cv_res$x
            t.uls = cv_res$ty
            y.uls = cv_res$y


            U.dot = 0
            U = 0
            for (i in 1:n_cv) {
              for (j in 1:length(ty[[i]])) {
                for (k in 1:length(tz[[i]])) {
                  temp = c(x[[i]][j, ], x[[i]][j, ] * (ty[[i]][j] - s))
                  temp1 = c(z[[i]][k, ], z[[i]][k, ] * (tz[[i]][k] - s))
                  k1 = K1(tz[[i]][k] - s, h2.temp) * K1(ty[[i]][j] - s, h1)
                  U = U + k1 * c(y[[i]][j] - temp %*% est.full) * temp1
                  U.dot = U.dot + k1 * temp1 %*% t(temp1)
                }
              }
            }
            est.g = solve(U.dot) %*% U


            rsd = fn = 0
            for (i in 1:n_cv) {
              for (j in 1:length(ty[[i]])) {
                for (k in 1:length(tz[[i]])) {
                  k1 = K1(ty[[i]][j] - s, h1) * K1(tz[[i]][k] - s, h2.temp)
                  rsd = rsd + k1 * c(y[[i]][j] - x[[i]][j, ] %*% est.full.b -
                                       z[[i]][k, ] %*% est.g[c(1:dim_z)]) ^ 2
                  fn = fn + k1
                }
              }
            }


            aspe = aspe + rsd / fn
          }
          mse[hh2] = aspe
        }
        Mse = Mse + mse
      }
      h2 = h2.can[which(Mse == min(Mse))]
    }


  } else{
    h1 = bd[1]
    h2 = bd[2]
  }



  x = y = z = ty = tz = list()
  ny = nz = rep(0, N)
  for (i in 1:N) {
    x[[i]] = matrix(data_res$x[which(data_res$id_y == i), ], ncol = dim_x)
    y[[i]] = data_res$y[which(data_res$id_y == i)]
    z[[i]] = matrix(data_cov$z[which(data_cov$id_z == i), ], ncol = dim_z)
    ty[[i]] = data_res$ty[which(data_res$id_y == i)]
    tz[[i]] = data_cov$tz[which(data_cov$id_z == i)]
    ny[i] = length(which(data_res$id_y == i))
    nz[i] = length(which(data_cov$id_z == i))
  }

  est.b = est.g = se.b = se.g = NULL


  for (kk in 1:nt) {
    s = time[kk]

    if (method == 1) {
      U.dot = 0
      U = 0
      for (i in 1:N) {
        for (j in 1:ny[i]) {
          for (k in 1:nz[i]) {
            temp = c(x[[i]][j, ], x[[i]][j, ] * (ty[[i]][j] - s), z[[i]][k, ], z[[i]][k, ] *
                       (tz[[i]][k] - s))
            k1 = K1(ty[[i]][j] - s, h1) * K1(tz[[i]][k] - s, h2)
            U = U + k1 * y[[i]][j] * temp
            U.dot = U.dot + k1 * temp %*% t(temp)
          }
        }
      }
      est.os = solve(U.dot) %*% U
      b.hat.os = est.os[(1:dim_x)]
      g.hat.os = est.os[(2 * dim_x + 1):(2 * dim_x + dim_z)]


      U.sq = 0
      for (i in 1:N) {
        temp.U = 0
        for (j in 1:length(ty[[i]])) {
          for (k in 1:length(tz[[i]])) {
            temp = c(x[[i]][j, ], x[[i]][j, ] * (ty[[i]][j] - s), z[[i]][k, ], z[[i]][k, ] *
                       (tz[[i]][k] - s))
            k1 = K1(ty[[i]][j] - s, h1) * K1(tz[[i]][k] - s, h2)
            temp.U = temp.U + k1 * temp * as.numeric(y[[i]][j] - temp %*% est.os)
          }
        }
        U.sq = U.sq + temp.U %*% t(temp.U)
      }
      se.os = diag(solve(U.dot) %*% U.sq %*% solve(U.dot)) ^ .5
      b.se.os = se.os[(1:dim_x)]
      g.se.os = se.os[(2 * dim_x + 1):(2 * dim_x + dim_z)]


      rsd_x = rsd_z = fn_x = fn_z = 0
      if (length(scb) == 2) {
        mtp = matrix(ifelse(runif(N * B) > 0.5, -1, 1), B)
        for (i in 1:N) {
          for (j in 1:ny[i]) {
            for (k in 1:nz[i]) {
              k1 = K1(tz[[i]][k] - s, h1) * K1(ty[[i]][j] - s, h2)
              temp = c(x[[i]][j, ], x[[i]][j, ] * (ty[[i]][j] - s), z[[i]][k, ], z[[i]][k, ] *
                         (tz[[i]][k] - s))
              rsd_x = rsd_x + matrix(k1 * x[[i]][j, ] * as.numeric(y[[i]][j] -
                                                                     temp %*% est.os)) %*% mtp[, i]
              rsd_z = rsd_z + matrix(k1 * z[[i]][k, ] * as.numeric(y[[i]][j] -
                                                                     temp %*% est.os)) %*% mtp[, i]
              fn_x = fn_x + k1 * x[[i]][j, ] %*% t(x[[i]][j, ])
              fn_z = fn_z + k1 * z[[i]][k, ] %*% t(z[[i]][k, ])
            }
          }
        }
        Qn.x = abs(solve(fn_x) %*% rsd_x)
        Qn.z = abs(solve(fn_z) %*% rsd_z)
      }

      est.b = cbind(est.b, b.hat.os)
      est.g = cbind(est.g, g.hat.os)
      se.b = cbind(se.b, b.se.os)
      se.g = cbind(se.g, g.se.os)
    }

    ############### two-step ###########################


    if (method > 1) {
      if (method == 2) {
        x.center = x
        y.center = y

        x.uls = data_res$x
        t.uls = data_res$ty
        y.uls = data_res$y
        Ex = list()
        Ey = list()
        for (i in 1:N) {
          Ex[[i]] = x[[i]]
          Ey[[i]] = y[[i]]

          for (k in 1:ny[i]) {
            temp = K1(ty[[i]][k] - t.uls, h1) / sum(K1(ty[[i]][k] - t.uls, h1))
            Ex[[i]][k, ] = t(x.uls) %*% temp
            Ey[[i]][k] = t(y.uls) %*% temp

          }
          x.center[[i]] = x[[i]] - Ex[[i]]
          y.center[[i]] = y[[i]] - Ey[[i]]
        }


        X = NULL
        U.dot = U = 0
        for (i in 1:N) {
          for (j in 1:length(ty[[i]])) {
            temp = c(x.center[[i]][j, ], x.center[[i]][j, ] * (ty[[i]][j] - s))
            k1 = K1(ty[[i]][j] - s, h1)
            U = U + k1 * temp * y.center[[i]][j]
            U.dot = U.dot + k1 * temp %*% t(temp)
          }
        }
        est.1 = solve(U.dot) %*% U
        b.hat.ts = est.1[c(1:dim_x)]

        U.ts.b = 0
        for (i in 1:N) {
          sum.b = 0
          for (j in 1:length(ty[[i]])) {
            k1 = K1(ty[[i]][j] - s, h1)
            temp = c(x.center[[i]][j, ], x.center[[i]][j, ] * (ty[[i]][j] -
                                                                 s))
            sum.b = sum.b + k1 * (y.center[[i]][j] - temp %*% est.1) %*% temp
          }
          U.ts.b = U.ts.b + t(sum.b) %*% sum.b
        }
        sd.hat = (diag(solve(U.dot) %*% U.ts.b %*% solve(U.dot))) ^ .5
        b.sd.ts = sd.hat[c(1:dim_x)]


        est.b = cbind(est.b, b.hat.ts)
        se.b = cbind(se.b, b.sd.ts)
      }

      if (method == 3) {
        X = NULL
        U.dot = U = 0
        for (i in 1:N) {
          for (j in 1:length(ty[[i]])) {
            temp = c(1, (ty[[i]][j] - s), x[[i]][j, ], x[[i]][j, ] * (ty[[i]][j] - s))
            k1 = K1(ty[[i]][j] - s, h1)
            U = U + k1 * temp * y[[i]][j]
            U.dot = U.dot + k1 * temp %*% t(temp)

          }
        }
        est.1 = solve(U.dot) %*% U
        b.hat.ts = est.1[c(3:(dim_x + 2))]

        U.ts.b = ifelse(length(beta) == 1, 0, diag(rep(0, length(beta))))
        for (i in 1:N) {
          UU = 0
          for (j in 1:length(ty[[i]])) {
            k1 = K1(ty[[i]][j] - s, h1)
            temp = c(1, (ty[[i]][j] - s), x[[i]][j, ], x[[i]][j, ] * (ty[[i]][j] -
                                                                        s))
            UU = UU +  k1 * c(y[[i]][j] - temp %*% est.1) %*% temp
          }
          U.ts.b = U.ts.b + t(UU) %*% UU
        }
        se.hat = diag(solve(U.dot) %*% U.ts.b %*% solve(U.dot)) ^ 0.5
        b.sd.ts = se.hat[c(3:(dim_x + 2))]


        est.b = cbind(est.b, b.hat.ts)
        se.b = cbind(se.b, b.sd.ts)
      }



      ############estimate gamma#########################

      U = U.dot = 0
      for (i in 1:N) {
        for (j in 1:length(ty[[i]])) {
          for (k in 1:length(tz[[i]])) {
            temp = c(z[[i]][k, ], z[[i]][k, ] * (tz[[i]][k] - s))
            k1 = K1(tz[[i]][k] - s, h2) * K1(ty[[i]][j] - s, h1)
            U = U + k1 * temp %*% (y[[i]][j] - x[[i]][j, ] %*% c(b.hat.ts))
            U.dot = U.dot + k1 * temp %*% t(temp)
          }
        }
      }
      est.2 = solve(U.dot) %*% U
      g.hat.ts = est.2[c(1:dim_z)]

      U.ts.g = 0
      for (i in 1:N) {
        UU = 0
        for (j in 1:length(ty[[i]])) {
          for (k in 1:length(tz[[i]])) {
            temp = c(z[[i]][k, ], z[[i]][k, ] * (tz[[i]][k] - s))
            k1 = K1(tz[[i]][k] - s, h2) %*% t(K1(ty[[i]][j] - s, h1))
            UU = UU + temp %*% (k1 * (y[[i]][j] - temp %*% c(est.2) - x[[i]][j, ] %*%
                                        c(b.hat.ts)))
          }
        }
        U.ts.g = U.ts.g + UU %*% t(UU)
      }

      g.se.2sg = (diag(solve(U.dot) %*% U.ts.g %*% solve(U.dot))) ^ .5
      g.se.2sg = g.se.2sg[c(1:dim_z)]
      est.g = cbind(est.g, g.hat.ts)
      se.g = cbind(se.g, g.se.2sg)

      est = c(b.hat.ts, est.2)

      if ((method == 2) & (length(scb) == 2)) {
        rsd_x = rsd_z = fn_x = fn_z = 0
        mtp = matrix(ifelse(runif(N * B) > 0.5, -1, 1), B)
        for (i in 1:N) {
          for (j in 1:ny[i]) {
            for (k in 1:nz[i]) {
              k1 = K1(tz[[i]][k] - s, h1) * K1(ty[[i]][j] - s, h2)
              temp = c(x.center[[i]][j, ], x.center[[i]][j, ] * (ty[[i]][j] -
                                                                   s))
              temp1 = c(x[[i]][j, ], z[[i]][k, ], z[[i]][k, ] * (tz[[i]][k] -
                                                                   s))
              rsd_x = rsd_x + matrix(k1 * x.center[[i]][j, ] * as.numeric(y.center[[i]][j] -
                                                                            temp %*% est.1)) %*% mtp[, i]
              rsd_z = rsd_z + matrix(k1 * z[[i]][k, ] * as.numeric(y[[i]][j] -
                                                                     temp1 %*% est)) %*% mtp[, i]
              fn_x = fn_x + k1 * x.center[[i]][j, ] %*% t(x.center[[i]][j, ])
              fn_z = fn_z + k1 * z[[i]][k, ] %*% t(z[[i]][k, ])
            }
          }
        }
        Qn.x = abs(solve(fn_x) %*% rsd_x)
        Qn.z = abs(solve(fn_z) %*% rsd_z)

      }



      if ((method == 3) & (length(scb) == 2)) {
        rsd_x = rsd_z = fn_x = fn_z = 0
        mtp = matrix(ifelse(runif(N * B) > 0.5, -1, 1), B)
        for (i in 1:N) {
          for (j in 1:ny[i]) {
            for (k in 1:nz[i]) {
              k1 = K1(tz[[i]][k] - s, h1) * K1(ty[[i]][j] - s, h2)
              temp = c(1, (ty[[i]][j] - s), x[[i]][j, ], x[[i]][j, ] * (ty[[i]][j] -
                                                                          s))
              temp1 = c(x[[i]][j, ], z[[i]][k, ], z[[i]][k, ] * (tz[[i]][k] -
                                                                   s))
              rsd_x = rsd_x + matrix(k1 * x[[i]][j, ] * as.numeric(y[[i]][j] -
                                                                     temp %*% est.1)) %*% mtp[, i]
              rsd_z = rsd_z + matrix(k1 * z[[i]][k, ] * as.numeric(y[[i]][j] -
                                                                     temp1 %*% est)) %*% mtp[, i]
              fn_x = fn_x + k1 * x[[i]][j, ] %*% t(x[[i]][j, ])
              fn_z = fn_z + k1 * z[[i]][k, ] %*% t(z[[i]][k, ])
            }
          }
        }
        Qn.x = abs(solve(fn_x) %*% rsd_x)
        Qn.z = abs(solve(fn_z) %*% rsd_z)

      }

    }

    c_alpha.x = ifelse(Qn.x < c_alpha.x, c_alpha.x, Qn.x)
    c_alpha.z = ifelse(Qn.z < c_alpha.z, c_alpha.z, Qn.z)
  }

  c_x = c_z = 0
  if (length(scb) == 2) {
    c_x = apply(c_alpha.x, 1, quantile, probs = alpha)
    c_z = apply(c_alpha.z, 1, quantile, probs = alpha)
  }

  return(
    list(
      est.b = est.b,
      est.g = est.g,
      se.b = se.b,
      se.g = se.g ,
      c_alpha_x = c_x,
      c_alpha_z = c_z
    )
  )



}
