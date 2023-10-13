#' @title Longitudinal data analysis
#'
#' @description This function provide regression analysis of mixed sparse synchronous and asynchronous longitudinal covariates.
#'
#' @param  data_res  An object of class tibble. The structure of the tibble must be: tibble(id_y=ID, ty=measurement time for response, y=observation for response, x=matrix(observation for synchronous covariates), x_add=matrix(observation for uninterested synchronous covariates)).
#'
#' @param  data_cov  An object of class tibble. The structure of the tibble must be: tibble(id_z=ID, tz=measurement time for response, z=matrix(observation for asynchronous covariates)).
#'
#' @param  N An object of class integer. The sample size.
#'
#' @param  bd An object of class vector. If use auto bandwidth selection, the structure of the vector must be: d=c(the maximum bandwidth, the minimum bandwidth, the fold of cross-validation, the number of bandwidth divided). If use fixed bandwidth, bd=c(the chosen bandwidth).
#'
#' @param  omit An object of class integer indicating the method used to do estimation for synchronous covariates. If use plm method, omit=1; if use centering method, omit=2; if use additional covariates information, omit=3.
#'
#' @param  method An object of class integer indicating the method used to do estimation for asynchronous covariates. If only deal with omit variable, method=0; if use two-stage method, method=1; if use kernel smoothing, method=2.
#'
#' @return a list with the following elements:
#' \item{est}{The estimation for the corresponding parameters.}
#' \item{se}{The estimation of standard error for the estimated parameters.}
#'
#' @import  dplyr  tibble
#' @importFrom MASS mvrnorm
#' @importFrom  stats quantile runif
#' @importFrom  dlm bdiag
#' @examples
#'
#'library(MASS)
#'library(tibble)
#'library(dplyr)
#' N=100
#' ty=tz=y=x=z=id_y=id_z=list()
#' a=b=g=1
#' ny=rpois(N,5)+1
#' nz=rpois(N,5)+1
#' for(i in 1:N){
#'   ty[[i]]=as.matrix(runif(ny[i]))
#'   tz[[i]]=as.matrix(runif(nz[i]))
#'   t.temp=rbind(tz[[i]],ty[[i]])
#'   n.temp=nz[i]+ny[i]
#'   corr=exp(-abs(rep(1,n.temp)%*%t(t.temp)-t.temp%*%t(rep(1,n.temp))))
#'   corr.e=2^(-abs(rep(1,n.temp)%*%t(t.temp)-t.temp%*%t(rep(1,n.temp))))
#'   MX=t.temp^.5
#'   MZ=rep(0, n.temp)
#'   x.temp=mvrnorm(1,MX,corr)
#'   z.temp=mvrnorm(1,MZ, corr)
#'   z[[i]]=as.matrix(z.temp[1:nz[i]])
#'   x[[i]]=as.matrix(x.temp[-(1:nz[i])])
#'   id_z[[i]]=rep(i,nz[i])
#'   id_y[[i]]=rep(i,ny[i])
#'   y.temp=a+g*z.temp+x.temp*b+as.matrix(mvrnorm(1,rep(0,n.temp),corr.e))
#'   y[[i]]=as.matrix(y.temp[-(1:nz[i])])
#' }
#' data_cov=tibble(id_z=unlist(id_z),tz=unlist(tz),z=matrix(unlist(z),length(unlist(z))))
#' data_res=tibble(id_y=unlist(id_y),ty=unlist(ty),x=matrix(unlist(x),length(unlist(x))),y=unlist(y))
#' bd=0.1
#' omit=1
#' method=1
#' lda(data_res,data_cov,N,bd,omit,method)
#' @export




lda <- function(data_res, data_cov, N, bd, omit, method) {
  K1 = function(x, h) {
    0.75 * (1 - (x / h) ^ 2) / h * (abs(x) < h)
  }

  u = function(tx, ty, x, y, h, n) {
    U.dot = 0
    U = 0
    for (i in 1:n) {
      for (j in 1:length(ty[[i]])) {
        for (k in 1:length(tx[[i]])) {
          k1 = K1(tx[[i]][k] - ty[[i]][j], h)
          U = U + k1 * x[[i]][k, ] * y[[i]][j]
          U.dot = U.dot + k1 * x[[i]][k, ] %*% t(x[[i]][k, ])
        }
      }
    }
    return(list(U, U.dot))
  }

  dim_x = dim(data_res$x)[2]
  dim_z = dim(data_cov$z)[2]
  id_y=data_res$id_y
  id_z=data_cov$id_z


  if (length(bd) == 4) {
    cv_num = bd[3]
    bd_num = bd[4]
    h.can = seq(bd[1], bd[2], length = bd_num)
    gr = sample(1:cv_num, N, replace = T)
    mse = matrix(0, nrow = cv_num, ncol = bd_num)


    for (hh in 1:bd_num) {
      h.temp = h.can[hh]
      for (cv in 1:cv_num) {
        cv_class = which(gr != cv)
        cv_res = data_res %>% group_by(id_y) %>% filter(id_y %in% cv_class)
        cv_cov = data_cov %>% group_by(id_z) %>% filter(id_z %in% cv_class)
        x.uls = cv_res$x
        t.uls = cv_res$ty
        y.uls = cv_res$y
        # x0=x[gr!=cv]; y0=y[gr!=cv]; ty0=ty[gr!=cv]; z0=z[gr!=cv]; tz0=tz[gr!=cv]; ny0=ny[gr!=cv]
        x.center = list()
        y.center = list()


        x0 = y0 = ty0 = tz0 = z0 = list()
        ny0 = NULL
        for (i in 1:length(cv_class)) {
          x0[[i]] = matrix(cv_res$x[which(cv_res$id_y == cv_class[i]), ], ncol = dim_x)
          y0[[i]] = cv_res$y[which(cv_res$id_y == cv_class[i])]
          z0[[i]] = matrix(cv_cov$z[which(cv_cov$id_z == cv_class[i]), ], ncol =
                             dim_z)
          ty0[[i]] = cv_res$ty[which(cv_res$id_y == cv_class[i])]
          tz0[[i]] = cv_cov$tz[which(cv_cov$id_z == cv_class[i])]
          ny0[i] = length(y0[[i]])
        }

        Ex = list()
        Ey = list()
        for (i in 1:length(cv_class)) {
          Ex[[i]] = x0[[i]]
          Ey[[i]] = y0[[i]]
          for (k in 1:ny0[i]) {
            temp = K1(ty0[[i]][k] - t.uls, h.temp) / sum(K1(ty0[[i]][k] - t.uls, h.temp))
            Ex[[i]][k, ] = t(x.uls) %*% temp
            Ey[[i]][k] = t(y.uls) %*% temp


          }
          x.center[[i]] = x0[[i]] - Ex[[i]]
          y.center[[i]] = y0[[i]] - Ey[[i]]
        }


        Y = unlist(y.center)
        X = NULL
        for (i in 1:length(cv_class)) {
          X = rbind(X, x.center[[i]])
        }

        b0 = solve(t(X) %*% X) %*% t(X) %*% Y

        rsd = list()
        Z = list()
        for (i in 1:sum(gr != cv)) {
          rsd[[i]] = y0[[i]] - x0[[i]] %*% b0
          Z[[i]] = cbind(1, z0[[i]])
        }
        U = u(tz0, ty0, Z, rsd, h.temp, sum(gr != cv))
        est0 = solve(U[[2]]) %*% U[[1]]


        cv_class1 = which(gr == cv)
        cv_res1 = data_res %>% group_by(id_y) %>% filter(id_y %in% cv_class1)
        cv_cov1 = data_cov %>% group_by(id_z) %>% filter(id_z %in% cv_class1)

        x1 = y1 = ty1 = tz1 = z1 = list()
        ny0 = NULL
        for (i in 1:length(cv_class1)) {
          x1[[i]] = matrix(cv_res1$x[which(cv_res1$id_y == cv_class1[i]), ], ncol =
                             dim_x)
          y1[[i]] = cv_res1$y[which(cv_res1$id_y == cv_class1[i])]
          z1[[i]] = matrix(cv_cov1$z[which(cv_cov1$id_z == cv_class1[i]), ], ncol =
                             dim_z)
          ty1[[i]] = cv_res1$ty[which(cv_res1$id_y == cv_class1[i])]
          tz1[[i]] = cv_cov1$tz[which(cv_cov1$id_z == cv_class1[i])]
        }

        mse.temp = 0
        KK = 0
        temp = 0
        for (i in 1:sum(gr == cv)) {
          for (j in 1:length(ty1[[i]])) {
            for (k in 1:length(tz1[[i]])) {
              KK = KK + K1(tz1[[i]][k] - ty1[[i]][j], h.temp)
              temp = temp + K1(tz1[[i]][k] - ty1[[i]][j], h.temp) * (y1[[i]][j] -
                                                                       est0[1] - est0[-1] %*% z1[[i]][k, ] - x1[[i]][j, ] %*% b0) ^ 2
            }
          }
        }
        if (KK != 0) {
          mse.temp = mse.temp + temp / KK
        }
        mse[cv, hh] = mse.temp
      }
    }
    MSE = apply(mse, 2, mean)
    h = h.can[which.min(MSE)]

  } else{
    h = bd
  }



  S.matrix = function(ty, h) {
    ty = as.matrix(unlist(ty))
    n = length(ty)
    K = matrix(NA, ncol = 1, nrow = n)
    S = matrix(NA, nrow = n, ncol = n)
    for (i in 1:n) {
      K = K1(ty[i] - ty, h)
      #w1=as.numeric(t(ty-ty[i])%*%K); w2=as.numeric(t((ty-ty[i])^2)%*%K); w=as.numeric(t(ty-ty[i]-w2/w1)%*%K); S[i,]=t((ty-ty[i]-w2/w1)*K)/w
      s1 = as.numeric(t(ty[i] - ty) %*% K)
      s2 = as.numeric(t((ty[i] - ty) ^ 2) %*% K)
      w = as.numeric(K * (s2 - (ty[i] - ty) * s1))
      S[i, ] = t(w) / sum(w)
    }
    return(S)
  }


  if (omit == 1 & method == 0) {
    x = y = ty = list()
    ny = rep(0, N)
    for (i in 1:N) {
      x[[i]] = matrix(data_res$x[which(data_res$id_y == i), ], ncol = dim_x)
      y[[i]] = data_res$y[which(data_res$id_y == i)]
      ty[[i]] = data_res$ty[which(data_res$id_y == i)]
      ny[i] = length(which(data_res$id_y == i))
    }



    x.uls = as.matrix(data_res$x)
    y.uls = as.matrix(data_res$y)
    t.uls = as.matrix(data_res$ty)
    n = dim(x.uls)[1]
    S = S.matrix(t.uls, h)
    beta = solve(t(x.uls) %*% t(diag(1, n) - S) %*% (diag(1, n) - S) %*% x.uls) %*%
      t(x.uls) %*% t(diag(1, n) - S) %*% (diag(1, n) - S) %*% y.uls
    b.hat.plm = beta
    a.hat.plm = S %*% (y.uls - x.uls %*% beta)
    #sum(is.na(S))


    e = list(NA)
    j = 0
    for (i in 1:N) {
      ind = c(j + 1, j + length(y[[i]]))
      e[[i]] = y[[i]] - (a.hat.plm[ind[1]:ind[2]] + x[[i]] %*% beta)
      e[[i]] = e[[i]] %*% t(e[[i]])
      j = j + length(y[[i]])
    }

    D = t(x.uls) %*% t(diag(1, n) - S) %*% (diag(1, n) - S) %*% x.uls
    C = bdiag(e)
    V = t(x.uls) %*% t(diag(1, n) - S) %*% C %*% (diag(1, n) - S) %*% x.uls
    se.b.plm = sqrt(diag(solve(D) %*% V %*% solve(D)))

    return(list(est = b.hat.plm, se = se.b.plm))
  }



  if (omit == 2 & method == 0) {
    x = y = ty = list()
    ny = rep(0, N)
    for (i in 1:N) {
      x[[i]] = matrix(data_res$x[which(data_res$id_y == i), ], ncol = dim_x)
      y[[i]] = data_res$y[which(data_res$id_y == i)]
      ty[[i]] = data_res$ty[which(data_res$id_y == i)]
      ny[i] = length(which(data_res$id_y == i))
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
        temp = K1(ty[[i]][k] - t.uls, h) / sum(K1(ty[[i]][k] - t.uls, h))
        if (ny[i] == 1) {
          Ex[[i]] = matrix(t(x.uls) %*% temp, 1)
          Ey[[i]] = t(y.uls) %*% temp
        } else{
          Ex[[i]][k, ] = t(x.uls) %*% temp
          Ey[[i]][k] = t(y.uls) %*% temp
        }

      }
      x.center[[i]] = x[[i]] - Ex[[i]]
      y.center[[i]] = y[[i]] - Ey[[i]]
    }

    X = NULL
    for (i in 1:N) {
      X = rbind(X, x.center[[i]])
    }
    Y = unlist(y.center)

    beta = solve(t(X) %*% X) %*% t(X) %*% Y

    U.sq.g = U.dot.g = ifelse(length(beta) == 1, 0, diag(rep(0, length(beta))))
    for (i in 1:N) {
      # U.sq.g=U.sq.g+(sum(x.center[[i]]*(y.center[[i]]-x.center[[i]]%*%c(beta))))^2
      U.sq.g = U.sq.g + t(c(y.center[[i]] - x.center[[i]] %*% c(beta)) %*%
                            x.center[[i]]) %*% (c(y.center[[i]] - x.center[[i]] %*% c(beta)) %*% x.center[[i]])
      U.dot.g = U.dot.g + t(x.center[[i]]) %*% x.center[[i]]

    }
    b.se.center = (diag(solve(U.dot.g) %*% U.sq.g %*% solve(U.dot.g))) ^
      .5
    return(list(est = beta, se = b.se.center))
  }



  if (omit == 3 & method == 0) {
    x1.uls = data_res$x
    x_add.uls = data_res$x_add
    t.uls = data_res$ty
    y.uls = data_res$y
    Ex1 = Ex2 = list()

    x2.matrix = matrix(c(rep(1, length(x_add.uls)), x_add.uls), length(x_add.uls))

    res_y = y.uls - x2.matrix %*% solve(t(x2.matrix) %*% x2.matrix) %*% t(x2.matrix) %*%
      y.uls
    res_x1 = x1.uls - x2.matrix %*% solve(t(x2.matrix) %*% x2.matrix) %*%
      t(x2.matrix) %*% x1.uls

    ####################### PLM  #######################

    x.uls = res_x1
    y.uls = res_y
    n = length(y.uls)

    S = S.matrix(t.uls, h)
    beta = solve(t(x.uls) %*% t(diag(1, n) - S) %*% (diag(1, n) - S) %*% x.uls) %*%
      t(x.uls) %*% t(diag(1, n) - S) %*% (diag(1, n) - S) %*% y.uls
    b.hat.po = beta
    a.hat.po = S %*% (y.uls - x.uls %*% beta)


    e = x = y = list()
    t = list()
    j = 0
    id = data_res$id_y
    x[[1]] = x.uls[which(id == 1), ]
    y[[1]] = y.uls[which(id == 1)]
    #ind=c(j+1,j+length(x[[1]]))
    e[[1]] = y[[1]] - (a.hat.po[which(id == 1)] + x[[1]] %*% beta)
    e[[1]] = e[[1]] %*% t(e[[1]])
    for (i in 2:N) {
      x[[i]] =  (x.uls[which(id == i), ])
      y[[i]] =  (y.uls[which(id == i)])
      t[[i]] =  (t.uls[which(id == i)])
      e[[i]] = y[[i]] - (a.hat.po[which(id == i)] + x[[i]] %*% beta)
      e[[i]] = e[[i]] %*% t(e[[i]])
      j = j + length(x[[i]])
    }

    D = t(x.uls) %*% t(diag(1, n) - S) %*% (diag(1, n) - S) %*% x.uls
    C = bdiag(e)
    V = t(x.uls) %*% t(diag(1, n) - S) %*% C %*% (diag(1, n) - S) %*% x.uls
    se.b.po = sqrt(as.numeric(diag(solve(D) %*% V %*% solve(D))))
    return(list(est = b.hat.po, se = se.b.po))
  }



  if (method == 1) {
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
        temp = K1(ty[[i]][k] - t.uls, h) / sum(K1(ty[[i]][k] - t.uls, h))
        if (ny[i] == 1) {
          Ex[[i]] = matrix(t(x.uls) %*% temp, 1)
          Ey[[i]] = t(y.uls) %*% temp
        } else{
          Ex[[i]][k, ] = t(x.uls) %*% temp
          Ey[[i]][k] = t(y.uls) %*% temp
        }

      }
      x.center[[i]] = x[[i]] - Ex[[i]]
      y.center[[i]] = y[[i]] - Ey[[i]]
    }

    X = NULL
    for (i in 1:N) {
      X = rbind(X, x.center[[i]])
    }
    Y = unlist(y.center)

    beta = solve(t(X) %*% X) %*% t(X) %*% Y


    U.sq.g = U.dot.g = ifelse(length(beta) == 1, 0, diag(rep(0, length(beta))))
    for (i in 1:N) {
      # U.sq.g=U.sq.g+(sum(x.center[[i]]*(y.center[[i]]-x.center[[i]]%*%c(beta))))^2
      U.sq.g = U.sq.g + t(c(y.center[[i]] - x.center[[i]] %*% c(beta)) %*%
                            x.center[[i]]) %*% c(y.center[[i]] - x.center[[i]] %*% c(beta)) %*% x.center[[i]]
      U.dot.g = U.dot.g + t(x.center[[i]]) %*% x.center[[i]]

    }
    b.se.center = (diag(solve(U.dot.g) %*% U.sq.g %*% solve(U.dot.g))) ^
      .5

    rsd = list()
    Z = list()


    for (i in 1:N) {
      rsd[[i]] = y[[i]] - x[[i]] %*% c(beta)
      Z[[i]] = cbind(1, z[[i]])
    }
    U = u(tz, ty, Z, rsd, h, N)
    est = solve(U[[2]]) %*% U[[1]]
    g.hat.2sg = est[-1]
    a.hat.2sg = est[1]

    U.sq = 0
    for (i in 1:N) {
      temp.U = 0
      for (j in 1:length(ty[[i]])) {
        for (k in 1:length(tz[[i]])) {
          temp = Z[[i]][k, ]
          k1 = K1(tz[[i]][k] - ty[[i]][j], h)
          temp.U = temp.U + k1 * temp * as.numeric(rsd[[i]][j] - temp %*%
                                                     est)
        }
      }
      U.sq = U.sq + temp.U %*% t(temp.U)
    }
    se_est = diag(solve(U[[2]]) %*% U.sq %*% solve(U[[2]])) ^ .5
    g.se.2sg = se_est[-1]
    a.se.2sg = se_est[1]

    return(list(
      est = c(a.hat.2sg, beta, g.hat.2sg),
      se = c(a.se.2sg, b.se.center, g.se.2sg)
    ))
  }


  ############# kernel smoothing #######################################


  if (method == 2) {
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


    U.dot = 0
    U = 0
    for (i in 1:N) {
      for (j in 1:length(ty[[i]])) {
        for (k in 1:length(tz[[i]])) {
          temp = c(1, x[[i]][j, ], z[[i]][k, ])
          k1 = K1(tz[[i]][k] - ty[[i]][j], h)
          U = U + k1 * temp * y[[i]][j]
          U.dot = U.dot + k1 * temp %*% t(temp)
        }
      }
    }
    est = solve(U.dot) %*% U
    a.hat.ks = est[1]
    b.hat.ks = est[2:(dim_x + 1)]
    g.hat.ks = est[-(1:(dim_x + 1))]


    U.sq = 0
    for (i in 1:N) {
      temp.U = 0
      for (j in 1:length(ty[[i]])) {
        for (k in 1:length(tz[[i]])) {
          temp = c(1, x[[i]][j, ], z[[i]][k, ])
          k1 = K1(tz[[i]][k] - ty[[i]][j], h)
          temp.U = temp.U + k1 * temp * as.numeric(y[[i]][j] - temp %*%
                                                     est)
        }
      }
      U.sq = U.sq + temp.U %*% t(temp.U)
    }
    se = diag(solve(U.dot) %*% U.sq %*% solve(U.dot)) ^ .5
    a.se.ks = se[1]
    b.se.ks = se[2:(dim_x + 1)]
    g.se.ks = se[-(1:(dim_x + 1))]

    return(list(
      est = c(a.hat.ks, b.hat.ks, g.hat.ks),
      se = c(a.se.ks, b.se.ks, g.se.ks)
    ))

  }





}
