library(reshape2)
library(drc)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(scales)
library(princurve)
library(stringr)
library(tidyverse)

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

profile_grid <- function(object, perc=50, level = 0.95, type="relative",
                         display.contour=TRUE, precision=50)
{  
  iVec <- seq(confint(object)["b:(Intercept)",1],confint(object)["b:(Intercept)",2], length=precision)
  jVec <- seq(confint(object)["e:(Intercept)",1],confint(object)["e:(Intercept)",2], length=precision)
  
  df = data.frame(expand.grid(iVec,jVec))
  names(df) = c("iVec","jVec")
  
  fcn = function(iVec, jVec) {
    # fit = drm(object$call$formula, data = object$data,
    #           fct = LL.4(fixed = c(iVec, NA, NA, jVec)))
    theta = rep(NA, nrow(object$parmMat))
    theta[1] = iVec
    theta[length(theta)] = jVec
    fit = update(object, fct = do.call(object$fct$name, list(theta)))
    -2*(loglikelihood(fit) - loglikelihood(object)) - qchisq(level, 1)
  }
  cbind(df, data.frame(ll=mapply(fcn, df$iVec, df$jVec)))
}


loglikelihood = function(object, theta=NULL) {
  if (is.null(theta)) theta = object$parmMat
  f = object$fct$fct
  df = object$fct$deriv1
  dose = object$data[[object$dataList$names$dName]]
  resp = object$data[[object$dataList$names$orName]]
  weights = object$data[["weights"]]
  N = length(object$data[[object$dataList$names$orName]])
  if (object$type == "continuous") {
    sigma2 = sum((f(dose, t(theta))-resp)^2)/N
    ll=sum(log(1 / (sqrt(2*pi*sigma2))) - (1/(2 * sigma2)) * 
             (resp - f(dose, t(theta)))^2)
  } else if (object$type == "Poisson") { 
    ll=sum(-f(dose, t(theta)) + resp * log(f(dose, t(theta))) - log(factorial(resp)))
  } else if (object$type == "binomial") {
    total = (object$data)[, 5]
    success = total * (object$data)[, 2]
    prob = f(dose, t(theta))
    ll = sum((log(choose(total, success)) + 
                resp * log(f(dose, t(theta))) + (1-resp) * log(1-f(dose, t(theta)))))
    ll = sum(log(choose(total, success)) + 
               resp * weights * log(f(dose, t(theta))/(1 - f(dose, t(theta)))) + 
               weights * log(1 - f(dose, t(theta))))
  }
  ll
}

score = function(object, theta=NULL) {
  if (is.null(theta)) theta = object$parmMat
  f = object$fct$fct
  df = object$fct$deriv1
  dose = object$data[[object$dataList$names$dName]]
  resp = object$data[[object$dataList$names$orName]]
  weights = object$data[["weights"]]
  N = length(object$data[[object$dataList$names$orName]])
  if (object$type == "continuous") {
    sigma2 = sum((f(dose, t(theta))-resp)^2)/N
    grad = (1/sigma2) * colSums((resp - f(dose, t(theta))) * df(dose, t(theta)))
  } else if (object$type == "Poisson") {
    grad = colSums(((resp / f(dose, t(theta))) - 1) * df(dose, t(theta)))
  } else if (object$type == "binomial") {
    grad = colSums(weights * df(dose, t(theta)) * 
                     (resp/f(dose, t(theta)) + resp/(1 - f(dose, t(theta))) - 1/(1 - f(dose, t(theta)))))
  }
  grad = as.matrix(grad, byrow=TRUE)
  grad[c(1,nrow(grad))]
}


confidence_region <- function(object, level = 0.95, lambda=NULL, epsilon=1e-3, maxiter=2000)
{ 
  profile_theta = function(object, point) {
    fixed[c(1,4)] = point
    fct = fct_gen(fixed)
    if (sum(is.na(fixed)) == 0) {
      theta = point
    } else {
      temp = update(object, fct = fct)
      theta = fixed[1:4]
      theta[is.na(theta)] = coef(temp)
      theta = theta[is.na(original_fixed)]
    }
    return(theta)
  }
  
  descent = function(step=delta, point) {
    beginning = point
    theta = profile_theta(object, point)
    new_loglik = as.vector(-2*(loglikelihood(object, theta) - loglikelihood(object)) - qchisq(level, 1))
    dff <<- rbind(dff, data.frame(b=point[1], e=point[2],  ll=abs(new_loglik)))
    
    descent_counter = 0
    while (!between(new_loglik, -epsilon, epsilon)) {
      grad = score(object, matrix(theta, byrow = TRUE))
      grad = (new_loglik * grad * step) / (norm(grad, type = "2"))
      point = point + grad 
      theta = profile_theta(object, point)
      new_loglik = as.vector(-2*(loglikelihood(object, theta) - loglikelihood(object)) - qchisq(level, 1))
      dff <<- rbind(dff, data.frame(b=point[1], e=point[2],  ll=abs(new_loglik)))
      descent_counter = descent_counter + 1
      candidate_counter = candidate_counter + 1 
      # print(c(descent_counter, delta))
      if ((descent_counter == 50) & is.null(lambda)) {
        step = step / 2
        delta <<- delta / 2
        point = beginning
        descent_counter = 0}
      if (descent_counter == maxiter) {
        stop("convergence error")
      }
    }
    
    dff <<- rbind(dff, data.frame(b=point[1], e=point[2],  ll=abs(new_loglik)))
    df <<- rbind(df, data.frame(b=point[1], e=point[2],  ll=abs(new_loglik)))
    point
  }
  
  generate_rotation_matrix = function(angle) {
    matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), nrow=2, byrow = TRUE)
  }
  
  find_contour = function (rotation) {
    point = beginning
    theta = profile_theta(object, point)
    grad = score(object, matrix(theta, byrow = TRUE))
    grad = (grad * delta) / (norm(grad, type = "2"))
    point = point + rotation %*% grad
    counter = 1
    while (norm(beginning-point, type="2")>=0.5*delta) {
      theta = profile_theta(object, point)
      grad = score(object, matrix(theta, byrow = TRUE))
      grad = (grad * delta) / (norm(grad, type = "2"))
      point = point + rotation %*% grad
      point = descent(step=delta, point)
      # print(counter)
      # print(point)
      if (counter %% 10 == 0) {
        print(ggplot(df, aes(x=b, y=e)) + geom_point() + theme_bw())
      }
      counter = counter + 1
      if (counter == ceiling(maxiter / 2)) {
        break
      } 
    }
    counter
  }
  
  start_time = Sys.time()
  
  original_object = object
  
  # Scaling the x-axis for computational purposes
  scaling_coef = abs(object$coefficients[["b:(Intercept)"]] / object$coefficients[["e:(Intercept)"]])
  object$origData[[object$call[[2]][[3]]]] =
    object$origData[[object$call[[2]][[3]]]] * scaling_coef
  
  if (is(object$fct, "llogistic")) {
    fct_gen = llogistic 
  } else if (is(object$fct, "Weibull-1")) {
    fct_gen = weibull1
  } else {
    stop(paste(object$fct$name, "function not supported",sep=" "))
  }
  
  fixed = object$fct$fixed
  original_fixed = fixed
  
  object = update(object, data = object$origData)
  
  
  # Setting up tuning parameter
  if (is.null(lambda)) {
    delta = 1 * min(unname(coef(object)[c("b:(Intercept)","e:(Intercept)")]))
    delta = level * delta / sqrt(length(object$data[[object$dataList$names$orName]]))
  } else {
    delta = lambda * min(unname(coef(object)[c("b:(Intercept)","e:(Intercept)")]))
    delta = level * delta / sqrt(length(object$data[[object$dataList$names$orName]]))
  }
  
  df <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(df) <- c("b", "e", "ll")
  dff <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(dff) <- c("b", "e", 'll')
  candidate_counter = 0
  
  # Searching for a point on the boundary
  beginning = unname(coef(object)[c("b:(Intercept)","e:(Intercept)")])
  
  if (is.null(lambda)) {
    repeat {
      point = try(descent(step=delta, beginning), silent = TRUE)
      if (inherits(point,"try-error")) {
        delta = delta / 2
        # print(delta)
      } else break
    }
  } else {
    point = descent(step=delta, beginning)
  }
  
  ggplot(dff, aes(x=b, y=e, color=ll)) + geom_point() +scale_color_gradient(low="red", high="blue")
  
  # Moving alongside the boundary
  beginning = point
  runs = find_contour(rotation = generate_rotation_matrix(pi/2))
  if (runs == ceiling(maxiter / 2)) {find_contour(rotation = generate_rotation_matrix(3*pi/2))}
  
  if (runs == ceiling(maxiter / 2))  warning("max_iter reached")
  
  df$e = df$e / scaling_coef
  
  list(CR=df, name=object$fct$name, candidates=dff, exec_time=Sys.time()-start_time, diagnostic_plot = ggplot(dff, aes(x=b, y=e, color=ll)) + geom_point() +scale_color_gradient(low="red", high="blue"))
}

EDfct = function(df, object, perc) {
  if (is(object$fct, "llogistic")) {
    return(exp((1/df$b)*log((1-(1 - perc/100))/(1 - perc/100))+log(df$e))) 
  } else if (is(object$fct, "Weibull-1")) {
    return(exp((1/df$b) * log(log(1/(1 - perc/100))) + log(df$e)))
  } else {
    stop(paste(object$fct$name, "function not supported",sep=" "))
  }
}

EDall = function(CR, object, perc, level = 0.95, lambda=NULL, epsilon=1e-2) {
  for (p in perc) {
    CR[[paste("ED",p,sep="")]] = EDfct(CR, object, perc = p)
  }
  CR
}

diag_plot = function(CR) {
  col_names = names(CR$CR)[str_detect(names(CR$CR), "ED")]
  figs = list()
  i = 1
  for (col in col_names) {
    ed_min = CR$CR[which.min(CR$CR[[col]]),]
    ed_max = CR$CR[which.max(CR$CR[[col]]),]
    
    figs[[i]] = ggplot(CR$CR, aes_string(x="b", y="e", color=col)) + theme_bw()+ geom_point() +
      scale_color_gradient(low="coral3", high="chartreuse3") +
      geom_point(data = ed_min, 
                 aes_string(x="b", y="e", color=col), size=3,stroke=3, 
                 fill="white", shape=21) +
      geom_point(data = ed_max, 
                 aes_string(x="b", y="e", color=col), size=3,stroke=3, 
                 fill="white", shape=21) + 
      ggtitle(paste("Diagnostic plot for ", col, sep="")) +
      theme(plot.title = element_text(hjust = 0.5))
    i = i + 1
    
  }
  # grid.arrange(figs)
  do.call(grid.arrange, c(figs, list(ncol=length(figs))))
}

EDx = function(object, perc=c(50), level=0.95, lambda=NULL, epsilon=1e-2, maxiter=2000) {
  CR = confidence_region(object, level, lambda, epsilon, maxiter = maxiter)
  CR$CR = EDall(CR$CR, object, perc, level, lambda, epsilon)
  CI = data.frame()
  lower = paste((1-level) * 100 / 2,"%")
  upper = paste(100 - (1-level) * 100 / 2, "%")
  for (p in perc) {
    row = data.frame(paste("ED",p,sep=""), 
                     quiet(ED(object, respLev = p)[,"Estimate"]),
                     min(CR$CR[[paste("ED",p,sep="")]]),
                     max(CR$CR[[paste("ED",p,sep="")]]))
    colnames(row) = c("ED", "Estimate", lower, upper)
    CI = rbind(CI, row)
  }
  rownames(CI) <- CI[,1]
  CI[,1] <- NULL
  CR$CI = CI
  CR$ED_diag = diag_plot(CR)
  class(CR) = "profileCI"
  CR
}


print.profileCI = function(CR) {
  # print(CR$ED_diag)
  print(CR$CI)
}

