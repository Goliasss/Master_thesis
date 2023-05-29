source("profileCI_source.R")

# Fit a model to the ryegrass example data
m1 = drm(rootl ~ conc, data = ryegrass, fct = LL.4())

# Run the default EDx() function -- returns confidence interval for the ED50
EDx(m1)

# Specify custom EDs
EDx(m1, c(1,5,10))

# Change the lambda parameter
EDx(m1, lambda = 0.01)
