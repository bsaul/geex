library(Ryacas)

### an example of a function in one variable

# define symbols that you will use in your functions
x <- Sym("x")

# define the function
my_func <- function(x) {
  return(x/(x^2 + 3))
}

# take the derivative
yacas_deriv <- yacas(deriv(my_func(x), x))

# convert the derivative to a function
my_deriv <- function(x) {
  eval(parse(text = yacas_deriv$text))
}

# or do it manaully
my_deriv2 = function(x) {
  eval(expression((x^2 + 3 - 2 * x^2)/(x^2 + 3)^2)) # copy yacas_deriv output
}

my_deriv(0)
my_deriv2(0)

# test it out...
xvals=seq(-10,10,.01) 
plot(xvals,sapply(xvals,my_func),type='l',xlab="x",ylab="f(x)")
points(0,my_func(0))
abline(a=my_func(0),b=my_deriv(0),lty=2) # deriv at x=0
points(-5,my_func(-5))
abline(a=my_func(-5)+my_deriv(-5)*(5),b=my_deriv(-5),lty=2) # deriv at x=-5
points(5,my_func(5))
abline(a=my_func(5)+my_deriv(5)*(-5),b=my_deriv(5),lty=2) # deriv at x=5


### an example with a constant vector

# define theta 
theta <- Sym("theta")

# a function that depends on theta
func <- function(x) {
    return(theta[1] * x + theta[2])
}

# Let's integrate this
Func <- yacas(Integrate(func(x), x))
# returns (x^2*theta)/2+NA*x; which is not quite what we want...

# To work around this problem, define your functions like this:
a <- Sym("a")
b <- Sym("b")
func2 <- function(x) {
  return(a * x + b)
}

# Now integrate 
Func2 <- yacas(Integrate(func2(x), x))
# returns (x^2*a)/2+b*x;

# Now replace a and b by the thetas:
Func2 <- gsub("a", "theta[1]", Func2$text)
Func2 <- gsub("b", "theta[2]", Func2)

my_integral <- function(x,theta) {
  # eval(parse(text = Func2)) # auto way
  eval(x^2 * theta[1]/2 + theta[2] * x) # manual way
}
my_integral(2,c(1,1))

### an example with a multivariate function

p.0000 = Sym("p0000") # note: cannot put use "." in symbol
p.0001 = Sym("p0001")

func = function(p.0000,p.0001) {
  return(p.0001+p.0000+2*p.0001*p.0000)
}

deriv.0000_expression = yacas(deriv(func(p.0000,p.0001), p.0000))
deriv.0000_expression

deriv.0001_expression = yacas(deriv(func(p.0000,p.0001), p.0001))
deriv.0001_expression

### example of obtaining the score vector and hessian matrix of a multivariate function

# define symbols that you will use in your functions
x <- Sym("x")
y <- Sym("y")
a <- Sym("a")
b <- Sym("b")

score.expression = yacas("D({x,y})
    Exp(x)/(1+Exp(x))+b*Exp(y)") 
score.expression

score.expression = yacas("D({x,y})
    a*Ln(x)+Ln(y)+b*(x^2)*(y^2)") 
score.expression

score.expression = yacas("D({x,y})
  a*x^3+y^3+b*x^2*y^2")
temp = list(a * (3 * x^2) + b * (2 * x) * y^2, 3 * y^2 + b * x^2 * (2 * y))
score = unlist(temp)
score = gsub(" ","",score)
for(i in 1:2) {cat(score[i],",",sep="")}
score.vector = c(
    ((a*(3*(x^2)))+((b*(2*x))*(y^2))),((3*(y^2))+((b*(x^2))*(2*y)))
  )

hessian.expression = yacas("HessianMatrix(
    x^3+y^3+(x^2)*(y^2), 
  {x,y})") 
temp = list(list(6 * x + 2 * y^2, 2 * x * (2 * y)), list(2 * x * (2 * y), 6 * y + 2 * x^2))
hessian = matrix(unlist(temp),nrow=2,ncol=2,byrow=T)
hessian = gsub(" ","",hessian)
for(i in 1:2) for (j in 1:2) {cat(hessian[i,j],",",sep="")}
hessian.matrix = matrix(c( 
    ((6*x)+(2*(y^2))),((2*x)*(2*y)),((2*x)*(2*y)),((6*y)+(2*(x^2)))
  ),nrow=2,ncol=2,byrow=T)

