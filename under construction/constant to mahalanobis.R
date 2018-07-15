#This function shows that 

# If modifying the variance, then const = (1/MH^2)...
#     MH = 1/sqrt(const). Thus,
#     sqrt(const) = 1/MH --> const = (1/MH^2)

# If modifying the means, then const = MH

S = diag(3)
x = c(0,0,0)
y = c(1,0,0)

N = 100
const_vec = seq(1/8, 8, len = N)
MH_vec = rep(NA, N)

for (i in 1:N){

  c = const_vec[i]
  Si_c = solve(c*S)
  MH_c= as.numeric(sqrt(t(x-y)%*%Si_c%*%(x-y)))
  
  MH_vec[i] = MH_c
  
}

plot(const_vec, MH_vec)
lines(const_vec, 1/sqrt(const_vec))


MH_vec2 = rep(NA,N)
Si = solve(S)
for (i in 1:N){
  
  c = const_vec[i]
  MH_c= as.numeric(sqrt(t(c*x-c*y)%*%Si%*%(c*x-c*y)))
  
  MH_vec2[i] = MH_c
  
}

plot(const_vec, MH_vec2)
lines(const_vec, const_vec)


