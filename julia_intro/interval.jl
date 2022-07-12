import Base.+
import Base.-
import Base.*
import Base./
import Base.zero
import Base.zeros
struct Interval{T}
    lb::T
    ub::T
end

+(a::Interval, b::Interval) = Interval(a.lb + b.lb, a.ub + b.ub)
+(a::Interval, b::Number) = Interval(a.lb+b, a.ub+b)
+(a::Number, b::Interval) = b+a
-(a::Interval, b::Interval) = Interval(a.lb - b.ub, a.ub - b.lb)
-(a::Interval, b::Number) = Interval(a.lb-b, a.ub-b)
-(a::Number, b::Interval) = -1*(b-a)
*(a::Interval, b::Interval) = begin 
    vertices = [a.lb*b.ub, a.lb*b.lb, a.ub*b.lb, a.ub*b.ub] 
    return Interval(minimum(vertices), maximum(vertices))
end
*(a::Number, b::Interval) = Interval(a*b.lb, a*b.ub)
*(a::Interval, b::Number) = b*a
*(A::Matrix, b::Vector{<:Interval}) = begin
    n,m = size(A)
    res = [Interval(0.0,0.0) for i in 1:n]
    for i in 1:n
        res[i] = sum(A[i,j]*b[j] for j in 1:m)
    end
    return res
end

# discrete dynamical system
# x_{k+1} = x_k + x_k*(c[1] - c[2] y_k) 
# y_{k+1} = y_k + y_k*(-c[3] + c[4] x_k)

LV_recursion(x,c) = [x[1]*(1+c[1]) - c[2]*x[1]*x[2];
                     x[2]*(1-c[3]) + c[4]*x[1]*x[2]]
c = 0.01*rand(4)

# discrete chemical system
# A -> B <-> 2C
k = [1.0, 3.0, 0.4]
S = [-1 0 0;
      1 -1 1
      0 2 -2]
rxns(x,k) = [k[1]*x[1],k[2]*x[2],k[3]*x[3]*x[3]]
f(x,S,k) = S*rxns(x,k)
reactor_recursion(x,S,k,dt) = x + dt*f(x,S,k)

N = 2000
dt = 0.001
x0 = [1.0, 0.0, 0.0]
xs = [x0]
for i in 1:N
    push!(xs, reactor_recursion(xs[end], S, k, dt))
end

plot([x[1] for x in xs])
plot!([x[2] for x in xs])
plot!([x[3] for x in xs])
# interval arithmetic

I_x0 = [Interval(1.0, 1.1), Interval(0.0,0.0), Interval(0.0,0.0)]
I_xs = [I_x0]
for i in 1:N
    push!(I_xs, reactor_recursion(I_xs[end], S, k, dt))
end
plot([x[1].lb for x in I_xs], linestyle = :dash, color = :white, linewidth = 2)
plot!([x[1].ub for x in I_xs], linestyle = :dash, color = :white, linewidth = 2)
plot!([x[2].lb for x in I_xs], linestyle = :dash, color = :red, linewidth = 2)
plot!([x[2].ub for x in I_xs], linestyle = :dash, color = :red, linewidth = 2)
plot!([x[3].lb for x in I_xs], linestyle = :dash, color = :blue, linewidth = 2)
plot!([x[3].ub for x in I_xs], linestyle = :dash, color = :blue, linewidth = 2)

