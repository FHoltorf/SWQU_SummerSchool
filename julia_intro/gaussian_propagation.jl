using LinearAlgebra, ForwardDiff

struct LinearGaussianApproximation
    μ
    Σ
end

const LGA = LinearGaussianApproximation
import Base.+
import Base.-
import Base.*

function *(a::LGA, b::Number)
    return LGA(a.μ*b, b^2*b.Σ)
end
function *(a::Number, b::LGA)
    return b*a
end
function +(a::LGA, b::LGA) 
    return LGA(a.μ + b.μ, a.Σ + b.Σ)
end
function +(a::LGA, b::Number) 
    return LGA(a.μ + b, a.Σ)
end
function +(a::Number, b::LGA) 
    return b+a
end


