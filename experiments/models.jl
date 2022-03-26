using JLD2

function oscillator()
    A = [0.0 1; -4π 0]
    X0 = BallInf([0.0, 10.0], 0.1)
    X = Universe(2)
    P = @ivp(x' = A*x, x ∈ X, x(0) ∈ X0)
    return P
end

function freedom()
    m1 = 1.0
    m2 = 1.0
    k1 = 1e4
    k2 = 1.0

    M = [m1 0; 0 m2]
    K = [(k1+k2) -k2; -k2 k2]

    sys = SecondOrderLinearContinuousSystem(M, zeros(size(M)), K)

    n = 2
    A = [zeros(n, n) Matrix(1.0I, n, n)
         -inv(M)*K    zeros(n, n)       ]

    X0 = Hyperrectangle([1.0, 10.0, 0.0, 0.0], [0.1, 0.5, 0.5, 0.5])
    P = @ivp(x' = Ax, x(0) ∈ X0)
    P = normalize(P)
    return P
end

function iss()
    ISS_path = joinpath(@__DIR__, "ISS.jld2")
    @load ISS_path A B

    A = Matrix(A)  # some methods cannot deal with sparse matrices
    U = Hyperrectangle(low=[0.0, 0.8, 0.9], high=[0.1, 1., 1.])
    U = ConstantInput(linear_map(B, U))
    X = Universe(270)
    X0 = BallInf(zeros(size(A, 1)), 0.0001)  # -0.0001 <= xi <= 0.0001 for all i
    P = @ivp(x' = A*x + u, x(0) ∈ X0, u ∈ U, x ∈ X)
    return P
end
