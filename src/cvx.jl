#%% Make the Convex.jl module available
using Convex
using LinearAlgebra
using SCS

#%% Generate random problem data
Λ = 30
L = 4
N = 32

T_Λ = randn(N, (Λ+1)^2)

#%% Create variable for encoding matrix
A = Variable((L+1)^2, N)

#%% Create identity matrix
C = hcat( I((L+1)^2), zeros(((L+1)^2), (Λ+1)^2 - (L+1)^2) )

#%% weighting of regularization parameter
gamma = 0.2

#%% The problem is to minimize ||Ax - b||^2 with L2 regularization
problem = minimize(sumsquares(A * T_Λ - C) + gamma * sumsquares(A))

#%% The problem is to minimize ||Ax - b||^2 with L1 regularization
problem = minimize(norm(vec(A * T_Λ - C), 2) + gamma * norm(vec(A), 1))

#%% Solve the problem by calling solve!
@time solve!(problem, () -> SCS.Optimizer())

#%% Check the status of the problem
problem.status

#%% Get the optimum value
problem.optval

#%% Print the solution
E_L = A.value

#%% Get closed form solution for E_L
T_L = @view T_Λ[:,1:(L+1)^2]
E_closed = T_L' * inv((T_Λ * T_Λ') + (gamma * I))

#%% Get rms error
rms_error = norm(E_L - E_closed, 2)
