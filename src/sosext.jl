export sosextractsequence

# Atom extraction
function sosextractsequence(s::AbstractSwitchedSystem, d::Integer; solver::AbstractMathProgSolver=JuMP.UnsetSolver(), tol=1e-5)
    lyap = getlyap(s, d; solver=solver, tol=tol)
    for σ in 1:2
        μ = lyap.dual[σ]
        X = monomials(variables(μ), d)
        ν = matmeasure(μ, X)
        @show extractatoms(ν, 1e-4).support
    end
end
