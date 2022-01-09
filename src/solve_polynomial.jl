
"""
    Converts system of symbolic polynoms into a system of AA polynoms
"""
function polynormalize(system::AbstractArray)
    # find all vars present in system
    vars = union!(map(get_variables, system)...)
    ring, polyvars = AbstractAlgebra.PolynomialRing(AbstractAlgebra.QQ, map(string, vars))

    # sym => monom
    var2poly = Dict(map(=>, vars, polyvars))

    polys = Vector{AbstractAlgebra.elem_type(ring)}(undef, length(system))
    # convert each symbolic poly into an AA poly
    for (i, expr) in enumerate(system)
        polys[i] = substitute(expr, var2poly)
    end

    polys
end

"""
    Solve polynomial `system`.

    rhs are assumed to equal zero
"""
function solve_poly(system::AbstractArray)
    # solves system .= 0

    # convert inputs into polynomials from AA
    polys = polynormalize(system)

    # compute the Groebner basis
    # basis = GroebnerBasis.groebner(polys)
    #
    # GroebnerBasis has bo suitable api yet,
    # a placeholder for now
    basis = polys

    #=
    Now basis is a system that admits shape form,
    i.e it looks like
        { f₁(xn), f₂(xn,xn-1) ... fₙ(xn,..,x1) }

    Such system can be solved iteratively
    by finding roots of univariate poly f₁(xn),
    substituting them into f₂, solving f₂, and so on
    =#

    # solve system in shape position as above
    solve_shape(basis)
end

# solve one univariate polynomial equation
function solve_shape(poly)
    println(poly)
    coeffs   = reverse(AbstractAlgebra.coefficients_of_univariate(poly))
    unipoly  = Polynomials.Polynomial(coeffs)
    Polynomials.roots(unipoly)
end

"""
    Solve polynomial `system` in shape position
"""
function solve_shape(system::AbstractArray)
    length(system) == 1 && return solve_shape(first(system))

    ring   = AbstractAlgebra.parent(first(system))
    ground = AbstractAlgebra.base_ring(ring)
    vars   = AbstractAlgebra.gens(ring)

    roots  = Vector{Vector{AbstractAlgebra.elem_type(ground)}}(undef, 0)

    poly = first(system)
    # select a `poly` from `system`, find its roots,
    # and for each `root`..
    for root in solve_shape(poly)
        vars[length(vars) - length(system) + 1] = ring(root)
        # ..substiture all other polynomials in system with this `root`..
        specialized = map(f -> AbstractAlgebra.evaluate(f, vars), system[2:end])
        # ..and collect solutions of a specialized system
        for specroots in solve_shape(specialized)
            push!(roots, [root, specroots...])
        end
    end

    roots
end

##### EXAMPLE ######

@variables x y z

# for now, already in shape position
system = [
        x - 1,
        x*y + 1,
        x*z*y - 2
]

println( solve_poly(system) )

####################
