module FFT

function twiddle!(wa::Vector{Complex{T}}) where T<:AbstractFloat
    Nby2 = length(wa)
    sinθ, cosθ = sincospi(inv(Nby2)) # sin(θ), cos(θ)

    cosΘ = 1.0 # cos(kθ)
    sinΘ = 0.0 # sin(kθ)

    @inbounds wa[1] = complex(cosΘ,  sinΘ)
    isone(Nby2) && return wa

    Nby4 = Nby2 >> 1
    @inbounds wa[Nby4+1] = complex(sinΘ, -cosΘ)
    isone(Nby4) && return wa

    Nby8 = Nby2 >> 2
    if Nby8 > 1
        jj = Nby4
        kk = Nby4 + 2
        ll = Nby2
        for ii in 2:Nby8
            cosϕ = cosΘ * cosθ + sinΘ * sinθ
            sinϕ = sinΘ * cosθ - cosΘ * sinθ
            cosΘ = cosϕ
            sinΘ = sinϕ
            @inbounds begin
                wa[ii] = complex( cosΘ,  sinΘ)
                wa[jj] = complex(-sinΘ, -cosΘ)
                wa[kk] = complex( sinΘ, -cosΘ)
                wa[ll] = complex(-cosΘ,  sinΘ)
            end
            jj -= 1
            kk += 1
            ll -= 1
        end
    end

    if Nby2 > 2
        ii = Nby8 + 1
        @inbounds wa[ii] = complex( 0.7071067811865476, -0.7071067811865476)
        ii += Nby4
        @inbounds wa[ii] = complex(-0.7071067811865476, -0.7071067811865476)
    end
    return wa
end

"""
Cooley-Tukey butterfly computation.

    ctb!(ya::VecI{Complex{T}},
         xa::VecI{Complex{T}},
         wa::VecI{Complex{T}},
         si::Int, hs::Int, ns::Int,
         ss::Int, pd::Int) where T<:AbstractFloat
-----------------------------------
1. `ya`: destination array
2. `xa`: source array
3. `wa`: twiddle factors array
4. `si`: starting index
5. `hs`: half of FFT size
6. `ns`: number of steps
7. `ss`: step size of the current subproblem
8. `pd`: index difference of butterfly computation pair
"""
function ctb!(
        ya::VecI{Complex{T}},
        xa::VecI{Complex{T}},
        wa::VecI{Complex{T}},
        si::Int,
        hs::Int,
        ns::Int,
        ss::Int,
        pd::Int
    ) where T<:AbstractFloat
    wi = 1
    yi = xi = si
    for _ in 1:ns
        xj = xi + hs
        @inbounds begin
            ya[yi]    = (xa[xi] + xa[xj])
            ya[yi+pd] = (xa[xi] - xa[xj]) * wa[wi]
        end
        yi += ss
        xi += pd
        wi += pd
    end
    return nothing
end

end
