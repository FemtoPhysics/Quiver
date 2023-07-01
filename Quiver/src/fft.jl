module FFT

const VecI = AbstractVector

## power of radix 2 (flooring mode)
function pwr2(x::Int64)
    x > 0 || error("pwr2(x): x should be positive!")
    r = 0
    x > 4294967295 && (x >>= 32; r += 32)
    x >      65535 && (x >>= 16; r += 16)
    x >        255 && (x >>=  8; r +=  8)
    x >         15 && (x >>=  4; r +=  4)
    x >          3 && (x >>=  2; r +=  2)
    x >          1 && (x >>=  1; r +=  1)

    return r
end

function twiddle!(wa::VecI{Complex{T}}) where T<:AbstractFloat
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

"""
Decimation-in-time FFT with a naturally ordered input-output.

    ditnn!(sa::VecI{Complex{T}},
           ba::VecI{Complex{T}},
           wa::VecI{Complex{T}},
           hs::Int) where T<:AbstractFloat
-------------------------------------------------------------
1. `sa`: signal array
2. `ba`: buffer array
3. `wa`: twiddle factors array
4. `sf`: switch flag
"""
function ditnn!(sa::VecI{Complex{T}}, ba::VecI{Complex{T}},
                wa::VecI{Complex{T}}, hs::Int) where T<:AbstractFloat
    ns = hs
    pd = 1
    ss = 2
    sf = false

    while ns > 0
        if sf
            for si in 1:pd
                ctb!(sa, ba, wa, si, hs, ns, ss, pd)
            end
        else
            for si in 1:pd
                ctb!(ba, sa, wa, si, hs, ns, ss, pd)
            end
        end

        ns >>= 1
        pd <<= 1
        ss <<= 1
        sf = !sf
    end

    return nothing
end

"""
Radix2 Kernel for 1D FFT

Initialize and return a kernel with `fftsize` and specific type `T<:AbstractFloat`.

    FFTKernel{T}(fftsize::Integer) where T<:AbstractFloat

Initialize and return a kernel with `fftsize` and default type `Float64`.

    FFTKernel(fftsize::Integer)
"""
struct FFTKernel{T<:AbstractFloat}
    cache  ::Vector{Complex{T}}
    twiddle::Vector{Complex{T}}
    fftsize::Int
    ifswap ::Bool

    function FFTKernel{T}(fftsize::Integer) where T<:AbstractFloat # type-stability ✓
        # fftsize should be power of 2
        cache   = Vector{Complex{T}}(undef, fftsize)
        twiddle = Vector{Complex{T}}(undef, fftsize >> 1)

        return new{T}(cache, twiddle!(twiddle), fftsize, isone(pwr2(fftsize) & 1))
    end

    FFTKernel(fftsize::Integer) = FFTKernel{Float64}(fftsize)
end

"""
    fft!(x::AbstractVector{Complex{T}}, f::FFTKernel{T}) where T<:AbstractFloat

Perform Fast-Fourier Transform (FFT) in place on the signal sequence `x` using an FFT kernel `f`. The input sequence `x` will be overwritten by the FFT result.
"""
function fft!(x::VecI{Complex{T}}, f::FFTKernel{T}) where T<:AbstractFloat
    cache = f.cache
    if f.ifswap
        @simd for i in eachindex(cache)
            @inbounds cache[i] = x[i]
        end
        ditnn!(cache, x, f.twiddle, f.fftsize >> 1)
    else
        ditnn!(x, cache, f.twiddle, f.fftsize >> 1)
    end

    return x
end

"""
    fft(x::AbstractVector{Complex{T}}, f::FFTKernel{T}) where T<:AbstractFloat
    fft(x::AbstractVector{T},          f::FFTKernel{T}) where T<:AbstractFloat

Similar to `fft!`, the function will pass the **copy** of the given signal sequence to the `fft!` routine.
"""
fft(x::VecI{Complex{T}}, f::FFTKernel{T}) where T<:AbstractFloat = fft!(copy(x), f)

function fft(x::VecI{T}, f::FFTKernel{T}) where T<:AbstractFloat
    cx = similar(x, Complex{T})
    @simd for i in eachindex(cx)
        @inbounds cx[i] = x[i]
    end
    return fft!(cx, f)
end

end
