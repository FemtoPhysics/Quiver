import Quiver: FFT

@testset "Test Twiddle Factors" begin
    @test FFT.twiddle!(Vector{ComplexF64}(undef, 1)) == ComplexF64[1.0 + 0.0im]
    @test FFT.twiddle!(Vector{ComplexF64}(undef, 2)) == ComplexF64[1.0 + 0.0im, 0.0 - 1.0im]
    let tmp = 0.5 * sqrt(2.0)
        @test FFT.twiddle!(Vector{ComplexF64}(undef, 4)) == ComplexF64[
            1.0 + 0.0im, complex(tmp, -tmp), 0.0 - 1.0im, complex(-tmp, -tmp)
        ]
    end
    @test FFT.twiddle!(Vector{ComplexF64}(undef, 8)) ≈ ComplexF64[
        1.0 + 0.0im,
        0.9238795325112867 - 0.3826834323650898im,
        0.7071067811865476 - 0.7071067811865476im,
        0.3826834323650898 - 0.9238795325112867im,
        0.0 - 1.0im,
        -0.3826834323650898 - 0.9238795325112867im,
        -0.7071067811865476 - 0.7071067811865476im,
        -0.9238795325112867 - 0.3826834323650898im
    ]
end

@testset "Radix-2 FFT" begin
    dat = [1.0 + 0.0im,  2.0 - 1.0im, 0.0 - 1.0im, -1.0 + 2.0im]
    ans = [2.0 + 0.0im, -2.0 - 2.0im, 0.0 - 2.0im,  4.0 + 4.0im]

    @test ans == FFT.fft(dat, FFT.FFTKernel(4))
end
