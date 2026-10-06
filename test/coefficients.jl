@testset "SC representation" begin

    coeffs = allocate_coefficients(
        Float64,
        3,
        SC(),
    )

    for n in 0:3
        for m in 0:n
            setC!(coeffs, n, m, 100n + m)

            if m > 0
                setS!(coeffs, n, m, -(100n + m))
            end
        end
    end

    for n in 0:3
        for m in 0:n
            @test C(coeffs, n, m) == 100n + m

            if m == 0
                @test_throws ArgumentError S(coeffs, n, m)
            else
                @test S(coeffs, n, m) == -(100n + m)
            end
        end
    end
    
    expected = [NaN   NaN   NaN     0.0   NaN   NaN   NaN
                NaN   NaN  -101.0 100.0 101.0   NaN   NaN
                NaN  -202.0 -201.0 200.0 201.0 202.0   NaN
            -303.0 -302.0 -301.0 300.0 301.0 302.0 303.0]
    @test isequal(coeffs.data, expected)
end

@testset "CS representation" begin

    coeffs = allocate_coefficients(
        Float64,
        3,
        CS(),
    )

    for n in 0:3
        for m in 0:n
            setC!(coeffs, n, m, 100n + m)

            if m > 0
                setS!(coeffs, n, m, -(100n + m))
            end
        end
    end

    for n in 0:3
        for m in 0:n
            @test C(coeffs, n, m) == 100n + m

            if m == 0
                @test_throws ArgumentError S(coeffs, n, m)
            else
                @test S(coeffs, n, m) == -(100n + m)
            end
        end
    end
    
    expected = [0.0  -101.0  -201.0  -301.0
                100.0   101.0  -202.0  -302.0
                200.0   201.0   202.0  -303.0
                300.0   301.0   302.0   303.0]
    @test isequal(coeffs.data, expected)
end

