using OpticalFibers
using OpticalFibers.ModeSolvers
using Test


@testset "ArrayFields" begin
     # Définir une grille et une fonction
    x = (-1:0.1:1)u"m"
    f(x, y) = VectorValue((x/1u"m")^2, (y/1u"m")^2)u"V"
    af = ArrayField([x, x], f.(x, x'))

    # Tester l'évaluation de la fonction
    point = Point(0.5u"m", 1.0u"m")
    result = af(point)
    @test result == VectorValue(0.25u"V", 1.0u"V")

    # Tester une opération scalaire
    scaled_af = 2 * af
    scaled_result = scaled_af(point)
    @test scaled_result == VectorValue(0.5u"V", 2.0u"V")

    # Tester une opération vectorielle (dot product)
    v = VectorValue(1.0u"V", 1.0u"V")
    dot_af = dot(af, v)
    dot_result = dot_af(point)
    @test dot_result == 1.25u"V^2"

    # Tester une multiplication entre deux ArrayFields
    g(x, y) = VectorValue(x/1u"m",y/1u"m")u"A"
    af2 = ArrayField([x, x], g.(x, x'))
    product_af = af⋅af2
    product_result = product_af(point)
    @test product_result == 1.125u"W"
end

@testset "FunctionField" begin
    # Définir une fonction simple
    f(x) = VectorValue(x[1]^2/1u"m^2", x[2]^2/1u"m^2")u"V"
    ff = FunctionField(2, f)

    # Tester l'évaluation de la fonction
    point = VectorValue(2.0u"m", 3.0u"m")
    result = ff(point)
    @test result == VectorValue(4.0u"V", 9.0u"V")

    # Tester une opération scalaire
    scaled_ff = 2 * ff
    scaled_result = scaled_ff(point)
    @test scaled_result == VectorValue(8.0u"V", 18.0u"V")

    # Tester une opération vectorielle (dot product)
    v = VectorValue(1.0u"V", 1.0u"V")
    dot_ff = dot(ff, v)
    dot_result = dot_ff(point)
    @test dot_result == 13.0u"V^2"

    # Tester une multiplication entre deux FunctionFields
    f2(x) = VectorValue(x[1]/1u"m",x[2]/1u"m")u"A"
    ff2 = FunctionField(2, f2)
    product_ff = ff⋅ff2
    result = product_ff(point)
    @test result == 35u"W"
end

@testset "FEMField" begin
    domain = (-1, 1, -1, 1)
    cell_nums = (101, 101)
    model = CartesianDiscreteModel(domain, cell_nums) |> simplexify
    Ω = Triangulation(model)
    f(x) = VectorValue(x[1]^2, x[2]^2)
    cf=CellField(f,Ω)
    dΩ = Measure(Ω,4)
    ff=FEMField(cf, dΩ, u"m",u"V")
    # Test the evaluation of the function
    point = Point(0.5u"m", 1.0u"m")
    result = ff(point)
    @test result == VectorValue(0.25u"V", 1.0u"V")
   

    # Tester une opération scalaire
    scaled_ff = 2 * ff
    scaled_result = scaled_ff(point)
    @test scaled_result == VectorValue(0.5u"V", 2.0u"V")

    # Tester une opération vectorielle (dot product)
    v = VectorValue(1.0u"V", 1.0u"V")
    dot_ff = dot(ff, v)
    dot_result = dot_ff(point)
    @test dot_result == 1.25u"V^2"

    # Tester une multiplication entre deux FunctionFields
    f2(x) = VectorValue(x[1], x[2])
    cf2=CellField(f2,Ω)
    ff2=FEMField(cf2, dΩ, u"m",u"A")
    product_ff = ff⋅ff2
    result = product_ff(point)
    @test result == 1.125u"W"
end

@testset "Multistep" begin
    # Bures p 133
    mode0=multi_step_fiber_modes(0.6328u"µm",0,8.335u"µm",[1.462420,1.457420]);
    neff0=getproperty.(mode0,:neff);
    U0=2*pi*8.335/0.6328*sqrt.(1.462420^2 .-neff0.^2);
    @test isapprox(U0,[2.1845,4.9966,7.7642],atol=1E-4);

    mode3=multi_step_fiber_modes(0.6328u"µm",3,8.335u"µm",[1.462420,1.457420]);
    neff3=getproperty.(mode3,:neff);
    U3=2*pi*8.335/0.6328*sqrt.(1.462420^2 .-neff3.^2);
    @test isapprox(U3,[5.7740,8.7290],atol=1E-4);

end

@testset "FEM1D" begin
    modeFD=FEM1D(0.6328u"µm",0,x->1.462420.-(1.462420+1.457420)*(x[1]>=8.335u"µm"),CartesianDiscreteModel((0,50),20000)*u"µm",neigs=3,order=3);
    neffFD=getproperty.(modeFD,:neff);
    UFD=2*pi*8.335/0.6328*sqrt.(1.462420^2 .-neffFD.^2);
    @test isapprox(UFD,[2.1845,4.9966,7.7642],atol=1E-4);
end


