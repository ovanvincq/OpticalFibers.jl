using OpticalFibers
using OpticalFibers.ModeSolvers
using Test

@testset "Multistep" begin
    # Bures p 133
    mode0=multi_step_fiber_modes(0.6328,0,8.335,[1.462420,1.457420]);
    neff0=getproperty.(mode0,:neff);
    U0=2*pi*8.335/0.6328*sqrt.(1.462420^2 .-neff0.^2);
    @test isapprox(U0,[2.1845,4.9966,7.7642],atol=1E-4);

    mode3=multi_step_fiber_modes(0.6328,3,8.335,[1.462420,1.457420]);
    neff3=getproperty.(mode3,:neff);
    U3=2*pi*8.335/0.6328*sqrt.(1.462420^2 .-neff3.^2);
    @test isapprox(U3,[5.7740,8.7290],atol=1E-4);

end

@testset "ScalarFD" begin
    modeFD=FD(0.6328,0,3,x->1.462420.-(1.462420+1.457420)*(x.>=8.335),10000,200,order=200);
    neffFD=getproperty.(modeFD,:neff);
    UFD=2*pi*8.335/0.6328*sqrt.(1.462420^2 .-neffFD.^2);
    @test isapprox(UFD,[2.1845,4.9966,7.7642],atol=1E-4);
end


