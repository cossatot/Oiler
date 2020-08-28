using Revise
using Test

using Oiler

function test_TD_1()
    # all functions checked against original matlab
    X = [0.75 0.75];
    Y = [1. 1.];
    Z = [-0.01 -0.01];

    P1 = [0., 1., -0.01];
    P2 = [1., 0., -1.];
    P3 = [1., 2., -1.];

    Ss = 0.;
    Ds = 1.;
    Ts = 0.;
    nu = 0.25;

    function test_TDdisp_FS_1()
        ueMS, unMS, uvMS = Oiler.TD.TDdispFS(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, nu)
        @test ueMS == [-0.03332290083010279 -0.03332290083010279] 
        @test unMS == [-2.0816681711721685e-17 -2.0816681711721685e-17]
        @test uvMS == [0.12911857403885424 0.12911857403885424] 
    end

    function test_TDdisp_HarFunc_1()
        X = X[:]
        Y = Y[:]
        Z = Z[:]

        P1 = P1[:]
        P2 = P2[:]
        P3 = P3[:]

        ueFSC, unFSC, uvFSC = Oiler.TD.TDdisp_HarFunc(X, Y, Z, P1, P2, P3, Ss, 
                                                      Ds, Ts, nu)

        @test ueFSC == [0.039420629630340716 0.039420629630340716]
        @test unFSC == [-2.0816681711721685e-17 -2.0816681711721685e-17]
        @test uvFSC == [0.27688544590710584 0.27688544590710584]
    end

    function test_CoordTrans_1()
        x1 = [-0.069091336142293 -0.069091336142293]
        x2 = [0.208166817117217e-16 0.208166817117217e-16]
        x3 = [-0.110256667352425 -0.110256667352425]

        A = [-0.703544597920105 0   0.710651109010207;
             0                  -1.000000000000000  0;
             0.710651109010207  0   0.703544597920105]

        X1, X2, X3 = Oiler.TD.CoordTrans(x1, x2, x3, A)
        @test X1 == [-0.029745186623777956 -0.029745186623777956]
        @test X2 == [-0.208166817117217e-16 -0.208166817117217e-16]
        @test X3 == [-0.12667031735299014 -0.12667031735299014]
    end

    function test_TDdispHS_1()
        ue, un, uv = Oiler.TD.TDdispHS(X, Y, Z, P1, P2, P3, Ss, Ds; Ts=Ts, nu=nu)

        @test ue == [-0.023647457823539948 -0.023647457823539948]
        @test un == [-6.938893903907228e-17 -6.938893903907228e-17]
        @test uv == [0.2793337025929701 0.2793337025929701]
    end

    @testset begin
        test_TDdisp_FS_1()
        test_TDdisp_HarFunc_1()
        test_CoordTrans_1()
        test_TDdispHS_1()

    end

end

@testset "test TD unit tests" begin
    test_TD_1()
end