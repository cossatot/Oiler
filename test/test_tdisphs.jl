using Test

using Oiler

atol = 1e-10

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
    
    function test_TDdispHS_1()
        ue, un, uv = Oiler.TD.TDdispHS(X, Y, Z, P1, P2, P3, Ss, Ds; Ts=Ts, nu=nu)

        @test isapprox(ue, [-0.023647457823539948 -0.023647457823539948])
        @test isapprox(un, [-6.938893903907228e-17 -6.938893903907228e-17];atol=atol)
        @test isapprox(uv, [0.2793337025929701 0.2793337025929701])
    end


    function test_TDdisp_FS_1()
        ueMS, unMS, uvMS = Oiler.TD.TDdispFS(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, nu)
        @test isapprox(ueMS, [-0.03332290083010279 -0.03332290083010279])
        @test isapprox(unMS, [-2.0816681711721685e-17 -2.0816681711721685e-17];atol=atol)
        @test isapprox(uvMS, [0.12911857403885424 0.12911857403885424])
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

        @test isapprox(ueFSC, [0.039420629630340716 0.039420629630340716])
        @test isapprox(unFSC, [-2.0816681711721685e-17 -2.0816681711721685e-17])
        @test isapprox(uvFSC, [0.27688544590710584 0.27688544590710584])
    end

    function test_CoordTrans_arr_1()
        x1 = [-0.069091336142293 -0.069091336142293]
        x2 = [0.208166817117217e-16 0.208166817117217e-16]
        x3 = [-0.110256667352425 -0.110256667352425]

        A = [-0.703544597920105 0   0.710651109010207;
             0                  -1.000000000000000  0;
             0.710651109010207  0   0.703544597920105]

        X1, X2, X3 = Oiler.TD.CoordTrans(x1, x2, x3, A)
        @test isapprox(X1, [-0.029745186623777956 -0.029745186623777956])
        @test isapprox(X2, [-0.208166817117217e-16 -0.208166817117217e-16])
        @test isapprox(X3, [-0.12667031735299014 -0.12667031735299014])
    end
    
    function test_CoordTrans_scal_1()
        x1 = -0.069091336142293 
        x2 = 0.208166817117217e-16
        x3 = -0.110256667352425

        A = [-0.703544597920105 0   0.710651109010207;
             0                  -1.000000000000000  0;
             0.710651109010207  0   0.703544597920105]

        X1, X2, X3 = Oiler.TD.CoordTrans(x1, x2, x3, A)
        @test isapprox(X1, -0.02974518662377796)
        @test isapprox(X2, -0.208166817117217e-16)
        @test isapprox(X3, -0.126670317352990145)
    end

    function test_trimodefinder_1()
        x = [1., 1.]
        y = [0.874171929193456, 0.874171929193456]
        z = [0.527658448440079 0.527658448440079]
        p1 = [1., 1.40716026095111]
        p2 = [0., 0.]
        p3 = [2., 0.]

        tm = Oiler.TD.trimodefinder(x, y, z, p1, p2, p3)

        @test isapprox(tm, [1; 1])
    end

    function test_TDsetupD()
        x = [0.5276584484400789, 0.5276584484400789]
        y = [1.0, 1.0]
        z = [0.874171929193456, 0.874171929193456]
        alpha = 1.235677296317091
        bx = 0.
        by = 0.
        bz = 1.
        TriVertex = [0.0, 1.0, 1.4071602609511116]
        SideVec = [-0.0, -0.5792747270704139, 0.8151323761067877]

        u, v, w = Oiler.TD.TDSetupD(x, y, z, alpha, bx, by, bz, nu, TriVertex, 
            SideVec)

        @test isapprox(u, [0.04684803698063782, 0.04684803698063782])
        @test isapprox(v, [0.02183825650620384, 0.02183825650620384])
        @test isapprox(w, [0.026257913228717546, 0.026257913228717546])
    end

    function test_AngSetupFSC_1()
        bX = -0.7106511090102073
        bY = 0.
        bZ = 0.7035445979201053
        PA = [0.0, 1.0, -0.01]
        PB = [1.0, 0.0, -1.0]

        ue, un, uv = Oiler.TD.AngSetupFSC(X, Y, Z, bX, bY, bZ, PA, PB, nu)
        @test isapprox(ue, [-0.014164198784400879 -0.014164198784400879])
        @test isapprox(un, [-0.023154768178915285 -0.023154768178915285])
        @test isapprox(uv, [0.12756250832391175 0.12756250832391175])
    end

    function test_AngDisDisp_1()
        x = [0.5276584484400789, 0.5276584484400789]
        y = [-0.3087466704106312, -0.3087466704106312]
        z = [-0.4344560453028106, -0.4344560453028106]
        alpha = -1.9059153572727021
        bx = 0.
        by = [0.5792747270704139]
        bz = [0.8151323761067877]

        u, v, w = Oiler.TD.AngDisDisp(x, y, z, alpha, bx, by, bz, nu)

        @test isapprox(u, [0.04684803698063782, 0.04684803698063782])
        @test isapprox(v, [0.033011615434935415, 0.033011615434935415])
        @test isapprox(w, [0.008753325124405468, 0.008753325124405468])
    end

    function test_AngDisDispFSC_1a()
        y1 = [0.5303300858899106, 0.5303300858899106]
        y2 = [-0.5303300858899106, -0.5303300858899106]
        y3 = [0.0, 0.0]
        beta = -2.1815462594897985
        b1 = -0.502506218238858
        b2 = 0.502506218238858
        b3 = -0.7035445979201053
        a = 0.01

        v1, v2, v3 = Oiler.TD.AngDisDispFSC(y1, y2, y3, beta, b1, b2, b3, nu, a)

        @test isapprox(v1, [0.037872132405226674, 0.037872132405226674])
        @test isapprox(v2, [0.05849443652903753, 0.05849443652903753])
        @test isapprox(v3, [-0.11685100052124828, -0.11685100052124828])
    end

    function test_AngDisDispFSC_1b()
        y1 = Float64[]
        y2 = Float64[]
        y3 = Float64[]
        beta = 0.9600463940999948
        b1 = -0.502506218238858
        b2 = 0.502506218238858
        b3 = -0.7035445979201053
        a = 1.0

        v1, v2, v3 = Oiler.TD.AngDisDispFSC(y1, y2, y3, beta, b1, b2, b3, nu, a)
        
        @test isapprox(v1, Float64[])
        @test isapprox(v2, Float64[])
        @test isapprox(v3, Float64[])
    end


    @testset "test TDs against matlab" begin
        test_TDdispHS_1()
        test_TDdisp_FS_1()
        test_TDdisp_HarFunc_1()
        test_CoordTrans_arr_1()
        test_CoordTrans_scal_1()
        test_trimodefinder_1()
        test_TDsetupD()
        test_AngDisDisp_1()
        test_AngSetupFSC_1()
        test_AngDisDispFSC_1a()
        test_AngDisDispFSC_1b()

    end

end


function test_TD_point_ordering()
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

    ue, un, uv = Oiler.TD.TDdispHS(X, Y, Z, P1, P2, P3, Ss, Ds; Ts=Ts, nu=nu)

    ue231, un231, uv231 = Oiler.TD.TDdispHS(X, Y, Z, P2, P3, P1, Ss, Ds; Ts=Ts, nu=nu)
    ue312, un312, uv312 = Oiler.TD.TDdispHS(X, Y, Z, P3, P1, P2, Ss, Ds; Ts=Ts, nu=nu)
    
    ue321, un312, uv321 = Oiler.TD.TDdispHS(X, Y, Z, P3, P2, P1, Ss, Ds; Ts=Ts, nu=nu)
    
    ue213, un213, uv213 = Oiler.TD.TDdispHS(X, Y, Z, P2, P1, P3, Ss, Ds; Ts=Ts, nu=nu)
    ue132, un132, uv132 = Oiler.TD.TDdispHS(X, Y, Z, P1, P3, P2, Ss, Ds; Ts=Ts, nu=nu)
    
    @testset "test TD reorder" begin
        @test isapprox(ue, ue231)
        @test isapprox(ue, ue312)
    end
    
    @testset "test TD reverse" begin
        @test isapprox(ue, ue321)
    end

    @testset "test TD scramble" begin
        @test isapprox(ue, ue213)
        @test isapprox(ue, ue132)
    end
end
    









@testset "test TD unit tests" begin
    test_TD_1()
    test_TD_point_ordering()
end