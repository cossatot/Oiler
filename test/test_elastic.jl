using Test
using PyPlot

using Oiler

# set up
simple_fault_trace = [1. 1.;
                      1.5 1.5];

simple_fault = Oiler.Fault(trace = simple_fault_trace, dip = 45., 
                           dip_dir = "SE", usd = 1., lsd = 20.,
                           hw = "a", fw = "b")

gnss_1 = Oiler.VelocityVectorSphere(lond = 0., latd = 0., ve = 1., vn = 1., 
    vel_type = "gnss", fix = "a", mov = "r")
gnss_2 = Oiler.VelocityVectorSphere(lond = 1.25, latd = 1.251, ve = 1., vn = 1.,
    vel_type = "gnss", fix = "a", mov = "r")
gnss_3 = Oiler.VelocityVectorSphere(lond = 1.25, latd = 1.249, ve = 1., vn = 1.,
    vel_type = "gnss", fix = "a", mov = "r")

gnss_lons, gnss_lats = Oiler.get_coords_from_vel_array([gnss_1; gnss_2; gnss_3])
    

function test_fault_to_okada()

    x, y = Oiler.Faults.fault_oblique_merc(simple_fault, gnss_lons, gnss_lats)
    gnss_x, gnss_y = x[1:3], y[1:3]

    sx1, sy1, sx2, sy2 = x[4], y[4], x[5], y[5]

    # testing setup and oblique mercator projection, not fault_to_okada
    # values from ProjectSegCoords.m from Meade and Loveless
    @test isapprox([sx1 sy1; sx2 sy2], 
        [-9.850338894354840e+03 -1.138348096341169e-12; 
         -9.771721687583831e+03 -1.458853871039167e-12])

    @test isapprox(gnss_x, 
        1e4 .* [-1.000758826921680; -0.981095071914695; -0.981110799198045])
    
    @test isapprox(gnss_y, 
        [0.044903674262838; 0.075810489054318; -0.081423421703794])

    # now test fault_to_okada vs. fault_params_to_okada_form.m by M+L
    D = Oiler.fault_to_okada(simple_fault, sx1, sy1, sx2, sy2)

    @test isapprox(D["strike"], 6.283185307179583)
    @test isapprox(D["L"], 78.617206771008568)
    @test isapprox(D["W"], 26.870057685088810)
    @test isapprox(D["ofx"], -9.850338894354840e+03)
    @test isapprox(D["ofy"], -20.000000000001140)
    @test isapprox(D["ofxe"], -9.771721687583831e+03)
    @test isapprox(D["ofye"], -20.000000000001464)
    @test isapprox(D["tfx"], -9.850338894354840e+03)
    @test isapprox(D["tfy"], -1.000000000001139)
    @test isapprox(D["tfxe"], -9.771721687583831e+03)
    @test isapprox(D["tfye"], -1.000000000001459)
end

test_fault_to_okada()


function test_okada_dip_slip()
    x, y = Oiler.Faults.fault_oblique_merc(simple_fault, gnss_lons, gnss_lats)
    gnss_x, gnss_y = x[1:3], y[1:3]

    sx1, sy1, sx2, sy2 = x[4], y[4], x[5], y[5]

    D = Oiler.fault_to_okada(simple_fault, sx1, sy1, sx2, sy2)
    ves, vns, vus, ved, vnd, vud, vet, vnt, vut = Oiler.okada(D, 0., 1., 0.,
        gnss_x, gnss_y)
    
    # test values from using okada_partials.m by Meade and Loveless
    @test isapprox(ves, [0; 0; 0])
    @test isapprox(vns, [0; 0; 0])
    @test isapprox(vus, [0; 0; 0])
    @test isapprox(ved,    
        1.0e-03 .* [-0.243919705202018; 0.035030039954510;-0.034934183707360])
    @test isapprox(vnd,    
        [-0.000199992219480; -0.009262479825412; 0.024227041358989])
    @test isapprox(vud,    
        [-0.001992787293162; 0.023204884580011; 0.057062430886315])
    @test isapprox(vet, [0; 0; 0])
    @test isapprox(vnt, [0; 0; 0])
    @test isapprox(vut, [0; 0; 0])
end

test_okada_dip_slip()


function test_okada_strike_slip()
    x, y = Oiler.Faults.fault_oblique_merc(simple_fault, gnss_lons, gnss_lats)
    gnss_x, gnss_y = x[1:3], y[1:3]

    sx1, sy1, sx2, sy2 = x[4], y[4], x[5], y[5]

    D = Oiler.fault_to_okada(simple_fault, sx1, sy1, sx2, sy2)
    ves, vns, vus, ved, vnd, vud, vet, vnt, vut = Oiler.okada(D, 1., 0., 0.,
        gnss_x, gnss_y)
    
    # test values from using okada_partials.m by Meade and Loveless
    @test isapprox(ves,    
        [0.000278311207785; -0.019970149820352; 0.003279330179775])
    @test isapprox(vns,    
        [0.002821784569174; -0.000074894257506; 0.000072978696014])
    @test isapprox(vus,    
        1.0e-3 .* [-0.316155394920109; 0.038678751237683; -0.037786443714757])
    @test isapprox(ved, [0; 0; 0])
    @test isapprox(vnd, [0; 0; 0])
    @test isapprox(vud, [0; 0; 0])
    @test isapprox(vet, [0; 0; 0])
    @test isapprox(vnt, [0; 0; 0])
    @test isapprox(vut, [0; 0; 0])
end

test_okada_strike_slip()


function test_okada_tensile()
    x, y = Oiler.Faults.fault_oblique_merc(simple_fault, gnss_lons, gnss_lats)
    gnss_x, gnss_y = x[1:3], y[1:3]

    sx1, sy1, sx2, sy2 = x[4], y[4], x[5], y[5]

    D = Oiler.fault_to_okada(simple_fault, sx1, sy1, sx2, sy2)
    ves, vns, vus, ved, vnd, vud, vet, vnt, vut = Oiler.okada(D, 0., 0., 1.,
        gnss_x, gnss_y)
    
    # test values from using okada_partials.m by Meade and Loveless
    @test isapprox(ves, [0; 0; 0])
    @test isapprox(vns, [0; 0; 0])
    @test isapprox(vus, [0; 0; 0])
    @test isapprox(ved, [0; 0; 0])
    @test isapprox(vnd, [0; 0; 0])
    @test isapprox(vud, [0; 0; 0])
    @test isapprox(vet,    
        1.0e-03 .* [0.244273925530614; -0.035436156199882; 0.034506189911012])
    @test isapprox(vnt,    
        [0.000199966556297; -0.006285547974839; -0.006095019820654])
    @test isapprox(vut,    
        [0.001992761716366; -0.038169618939615; -0.038100199139836])
end

test_okada_tensile()

function test_okada_partials()
    x, y = Oiler.Faults.fault_oblique_merc(simple_fault, gnss_lons, gnss_lats)
    gnss_x, gnss_y = x[1:3], y[1:3]

    sx1, sy1, sx2, sy2 = x[4], y[4], x[5], y[5]

    D = Oiler.fault_to_okada(simple_fault, sx1, sy1, sx2, sy2)
    ves, vns, vus, ved, vnd, vud, vet, vnt, vut = Oiler.okada(D, 1., 1., 1.,
        gnss_x, gnss_y)
    
    # test values from using okada_partials.m by Meade and Loveless
    @test isapprox(ves,    
        [0.000278311207785; -0.019970149820352; 0.003279330179775])
    @test isapprox(vns,    
        [0.002821784569174; -0.000074894257506; 0.000072978696014])
    @test isapprox(vus,    
        1.0e-3 .* [-0.316155394920109; 0.038678751237683; -0.037786443714757])
    @test isapprox(ved,    
        1.0e-03 .* [-0.243919705202018; 0.035030039954510;-0.034934183707360])
    @test isapprox(vnd,    
        [-0.000199992219480; -0.009262479825412; 0.024227041358989])
    @test isapprox(vud,    
        [-0.001992787293162; 0.023204884580011; 0.057062430886315])
    @test isapprox(vet,    
        1.0e-03 .* [0.244273925530614; -0.035436156199882; 0.034506189911012])
    @test isapprox(vnt,    
        [0.000199966556297; -0.006285547974839; -0.006095019820654])
    @test isapprox(vut,    
        [0.001992761716366; -0.038169618939615; -0.038100199139836])
end

test_okada_partials()