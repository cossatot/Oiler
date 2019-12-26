using CSV

using Oiler

function check_closure(pole1::Oiler.EulerPoleSphere,
    pole2::Oiler.EulerPoleSphere, pole3::Oiler.EulerPoleSphere;
    tol::Float64=1e-10)
    

    p12 = pole2 - pole1
    p23 = pole3 - pole2
    p31 = pole1 - pole3

    p_sum = p12 + p23 + p31

    if p_sum.rotrate > tol
        p1 = pole1.mov
        p2 = pole2.mov
        p3 = pole3.mov

        println("problem with ($p1, $p2, $p3): $p_sum")
    end
end

pole_df = CSV.read("./test_data/poles.PA.txt")

function row_to_pole(row)
    Oiler.EulerPoleSphere(latd=row.lat, lond=row.lon, rotrate=row.rotrate,
    fix="PA", mov=row.plate)
end

poles = [row_to_pole(row) for row in eachrow(pole_df)]

pole_triples = []

np = length(poles)


for p1 in 1:np
    for p2 in p1+1:np
        for p3 in p2+1:np
            tri = Set([p1, p2, p3])
            if !(tri in pole_triples)
                push!(pole_triples, tri)
            end
        end
    end
    println(p1)
end

pole_triples = [collect(p) for p in pole_triples]



for (i, pt) in enumerate(pole_triples)
    check_closure(poles[pt[1]], poles[pt[2]], poles[pt[3]])

    if i % 1000 == 0
        println(i)
    end
end