@testset "Jaws" begin

    @testset "Constructors" begin

        x, y = rand(2), rand(2)
        jaws = Jaws(x, y)

        @test jaws.x == x
        @test jaws.y == y
        
        @test jaws == Jaws(x[1], x[2], y[1], y[2])
        
        fieldsize = rand()
        jaws = Jaws(fieldsize)
        @test getx(jaws) == 0.5*fieldsize*[-1, 1]
        @test gety(jaws) == 0.5*fieldsize*[-1, 1]
    end

    @testset "Methods" begin
        x, y = rand(2), rand(2)
        jaws = Jaws(x, y) 

        @test getx(jaws) == x
        @test gety(jaws) == y
        @test getpositions(jaws) == (x, y)

        @test getarea(jaws) == (x[2]-x[1])*(y[2]-y[1])
    end

    @testset "Intersection (∩)" begin
        fieldsize = 20.
        jaws = Jaws(fieldsize)

        x2, y2 =[-9., 12.], [-6., 13.]
        jaws2 = Jaws(x2, y2) ∩ jaws
        @test jaws2 == Jaws(x2[1], x1[2], y2[1], y1[2])

        jaws2 = Jaws(0.9*fieldsize)
        @test (jaws2 ∩ jaws) == jaws2

        jaws2 = Jaws(1.1*fieldsize)
        @test (jaws2 ∩ jaws) == jaws

    end

end
