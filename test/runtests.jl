using BasePlateFEA, Test, LinearAlgebra
const B = BasePlateFEA

@testset "BasePlateFEA" begin
    @testset "section" begin
        A, yCG = B.cs_data(2.8984, 2.8984, 0.6992, 0.1016)        # values of cs_data.m (Octave)
        @test A ≈ 1.02550976 rtol = 1e-8
        @test yCG ≈ 1.2338358088293568 rtol = 1e-8
    end

    @testset "spring law with unloading" begin
        law = B.elastic_hardening_spring(100.0, 10.0, 2.0)          # dy = 0.1
        @test B.spring_response(law, -1.0, 0.0) == (0.0, 0.0)
        @test B.spring_response(law, 0.05, 0.0)[1] ≈ 5.0
        f, k = B.spring_response(law, 0.6, 0.0)                     # hardening branch: 10 kip over 10 in
        @test f ≈ 10.5 && k ≈ 1.0
        f, k = B.spring_response(law, 0.55, 0.6)                    # elastic unloading from d = 0.6
        @test f ≈ 10.5 - 100 * 0.05 && k ≈ 100.0
        @test B.spring_response(law, 0.4, 0.6) == (0.0, 0.0)        # fully unloaded
        @test B.spring_response(law, 2000.0, 0.0)[1] ≈ 20.0         # constant beyond the last point
    end

    inp = BasePlateInput(geometry = barcelona_geometry(:C), fe_size = 0.3)
    mesh = generate_mesh(inp)

    @testset "mesh" begin
        g = inp.geometry
        @test sum(mesh.bp_area) ≈ g.L * g.W rtol = 1e-10
        @test length(mesh.bolt_nodes) == 2 && all(length(h) == g.hole_sectors for h in mesh.hole_nodes)
        @test count(mesh.hole) > 0
        @test all(abs(mesh.x[i, 3]) < 1e-12 for i in mesh.bp_nodes)
        # every element has positive area and the wall elements are 4-node
        @test all(length(mesh.conn[e]) == (mesh.part[e] == 1 ? 3 : 4) for e in 1:B.nelements(mesh))
    end

    # global equilibrium of a converged state: vertical force and moment about the plate x axis
    function check_equilibrium(res, N)
        m = res.model.mesh; f = res.frames[end]
        contact = f.contact_stress .* m.bp_area
        @test sum(contact) - sum(f.bolt_force) ≈ N rtol = 1e-5
        yc = inp.geometry.W / 2
        Mx = sum(contact[i] * (m.x[n, 2] - yc) for (i, n) in enumerate(m.bp_nodes)) -
             sum(f.bolt_force[i] * (m.x[n, 2] - yc) for (i, n) in enumerate(m.bolt_nodes))
        @test abs(Mx) ≈ f.M rtol = 0.03          # the in-plane springs and the penalty drilling stiffness carry a little
    end

    @testset "elastic analysis" begin
        res = analyze(inp; mesh, N = 10.0, Mmax = 16.0, verbose = false)
        @test res.status == :completed
        θ, M = moment_rotation(res)
        @test issorted(M) && issorted(θ)
        @test 5000 < M[2] / θ[2] < 10000                           # initial stiffness, kip-in/rad
        @test M[end] / θ[end] < 0.5 * M[2] / θ[2]                  # softening after uplift
        check_equilibrium(res, 10.0)
    end

    @testset "elastic-plastic analysis" begin
        inp_p = BasePlateInput(geometry = barcelona_geometry(:C), fe_size = 0.3,
                               plate_steel = bilinear_steel(53.2, 64.3, 0.15), upright_steel = bilinear_steel(68.9, 77.8, 0.10))
        res_e = analyze(inp; mesh, N = 10.0, Mmax = 100.0, θmax = 0.02, verbose = false)
        res_p = analyze(inp_p; mesh, N = 10.0, Mmax = 100.0, θmax = 0.02, verbose = false)
        @test res_p.status == :max_rotation
        @test maximum(x -> isnan(x) ? 0.0 : x, res_p.frames[end].fields[:peeq]) > 1e-3
        @test res_p.frames[end].M < res_e.frames[end].M              # yielding lowers the moment at the same rotation
        check_equilibrium(res_p, 10.0)
        d = viewer_data(res_p)
        @test length(d[:frames]) == length(res_p.frames) || length(d[:frames]) <= 43
        path = write_viewer(joinpath(mktempdir(), "v.html"), res_p)
        @test filesize(path) > 100_000
    end
    @testset "tall upright, steel curve, buckling modes, imperfection" begin
        st = yun_gardner_steel(64.45, 71.8)
        @test st.εp[1] == 0 && issorted(st.εp) && issorted(st.σ)
        @test st.σ[end] ≈ 71.8 * 1.06 rtol = 5e-3                    # true stress at εu = 0.06
        inp_t = BasePlateInput(geometry = barcelona_geometry(:C), fe_size = 0.3, H = 12.0, H_fine = 4.5, fe_size_upper = 0.6,
                               upright_plastic_height = 4.5, top_restraint = :planar, Hm = 85 / 25.4,
                               upright_steel = st)
        mesh_t = generate_mesh(inp_t)
        @test maximum(mesh_t.x[:, 3]) ≈ 12.0
        @test mesh_t.x[mesh_t.rot_nodes[1], 3] ≈ 3.3
        model = Model(inp_t; mesh = mesh_t)
        λ, Φ = buckling_modes(model; N = 10.0, M = 20.0, nmodes = 2)
        @test 1 < λ[1] <= λ[2] < 20
        mesh_i = apply_imperfection(mesh_t, Φ[:, 1], 0.03)
        @test maximum(sqrt.(sum(abs2, mesh_i.x .- mesh_t.x, dims = 2))) ≈ 0.03
        @test mesh_i.x[mesh_t.bp_nodes, :] == mesh_t.x[mesh_t.bp_nodes, :]
        res = analyze(Model(inp_t; mesh = mesh_i); N = 10.0, Mmax = 100.0, θmax = 0.004, verbose = false)
        @test res.status == :max_rotation
        f = res.frames[end]
        @test base_moment(res.model, f) > f.M                          # the axial force adds N·δ
        @test moment_rotation(res; second_order = true)[2][end] ≈ base_moment(res.model, f)
    end
end
