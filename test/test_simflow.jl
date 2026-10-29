using LegendHDF5IO
using Refit
using StatsBase
using Test

@testset "SimProd" begin
    # a fake production cycle: two simulations, two LAr and two A/E categories,
    # the HPGe groups "all" and "icpc" and the 2D histogram mul2 in the group "all"
    dir = mktempdir()
    write(
        joinpath(dir, "simflow-config.yaml"),
        """
        experiment: l1000test
        paths:
          tier:
            pdf: \$_/generated/tier/pdf
        """,
    )
    pdf_dir = mkpath(joinpath(dir, "generated", "tier", "pdf"))
    touch(joinpath(pdf_dir, "other-file.lh5"))

    binning = (binedges = (first = 0.0, last = 3.0, step = 1.0), closedleft = true)
    # the counts of each category, the "icpc" group has half of them
    counts = Dict(
        "is_spms_coincident/is_ss_aoe" => [1.0, 0.0, 0.0],
        "is_spms_coincident/is_low_aoe" => [0.0, 2.0, 0.0],
        "not_spms_coincident/is_ss_aoe" => [0.0, 0.0, 4.0],
        "not_spms_coincident/is_low_aoe" => [8.0, 0.0, 0.0],
    )
    for (simid, n_prim) in ("bulk_2vbb" => 10, "lar_Ar39" => 100)
        lh5open(joinpath(pdf_dir, "l1000test-$simid-tier_pdf.lh5"), "w") do f
            for (path, w) in counts, (group, scale) in ("all" => 1, "icpc" => 0.5)
                f["pdf/mul1/$path/$group/binning/axis_0"] = binning
                f["pdf/mul1/$path/$group/weights"] = scale .* w
            end
            f["pdf/mul2/all"] = Histogram((0.0:1.0:3.0, 0.0:1.0:2.0), [1.0 2; 3 4; 5 6])
            f["nr_sim_events"] = n_prim
        end
    end

    prod = SimProd(dir)
    @test prod.experiment == "l1000test"
    @test prod.pdf_dir == pdf_dir
    @test simids(prod) == [:bulk_2vbb, :lar_Ar39]
    @test pdf_file(prod, :bulk_2vbb) ==
          joinpath(pdf_dir, "l1000test-bulk_2vbb-tier_pdf.lh5")
    @test_throws ArgumentError pdf_file(prod, :missing)
    @test_throws ArgumentError SimProd(pdf_dir)

    # the histograms of the first file
    @test length(pdf_groups(prod)) == 9
    @test "mul1/not_spms_coincident/is_ss_aoe/icpc" in pdf_groups(prod)
    @test "mul2/all" in pdf_groups(prod)
    @test read_hist("pdf/" * pdf_groups(prod)[1], pdf_file(prod, :bulk_2vbb)) isa Histogram

    @test sort(cut_groups(prod)) == [
        [:is_low_aoe, :is_spms_coincident, :mul1],
        [:is_low_aoe, :mul1, :not_spms_coincident],
        [:is_spms_coincident, :is_ss_aoe, :mul1],
        [:is_ss_aoe, :mul1, :not_spms_coincident],
        [:mul2],
    ]
    @test sort(hpge_groups(prod)) == [:all, :icpc]
    @test sprint(show, prod) == "SimProd($(repr(dir)))"
    @test eval(Meta.parse(sprint(show, prod))).pdf_dir == pdf_dir
    @test sprint(show, MIME"text/plain"(), prod) == """
          SimProd at $dir
          ├─ experiment: l1000test
          ├─ 2 simulations:
          │    bulk_2vbb
          │    lar_Ar39
          ├─ cut groups:
          │    is_low_aoe ∧ is_spms_coincident ∧ mul1
          │    is_spms_coincident ∧ is_ss_aoe ∧ mul1
          │    is_low_aoe ∧ mul1 ∧ not_spms_coincident
          │    is_ss_aoe ∧ mul1 ∧ not_spms_coincident
          │    mul2
          └─ HPGe groups: all, icpc"""

    # all multiplicity-1 events, divided by the number of primaries
    pdf = read_pdf(prod, :bulk_2vbb)
    @test pdf isa Histogram
    @test pdf.edges == (0.0:1.0:3.0,)
    @test pdf.weights ≈ [9, 2, 4] ./ 10
    @test read_pdf(prod, :bulk_2vbb; hpge_group = :icpc).weights ≈ [4.5, 1, 2] ./ 10

    # a list of categories in any order, any category of the other cuts
    sel(groups) = read_pdf(prod, :bulk_2vbb; groups).weights .* 10
    @test sel(:mul1) ≈ [9, 2, 4]
    @test sel(:not_spms_coincident) ≈ [8, 0, 4]
    @test sel(:is_ss_aoe) ≈ [1, 0, 4]
    @test sel([:is_ss_aoe]) ≈ [1, 0, 4]
    @test sel([:not_spms_coincident, :is_ss_aoe]) ≈ [0, 0, 4]
    @test sel([:is_ss_aoe, :not_spms_coincident]) ≈ [0, 0, 4]
    @test sel([:is_ss_aoe, :mul1, :not_spms_coincident]) ≈ [0, 0, 4]
    # any of several lists of categories
    @test sel([[:is_ss_aoe], [:is_spms_coincident]]) ≈ [1, 2, 4]
    @test sel([[:is_ss_aoe, :is_spms_coincident], [:not_spms_coincident]]) ≈ [9, 0, 4]
    # a histogram matched by two lists is counted once
    @test sel([[:mul1], [:is_ss_aoe]]) ≈ [9, 2, 4]
    # an empty list has no cuts
    icpc(groups) = read_pdf(prod, :bulk_2vbb; groups, hpge_group = :icpc).weights
    @test icpc(Symbol[]) ≈ icpc(:mul1)
    @test icpc([]) ≈ icpc(:mul1)
    # in the group all, no cuts also selects the 2D histogram mul2, which cannot be added
    @test_throws ErrorException sel([])

    @test_throws ArgumentError read_pdf(prod, :bulk_2vbb; groups = :is_high_aoe)
    @test_throws ArgumentError read_pdf(prod, :bulk_2vbb; hpge_group = :bege)
    # categories that exclude each other select no event
    @test_throws ArgumentError read_pdf(
        prod,
        :bulk_2vbb;
        groups = [:is_ss_aoe, :is_low_aoe],
    )
    # strings are not categories
    @test_throws TypeError read_pdf(prod, :bulk_2vbb; groups = "is_ss_aoe")

    # several HPGe groups are summed
    @test read_pdf(prod, :bulk_2vbb; hpge_group = [:all, :icpc]).weights ≈
          [13.5, 3, 6] ./ 10
    @test read_pdf(prod, :bulk_2vbb; hpge_group = [:icpc]).weights ≈ [4.5, 1, 2] ./ 10
    @test read_pdf(prod, :bulk_2vbb; hpge_group = [:icpc, :bege]).weights ≈
          [4.5, 1, 2] ./ 10
    # has_pdf needs every HPGe group, read_pdf sums those it finds
    @test has_pdf(prod, :bulk_2vbb; hpge_group = [:icpc, :all])
    @test !has_pdf(prod, :bulk_2vbb; hpge_group = [:icpc, :bege])
    @test !has_pdf(prod, :bulk_2vbb; hpge_group = [:bege, :coax])
    # empty vectors: no HPGe group to check, no cuts
    @test has_pdf(prod, :bulk_2vbb; hpge_group = [])
    @test has_pdf(prod, :bulk_2vbb; hpge_group = Symbol[], groups = :is_high_aoe)
    @test !has_pdf(prod, :missing; hpge_group = [])
    @test has_pdf(prod, :bulk_2vbb; groups = [])
    @test has_pdf(prod, :bulk_2vbb; groups = Symbol[], hpge_group = :icpc)

    # check without reading
    @test has_pdf(prod, :bulk_2vbb)
    @test !has_pdf(prod, :missing)
    @test !has_pdf(prod, :missing; groups = :mul1)
    @test has_pdf(prod, :bulk_2vbb; groups = [:is_ss_aoe, :not_spms_coincident])
    @test has_pdf(prod, :bulk_2vbb; groups = [:not_spms_coincident, :is_ss_aoe])
    @test has_pdf(prod, :bulk_2vbb; groups = :is_ss_aoe, hpge_group = :icpc)
    @test has_pdf(prod, :bulk_2vbb; hpge_group = :icpc)
    @test has_pdf(prod, :bulk_2vbb; groups = [[:is_high_aoe], [:is_ss_aoe]])
    @test !has_pdf(prod, :bulk_2vbb; groups = [:is_ss_aoe, :is_low_aoe])
    @test !has_pdf(prod, :bulk_2vbb; groups = :is_high_aoe)
    @test !has_pdf(prod, :bulk_2vbb; hpge_group = :bege)
    # the 2D histogram mul2
    @test has_pdf(prod, :bulk_2vbb; groups = :mul2)
    @test has_pdf(prod, :bulk_2vbb; groups = :mul2, hpge_group = :all)
    @test !has_pdf(prod, :bulk_2vbb; groups = :mul2, hpge_group = :icpc)
    mul2 = read_pdf(prod, :bulk_2vbb; groups = :mul2)
    @test mul2.edges == (0.0:1.0:3.0, 0.0:1.0:2.0)
    @test mul2.weights ≈ [1 2; 3 4; 5 6] ./ 10
    @test_throws ArgumentError read_pdf(
        prod,
        :bulk_2vbb;
        groups = :mul2,
        hpge_group = :icpc,
    )
    # mul1 does not select mul2
    @test ndims(read_pdf(prod, :bulk_2vbb; groups = :mul1).weights) == 1
    @test has_pdf(prod, :bulk_2vbb, "mul1/not_spms_coincident/is_ss_aoe/icpc")
    @test has_pdf(prod, :bulk_2vbb, "/mul2/all")
    @test !has_pdf(prod, :bulk_2vbb, "mul1/not_spms_coincident/is_ss_aoe")
    @test !has_pdf(prod, :missing, "mul2/all")
    # has_pdf and read_pdf agree
    for groups in (:mul1, :is_ss_aoe, :is_high_aoe), hpge_group in (:all, :bege)
        @test has_pdf(prod, :bulk_2vbb; groups, hpge_group) == try
            read_pdf(prod, :bulk_2vbb; groups, hpge_group)
            true
        catch e
            e isa ArgumentError || rethrow()
            false
        end
    end

    # several simulations, in the order given
    pdfs = read_pdf(prod, [:lar_Ar39, :bulk_2vbb]; groups = :is_ss_aoe)
    @test collect(keys(pdfs)) == [:lar_Ar39, :bulk_2vbb]
    @test pdfs[:lar_Ar39].weights ≈ [1, 0, 4] ./ 100
    @test pdfs[:bulk_2vbb].weights ≈ [1, 0, 4] ./ 10
    @test build_forward_model(pdfs; e_low = 0, e_high = 3, binning = 1) isa Tuple

    # the lists of simulations and histograms are kept: files added later are not seen,
    # also when the first file changes, until a new SimProd is made
    foreach(i -> touch(joinpath(pdf_dir, "l1000test-sim$i-tier_pdf.lh5")), 10:19)
    @test simids(prod) == [:bulk_2vbb, :lar_Ar39]
    @test !has_pdf(prod, :sim10)
    @test_throws ArgumentError pdf_file(prod, :sim10)
    lh5open(joinpath(pdf_dir, "l1000test-bulk_2vbb-tier_pdf.lh5"), "cw") do f
        f["pdf/mul3/all"] = Histogram(0.0:1.0:3.0, [1.0, 2, 3])
    end
    @test !has_pdf(prod, :bulk_2vbb; groups = :mul3)
    # the kept lists cannot be changed from outside
    push!(simids(prod), :sim10)
    @test !has_pdf(prod, :sim10)
    prod = SimProd(dir)
    @test has_pdf(prod, :sim10)
    @test has_pdf(prod, :bulk_2vbb; groups = :mul3)

    # a long list of simulations is cut in the REPL, but not elsewhere
    @test occursin("sim15\n", sprint(show, MIME"text/plain"(), prod))
    limited = sprint(show, MIME"text/plain"(), prod; context = :limit => true)
    @test occursin("12 simulations:", limited)
    @test occursin("│    sim10\n│    ⋮\n│    sim19\n", limited)
    @test !occursin("sim11", limited)
end
