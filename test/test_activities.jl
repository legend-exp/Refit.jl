using Distributions
using IntervalSets
using LegendHDF5IO
using Refit
using StatsBase
using Test
using Unitful

@testset "activity_prior" begin
    d = activity_prior(; val = 83u"μBq", unc = 113u"μBq")
    @test d == truncated(Normal(83, 113); lower = 0)
    @test minimum(d) == 0

    # the quantile cl of the half-normal is the limit
    d = activity_prior(; limit = 12u"μBq", cl = 0.9)
    @test quantile(d, 0.9) ≈ 12
    @test median(d) < 12 / 2
    @test minimum(d) == 0

    # the distribution of the activity in the unit asked for, from any unit of activity
    @test activity_prior(; val = 2u"mBq", unc = 1u"mBq") ==
          truncated(Normal(2000, 1000); lower = 0)
    @test activity_prior(; val = 2u"mBq", unc = 1u"mBq", unit = u"Bq") ≈
          truncated(Normal(0.002, 0.001); lower = 0)
    @test quantile(activity_prior(; limit = 1u"Bq", cl = 0.9, unit = u"mBq"), 0.9) ≈ 1000

    # activities with units
    @test_throws TypeError activity_prior(; val = 83, unc = 113)
    @test_throws Unitful.DimensionError activity_prior(; val = 1u"kg", unc = 1u"kg")
    @test_throws ArgumentError activity_prior(; val = 1u"μBq")
    @test_throws ArgumentError activity_prior(;
        val = 1u"μBq",
        unc = 1u"μBq",
        limit = 2u"μBq",
        cl = 0.9,
    )
    @test_throws ArgumentError activity_prior()
    @test_throws ArgumentError activity_prior(; limit = 2u"μBq", cl = 1.5)
end

@testset "read_activities" begin
    # a part made of two materials, stood for by the same simulation, a limit and
    # quantities of each kind
    file = tempname() * ".yaml"
    write(
        file,
        """
        clamp_ultem:
          description: clamp, Ultem part
          mass_in_kg: 0.5
          specific_activity_in_uBq_per_kg:
            clamps_Pb214: { val: 80, unc: 100 }
        clamp_bronze:
          mass_in_kg: 2
          specific_activity_in_uBq_per_kg:
            clamps_Pb214: { val: 10, unc: 1 }
            clamps_Pb212: { limit: 4, cl: 0.9 }
        reflector:
          area_in_m2: 3
          specific_activity_in_uBq_per_m2:
            reflector_Pb214: { val: 2, unc: 1 }
        sipms:
          number_of_pieces: 100
          activity_per_piece_in_uBq:
            sipms_Pb214: { limit: 0.5, cl: 0.95 }
        """,
    )
    activities = read_activities(file)
    @test Set(keys(activities)) == Set([:clamp_ultem, :clamp_bronze, :reflector, :sipms])
    @test Set(keys(activities[:clamp_bronze])) == Set([:clamps_Pb214, :clamps_Pb212])

    # in another unit
    @test read_activities(file; unit = u"mBq")[:clamp_bronze][:clamps_Pb214] ≈
          truncated(Normal(0.02, 0.002); lower = 0)

    # the quantity times the activity per unit, in µBq
    @test activities[:clamp_ultem][:clamps_Pb214] == truncated(Normal(40, 50); lower = 0)
    @test activities[:clamp_bronze][:clamps_Pb214] == truncated(Normal(20, 2); lower = 0)
    @test activities[:reflector][:reflector_Pb214] == truncated(Normal(6, 3); lower = 0)
    @test quantile(activities[:clamp_bronze][:clamps_Pb212], 0.9) ≈ 8
    @test quantile(activities[:sipms][:sipms_Pb214], 0.95) ≈ 50

    # one quantity with the matching activity
    for bad in (
        "a:\n  mass_in_kg: 1\n  specific_activity_in_uBq_per_m2:\n    s: { val: 1, unc: 1 }\n",
        "a:\n  mass_in_kg: 1\n  area_in_m2: 1\n  specific_activity_in_uBq_per_kg:\n    s: { val: 1, unc: 1 }\n",
    )
        write(file, bad)
        @test_throws ArgumentError read_activities(file)
    end
end

@testset "activity_components and background_index" begin
    # flat spectra: 1e-3 counts per keV and per simulated primary
    flat = Histogram(0.0:1.0:4000.0, fill(1e-3, 4000))
    double = Histogram(0.0:1.0:4000.0, fill(2e-3, 4000))
    pdfs = Dict(:clamps_Pb214 => flat, :reflector_Pb214 => double, :unused => flat)
    activities = Dict(
        :clamp_ultem =>
            Dict(:clamps_Pb214 => activity_prior(; val = 40u"μBq", unc = 50u"μBq")),
        :clamp_bronze =>
            Dict(:clamps_Pb214 => activity_prior(; val = 20u"μBq", unc = 2u"μBq")),
        :reflector =>
            Dict(:reflector_Pb214 => activity_prior(; val = 6u"μBq", unc = 3u"μBq")),
    )

    # one source for each component and simid, the same spectrum for both parts
    components, priors = activity_components(pdfs, activities)
    names =
        [:clamp_ultem_clamps_Pb214, :clamp_bronze_clamps_Pb214, :reflector_reflector_Pb214]
    @test Set(keys(components)) == Set(names) == Set(keys(priors))
    @test components[:clamp_ultem_clamps_Pb214] === pdfs[:clamps_Pb214]
    @test components[:clamp_bronze_clamps_Pb214] === pdfs[:clamps_Pb214]
    @test priors[:clamp_ultem_clamps_Pb214] == activities[:clamp_ultem][:clamps_Pb214]
    @test_throws ArgumentError activity_components(Dict(:clamps_Pb214 => flat), activities)

    # one sample: the activities in µBq. A flat spectrum gives the same index in any
    # window: 1e-6 * a * 1e-3 counts per second and per keV, for the mass
    livetime, mass = 1u"yr", 1000u"kg"
    model, bin_edges = build_forward_model(
        components;
        norm = ustrip(NoUnits, u"μBq" * livetime),
        e_low = 1000,
        e_high = 3000,
    )
    params = (
        clamp_ultem_clamps_Pb214 = 40.0,
        clamp_bronze_clamps_Pb214 = 20.0,
        reflector_reflector_Pb214 = 6.0,
    )
    samples = (v = [params], weight = [1])
    expected(a) = uconvert(u"keV^-1*kg^-1*yr^-1", a * 1e-9u"s^-1" / (mass * u"keV"))
    bi = background_index(samples, model, bin_edges; livetime, mass, n = 10)
    @test collect(keys(bi)) == [:total; sort(names)]
    @test all(length.(values(bi)) .== 10)
    @test unit(bi[:total][1]) == u"keV^-1*kg^-1*yr^-1"
    @test all(bi[:clamp_ultem_clamps_Pb214] .≈ expected(40))
    @test all(bi[:reflector_reflector_Pb214] .≈ expected(2 * 6))
    @test all(bi[:total] .≈ expected(40 + 20 + 12))
    # half the mass, twice the index; the same exposure in other units; a window of 10 keV
    @test background_index(samples, model, bin_edges; livetime, mass = 500u"kg", n = 1)[:total][1] ≈
          2 * expected(72)
    @test background_index(
        samples,
        model,
        bin_edges;
        livetime = 365.25u"d",
        mass = 1e6u"g",
        n = 1,
    )[:total][1] ≈ expected(72)
    @test background_index(
        samples,
        model,
        bin_edges;
        livetime,
        mass,
        window = 2000 .. 2010,
        n = 1,
    )[:total][1] ≈ expected(72)
    # the live time and the mass need units
    @test_throws TypeError background_index(
        samples,
        model,
        bin_edges;
        livetime = 3.15e7,
        mass,
        n = 1,
    )
    @test_throws TypeError background_index(
        samples,
        model,
        bin_edges;
        livetime,
        mass = 1000,
        n = 1,
    )

    # the counts in energy ranges, which are added: 1e-9 * a counts per second and per keV
    per_kev(a) = a * 1e-9 * ustrip(u"s", livetime)
    counts = get_predicted_counts(samples, model, bin_edges, 2000 .. 2010; n = 1)
    @test counts[:total][1] ≈ 10 * per_kev(72)
    @test counts[:clamp_bronze_clamps_Pb214][1] ≈ 10 * per_kev(20)
    two =
        get_predicted_counts(samples, model, bin_edges, [2000 .. 2010, 2100 .. 2130]; n = 1)
    @test two[:total][1] ≈ 40 * per_kev(72)
    @test_throws ArgumentError get_predicted_counts(
        samples,
        model,
        bin_edges,
        2000.5 .. 2010;
        n = 1,
    )
    @test_throws ArgumentError get_predicted_counts(
        samples,
        model,
        bin_edges,
        2000 .. 2010;
        n = 0,
    )
    # the bins must cover the ranges, also without regions left out of the model
    blinded, blinded_edges = build_forward_model(
        components;
        norm = ustrip(NoUnits, u"μBq" * livetime),
        binning = [1000:1:2050, 2060:1:3000],
        e_low = 1000,
        e_high = 3000,
    )
    @test_throws ArgumentError get_predicted_counts(
        samples,
        blinded,
        blinded_edges,
        2040 .. 2070;
        n = 1,
    )
    @test get_predicted_counts(samples, blinded, blinded_edges, 2060 .. 2070; n = 1)[:total][1] ≈
          10 * per_kev(72)

    # the values of each sample are summed in groups
    grouped = group_histograms(bi, ["clamps" => r"clamps", "rest" => :rest])
    @test collect(keys(grouped)) == [:total, :clamps, :rest]
    @test grouped[:clamps] ≈ bi[:clamp_ultem_clamps_Pb214] .+ bi[:clamp_bronze_clamps_Pb214]
    @test grouped[:rest] == bi[:reflector_reflector_Pb214]

    # the samples are drawn with their weights
    inactive = map(_ -> 0.0, params)
    weighted = (v = [inactive, params], weight = [0, 3])
    @test all(
        background_index(weighted, model, bin_edges; livetime, mass, n = 100)[:total] .≈
        expected(72),
    )
end

@testset "activities of a production cycle" begin
    # a fake cycle with the activities in its configuration directory
    dir = mktempdir()
    write(
        joinpath(dir, "simflow-config.yaml"),
        """
        experiment: l1000test
        paths:
          config: \$_/inputs/simprod/config
          tier:
            pdf: \$_/generated/tier/pdf
        """,
    )
    pdf_dir = mkpath(joinpath(dir, "generated", "tier", "pdf"))
    for (simid, n_prim) in ("clamps_Pb214" => 10, "reflector_Pb214" => 100)
        lh5open(joinpath(pdf_dir, "l1000test-$simid-tier_pdf.lh5"), "w") do f
            for cut in ("is_ss_aoe", "is_low_aoe")
                f["pdf/mul1/$cut/all"] = Histogram(0.0:1.0:4.0, [1.0, 2.0, 3.0, 4.0])
            end
            f["nr_sim_events"] = n_prim
        end
    end
    config_dir =
        mkpath(joinpath(dir, "inputs", "simprod", "config", "bkgmodel", "l1000test"))
    write(
        joinpath(config_dir, "components.yaml"),
        """
        clamp:
          mass_in_kg: 2
          specific_activity_in_uBq_per_kg:
            clamps_Pb214: { val: 10, unc: 1 }
        reflector:
          area_in_m2: 3
          specific_activity_in_uBq_per_m2:
            reflector_Pb214: { val: 2, unc: 1 }
        """,
    )

    prod = SimProd(dir)
    activities = read_activities(prod)
    @test activities[:clamp][:clamps_Pb214] == truncated(Normal(20, 2); lower = 0)

    components, priors = activity_components(prod, activities; groups = :is_ss_aoe)
    @test Set(keys(components)) == Set([:clamp_clamps_Pb214, :reflector_reflector_Pb214])
    @test components[:clamp_clamps_Pb214].weights ≈ [1, 2, 3, 4] ./ 10
    @test components[:reflector_reflector_Pb214].weights ≈ [1, 2, 3, 4] ./ 100
    @test components[:clamp_clamps_Pb214].weights ≈
          activity_components(prod, activities; groups = :mul1)[1][:clamp_clamps_Pb214].weights ./
          2

    activities[:fibers] =
        Dict(:fibers_Pb214 => activity_prior(; val = 1u"μBq", unc = 1u"μBq"))
    @test_throws ArgumentError activity_components(prod, activities)
end
