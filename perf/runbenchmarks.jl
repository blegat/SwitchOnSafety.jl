using SwitchedSystems
using BenchmarkTools

function create_benchmark_suite()
    suite = BenchmarkGroup()
    mechanism = create_floating_atlas()
    remove_fixed_joints!(mechanism)

    let
        state = MechanismState(ScalarType, mechanism)
        result = DynamicsResult(ScalarType, mechanism)
        suite["mass_matrix"] = @benchmarkable mass_matrix!($(result.massMatrix), $state) setup = rand!($state)
    end

    let
        state = MechanismState(ScalarType, mechanism)
        result = DynamicsResult(ScalarType, mechanism)
        torques = Vector{ScalarType}(num_velocities(mechanism))
        suite["inverse_dynamics"] = @benchmarkable(
            inverse_dynamics!($torques, $(result.jointWrenches), $(result.accelerations), $state, v̇, externalWrenches),
            setup = (
                v̇ = rand(num_velocities($mechanism));
                externalWrenches = Dict(body => rand(Wrench{ScalarType}, root_frame($mechanism)) for body in non_root_bodies($mechanism));
                rand!($state)
            )
        )
    end

    let
        state = MechanismState(ScalarType, mechanism)
        result = DynamicsResult(ScalarType, mechanism)
        suite["dynamics"] = @benchmarkable(dynamics!($result, $state, externalWrenches),
            setup=(
                rand!($state);
                externalWrenches = Dict(body => rand(Wrench{ScalarType}, root_frame($mechanism)) for body in non_root_bodies($mechanism))
            )
        )
    end

    suite
end

function runbenchmarks()
    suite = create_benchmark_suite()
    tune!(suite)
    Profile.clear_malloc_data()
    results = run(suite, verbose = true)
    showall(results)
    println()
end

runbenchmarks()
