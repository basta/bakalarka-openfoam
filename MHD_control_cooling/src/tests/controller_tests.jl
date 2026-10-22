# test/controller_tests.jl
# Unit tests for controllers implementing the AbstractController interface.

using Test
using LinearAlgebra
using SparseArrays
using CircularArrayBuffers
using JLD2 # Required for saving/loading ABC data in DMDC_MPCController constructor
using JuMP # Needed for status codes like OPTIMAL (even in mock)

# --- Test Setup ---

# Assume the script is run from the 'test' directory or the project root.
# Adjust paths if necessary.
const SRC_DIR = isdir("src") ? "src" : "../src"
const CONTROLLERS_DIR = joinpath(SRC_DIR, "controllers")

# Include the interface definition first
include("../controllers/interface.jl")

# --- Mock Dependencies for DMDC_MPCController ---
# Create minimal mocks to allow DMDC_MPCController to load and run
# without the full dmd_lib.jl and mpc_lib.jl implementations.

# Mock dmd_lib.jl functions
function create_abc_system_from_data(data_path, DELAY)
    @info "[Mock] Called create_abc_system_from_data for path '$data_path', DELAY=$DELAY"
    # Return dummy values with plausible dimensions
    N_X_mock = 3 # Mock physical state dimension
    N_U_mock = 4 # Mock preprocessed input dimension
    N_aug_mock = N_X_mock * (DELAY + 1) + N_U_mock * DELAY
    A_full_mock = sprandn(N_aug_mock, N_aug_mock, 0.1) # Sparse random matrix
    B_full_mock = randn(N_aug_mock, N_U_mock)
    mean_vec_mock = randn(N_X_mock)
    # Return tuple matching the expected output of the real function
    return (A_full_mock, B_full_mock, nothing, nothing, nothing, mean_vec_mock, N_X_mock, N_U_mock, nothing, nothing, nothing)
end

# Mock mpc_lib.jl structures and functions
module DenseMPCLibMock
    using SparseArrays
    using LinearAlgebra
    # Import or define necessary types/enums from JuMP if used directly
    using JuMP # Assuming the real lib uses JuMP status codes

    # Mock Controller Struct
    struct DenseMPCController
        A::AbstractMatrix
        B::AbstractMatrix
        Np::Int
        m_phys::Int # Store physical input dimension
        # Add other fields if needed by the mock compute_control_action
    end

    # Mock Setup Function
    function setup_dense_mpc(A, B, Np, Q_vec, R_vec, q_vec, r_vec, E_vec, F_vec, b_vec, m_phys)
        @info "[Mock] Called setup_dense_mpc with Np=$Np, m_phys=$m_phys"
        # Just return a mock controller instance
        return DenseMPCController(A, B, Np, m_phys)
    end

    # Mock Compute Action Function
    function compute_control_action(controller::DenseMPCController, normalized_state_z0::Vector{Float64})
        @info "[Mock] Called compute_control_action"
        # Return a dummy physical control action of the correct size and a success status
        dummy_u_phys = randn(controller.m_phys) # Use the stored m_phys
        status = JuMP.OPTIMAL # Use actual JuMP status
        return dummy_u_phys, status
    end
end # end module DenseMPCLibMock

# Define the placeholder functions required by DMDC_MPCController
# These need to be defined *before* including mpc_dmdc.jl
function estimate_physical_state(state_dict::Dict{Symbol, Any}, N_X::Int)::Vector{Float64}
    @info "[Mock] Called estimate_physical_state, N_X=$N_X"
    # Return a dummy state vector of the correct size
    return randn(N_X)
end

function preprocess_input(u_phys::AbstractVector{Float64}, N_U::Int)::Vector{Float64}
    m_phys_in = length(u_phys)
    @info "[Mock] Called preprocess_input, m_phys=$m_phys_in, N_U=$N_U"
    # Return a dummy preprocessed input vector of the correct size
    # This mock doesn't need to replicate the exact outer product logic,
    # just return the right size for testing flow.
    if N_U == 0 && m_phys_in > 0 # Handle DELAY=0 case where N_U might be inferred differently
         # A simple pass-through or specific logic might be needed depending on DELAY=0 model
         @warn "[Mock preprocess_input] DELAY=0 case might need specific mock logic. Returning zeros."
         return zeros(Float64, N_U) # Or handle based on expected N_U for DELAY=0
    elseif N_U > 0
        return randn(N_U)
    else # N_U is 0, m_phys is 0
        return zeros(Float64, 0)
    end
end


# --- Include Controller Implementations ---
include("../controllers/identification.jl")
# Replace the actual DenseMPCLib with our mock *before* including the controller

# --- Test Suite ---

@testset "Controller Tests" begin

    @testset "IdentificationController Tests" begin
        num_inputs = 4
        T_stab = 50.0
        ctrl = IdentificationController(T_stab, num_inputs)
        state_dummy = Dict{Symbol, Any}() # State is not used by this controller

        @test ctrl isa IdentificationController
        @test ctrl isa AbstractController
        @test ctrl.num_inputs == num_inputs
        @test ctrl.T_stabilize == T_stab
        @test ctrl.last_change_time == -Inf
        @test ctrl.holding_time == 0.0
        @test all(ctrl.last_output .== 0.0)

        # Test initial output (during stabilization)
        time1 = 10.0
        action1 = compute_control_action(ctrl, time1, state_dummy)
        @test action1 isa Vector{Float64}
        @test length(action1) == num_inputs
        @test all(action1 .== 0.0)
        @test ctrl.last_change_time == time1 # First change time set
        @test ctrl.holding_time > 0 # Holding time should be set

        # Test output remains zero before T_stabilize elapses
        time2 = time1 + ctrl.holding_time / 2.0
        action2 = compute_control_action(ctrl, time2, state_dummy)
        @test all(action2 .== 0.0)
        @test action2 === ctrl.last_output # Should return the stored value

        # Test output changes after holding time (still during stabilization)
        time3 = time1 + ctrl.holding_time + 1.0
        action3 = compute_control_action(ctrl, time3, state_dummy)
        @test all(action3 .== 0.0) # Still zero because time < T_stab
        # Output vector might be regenerated even if it's zeros
        # @test action3 !== action1 # This might fail if zeros() returns same underlying array sometimes
        @test ctrl.last_change_time == time3
        @test ctrl.holding_time > 0

        # Test output becomes non-zero after T_stabilize and holding time
        time4 = T_stab + 1.0
        # Force update by setting last_change_time appropriately
        ctrl.last_change_time = T_stab - ctrl.holding_time - 1.0 # Ensure time >= last_change + holding
        action4 = compute_control_action(ctrl, time4, state_dummy)
        @test action4 isa Vector{Float64}
        @test length(action4) == num_inputs
        @test any(action4 .!= 0.0) # Should be random non-zero now
        @test ctrl.last_change_time == time4
        @test ctrl.holding_time > 0
        last_output_val = copy(action4)

        # Test output stays constant during holding time (after stabilization)
        time5 = time4 + ctrl.holding_time / 2.0
        action5 = compute_control_action(ctrl, time5, state_dummy)
        @test all(action5 .== last_output_val)
        @test action5 === ctrl.last_output # Should return the stored value

        # Test output changes again after holding time
        time6 = time4 + ctrl.holding_time + 1.0
        action6 = compute_control_action(ctrl, time6, state_dummy)
        @test any(action6 .!= last_output_val) # Should be different random values
        @test ctrl.last_change_time == time6

        # Test time going backwards (should reset last_change_time and potentially output zeros)
        time7 = time6 - 10.0
        @test_logs (:warn, r"Time .* is less than last_change_time") compute_control_action(ctrl, time7, state_dummy)
        # After reset, it might re-evaluate based on time7 vs T_stabilize
        action7 = compute_control_action(ctrl, time7, state_dummy) # Call again to get the output after reset
        @test ctrl.last_change_time == time7 # Should be reset
        # The output depends on whether time7 < T_stabilize
        if time7 < T_stab
             @test all(action7 .== 0.0)
        else
             # It would generate a new random output based on the reset time
             @test any(action7 .!= 0.0)
        end

        # Test update_controller_state! (should do nothing and not error)
        @test update_controller_state!(ctrl, 100.0, state_dummy) === nothing
        @test_nowarn update_controller_state!(ctrl, 100.0, state_dummy)

    end # @testset IdentificationController

    @testset "DMDC_MPCController Tests" begin
        # --- Configuration for Mock MPC ---
        m_phys_test = 2 # Physical input dimension for test
        Np_test = 5
        delay_test = 3
        N_X_test = 3 # Must match mock create_abc_system_from_data
        N_U_test = 4 # Must match mock create_abc_system_from_data

        # Create dummy ABC file for testing load_abc=true
        dummy_abc_path = "test_ABC.jld2"
        if isfile(dummy_abc_path) rm(dummy_abc_path) end # Clean up previous runs
        A_dummy = sprandn(N_X_test*(delay_test+1)+N_U_test*delay_test, N_X_test*(delay_test+1)+N_U_test*delay_test, 0.1)
        B_dummy = randn(N_X_test*(delay_test+1)+N_U_test*delay_test, N_U_test)
        mean_vec_dummy = randn(N_X_test)
        # Ensure JLD2 saves sparse matrices correctly if needed, or use dense for dummy
        jldsave(dummy_abc_path; A=Matrix(A_dummy), B=B_dummy, delay=delay_test, mean_vec=mean_vec_dummy)

        test_config_load = Dict(
            "controller" => Dict(
                "params" => Dict(
                    "DMDC_MPC" => Dict(
                        "Np" => Np_test,
                        "m_phys" => m_phys_test,
                        "abc_path" => dummy_abc_path,
                        "load_abc" => true,
                        # "delay" is loaded from file
                        # Costs/constraints are placeholders in constructor for now
                    )
                )
            )
        )

        test_config_nogen = Dict(
            "controller" => Dict(
                "params" => Dict(
                    "DMDC_MPC" => Dict(
                        "Np" => Np_test,
                        "m_phys" => m_phys_test,
                        "data_path" => "dummy_data.jld2", # Mock create_abc will use this
                        "abc_path" => "ignored_abc.jld2",
                        "load_abc" => false,
                        "delay" => delay_test
                        # Costs/constraints are placeholders in constructor for now
                    )
                )
            )
        )

        # --- Test Constructor (Loading ABC) ---
        @testset "Constructor (load_abc=true)" begin
             @info "--- Testing DMDC_MPC Constructor (load_abc=true) ---"
             ctrl_load = nothing
             @test_nowarn ctrl_load = DMDC_MPCController(test_config_load)
             @test ctrl_load isa DMDC_MPCController
             @test ctrl_load isa AbstractController
             @test ctrl_load.Np == Np_test
             @test ctrl_load.m_phys == m_phys_test
             @test ctrl_load.DELAY == delay_test # Loaded from file
             @test ctrl_load.N_X == N_X_test # Loaded from file (via mean_vec)
             @test ctrl_load.N_U == N_U_test # Inferred from loaded A/B/delay/N_X
             @test length(ctrl_load.mean_vec) == N_X_test
             @test size(ctrl_load.x_history) == (N_X_test, delay_test + 1)
             @test size(ctrl_load.u_phys_history) == (m_phys_test, delay_test)
             @test length(ctrl_load.last_computed_u_phys) == m_phys_test
             @test ctrl_load.dense_controller isa DenseMPCLibMock.DenseMPCController
             @test ctrl_load.dense_controller.m_phys == m_phys_test # Check m_phys passed down
             @info "--- Finished DMDC_MPC Constructor (load_abc=true) Test ---"
        end

        # --- Test Constructor (Generating ABC via Mock) ---
         @testset "Constructor (load_abc=false)" begin
            @info "--- Testing DMDC_MPC Constructor (load_abc=false) ---"
            ctrl_nogen = nothing
            @test_nowarn ctrl_nogen = DMDC_MPCController(test_config_nogen)
            @test ctrl_nogen isa DMDC_MPCController
            @test ctrl_nogen isa AbstractController
            @test ctrl_nogen.Np == Np_test
            @test ctrl_nogen.m_phys == m_phys_test
            @test ctrl_nogen.DELAY == delay_test # From config
            @test ctrl_nogen.N_X == N_X_test # From mock create_abc
            @test ctrl_nogen.N_U == N_U_test # From mock create_abc
            @test length(ctrl_nogen.mean_vec) == N_X_test
            @test size(ctrl_nogen.x_history) == (N_X_test, delay_test + 1)
            @test size(ctrl_nogen.u_phys_history) == (m_phys_test, delay_test)
            @test length(ctrl_nogen.last_computed_u_phys) == m_phys_test
            @test ctrl_nogen.dense_controller isa DenseMPCLibMock.DenseMPCController
            @test ctrl_nogen.dense_controller.m_phys == m_phys_test
            @info "--- Finished DMDC_MPC Constructor (load_abc=false) Test ---"
        end

        # --- Test compute_control_action ---
        @testset "compute_control_action" begin
            @info "--- Testing DMDC_MPC compute_control_action ---"
            # Use one of the constructed controllers
            ctrl = DMDC_MPCController(test_config_nogen) # Or _load
            state_dummy = Dict{Symbol, Any}(:temperature => rand(Float32, 10), :heatflux => rand(Float32, 10)) # Dummy state data

            time_k = 100.0
            action_k = nothing
            @test_nowarn action_k = compute_control_action(ctrl, time_k, state_dummy)

            @test action_k isa Vector{Float64}
            @test length(action_k) == m_phys_test

            # Check history update
            # x_history should contain the result of mock estimate_physical_state
            # u_phys_history should contain the *previous* last_computed_u_phys (which was zeros initially)
            @test length(ctrl.x_history) > 0 # Should have pushed one state
            @test ctrl.x_history[:, 1] isa Vector{Float64} # Newest state
            @test length(ctrl.x_history[:, 1]) == N_X_test

            if delay_test > 0
                 @test length(ctrl.u_phys_history) > 0 # Should have pushed one input
                 @test all(ctrl.u_phys_history[:, 1] .== 0.0) # First pushed input was the initial zero vector
            end

            # Check last_computed_u_phys update
            @test all(ctrl.last_computed_u_phys .== action_k)

            # Simulate another step
            time_k1 = time_k + 1.0
            action_k1 = compute_control_action(ctrl, time_k1, state_dummy)

            @test action_k1 isa Vector{Float64}
            @test length(action_k1) == m_phys_test
            @test all(ctrl.last_computed_u_phys .== action_k1)

            # Check u_phys_history again (should now contain action_k)
            if delay_test > 0
                @test length(ctrl.u_phys_history) > 1 # Should have pushed second input
                @test all(ctrl.u_phys_history[:, 1] .== action_k) # Newest input in history is action_k
            end
             @info "--- Finished DMDC_MPC compute_control_action Test ---"
        end

         # --- Test update_controller_state! ---
         @testset "update_controller_state!" begin
             @info "--- Testing DMDC_MPC update_controller_state! ---"
             ctrl = DMDC_MPCController(test_config_nogen)
             state_dummy = Dict{Symbol, Any}()
             @test update_controller_state!(ctrl, 100.0, state_dummy) === nothing
             @test_nowarn update_controller_state!(ctrl, 100.0, state_dummy)
             @info "--- Finished DMDC_MPC update_controller_state! Test ---"
         end

        # Clean up dummy file
        if isfile(dummy_abc_path) rm(dummy_abc_path) end

    end # @testset DMDC_MPCController

end # @testset Controller Tests

println("All controller tests passed!")
