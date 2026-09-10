# Solver regression tests

Enable `ONDSELSOLVER_BUILD_TESTS` in the normal solver build. There is one
backend: the OndselSolver library being built. No separate FreeCADMbD checkout,
CAD application, model download or backend-selection flag is needed by the
dynamics tests.

```sh
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DONDSELSOLVER_BUILD_TESTS=ON
cmake --build build-release --config Release --target test_run test_dynamics
ctest --test-dir build-release -C Release --output-on-failure
```

For a multi-configuration generator, set `CMAKE_CONFIGURATION_TYPES=Release`
when configuring if only Release artifacts are wanted. An installed Google Test
can be selected with `GTest_DIR`; otherwise the existing FetchContent fallback
is used. Windows DLL-based Google Test packages are copied beside the test
executables with CMake 3.21 or newer.

## Independent dynamics benchmarks

`TestDynamics.cpp` creates its models in memory. Run just these tests with:

```sh
ctest --test-dir build-release -C Release -L dynamics --output-on-failure
```

The tests check full sampled histories, finite results, motion, forces and energy,
not merely whether the solver returns without throwing:

- Free fall: exact trajectory, velocity, acceleration and energy.
- Fixed bodies retain their input translation and orientation, including
  half-turns and offset centres of mass, while another body falls independently.
- Physical pendulum: independent RK4 reference, tighter-tolerance convergence,
  energy drift and step-refined hinge-force balance. Output spacing is distinct
  from the maximum integration step.
- Constant-torque rotor: angular acceleration and work/energy balance.
- Spring/mass/damper: analytical underdamped response and decreasing energy.
- Prescribed motion with an independently falling body: mixed driven/free dynamics.
- Prescribed rotor and slider/crank: inverse-dynamics reaction checks.
- Metre versus millimetre inputs: consistent force, inertia and torque scaling.
- Load Jacobians: independent central differences at offset markers on two moving
  bodies, including nonzero angular and translational velocities.
- Typed-load serialization, invalid inputs, coincident spring endpoints, output
  at an off-grid end time, nonzero start time and replacement of results on rerun.
- Rigid translational/rotational limits, compliant stop rebound and damping,
  compliant distance-contact rebound, distance-limit serialization, and rigid
  cursor-driven contact pushing.
- Fifth-order differentiation at small time steps and rejection of singular time
  nodes, guarding against reuse of a stale integration matrix.

Models use kg, m and s unless explicitly scaled to kg, mm and s. Tolerances are
stated beside each assertion. These tests do not qualify every joint type,
long-duration stability, arbitrary-shape contact, friction or flexible bodies.

## Legacy tests

`test.cpp` retains the existing file-based kinematic and dragging smoke tests.
Many only check that execution completes, and the legacy kinematic wrapper can
swallow errors; they are not substitutes for the numerical assertions above.
Run them serially because several write the same output filenames.

The three `runPreDragBackhoe` cases cover assembly-level marker parsing as well
as dragging. They must remain enabled in Release builds, where parser bounds
errors are not hidden by debug assertions.
All three pass in the current suite. They are deliberately kept enabled so
changes to drag weighting or parser bounds cannot silently reintroduce the
Release-only crashes seen while bringing up the dynamics port.
