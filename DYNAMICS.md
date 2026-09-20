# Forward rigid-body dynamics

Forward dynamics is implemented in the bundled OndselSolver library. FreeCADMbD
is a source of ported code, **not** an additional solver dependency or runtime
backend. The existing assembly-solving and `runKINEMATIC()` entry points remain.

Kinematic simulations retain upstream's drag-only limit behavior. Forward
dynamics includes active stops in both initial-condition and integration
equations, so impacts and separating reactions are handled during the run.

This release changes C++ class layouts and virtual interfaces. Rebuild callers
against the new headers and library; it is not a binary-compatible replacement
for the version-1 shared library. The shared-library ABI version is now 2.

## Entry point and inputs

Construct an `ASMTAssembly` using the existing parts, mass markers, attachment
markers and joint classes, then call `assembly->runDYNAMIC()`.

Moving bodies need finite positive mass and three finite positive principal
inertias about their centre of mass. Place/orient the principal mass marker at
that centre in the body's coordinate system. Supply consistent length, mass and
time units throughout; there is no automatic conversion from FreeCAD's document
units in this solver API. Angles and angular velocities use radians.

Set `constantGravity` using the existing `ASMTConstantGravity` API, including an
explicit zero vector if gravity is unwanted. `ASMTSimulationParameters` supplies
the increasing time interval, output interval, positive minimum/maximum steps,
integration and corrector tolerances, maximum iterations and BDF order (1–5).
`seterrorTol()` sets both kinematic and dynamic tolerances; the individual fields
allow finer control. `iterMaxDyn` is separate from the kinematic iteration limit.

Current placements and velocities are the initial inputs to each run. As in the
existing solver API, producing results updates the ASMT objects. Explicitly
restore the desired initial state before repeating a simulation. A dynamics run
clears previous output histories. Exceptions propagate to the caller; a failed
run can leave a partial history and must not be presented as a completed result.

## Loads

`ASMTForceTorque::With()` creates a typed load. Set its name and the full names of
markers I and J, configure it, and attach it with `assembly->addForceTorque(load)`.

```cpp
auto load = MbD::ASMTForceTorque::With();
load->setName("DriveTorque");
load->setMarkerI(bodyMarker->fullName(""));
load->setMarkerJ(groundMarker->fullName(""));
load->setTorque3D(0.0, 0.0, 0.8);
assembly->addForceTorque(load);
assembly->runDYNAMIC();
```

`setForce3D(x,y,z)` and `setTorque3D(x,y,z)` select a constant, **world-resolved**
wrench acting at marker I. The balancing wrench acts at J. For separated markers,
the opposite force is transported to J with the balancing couple, so total force
and moment are both zero. To apply an external load to a body, use that body's
marker as I and a grounded marker as J. These are not body-following vectors.

Alternatively, `setSpringDamper(k,c,l0)` selects a line spring/damper. With
`d = rJ-rI`, `l = |d|`, `u = d/l`, and `v = vJ-vI`, its force on I is:

```
F_I = [k (l-l0) + c (u dot v)] u
F_J = -F_I
```

Stiffness, damping and rest length must be finite and nonnegative. Coincident
attachment points are rejected because the force direction is undefined. Load
markers must be fixed to their bodies. The force/torque setters switch back to
constant-wrench mode. Analytical position and velocity Jacobians include marker
offsets and quaternion derivatives. No symbolic load-expression parser is added.

The ASMT serialization uses explicit `WorldWrench` and `LinearSpringDamper` load
records. It does not interpret FreeCADMbD's expression-based load records.

## Outputs and scope

Existing ASMT histories contain sampled poses, velocities, accelerations and
joint/motion reactions. Load `fxs/fys/fzs` and `txs/tys/tzs` histories are world
components acting at marker I (torque is about that marker). The output includes
input and assembled-initial-condition samples at the start, regular output times,
and the exact requested end time, even when it is off the regular output grid.
Check convergence separately for motion and reactions: a satisfactory trajectory
does not by itself establish sufficiently accurate acceleration or joint forces.

This is the solver foundation; FreeCAD's Assembly module supplies the mass/unit
adapter, persistent simulation inputs, task panels, plotting and playback.
Translation and rotation limits support both rigid unilateral constraints and
compliant stiffness/damping behavior. `ASMTDistanceLimit` supplies the same
behavior for exact sphere-sphere contact in Assembly, including rigid drag
separation and compliant dynamic impact. Runtime-evaluated force callbacks let
FreeCAD's Assembly module supply transformed BRep contact forces without adding
an OpenCASCADE dependency to OndselSolver; Assembly also uses these callbacks
for regularized tangential contact friction. General solver-native friction, flexible bodies, FEM
coupling and cancellable background integration are not implemented or qualified
for forward dynamics.
Existing kinematic use of unsupported features is unchanged.

## Port provenance

The DAE/BDF integrator, corrector and associated dynamics hooks were adapted from
`aiksiongkoh/FreeCADMbD`, revision `ef7d572` (2025-12-09), which retains the
constraint architecture used by this bundled solver. Copyright/license headers
are retained. The port preserves C++17 and does not introduce Boost or a second
solver library. Later FreeCADMbD class/API restructuring was not imported.

Local adaptations include the typed load layer, tolerance wiring, initial
momentum derivatives, the quaternion-acceleration setter fix, safe column-vector
results, scaled Taylor-matrix inversion with propagated failures, startup
interpolation, final-time sampling, input validation and rerun history handling.
See `tests/README.md` and `tests/TestDynamics.cpp` for executable qualification.
