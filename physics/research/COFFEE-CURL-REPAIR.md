# Curl-preserving field candidate

The selected numerical field for the Vela-informed visualization is the
97³, three-rotation **curl-radiative / no-GLM** experiment. It is a finite,
force-free numerical model with an enlarged conducting star. It is not a
quantitatively matched Vela reconstruction. The 65³ counterpart measures
resolution sensitivity. No snapshot is automatically released as a validated
emission calculation.

## Source and magnetic constraint

The original analytic dipole was singular at the origin. Directly blending
its magnetic components into the evolved field also injected magnetic
divergence. The replacement initializes

\[
\mathbf A=(\boldsymbol\mu\times\mathbf r)f(r),\qquad
\mathbf B=D_4\times\mathbf A.
\]

Outside the numerical stellar radius, \(f=r^{-3}\). Inside, a finite Taylor
continuation matches the first four derivatives at the surface. The same
fourth-order centered derivative is used to initialize the curl and to evolve
Faraday's law, so \(D_4\cdot(D_4\times\mathbf A)=0\) to roundoff in the interior.

The star imposes only the conducting electric relation
\(\mathbf E=-(\boldsymbol\Omega\times\mathbf r)\times\mathbf B\).
No evolved magnetic component is reset to a prescribed dipole. The transition
is a compact, smooth layer two spatial cells wide, centered on the stated
stellar radius. Its width changes physically under grid refinement; this
remains a numerical source-size assumption. The source moment and its phase
are measured from each output rather than assumed to remain exact.

Uniform GLM cleaning passed its independent Fourier-mode check. However,
with the finite-domain guard treatment it transported boundary divergence
into the usable interior. The selected calculation disables that magnetic
gradient correction. Faraday's discrete curl and constant-coefficient KO
dissipation preserve the interior magnetic constraint; the actual interior
divergence is measured in every candidate output. The GLM scalar can still
be stored by upstream code but is dynamically decoupled from the fields.

## Outer boundary

Independent magnetic-component damping was removed. A separate electric
penalty layer relaxes tangential electric field toward

\[
\mathbf E_t=-\mathbf n\times(\mathbf B-\mathbf B_{\rm initial}),
\]

while preserving normal electric field. Here \(\mathbf n\) is the normalized,
conductivity-weighted combination of outward cube-face normals. The cubic
layer's physical peak rate is held fixed as resolution changes. The full-step
electric relaxation uses exponential splitting. The initial background is
subtracted to avoid driving a static dipole or constant guide field.

For a normally outgoing plane wave, this penalty vanishes. It damps the
incoming characteristic and never directly modifies magnetic components.
This is an approximate radiation boundary, **not a perfectly matched layer**.
An ordinary electric-conductivity absorber was also tested and rejected:
it returned about 6–20% of a normally incident pulse's energy over the rates
tested. Those results and the separate source tree remain as records.

The final compiled boundary tests use the same no-GLM mode as the candidate:

* Normally incident Gaussian guide-field pulse: returned incoming energy is
  below \(10^{-26}\) of incident energy at both 96 and 192 cells. Residual
  outgoing energy inside the measuring window is \(2.49\times10^{-4}\) and
  \(2.12\times10^{-7}\), respectively
* A 45-degree Gaussian pulse has a sampled, D4-solenoidal magnetic field
  because x/y spacings are equal. After the incident pulse clears the central
  window, the finite cube returns 3.0869% of the peak incident window energy
  and its central D4 divergence stays below \(2.85\times10^{-7}\). A wider
  domain supplies a causally isolated reference. This is one pulse test, not
  a general current-sheet or steady-wind absorption certificate

## Reproduction and provenance

Build with `run_coffee_build_radiative.sh`. Run the bounded model with
`run_coffee_radiative_experiment.py --no-cleaning`; the saved configuration
defines the actual physical parameters. `run_coffee_curl_audit.py` measures
source drift and native divergence. The separate independent source audit
tests moment fitting and a D4-compatible closed-surface flux identity.

`third_party/coffee-curl-radiative-provenance.json` records the pinned
upstream commit, complete patch chain, edited source hashes, binary hashes,
build entry point and candidate configuration. The native upstream project
and its modifications retain the upstream GPL-2.0 license. Unchanged sources
and dependencies are shared by symlink with the preceding local variants;
edited source files are private to the new variant.

The candidate directories are:

* `third_party/coffee-runs/pulsar-a75-n65-3rot-curl-radiative-noglm`
* `third_party/coffee-runs/pulsar-a75-n97-3rot-curl-radiative-noglm`

Retained experiments with GLM or with magnetic/electric component damping
are historical diagnostics. They are not interchangeable with this selected
field. Every result must identify its actual snapshot and source hashes.

## Limits retained in the visualization

The numerical star radius is 0.5 light-cylinder radii, much larger than Vela's
physical stellar radius. The domain and current sheet have finite resolution.
The conducting source diffuses and lags under evolution; its measured errors
must accompany comparisons. Interpolating native fields does not preserve
the native discrete divergence identity. The radiation population and
accelerating field are prescribed model inputs. Pulse and spectral fits do
not establish a unique reconstruction, and no curvature-only model accounts
for Vela's full observed broadband emission.
