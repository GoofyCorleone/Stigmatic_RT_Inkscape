# Aplanatism

<p align="center">
  <a href="../README.md">Overview &amp; Install</a> &nbsp;·&nbsp;
  <a href="basic-elements.md">Basic optical elements</a> &nbsp;·&nbsp;
  <a href="stigmatism.md">Stigmatism</a> &nbsp;·&nbsp;
  <a href="aplanatism.md"><b>Aplanatism</b></a>
</p>

---

## Why stigmatism is not enough

A stigmatic surface images **one** point perfectly. Move the object sideways by a
millimetre and the guarantee is gone — the off-axis image grows a comet-shaped flare, which
is exactly what the aberration called **coma** looks like.

To image an *extended* object you need the neighbourhood of the axial point to be imaged
sharply too. Abbe's answer (1873) is the **sine condition**: the system must keep

```
sin u₀ / sin u_N  =  (n_N / n₀) · M
```

*constant* for every ray leaving the axial object point, where `u₀` and `u_N` are the angles
the ray makes with the axis on entering and leaving, and `M` is the lateral magnification.
A system that is both stigmatic and satisfies the sine condition is **aplanatic**: free of
spherical aberration *and* coma, simultaneously.

## The generalized condition

Chapter 4 of Silva-Lora's thesis combines the sine condition with the rigorous-stigmatism
formulation and arrives at a condition written directly in the GOTS shape parameters. For a
system of `N` Cartesian surfaces:

```
𝓜(N, ρₖ) = ∏ₖ (n_{k+1}/n_k) · [ 1/(d_{k+1} − ζₖ) − 1/V̄C̄ₖ ] / [ 1/(dₖ − ζₖ) − 1/V̄C̄ₖ ] = 1
```

for every ray. `V̄C̄ₖ` is the distance from the vertex to the point where the surface normal
crosses the axis — for a sphere it is just the radius, and for a general Cartesian oval it
depends on where the ray strikes:

```
1/V̄C̄ = [ O + 2(2T − S·O)·ρ² / (1 + rad) ] / rad ,   rad = √(1 + (2S − O·OG)·ρ²)
```

> The form published in the thesis (Eq. 92) divides by `G` and by `O`, both of which vanish
> in two of the four families below. Using the identity `T = S²/OG`, which holds identically
> for the GOTS parameters, it can be rewritten as above — finite in every degenerate limit.
> That is the form `inversa_vc` implements.

Two consequences are worth keeping in mind:

- **𝓜 is a product over surfaces.** Chaining aplanatic elements gives an aplanatic system,
  which is why the three-lens relay further down works.
- **Aplanatism does not remove the other aberrations.** Astigmatism, field curvature and
  distortion survive. This matters for reading the figures — see
  [where the image really forms](#where-the-image-really-forms).

## The four rigorously aplanatic families

Setting `𝓜 ≡ 1` exactly, rather than minimising it numerically, yields four closed-form
solutions (§4.4 of the thesis):

| Type | Condition | Surfaces | Image | Extra |
|---|---|---|---|---|
| **Tipo-0** | `2Sₖ = GₖOₖ²` | spheres at the **Young points** | virtual | also free of astigmatism (**anastigmatic**) |
| **Tipo-1** | `d₁ → ∞` | identical conics, collimated interior | real (`M = −1`, biconvex) or virtual (`M = +1`, meniscus) | the classic biconic |
| **Tipo-2** | `Oₖ = 0` | flat-vertex ovoids | virtual | generalises the plane-parallel plate; ophthalmic-style |
| **Tipo-3** | `Gₖ = 0` | ovoids | virtual | positive meniscus |

All four are implemented in `disenar_aplanatica`, and
[`tests/unit/test_aplanatico.py`](../inkscape-raytracing/tests/unit/test_aplanatico.py)
checks them against Tables 6–9 of the thesis. The conjugates come out exactly as tabulated:
`d₁ = 26.666667, d₂ = −8` (Table 6), `d₁ = ∞, d₂ = 130` (Table 7), `d₂ = 4.444444`
(Table 8), `d₁ = 41.481481, d₂ = −22.4` (Table 9).

Note that in air all four have `|M| = 1`. That is not a limitation of the implementation but
a property of the theory: types 2 and 3 require `n₂ = n₀` for `𝓜 = 1` (Eq. 112), and a
type-0 air-to-air singlet has `M = (n₀/n₂)² = 1`. Magnification other than unity needs
immersion or a non-rigorous, optimised design.

## The extension

**Extensions → Optics → Lente Aplanática**, four tabs.

### Tab 1 — Diseño

![Aplanatic lens, design tab](img/ui-aplanatica-1.svg)

| Field | Meaning |
|---|---|
| `Modo` | *Tipo* uses one of the four closed-form families; *Superficies libres* lets you enter conjugates yourself (next tab) |
| `Tipo` | which family (0–3) |
| `Rama` | type-1 only: `real` gives the biconvex, real-image branch; `virtual` the meniscus |
| `n₀ / n₁ / n₂` | object space, glass, image space. Types 1–3 force `n₂ = n₀` and warn if you set otherwise |
| `ζ₀ / ζ₁` | the two vertices |
| `d₀` | object position |

The important difference from the stigmatic generators: **you do not enter the image
position**. `d₂` is *derived* from the aplanatism condition and drawn on the canvas — a
filled circle if the image is real, hollow if virtual.

### Tab 2 — Superficies libres (GOTS)

![Aplanatic lens, free surfaces tab](img/ui-aplanatica-2.svg)

The general mode. You give all three conjugates `d₀, d₁, d₂` and the extension derives the
GOTS parameters of both surfaces with `calcular_gots`. The singlet is rigorously stigmatic
for the axial pair but in general **not** aplanatic — the Análisis tab tells you how far off
it is, via the `(M − 1)_RMS` metric of Eq. 94.

This is the mode for exploring your own designs, and for reproducing the σ sweep of Fig. 16
of the thesis, where the best-form singlet is found by minimising that metric.

### Tab 3 — Apertura y Haz

![Aplanatic lens, aperture tab](img/ui-aplanatica-3.svg)

Aperture (`0` = automatic, clipped where the two surfaces would cross), the **physical
pupil** and the ray count. The pupil matters more than it looks: limiting the beam with a
real stop rather than a fixed cone angle is what makes different field points comparable,
because they all cross the same opening. Every figure in the thesis is built that way.

### Tab 4 — Análisis del aplanatismo

![Aplanatic lens, analysis tab](img/ui-aplanatica-4.svg)

This is where the extension answers the question "did I actually get an aplanatic lens, and
where does it form the image?".

| Field | Meaning |
|---|---|
| `Número de fuentes puntuales` | how many field points to trace |
| `Altura máxima del campo` | the largest object height |
| `Emitir los haces…` | write the seed rays so you can trace in Inkscape |
| `Dibujar superficie imagen y focos reales` | trace each fan, mark its **real** meridional focus with a hollow circle, join them with a purple curve, and draw the paraxial plane in grey for contrast |
| `Anotar la verificación en el dibujo` | write the Abbe sine ratio, `(M − 1)_RMS`, the spot size and the departure from the paraxial plane onto the canvas, with an explicit `APLANÁTICO en todo rigor` / `NO aplanático` verdict |

## Where the image really forms

This is the part that is easy to get wrong, so the figures are built to make it obvious.

An aplanatic lens images off-axis points **sharply** — that is the whole point. But those
sharp images do **not** lie on the paraxial image plane: the locus of the real foci is a
**curved surface**. The thesis says so explicitly in Fig. 17 ("cuya imagen se forma sobre
una superficie curva") and discusses it on p. 101, where it also separates the meridional
from the sagittal focal surface.

So each figure below traces **four point sources** at different field heights through a
physical pupil, finds the real meridional focus of every fan, and draws:

- **hollow circles** — the real foci,
- a **purple dashed curve** through them — the image surface,
- a **grey line** — the paraxial plane, for reference only,
- an **inset** magnifying the image region along `z` (×3 to ×122 depending on the type),
  because at true scale the departure is only tenths of a millimetre.

| | |
|---|---|
| ![Tipo-0](../ejemplo_aplanatica_tipo0.svg) | ![Tipo-1](../ejemplo_aplanatica_tipo1.svg) |
| ![Tipo-2](../ejemplo_aplanatica_tipo2.svg) | ![Tipo-3](../ejemplo_aplanatica_tipo3.svg) |

Reading the numbers on the figures:

| Type | Meridional spot at the real foci | Departure from the paraxial plane |
|---|---|---|
| Tipo-0 | 2.7 µm | 0.21 mm |
| Tipo-1 | 23 µm | 3.78 mm |
| Tipo-2 | 1.3 µm | 0.11 mm |
| Tipo-3 | 15 µm | 2.32 mm |

Sharp foci, curved surface. For the anastigmatic tipo-0 the residual matches the ≈60 µm the
thesis reports at full field (Fig. 22).

## A three-lens aplanatic system

Because `𝓜` is a product over surfaces, aplanatic elements chain. Three type-1 lenses form
a 1:1 relay of six Cartesian surfaces — the real image of each is the object of the next:

![Three-lens relay](../ejemplo_aplanatica_3lentes.svg)

| Quantity | Value |
|---|---|
| Total magnification | `−1.0000` |
| `(M − 1)_RMS` over the six surfaces | **4.9·10⁻¹⁶** |
| Abbe sine ratio of the whole train | `−1.000000`, constant to `4·10⁻¹⁵` |
| Conjugates | z = 0 → 130 → 260 → 390 mm |

[`ejemplo_aplanatica_3lentes_inkscape.svg`](../ejemplo_aplanatica_3lentes_inkscape.svg)
carries the same geometry with the rays as `optics:beam` seeds. Open it, press `Ctrl+A` and
run **Ray Tracing**: all 36 beams land on the analytically predicted foci to within
**0.03–3.8 µm**.

## Tutorial

**Design a lens.**

1. **Extensions → Optics → Lente Aplanática**.
2. *Diseño*: leave `Modo = Tipo`, pick `Tipo-1`, `Rama = biconvexa`. Keep the defaults
   `n₁ = 1.8`, `ζ₀ = 60`, `ζ₁ = 70`, `d₀ = 0` — they reproduce Table 7 of the thesis.
3. *Apertura y Haz*: pupil radius `8` at `z = 55`, 9 rays.
4. *Análisis*: 4 field points up to `9`, both checkboxes on.
5. **Apply.**

You get the lens, the pupil, four coloured fans, the curved image surface with its real
foci, and a block of text stating the sine ratio, `(M − 1)_RMS = 0` and the verdict
`APLANÁTICO en todo rigor`.

**See a design fail the test.** Switch *Modo* to `Superficies libres` and, on tab 2, enter
`d₁ = 599.172082`, `d₂ = 200`. Apply again: the annotation now reads
`(M − 1)_RMS ≈ 1.5·10⁻²` and `NO aplanático`, and the foci scatter. That design *is*
rigorously stigmatic on axis — it just is not aplanatic.

**Reproduce the figures without Inkscape.**

```bash
python generar_ejemplos_aplanaticos.py
```

writes the four type figures plus the three-lens relay and its traceable variant.

## How it is verified

The aplanatic layer of `gots_util.py`:

| Function | What it does |
|---|---|
| `disenar_aplanatica(tipo, n0, n1, ζ0, ζ1, d0, …)` | conjugates, per-surface GOTS and lateral magnification for each of the four families |
| `inversa_vc(sup, ρ)` | `1/V̄C̄ₖ` (Eq. 92) in the degeneracy-free form given above |
| `factor_aplanatismo` / `condicion_aplanatismo` | the per-surface factor and the product `𝓜(N, ρₖ)` (Eqs. 78/93) |
| `rms_aplanatismo(valores)` | `(M − 1)_RMS` (Eq. 94) — the error metric for approximate designs |
| `perfil_superficie_aplanatica(sup, r_max, N)` | profile sampled to a *requested* aperture (types 1–3 have an unbounded `ρ` domain, so the generic sampler would put nearly every node outside the lens) |
| `trazar_rayo(superficies, P, d, r_max)` | exact ray trace through the analytic surfaces, per-surface apertures allowed |
| `foco_meridional(origenes, direcciones)` | best-crossing point of a fan and its spot radius |
| `superficie_imagen(…)` | the curved image surface: real focus and spot per field height |
| `seno_abbe(…)` | the sine ratio on the axial pair and its dispersion |
| `alturas_estigmaticas(diseno, ρ₀)` | chains the exact stigmatic geometry to get `ρ₀ → ρ₁` without a tracer |

The module traces its own designs exactly (`trazar_rayo`, the method of chapter 3 of the
thesis), independently of the Bézier approximation written into the SVG. Against that
reference:

| Check | Result |
|---|---|
| Axial focus of all four types | exact to **< 1 nm**, spot **0 nm** |
| Abbe sine ratio | equals `(n_N/n₀)·M` to **10⁻¹⁶**, constant across the fan |
| Identity `sin u₀/sin u_N = (n_N/n₀)·M·𝓜(ρₖ)` (Eqs. 77/86/87) | holds to **10⁻¹⁵** on aplanatic *and* non-aplanatic designs |
| Through the real `render.py` pipeline | foci to **≤ 6·10⁻⁵ mm**, sine ratio constant to ~10⁻⁶ |

That last identity is the strongest check available: it ties `inversa_vc`,
`condicion_aplanatismo` and the exact tracer together, and fails loudly if any one of them
is wrong.

A geometric detail that makes the tracer clean: by the very property that defines `V̄C̄ₖ`
(Eq. 90), the normal at any point of a Cartesian surface crosses the axis at
`Cₖ = ζₖ + V̄C̄ₖ`. So the normal direction is `n̂ ∝ (q·(z − ζ) − 1, q·r)` with `q = 1/V̄C̄ₖ`
— which stays finite even for a flat surface, where `q = 0` and the normal is axial.

## Two traps worth knowing

**Pick the right branch of the surface.** A ray coming from the left hits the *far* side of
a sphere, or the back branch of an ovoid, long before it reaches the designed surface. The
intersection search must be restricted to the useful cap (`ρ ≤ ρ(r_max)`); without that
bound, types 0 and 3 silently produced completely wrong foci.

**Do not iterate an aperture estimate.** Sizing a lens by tracing, measuring the maximum ray
height and shrinking the aperture to match converges to a *false* fixed point: once the
aperture tightens, the rays that exceed it get vignetted and stop contributing to the
maximum. Measure once with a generous trial aperture, then clamp to the radius where the
two surfaces of the lens cross — beyond that the element has no thickness at all.

## A note on Table 5 of the thesis

The free-surface mode reproduces surface 0 of Table 5 exactly
(`G = −2.719771, O = 0.025006, T = −2.613020·10⁻⁸, S = −4.215636·10⁻⁵`) but not surface 1.
The reason is an inconsistency in the source: the table places the final image at `d₂ = 200`
while the surrounding text specifies `d₂ = 120 mm`. Sweeping `d₂` with `d₁ = 599.172082`
fixed puts the minimum of `(M − 1)_RMS` at **d₂ ≈ 120.17**, confirming the text rather than
the table.

## What is not implemented yet

- **Approximate (best-form) design** by sweeping `d₁` to minimise `(M − 1)_RMS` — the metric
  and the exact tracer are both in place, only the optimisation loop is missing (§4.3.7 and
  Fig. 16 of the thesis).
- **Primary aberration coefficients** (chapter 5) and **achromatism** (chapter 6).

---

Back to [Overview](../README.md).
