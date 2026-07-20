# Stigmatism

<p align="center">
  <a href="../README.md">Overview &amp; Install</a> &nbsp;·&nbsp;
  <a href="basic-elements.md">Basic optical elements</a> &nbsp;·&nbsp;
  <a href="stigmatism.md"><b>Stigmatism</b></a> &nbsp;·&nbsp;
  <a href="aplanatism.md">Aplanatism</a>
</p>

---

## The theory in one page

A surface is **stigmatic** for a pair of points `P₀` and `P₁` when *every* ray leaving `P₀`
arrives at `P₁` — not approximately, not for small angles, but exactly. Fermat's principle
says this happens when the optical path length is the same along all of them:

```
n₀ · |P − P₀|  +  n₁ · |P − P₁|  =  constant
```

for every point `P` on the surface. That equation defines a quartic surface of revolution
known since Descartes as the **Cartesian oval**. A sphere is the special case; the general
oval is what you need when the conjugates are not at the Young points.

Two details matter in practice:

**Virtual images flip a sign.** If the image is virtual the invariant is the *difference*
of optical paths, `n₀·l₀ − n₁·l₁ = k`, not the sum. Using the sum regardless produces a
surface that simply is not stigmatic. Internally the code carries signed indices
`m₀ = ±n₀`, `m₁ = ±n₁` so both cases come out of one formula.

**The shape parameters.** Silva-Lora's formulation (the GOTS parameters) writes the surface
profile as

```
τ(ρ) = (O + T·ρ²)·ρ² / [ 1 + S·ρ² + √(1 + (2S − O·OG)·ρ²) ]
z(ρ) = ζ + τ(ρ)
```

where `ρ` is the distance from the vertex to the point on the surface, and

| Symbol | Meaning | Units |
|---|---|---|
| `O` | paraxial curvature, `1/R`, signed | 1/length |
| `T` | fourth-order deformation | 1/length³ |
| `S` | conic-coupling term | 1/length² |
| `G` | axial shift; `OG = O·G` | length |
| `ζ` | vertex position on the axis | length |

`calcular_gots(n_k, n_k1, ζ, d_k, d_k1)` converts the physical description — the indices
and the two conjugate positions — into `(G, O, T, S, OG)` through Silva-Lora's Eqs. 10–13.
Everything else in the project is built on top of that one function.

### The mathematical core (`gots_util.py`)

A self-contained module, numpy only, shared by all three extensions and by the standalone
scripts:

| Function | What it does |
|---|---|
| `calcular_gots(n_k, n_k1, ζ, d_k, d_k1)` | physical parameters → GOTS (Eqs. 10–13), with degeneracy checks |
| `perfil_superficie(params, N, r_max)` | meridional profile `(r, z)` from GOTS (Eq. 16) |
| `perfil_ovoide_descartes(n0, n1, ζ, d_obj, d_img, N)` | closed Descartes oval straight from the Fermat invariant, signed for real/virtual images |
| `encontrar_apertura(p0, p1)` | radius at which two surface profiles intersect — the largest physical lens |
| `calcular_d1_sigma(σ, ζ0, ζ1, d0, d2, n0, n1, n2)` | intermediate image position realising a given shape factor |
| `puntos_a_bezier_path_str(pts, cerrar)` | points → SVG cubic Bézier with exact non-uniform (Bessel) tangents |

The aplanatic layer built on top of this is documented on the
[Aplanatism](aplanatism.md#how-it-is-verified) page.

Going the other way — GOTS back to physical parameters — is not implemented, but it is
straightforward algebra: `O` and `S` give two equations in the two unknowns `ξ = d_k − ζ`
and `η = d_{k+1} − ζ`.

## Superficie Cartesiana — a single surface

**Extensions → Optics → Superficie Cartesiana** draws one refracting surface. It has three
tabs.

### Tab 1 — Parámetros Físicos

![Superficie Cartesiana, physical tab](img/ui-superficie-1.svg)

This is the tab you normally use: describe the optics you want and let the extension work
out the shape.

| Field (Spanish) | Meaning |
|---|---|
| `n₁` | refractive index *before* the surface (object side) |
| `n₂` | refractive index *after* the surface |
| `d₀ (punto objeto)` | axial position of the object point |
| `d₁ (punto imagen)` | axial position of the image point |
| `ζ (vértice)` | axial position of the surface vertex |

All positions are absolute coordinates on the optical axis, so the sign convention takes
care of itself: an object to the left of the vertex has `d₀ < ζ`, and an image to the right
has `d₁ > ζ`. If `d₁ < ζ` the image is virtual and the extension automatically switches to
the difference invariant and draws a concave surface closed by a rectangular glass body.

### Tab 2 — Parámetros GOTS

![Superficie Cartesiana, GOTS tab](img/ui-superficie-2.svg)

Skip the physics and give the shape parameters directly. This is for teaching figures and
parameter sweeps: fix three of `(G, O, T, S)` and vary the fourth to see what each one does
to the profile.

### Tab 3 — Geometría

![Superficie Cartesiana, geometry tab](img/ui-superficie-3.svg)

Presentation only: aperture radius (`0` = pick automatically), drawing units, whether to
draw the optical axis and the object/image markers, and how dense the test ray fan should
be.

### Result

A surface designed for `n₁ = 1`, `n₂ = 1.5`, object at −100 mm, image at +200 mm, before
and after tracing:

| Before | After tracing |
|---|---|
| ![Cartesian oval](../ejemplo_cartesiana.svg) | ![Cartesian oval traced](../ejemplo_cartesiana_traced.svg) |

Every ray meets the green image point — that is what "rigorously stigmatic" buys you, at a
numerical aperture where a sphere would smear the focus badly.

## Lente Ovoide (LSOE) — a stigmatic singlet

A **LSOE** (*Lente Singlete Ovoide Estigmática*) is a lens whose **two** surfaces are
Cartesian ovals, arranged so the intermediate image of the first is the object of the
second. The pair is then stigmatic for the overall conjugates.

There is one degree of freedom left: where to put that intermediate image. Moving it
redistributes curvature between the two faces without changing the conjugates, which is
exactly the classical **shape factor σ**:

| σ | Shape |
|---|---|
| `−1` | plano-convex, flat front |
| `0` | symmetric biconvex |
| `+1` | convex-plano, flat back |

`calcular_d1_sigma` solves for the intermediate image position that realises a requested σ.

### Tab 1 — Diseño de la lente

![LSOE, design tab](img/ui-lsoe-1.svg)

| Field | Meaning |
|---|---|
| `n₀ / n₁ / n₂` | object space, lens glass, image space |
| `ζ₀ / ζ₁` | vertices of the front and back surfaces (their difference is the axial thickness) |
| `d₀ / d₂` | object and image positions |
| `σ` | shape factor |

### Tab 2 — Apertura y Haz

![LSOE, aperture tab](img/ui-lsoe-2.svg)

Aperture radius (`0` lets the extension use the radius at which the two ovals intersect,
which is the physically largest lens you can cut), number of rays, fan half-angle, and the
axis/marker switches.

### Result

A biconvex LSOE for object at 0 mm, image at 200 mm, `n₁ = 1.6`, σ = 0, with 11 rays out to
±7°:

| Before | After tracing |
|---|---|
| ![LSOE](../ejemplo_lsoe.svg) | ![LSOE traced](../ejemplo_lsoe_traced.svg) |

## Divergent configurations

Both generators also produce rigorously stigmatic **virtual** images. Here a biconcave LSOE
and a concave Cartesian surface, both with the object at −100 mm and the virtual image at
−40 mm — between the object and the element, as a diverging system requires:

| Divergent LSOE — before | after |
|---|---|
| ![LSOE div](../ejemplo_lsoe_divergente.svg) | ![LSOE div traced](../ejemplo_lsoe_divergente_traced.svg) |

| Divergent surface — before | after |
|---|---|
| ![Cartesian div](../ejemplo_cartesiana_divergente.svg) | ![Cartesian div traced](../ejemplo_cartesiana_divergente_traced.svg) |

The dashed back-extrapolations of the refracted rays meet exactly at the virtual image
point.

> **Why the divergent surface is drawn with a body.** The full closed Descartes oval is not
> traceable in this configuration: its back cap sits between the object and the designed
> surface and would intercept the rays first. The extension therefore draws only the
> designed concave branch and closes it with a rectangular glass body toward `+z`.

## Composite systems

Stigmatism survives chaining, so you can build real instruments. Two LSOEs turning a point
source into a plane wavefront — the second lens images the first one's focus to infinity:

| Before | After tracing |
|---|---|
| ![Collimator](../ejemplo_colimador.svg) | ![Collimator traced](../ejemplo_colimador_traced.svg) |

Collimation error: **0.009°**.

And the full rig — a divergent LSOE whose virtual image is the object of a convergent LSOE,
then a beam splitter feeding a glass slab on one arm and a fold mirror on the other:

| Before | After tracing |
|---|---|
| ![Complete](../ejemplo_sistema_completo.svg) | ![Complete traced](../ejemplo_sistema_completo_traced.svg) |

Traced through the *actual* Inkscape extension, all 18 trajectories (9 rays × 2 arms)
converge at both foci to within **0.011 mm**.

### Rule of thumb

| Configuration | Object `d₀` | Image | Shape |
|---|---|---|---|
| Convergent, real image | any | on the far side of the element, or further away than the object | convex / biconvex |
| Divergent, virtual image | real (before the vertex) | same side as the object and **closer** to it | concave / biconcave |

If the extension stops with *"parámetro degenerado"*, the configuration is geometrically
impossible — object or image sitting on the vertex, `n₁ = n₂`, or `κ = n₁·η − n₀·ξ = 0`.
Move `d₀`, `d₁` or `ζ`.

## Tutorial: build one yourself

1. Run **Extensions → Optics → Lente Ovoide (LSOE)**.
2. On *Diseño*, set `n₁ = 1.6`, `ζ₀ = 0`, `ζ₁ = 10`, `d₀ = −80`, `d₂ = 110`, `σ = 0`.
3. On *Apertura y Haz*, leave the aperture at `0`, set 11 rays and a 7° half-angle.
4. Press **Apply**. You get the lens, the axis, the object/image markers and the seed rays.
5. Press `Ctrl+A`, then **Extensions → Optics → Ray Tracing**.

The rays should meet at the green marker. Now go back and change σ to `−1` and re-run: the
lens becomes plano-convex, the shape changes completely — and the focus stays exact. That
is the point of the whole formulation.

To do the same without opening Inkscape:

```bash
python generar_lsoe_svg.py       # writes lsoe_raytracing.svg, ready to trace
python generar_ejemplos.py       # regenerates the before/after figures above
```

## Numerical accuracy

Getting exact optics through an SVG round-trip is harder than the maths. Three fixes make
the in-Inkscape result agree with an independent exact tracer to ≈0.001° and a few µm:

1. **Spline→Bézier handles.** The profile is converted with the exact non-uniform
   (Bessel) formula `P₁ = pᵢ + mᵢ·hᵢ/3`. Fixed-tension handles were twice too long and
   bulged the surface between nodes at order `h²`, which visibly walked the focus.
2. **inkex path precision.** `render.py` reads `d` attributes directly and raises inkex's
   number template from `{:.6g}` to `{:.10g}`. The default six significant digits quantise
   coordinates to ~10⁻³ px, corrupting the normals of fine Bézier segments — errors that
   got *worse* with denser sampling, which is how the bug was eventually found.
3. **Point decimation and closure.** Nearly coincident profile points are merged before
   emitting the path, and a closed contour never repeats its first point: a zero-length
   closing chord degenerates the Bessel tangent and puts a microscopic crease at the
   vertex, which deflects precisely the axial rays.

Traced beams are also drawn compensating the containing layer's transform, and beams
outside any layer fall back to a root-level `generated_beams` layer instead of being
dropped.

---

Next: [Aplanatism](aplanatism.md) — removing coma as well, and where the image *really*
forms.
