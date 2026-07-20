# Basic optical elements

<p align="center">
  <a href="../README.md">Overview &amp; Install</a> &nbsp;·&nbsp;
  <a href="basic-elements.md"><b>Basic optical elements</b></a> &nbsp;·&nbsp;
  <a href="stigmatism.md">Stigmatism</a> &nbsp;·&nbsp;
  <a href="aplanatism.md">Aplanatism</a>
</p>

---

Everything on this page comes from the **upstream**
[inkscape-raytracing](https://github.com/damienBloch/inkscape-raytracing) extension by
Damien Bloch. It is the foundation the rest of the project builds on, so it is worth
understanding before moving to the stigmatic and aplanatic generators.

## The idea

There is no separate "optics file format". You draw ordinary Inkscape shapes, and then you
**tag** them with what they are made of. The tag lives in the object's description field
(`Object → Object Properties → Description`) and looks like this:

```
optics:mirror
optics:glass:1.5168
optics:beam
optics:beam_dump
optics:beam_splitter
```

Then you select everything and run the tracer, which walks the rays through the scene.
Because the tag is just text in the SVG, anything that can write an SVG can build an
optical system — which is exactly how the generators in this project work.

## The five materials

| Tag | Behaviour |
|---|---|
| `optics:beam` | **not** an object: a seed ray. The tracer takes the *end point* of the path as the ray origin and the path's final tangent as its direction |
| `optics:mirror` | reflects |
| `optics:glass:n` | refracts by Snell's law with refractive index `n`; the shape must be closed, and the tracer treats its inside as the medium |
| `optics:beam_splitter` | splits each incoming ray into a transmitted and a reflected one |
| `optics:beam_dump` | absorbs — use it for apertures, stops and screens |

The document border acts as an implicit beam dump, so rays never escape to infinity.

## Tagging a shape: *Set material as*

Select a shape and run **Extensions → Optics → Set material as**. It writes the right tag
into the description for you:

![Set material dialog](img/ui-set-material.svg)

| Field | Meaning |
|---|---|
| **Select material** | which `optics:` tag to write |
| **Optical index** | only used when the material is *Glass*; `1.5168` is N-BK7 at the sodium line |

You can equally type the tag by hand into **Object Properties → Description** — the tracer
only reads the text.

## Tracing: *Ray Tracing*

Select the objects you want included (usually `Ctrl+A`) and run
**Extensions → Optics → Ray Tracing**. It has no parameters:

![Ray Tracing dialog](img/ui-raytracing.svg)

The traced rays are written into a new layer called `generated_beams`, inside the layer
that held the seed beam. Because it is a normal layer you can restyle, hide or delete the
result and trace again after changing the scene.

## A first system, step by step

1. **Draw the optics.** A closed path for a lens, a straight segment for a mirror.
2. **Tag them.** Select the lens → *Set material as* → *Glass*, index `1.5`. Select the
   mirror → *Set material as* → *Mirror*.
3. **Draw the seed rays.** One short straight segment per ray. Remember the direction rule:
   the ray starts at the **end** of the segment and continues along its direction, so draw
   each segment *pointing where the light should go*.
4. **Tag the rays** as *Beam*.
5. **Select all** and run *Ray Tracing*.

A worked example that already contains every tag is
[`lsoe_raytracing.svg`](../lsoe_raytracing.svg) — open it in Inkscape, press `Ctrl+A`, and
run the tracer.

## What a traced scene looks like

Here is a single-pass rig that uses **every** primitive at once — lens, beam splitter,
mirror, glass block and two beam dumps — before and after tracing:

| Before | After tracing |
|---|---|
| ![Experiment](../ejemplo_experimento.svg) | ![Experiment traced](../ejemplo_experimento_traced.svg) |

Each ray refracts at the lens, splits at the 45° beam splitter into a transmitted and a
reflected copy, and then either crosses the glass slab (transmitted arm, picking up the two
parallel shifts a plate produces) or bounces off the fold mirror (reflected arm), before
being absorbed by a dump.

## Practical notes

**Draw beams in the direction of propagation.** This trips up everyone once. The seed is
the path's *end point*, so a segment drawn from the lens back to the source will send light
the wrong way.

**Glass shapes must be closed.** An open path has no inside, and the tracer cannot decide
when a ray is entering or leaving the medium.

**Beams must live inside an Inkscape layer.** Upstream, a beam sitting at the document root
was silently ignored. This fork adds a fallback that puts the generated layer at the root
instead, so externally generated SVGs (like all the examples here) trace correctly.

**Beam splitters double the ray count.** Every splitter multiplies the number of
trajectories by two, so a scene with several of them and a recursive path can explode.
Terminate arms with beam dumps.

**Precision.** For finely sampled curves, this fork raises inkex's number template from
`{:.6g}` to `{:.10g}` and reads path `d` attributes directly. Without that, coordinates get
quantised to ~10⁻³ px, which corrupts the normals of small Bézier segments and produces
refraction errors of up to ~1°. See the note at the end of
[Stigmatism](stigmatism.md#numerical-accuracy) for the full story.

---

Next: [Stigmatism](stigmatism.md) — surfaces that image a point *perfectly*.
