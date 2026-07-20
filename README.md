# Stigmatic_RT_Inkscape

**Aberration-free optical design inside Inkscape.** A fork of
[Damien Bloch's inkscape-raytracing](https://github.com/damienBloch/inkscape-raytracing)
that adds generators for **rigorously stigmatic** and **rigorously aplanatic**
refracting elements built from Cartesian ovals, plus the ray tracer to check them.

<p align="center">
  <a href="README.md"><b>Overview &amp; Install</b></a> &nbsp;·&nbsp;
  <a href="docs/basic-elements.md">Basic optical elements</a> &nbsp;·&nbsp;
  <a href="docs/stigmatism.md">Stigmatism</a> &nbsp;·&nbsp;
  <a href="docs/aplanatism.md">Aplanatism</a>
</p>

---

> **Honest disclaimer.** This repository is **not** an original work. It is a fork of the
> excellent [inkscape-raytracing](https://github.com/damienBloch/inkscape-raytracing)
> extension by **Damien Bloch**. The ray-tracing engine, the rendering pipeline and the SVG
> optical-object conventions are his. What is added here — the three generator extensions,
> the `gots_util.py` mathematical core and the standalone example scripts — was written by
> **Claude** (Opus 4.6 / 4.8) under my supervision, following the theory of Alberto Luis
> Silva-Lora's doctoral thesis. Please cite and star the original project:
> https://github.com/damienBloch/inkscape-raytracing

---

## What this is

Ordinary lens design starts from spheres and then fights the aberrations they introduce.
This project goes the other way: it builds surfaces that are **exactly** free of a given
aberration by construction, using the Cartesian-oval (GOTS) formulation, and then verifies
them by tracing rays through the very SVG geometry that Inkscape writes.

Three generators are added to Inkscape's **Extensions → Optics** menu:

| Extension | What it builds | Guarantee |
|---|---|---|
| **Superficie Cartesiana** | one Cartesian-oval refracting surface | exact point-to-point imaging for one conjugate pair |
| **Lente Ovoide (LSOE)** | a stigmatic singlet, both surfaces Cartesian ovals | same, shaped by a form factor σ |
| **Lente Aplanática** | an aplanatic singlet in four closed-form families | exactly **no spherical aberration and no coma**, simultaneously |

Everything they draw is tagged so the upstream **Ray Tracing** command picks it up, and the
results agree with an independent exact tracer to within microns.

## A taste of what comes out

A divergent stigmatic singlet whose *virtual* image feeds a convergent one, then a 45° beam
splitter, a glass slab, a fold mirror and two beam dumps — every primitive the extension
knows, in one pass. All 18 traced trajectories land on both foci within **0.011 mm**:

![Complete system](ejemplo_sistema_completo_traced.svg)

Three aplanatic lenses chained into a 1:1 relay of six Cartesian surfaces. The aplanatism
condition is a *product* over surfaces, so chaining aplanatic elements keeps the system
aplanatic — `(M − 1)_RMS = 4.9·10⁻¹⁶` across all six:

![Three-lens relay](ejemplo_aplanatica_3lentes.svg)

Two stigmatic singlets turning a point source into a plane wavefront, collimated to
**0.009°**:

![Collimator](ejemplo_colimador_traced.svg)

Where to go next: [basic elements](docs/basic-elements.md) for the primitives and the ray
tracer, [stigmatism](docs/stigmatism.md) for the Cartesian-oval surfaces and LSOE lenses,
[aplanatism](docs/aplanatism.md) for the aberration-free-in-two-senses lenses.

---

## Installation

### 1. Install Inkscape 1.2 or newer

Download from https://inkscape.org/release/. The extensions use only `numpy`, `lxml` and
`inkex`, all of which ship with modern Inkscape — there is nothing else to install for
normal use.

### 2. Find your user extensions directory

In Inkscape open **Edit → Preferences → System** and read the **User extensions** field.
The usual locations are:

| OS | Path |
|---|---|
| **Linux** | `~/.config/inkscape/extensions/` |
| **macOS** | `~/Library/Application Support/org.inkscape.Inkscape/config/inkscape/extensions/` |
| **Windows** | `%APPDATA%\inkscape\extensions\` |

### 3. Copy the extension files

```bash
git clone https://github.com/GoofyCorleone/Stigmatic_RT_Inkscape.git
cd Stigmatic_RT_Inkscape
```

Copy the **contents of** `inkscape-raytracing/inkscape_raytracing/` into that directory:

**Linux**
```bash
cp -r inkscape-raytracing/inkscape_raytracing/* ~/.config/inkscape/extensions/
```

**macOS**
```bash
cp -r inkscape-raytracing/inkscape_raytracing/* \
      ~/Library/Application\ Support/org.inkscape.Inkscape/config/inkscape/extensions/
```

**Windows** (PowerShell)
```powershell
Copy-Item -Recurse -Force `
  .\inkscape-raytracing\inkscape_raytracing\* `
  "$env:APPDATA\inkscape\extensions\"
```

Restart Inkscape. Under **Extensions → Optics** you should now see:

- **Ray Tracing** — traces rays through the scene *(upstream)*
- **Set material as** — tags a shape as mirror, glass, beam… *(upstream)*
- **Superficie Cartesiana** — one stigmatic surface *(new)*
- **Lente Ovoide (LSOE)** — stigmatic singlet *(new)*
- **Lente Aplanática** — aplanatic singlet *(new)*

> **Troubleshooting.** If the new entries do not appear, check that the `.inx` files landed
> directly in the extensions directory (not in a sub-folder) and look at
> **Extensions → Previous Extension** or Inkscape's error dialog for a Python traceback.
> On Windows, make sure the copy did not create a nested
> `extensions\inkscape_raytracing\` folder.

### 4. Optional: running the scripts outside Inkscape

The example figures in this documentation are produced by standalone scripts that do not
need Inkscape at all. For those you want a virtual environment:

```bash
python -m venv venv
source venv/bin/activate          # Windows: venv\Scripts\activate
pip install numpy scipy matplotlib

python generar_ejemplos.py            # the basic before/after figures
python generar_ejemplos_avanzados.py  # collimator, multi-element rig, full system
python generar_ejemplos_aplanaticos.py  # the four aplanatic types + 3-lens relay
python generar_capturas_ui.py         # regenerate the dialog images in docs/img/
```

To run the test suite:

```bash
cd inkscape-raytracing
python -m pytest tests/unit -q
```

> `tests/integration/` currently fails to collect: those upstream tests import a module
> layout this fork no longer has. Run `tests/unit` only.

---

## Repository layout

```
README.md                          this page
docs/
  basic-elements.md                upstream primitives and the ray tracer
  stigmatism.md                    Cartesian surfaces and LSOE lenses
  aplanatism.md                    aplanatic lenses
  img/                             dialog reference images
generar_lsoe_svg.py                minimal standalone LSOE generator
generar_ejemplos.py                basic before/after example figures
generar_ejemplos_avanzados.py      collimator, multi-element rig, full system
generar_ejemplos_aplanaticos.py    aplanatic types and three-lens relay
generar_capturas_ui.py             renders docs/img/ from the .inx files
ejemplo_*.svg                      pre-generated example figures
inkscape-raytracing/
  inkscape_raytracing/
    gots_util.py                   ★ mathematical core (numpy only)
    superficie_cartesiana.py/.inx  ★ single stigmatic surface
    lente_ovoide.py/.inx           ★ stigmatic singlet
    lente_aplanatica.py/.inx       ★ aplanatic singlet
    render.py/.inx                 ray tracer (upstream, precision-patched)
    raytracing/                    ray-tracing engine (upstream)
  tests/unit/                      geometry, math and aplanatic designs
```

A note on language: the documentation is in English, but **the source code and the
extension dialogs are in Spanish** (`fuente`, `apertura`, `superficies`, …), which is the
convention of the wider repository this project lives in. The dialog images throughout the
docs therefore show Spanish labels; each page explains what every field means.

---

## Credits and license

- Upstream engine and Inkscape integration: **Damien Bloch**,
  [inkscape-raytracing](https://github.com/damienBloch/inkscape-raytracing).
- Optical theory: **Alberto Luis Silva Lora**, *Estudio de las aberraciones primarias a
  partir de la teoría del estigmatismo riguroso*, Universidad Industrial de Santander, 2024,
  and the Silva-Lora & Torres papers cited therein.
- This fork follows the upstream license.
