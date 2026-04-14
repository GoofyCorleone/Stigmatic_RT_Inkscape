# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Context

New project combining Inkscape SVG output with Cartesian oval (Descartes ovoid) ray tracing. Likely produces publication-quality 2D vector diagrams of optical systems built on the GOTS formulation.

Closely related to `../RayTracing/` — the `gots` package there is the canonical ray-tracing engine. See `../RayTracing/CLAUDE.md` and `../RayTracing/README.md` for full details on:

- GOTS parameter computation (`calcular_gots`)
- Surface geometry, ray–surface intersection (quartic solver), vectorial Snell's law
- `SistemaOptico` multi-surface system and LSOE singlet factory
- STL export and existing matplotlib visualizations

## Language

Code, comments, and variable names are in **Spanish** (e.g., `fuente`, `trazar_rayos`, `angulo_max`). Maintain this convention.

## Dependencies

Expected stack (mirrors the sibling project):

```
numpy matplotlib scipy
```

For Inkscape SVG generation add `svgwrite` or use `matplotlib` with an SVG backend.

## Setup

No shared virtual environment. Create one locally:

```bash
cd Python/Inkscape_CartesianSurfaces_RayTracing
python -m venv venv && source venv/bin/activate
pip install numpy matplotlib scipy
```
