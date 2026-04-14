"""Genera un SVG con una lente LSOE y un haz de rayos para ray-tracing.

Crea el archivo lsoe_raytracing.svg con:
- La lente biconvexa LSOE (superficies cartesianas exactas, marcada como glass)
- Un haz de rayos divergentes desde el punto objeto (marcados como beam)
- El eje óptico y los marcadores de objeto/imagen

El SVG está listo para ser procesado con la extensión inkscape-raytracing:
  Extensions > Optics > Ray Tracing
"""

import sys
import os
import numpy as np

# Añadir el directorio de la extensión al path para importar gots_util
EXT_DIR = os.path.join(
    os.path.dirname(__file__),
    "inkscape-raytracing", "inkscape_raytracing"
)
sys.path.insert(0, EXT_DIR)

from gots_util import (
    calcular_gots,
    perfil_superficie,
    encontrar_apertura,
    calcular_d1_sigma,
    puntos_a_bezier_path_str,
)

# ── Parámetros del sistema ──────────────────────────────────────────────────
# LSOE de referencia (Tabla 4 de Silva-Lora 2024)
N0    = 1.0          # índice del espacio objeto (aire)
N1    = 1.6          # índice del vidrio
N2    = 1.0          # índice del espacio imagen (aire)
ZETA0 = 80.0         # mm — vértice superficie frontal
ZETA1 = 90.0         # mm — vértice superficie trasera
D0    = 0.0          # mm — punto objeto
D2    = 200.0        # mm — punto imagen
SIGMA = 0.0          # biconvexo simétrico

N_RAYOS   = 11       # número de rayos en el haz
ANGULO    = 7.0      # grados — apertura angular del haz

# ── Escala y documento SVG ──────────────────────────────────────────────────
# 1 mm → SCALE px (SVG usa px como unidad interna; elegimos 2 px/mm para
# que el documento sea razonable a 96 DPI)
SCALE   = 2.0        # px por mm
MARGEN  = 30.0       # mm de margen en el SVG

# ── Calcular geometría ──────────────────────────────────────────────────────
d1 = calcular_d1_sigma(SIGMA, ZETA0, ZETA1, D0, D2, N0, N1, N2)
print(f"d₁ (sigma={SIGMA}): {d1:.4f} mm")

p0 = calcular_gots(N0, N1, ZETA0, D0, d1)
p1 = calcular_gots(N1, N2, ZETA1, d1, D2)

r_ap = encontrar_apertura(p0, p1)
print(f"Apertura: r = {r_ap:.4f} mm")

r0_arr, z0_arr = perfil_superficie(p0, N=600, r_max=r_ap)
r1_arr, z1_arr = perfil_superficie(p1, N=600, r_max=r_ap)

# Normalizar: vértice frontal en el origen (x = z − ZETA0, y = −r)
def px(val_mm):
    return val_mm * SCALE

def xcoord(z_mm):
    return px(z_mm - ZETA0)

def ycoord_upper(r_mm):
    return -px(r_mm)

def ycoord_lower(r_mm):
    return px(r_mm)

# ── Construir puntos del path de la lente ──────────────────────────────────
puntos = []
# vértice frontal
puntos.append((xcoord(z0_arr[0]), 0.0))
# sup frontal — mitad superior
for z, r in zip(z0_arr[1:], r0_arr[1:]):
    puntos.append((xcoord(z), ycoord_upper(r)))
# borde superior (rim)
puntos.append((xcoord(z1_arr[-1]), ycoord_upper(r1_arr[-1])))
# sup trasera — mitad superior invertida
for z, r in zip(reversed(list(z1_arr[:-1])), reversed(list(r1_arr[:-1]))):
    puntos.append((xcoord(z), ycoord_upper(r)))
# sup trasera — mitad inferior
for z, r in zip(z1_arr[1:], r1_arr[1:]):
    puntos.append((xcoord(z), ycoord_lower(r)))
# borde inferior
puntos.append((xcoord(z0_arr[-1]), ycoord_lower(r0_arr[-1])))
# sup frontal — mitad inferior invertida
for z, r in zip(reversed(list(z0_arr[:-1])), reversed(list(r0_arr[:-1]))):
    puntos.append((xcoord(z), ycoord_lower(r)))

lens_path_d = puntos_a_bezier_path_str(puntos, cerrar=True)

# ── Haces de rayos desde el punto objeto ───────────────────────────────────
x_fuente  = xcoord(D0)
# Indicadores: desde el margen izquierdo del documento hasta la fuente
# doc_x_min = xcoord(D0) - px(MARGEN), así que la distancia fuente→borde = px(MARGEN)
L_haz     = px(MARGEN * 0.85)   # ~85% del margen para no tocar el borde
angulos  = np.linspace(-np.radians(ANGULO), np.radians(ANGULO), N_RAYOS)

beam_paths = []
for theta in angulos:
    # Indicador dibujado desde la fuente hacia la lente (+x).  El ray-tracer
    # toma como origen el endpoint del path y dirección −tangente final, por
    # lo que para «M source L end» el rayo nace en ``end`` con dirección
    # (end − source).  La extrapolación inversa pasa por la fuente, así que
    # es equivalente a un rayo originado en el punto objeto.
    x_end = x_fuente + L_haz * np.cos(theta)
    y_end = L_haz * np.sin(theta)
    beam_paths.append(
        f"M {x_fuente:.3f},0 L {x_end:.3f},{y_end:.3f}"
    )

# ── Dimensiones del documento ───────────────────────────────────────────────
doc_x_min = xcoord(D0) - px(MARGEN)
doc_x_max = xcoord(D2) + px(MARGEN)
doc_y_min = ycoord_upper(r_ap) - px(MARGEN * 0.5)
doc_y_max = ycoord_lower(r_ap) + px(MARGEN * 0.5)
doc_w     = doc_x_max - doc_x_min
doc_h     = doc_y_max - doc_y_min

# Transformación para poner el origen SVG en (0,0) — desplazar todo
tx = -doc_x_min
ty = -doc_y_min

# ── Generar SVG ─────────────────────────────────────────────────────────────
def svg_move(d_str, tx, ty):
    """Desplaza un path en el SVG (traslación simple)."""
    return d_str  # Los paths ya están en coordenadas absolutas; usaremos transform en el grupo

OUT_FILE = os.path.join(os.path.dirname(__file__), "lsoe_raytracing.svg")

beam_elements = "\n".join(
    f'''    <path
       style="fill:none;stroke:#ff6600;stroke-width:1px"
       d="{bp}">
      <desc>optics:beam</desc>
    </path>'''
    for bp in beam_paths
)

# Posición del eje óptico
x_eje_izq = xcoord(D0 - MARGEN * 0.8)
x_eje_der = xcoord(D2 + MARGEN * 0.8)

# Círculos objeto e imagen
r_marker  = px(1.5)
x_obj_svg = xcoord(D0)
x_img_svg = xcoord(D2)

svg_content = f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg
   xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
   xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.0.dtd"
   xmlns="http://www.w3.org/2000/svg"
   xmlns:xlink="http://www.w3.org/1999/xlink"
   version="1.1"
   width="{doc_w:.2f}"
   height="{doc_h:.2f}"
   viewBox="0 0 {doc_w:.2f} {doc_h:.2f}"
   inkscape:version="1.4"
   sodipodi:docname="lsoe_raytracing.svg">
  <sodipodi:namedview
     inkscape:current-layer="layer1"
     inkscape:document-units="px"
     units="px"/>
  <g
     inkscape:label="Layer 1"
     inkscape:groupmode="layer"
     id="layer1"
     transform="translate({tx:.3f},{ty:.3f})">

    <!-- Eje óptico -->
    <path
       style="fill:none;stroke:#aaaaaa;stroke-width:0.5px;stroke-dasharray:6,3"
       d="M {x_eje_izq:.3f},0 L {x_eje_der:.3f},0"/>

    <!-- Lente LSOE (σ={SIGMA}, n₁={N1}) -->
    <path
       style="opacity:0.8;fill:#b7c2dd;stroke:#000000;stroke-width:0.8px;stroke-linejoin:round"
       d="{lens_path_d}">
      <desc>optics:glass:{N1:.4f}</desc>
    </path>

    <!-- Punto objeto (rojo) -->
    <circle cx="{x_obj_svg:.3f}" cy="0" r="{r_marker:.3f}"
            style="fill:#dd2200;stroke:none"/>

    <!-- Punto imagen (verde) -->
    <circle cx="{x_img_svg:.3f}" cy="0" r="{r_marker:.3f}"
            style="fill:#00aa33;stroke:none"/>

    <!-- Haces de rayos desde el punto objeto -->
{beam_elements}

  </g>
</svg>
"""

with open(OUT_FILE, "w", encoding="utf-8") as f:
    f.write(svg_content)

print(f"\nSVG generado: {OUT_FILE}")
print(f"Documento: {doc_w:.0f} × {doc_h:.0f} px  (escala {SCALE} px/mm)")
print(f"Lente: vértices en z={ZETA0}–{ZETA1} mm, apertura r={r_ap:.1f} mm")
print(f"Objeto en z={D0} mm, imagen en z={D2} mm")
print(f"\nPara ver el ray-tracing:")
print(f"  1. Abrir Inkscape: open {OUT_FILE}")
print(f"  2. Seleccionar todo (Ctrl+A)")
print(f"  3. Extensions > Optics > Ray Tracing")
