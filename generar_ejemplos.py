"""Genera SVGs de ejemplo para el README, incluyendo una versión con los
rayos trazados exactamente por el motor canónico (``gots/sistema_optico``)
para comprobar que la convergencia coincide con el punto imagen teórico.

Salidas:
  - ``ejemplo_lsoe.svg``           LSOE biconvexa σ=0 con haz divergente.
  - ``ejemplo_lsoe_traced.svg``    La misma lente con los rayos refractados.
  - ``ejemplo_cartesiana.svg``     Óvalo de Descartes completo (una superficie).
  - ``ejemplo_cartesiana_traced.svg``
"""

import os
import sys
import numpy as np

AQUI = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(AQUI, "inkscape-raytracing", "inkscape_raytracing"))
sys.path.insert(0, os.path.join(AQUI, "..", "RayTracing"))

from gots_util import (
    calcular_gots,
    perfil_superficie,
    perfil_ovoide_descartes,
    encontrar_apertura,
    calcular_d1_sigma,
    puntos_a_bezier_path_str,
)
from gots.sistema_optico import SistemaOptico
from gots.superficie_cartesiana import SuperficieCartesiana
from gots.rayo import Rayo


SCALE  = 2.0        # px/mm
MARGEN = 30.0       # mm


# ── utilidades SVG ─────────────────────────────────────────────────────────

def px(v):        return v * SCALE
def xc(z, zref):  return px(z - zref)
def yc(r):        return -px(r)        # r>0 arriba del eje óptico


def svg_header(xmin, xmax, ymin, ymax):
    w = xmax - xmin
    h = ymax - ymin
    return (
        f'<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
        f'<svg xmlns="http://www.w3.org/2000/svg" version="1.1" '
        f'width="{w:.2f}" height="{h:.2f}" '
        f'viewBox="0 0 {w:.2f} {h:.2f}">\n'
        f'  <g transform="translate({-xmin:.3f},{-ymin:.3f})">\n'
    )


def svg_footer():
    return "  </g>\n</svg>\n"


def path_lente(puntos_xy, fill="#b7c2dd", stroke="#000", sw=0.6, beam_tag=False):
    d   = puntos_a_bezier_path_str(puntos_xy, cerrar=True)
    tag = "<desc>optics:glass:1.6</desc>" if beam_tag else ""
    return (
        f'    <path d="{d}" '
        f'style="fill:{fill};fill-opacity:0.7;stroke:{stroke};'
        f'stroke-width:{sw}px;stroke-linejoin:round"/>\n'
    )


def line(x1, y1, x2, y2, stroke="#888", sw=0.3, dash=None):
    da = f";stroke-dasharray:{dash}" if dash else ""
    return (f'    <path d="M {x1:.3f},{y1:.3f} L {x2:.3f},{y2:.3f}" '
            f'style="fill:none;stroke:{stroke};stroke-width:{sw}px{da}"/>\n')


def circ(cx, cy, r, color):
    return (f'    <circle cx="{cx:.3f}" cy="{cy:.3f}" r="{r:.3f}" '
            f'style="fill:{color};stroke:none"/>\n')


# ── geometría LSOE ─────────────────────────────────────────────────────────

def generar_lsoe(out_path, traced=False):
    N0, N1, N2 = 1.0, 1.6, 1.0
    Z0, Z1     = 80.0, 90.0
    D0, D2     = 0.0, 200.0
    SIGMA      = 0.0
    ZREF       = Z0

    d1 = calcular_d1_sigma(SIGMA, Z0, Z1, D0, D2, N0, N1, N2)
    p0 = calcular_gots(N0, N1, Z0, D0, d1)
    p1 = calcular_gots(N1, N2, Z1, d1, D2)
    r_ap = encontrar_apertura(p0, p1)

    r0, z0 = perfil_superficie(p0, N=600, r_max=r_ap)
    r1, z1 = perfil_superficie(p1, N=600, r_max=r_ap)

    # contorno cerrado de la lente
    pts = []
    pts.append((xc(z0[0], ZREF), 0.0))
    for zz, rr in zip(z0[1:], r0[1:]):
        pts.append((xc(zz, ZREF), yc(rr)))
    pts.append((xc(z1[-1], ZREF), yc(r1[-1])))
    for zz, rr in zip(reversed(list(z1[:-1])), reversed(list(r1[:-1]))):
        pts.append((xc(zz, ZREF), yc(rr)))
    for zz, rr in zip(z1[1:], r1[1:]):
        pts.append((xc(zz, ZREF), px(rr)))
    pts.append((xc(z0[-1], ZREF), px(r0[-1])))
    for zz, rr in zip(reversed(list(z0[:-1])), reversed(list(r0[:-1]))):
        pts.append((xc(zz, ZREF), px(rr)))

    # encuadre
    xmin = xc(D0 - MARGEN, ZREF); xmax = xc(D2 + MARGEN, ZREF)
    ymin = yc(r_ap) - px(MARGEN * 0.5); ymax = px(r_ap) + px(MARGEN * 0.5)

    out  = svg_header(xmin, xmax, ymin, ymax)
    # eje óptico
    out += line(xc(D0 - MARGEN * 0.8, ZREF), 0, xc(D2 + MARGEN * 0.8, ZREF), 0,
                stroke="#aaaaaa", sw=0.5, dash="6,3")
    # lente
    out += path_lente(pts, fill="#b7c2dd", stroke="#000", sw=0.8)

    # marcadores objeto/imagen
    out += circ(xc(D0, ZREF), 0, px(1.5), "#dd2200")
    out += circ(xc(D2, ZREF), 0, px(1.5), "#00aa33")

    N_RAYOS, AMAX = 11, 7.0
    angs = np.linspace(-np.radians(AMAX), np.radians(AMAX), N_RAYOS)

    if not traced:
        # Indicadores (no trazados) desde el objeto hacia la lente
        x_src = xc(D0, ZREF)
        Lh    = px(MARGEN * 0.85)
        for th in angs:
            xe = x_src + Lh * np.cos(th)
            ye = Lh * np.sin(th)
            out += (f'    <path d="M {x_src:.3f},0 L {xe:.3f},{ye:.3f}" '
                    f'style="fill:none;stroke:#ff6600;stroke-width:1px"/>\n')
    else:
        # Trazar con el motor canónico
        sup0 = SuperficieCartesiana.desde_parametros_fisicos(N0, N1, Z0, D0, d1)
        sup1 = SuperficieCartesiana.desde_parametros_fisicos(N1, N2, Z1, d1, D2)
        sis  = SistemaOptico()
        sis.agregar_superficie(sup0); sis.agregar_superficie(sup1)
        # Rayos en el plano (x-y): origen en (D0, 0, 0), direcciones (cos,sen,0)
        for th in angs:
            # eje óptico = z; rayo meridional en plano y-z
            r = Rayo(origen=np.array([0.0, 0.0, D0]),
                     direccion=np.array([0.0, np.sin(th), np.cos(th)]))
            res = sis.trazar_rayo(r)
            pts_r = res.puntos
            dirs  = res.direcciones
            Pf, Df = pts_r[-1], dirs[-1]
            z_end = D2 + MARGEN * 0.85
            t_ext = (z_end - Pf[2]) / Df[2]
            Pend = Pf + t_ext * Df
            poly = list(pts_r) + [Pend]
            segs = " ".join(
                f"{'M' if i == 0 else 'L'} {xc(P[2], ZREF):.3f},{-px(P[1]):.3f}"
                for i, P in enumerate(poly)
            )
            out += (f'    <path d="{segs}" '
                    f'style="fill:none;stroke:#ff6600;stroke-width:0.7px"/>\n')

    out += svg_footer()
    with open(out_path, "w") as f:
        f.write(out)
    print(f"{out_path}  (traced={traced})  apertura r={r_ap:.2f} mm")


# ── geometría single-surface Cartesiana ────────────────────────────────────

def generar_cartesiana(out_path, traced=False):
    N1, N2 = 1.0, 1.5
    Z      = 0.0
    D0, D1 = -100.0, 200.0
    ZREF   = Z

    zs, rs = perfil_ovoide_descartes(N1, N2, Z, D0, D1, N=400)

    pts = [(xc(z, ZREF), yc(r)) for z, r in zip(zs, rs)]
    for z, r in zip(zs[-2:0:-1], rs[-2:0:-1]):
        pts.append((xc(z, ZREF), px(r)))

    r_max = rs.max()
    xmin = xc(D0 - MARGEN, ZREF); xmax = xc(D1 + MARGEN, ZREF)
    ymin = yc(r_max) - px(MARGEN * 0.4); ymax = px(r_max) + px(MARGEN * 0.4)

    out  = svg_header(xmin, xmax, ymin, ymax)
    out += line(xc(D0 - MARGEN * 0.8, ZREF), 0,
                xc(D1 + MARGEN * 0.8, ZREF), 0,
                stroke="#aaaaaa", sw=0.5, dash="6,3")
    out += path_lente(pts, fill="#b7c2dd", stroke="#000", sw=0.8)
    out += circ(xc(D0, ZREF), 0, px(1.5), "#dd2200")
    out += circ(xc(D1, ZREF), 0, px(1.5), "#00aa33")

    N_R, AMAX = 9, 6.0
    angs = np.linspace(-np.radians(AMAX), np.radians(AMAX), N_R)

    if not traced:
        x_src = xc(D0, ZREF); Lh = px(MARGEN * 0.85)
        for th in angs:
            xe = x_src + Lh * np.cos(th)
            ye = Lh * np.sin(th)
            out += (f'    <path d="M {x_src:.3f},0 L {xe:.3f},{ye:.3f}" '
                    f'style="fill:none;stroke:#ff6600;stroke-width:1px"/>\n')
    else:
        # Una sola superficie refractiva: n1 → n2 en el vidrio, después de z≈0.
        # Modelamos el vidrio como medio de índice n2 llenando todo el semiespacio
        # z > z(r). El motor canónico espera un par de superficies; usamos sólo
        # una y aceptamos que el segundo segmento sea el rayo dentro del vidrio.
        sup = SuperficieCartesiana.desde_parametros_fisicos(N1, N2, Z, D0, D1)
        sis = SistemaOptico(); sis.agregar_superficie(sup)
        for th in angs:
            r = Rayo(origen=np.array([0.0, 0.0, D0]),
                     direccion=np.array([0.0, np.sin(th), np.cos(th)]))
            res = sis.trazar_rayo(r)
            pts_r = res.puntos; dirs = res.direcciones
            Pf, Df = pts_r[-1], dirs[-1]
            z_end = D1 + MARGEN * 0.85
            t_ext = (z_end - Pf[2]) / Df[2]
            Pend = Pf + t_ext * Df
            poly = list(pts_r) + [Pend]
            segs = " ".join(
                f"{'M' if i == 0 else 'L'} {xc(P[2], ZREF):.3f},{-px(P[1]):.3f}"
                for i, P in enumerate(poly)
            )
            out += (f'    <path d="{segs}" '
                    f'style="fill:none;stroke:#ff6600;stroke-width:0.7px"/>\n')

    out += svg_footer()
    with open(out_path, "w") as f:
        f.write(out)
    print(f"{out_path}  (traced={traced})  r_max={r_max:.2f} mm")


if __name__ == "__main__":
    generar_lsoe(os.path.join(AQUI, "ejemplo_lsoe.svg"),         traced=False)
    generar_lsoe(os.path.join(AQUI, "ejemplo_lsoe_traced.svg"),  traced=True)
    generar_cartesiana(os.path.join(AQUI, "ejemplo_cartesiana.svg"),        traced=False)
    generar_cartesiana(os.path.join(AQUI, "ejemplo_cartesiana_traced.svg"), traced=True)
