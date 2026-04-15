"""Ejemplos avanzados para el README.

1.  Colimador de dos lentes: fuente puntual → frente de onda esférico →
    dos LSOE estigmáticas → frente de onda plano.

2.  Experimento multielemento: combina lente LSOE, espejos, divisor de
    haz, bloque de vidrio refractivo y absorbente (beam dump) usando el
    motor canónico de inkscape-raytracing (``raytracing.World``).

Ambos escenarios se emiten como SVG listos para abrir en Inkscape (con
las etiquetas ``<desc>optics:...</desc>`` adecuadas) y en variante
«trazada» con los rayos ya propagados.
"""

import os
import sys
from math import cos, sin, pi, radians

import numpy as np

AQUI = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(AQUI, "inkscape-raytracing", "inkscape_raytracing"))
sys.path.insert(0, os.path.join(AQUI, "..", "RayTracing"))

from gots_util import (
    calcular_gots, perfil_superficie, encontrar_apertura,
    calcular_d1_sigma, puntos_a_bezier_path_str,
)
from gots.sistema_optico import SistemaOptico
from gots.superficie_cartesiana import SuperficieCartesiana
from gots.rayo import Rayo

# Motor canónico 2D
from raytracing import World, OpticalObject, Ray
from raytracing.vector import Vector, UnitVector
from raytracing.geometry import CubicBezier, CompoundGeometricObject
from raytracing.material import Mirror, BeamSplitter, BeamDump, Glass


# ──────────────────────────── utilidades SVG ────────────────────────────────

SCALE  = 2.0         # px / mm
MARGEN = 25.0        # mm


def px(v):
    return v * SCALE


def svg_header(xmin, xmax, ymin, ymax):
    w = xmax - xmin
    h = ymax - ymin
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
        f'<svg xmlns="http://www.w3.org/2000/svg" version="1.1" '
        f'width="{w:.2f}" height="{h:.2f}" '
        f'viewBox="0 0 {w:.2f} {h:.2f}">\n'
        f'  <g transform="translate({-xmin:.3f},{-ymin:.3f})">\n'
    )


def svg_footer():
    return "  </g>\n</svg>\n"


def path_str(pts, cerrar=True):
    return puntos_a_bezier_path_str(pts, cerrar=cerrar)


def linea(x1, y1, x2, y2, color="#888", sw=0.5, dash=None, op=1.0):
    da = f";stroke-dasharray:{dash}" if dash else ""
    return (
        f'    <path d="M {x1:.3f},{y1:.3f} L {x2:.3f},{y2:.3f}" '
        f'style="fill:none;stroke:{color};stroke-width:{sw}px;'
        f'opacity:{op}{da}"/>\n'
    )


def circulo(cx, cy, r, fill="#d22", stroke="none", sw=0.3):
    return (
        f'    <circle cx="{cx:.3f}" cy="{cy:.3f}" r="{r:.3f}" '
        f'style="fill:{fill};stroke:{stroke};stroke-width:{sw}px"/>\n'
    )


def texto(x, y, s, size=8, anchor="middle", color="#333"):
    return (
        f'    <text x="{x:.2f}" y="{y:.2f}" '
        f'style="font:{size}px sans-serif;fill:{color};'
        f'text-anchor:{anchor}">{s}</text>\n'
    )


# ─────────────────────── Ejemplo 1 — Colimador 2 lentes ─────────────────────

def _contorno_lsoe(p0, p1, r_ap, xshift):
    """Construye el contorno cerrado (en coords SVG) de una LSOE en Inkscape
    con el vértice frontal desplazado ``xshift`` mm hacia la derecha."""
    r0, z0 = perfil_superficie(p0, N=600, r_max=r_ap)
    r1, z1 = perfil_superficie(p1, N=600, r_max=r_ap)

    def X(z):  return px(z + xshift)
    def Yu(r): return -px(r)
    def Yl(r): return  px(r)

    pts = [(X(z0[0]), 0.0)]
    for z, r in zip(z0[1:], r0[1:]):
        pts.append((X(z), Yu(r)))
    pts.append((X(z1[-1]), Yu(r1[-1])))
    for z, r in zip(reversed(list(z1[:-1])), reversed(list(r1[:-1]))):
        pts.append((X(z), Yu(r)))
    for z, r in zip(z1[1:], r1[1:]):
        pts.append((X(z), Yl(r)))
    pts.append((X(z0[-1]), Yl(r0[-1])))
    for z, r in zip(reversed(list(z0[:-1])), reversed(list(r0[:-1]))):
        pts.append((X(z), Yl(r)))
    return pts


def generar_colimador(out_path, traced=False):
    """Fuente puntual → LSOE1 (converge a F intermedio) → LSOE2 (colima).

    El objeto de LSOE2 se coloca en F₁ (la imagen de LSOE1) y su imagen
    se envía a una distancia muy grande (≈10⁵ mm), produciendo un haz
    prácticamente colimado a la salida.
    """
    # LSOE 1 — foco corto
    N0, N1, N2 = 1.0, 1.6, 1.0
    Z0a, Z1a   = 0.0, 10.0
    D0a, D2a   = -60.0, 80.0
    SIGMA      = 0.0
    d1a        = calcular_d1_sigma(SIGMA, Z0a, Z1a, D0a, D2a, N0, N1, N2)
    p0a        = calcular_gots(N0, N1, Z0a, D0a, d1a)
    p1a        = calcular_gots(N1, N2, Z1a, d1a, D2a)
    r_apa      = encontrar_apertura(p0a, p1a) * 0.75

    # LSOE 2 — colimadora (imagen casi al infinito)
    SEP        = 80.0                 # foco intermedio coincide con d0 de lens2
    D_INF      = 1.0e5                # imagen «al infinito»
    Z0b        = Z1a + SEP            # = 90
    Z1b        = Z0b + 12.0
    D0b        = D2a                  # foco intermedio (en coordenadas globales)
    D2b        = D2a + D_INF
    d1b        = calcular_d1_sigma(SIGMA, Z0b, Z1b, D0b, D2b, N0, N1, N2)
    p0b        = calcular_gots(N0, N1, Z0b, D0b, d1b)
    p1b        = calcular_gots(N1, N2, Z1b, d1b, D2b)
    r_apb      = encontrar_apertura(p0b, p1b) * 0.75

    # Contornos (coordenadas SVG: xshift = 0, porque ya manejamos z global)
    pts_a = _contorno_lsoe(p0a, p1a, r_apa, xshift=0.0)
    pts_b = _contorno_lsoe(p0b, p1b, r_apb, xshift=0.0)

    # Encuadre
    x_src    = D0a
    x_exit   = Z1b + 40.0
    r_doc    = max(r_apa, r_apb) + 8.0
    xmin = px(x_src) - px(MARGEN * 0.4)
    xmax = px(x_exit) + px(MARGEN * 0.3)
    ymin = -px(r_doc)
    ymax =  px(r_doc)

    out = svg_header(xmin, xmax, ymin, ymax)
    # Eje óptico
    out += linea(xmin, 0, xmax, 0, color="#aaaaaa", sw=0.5, dash="6,3")

    # Lentes
    out += (f'    <path d="{path_str(pts_a)}" '
            f'style="fill:#b7c2dd;fill-opacity:0.7;stroke:#000;stroke-width:0.7px;'
            f'stroke-linejoin:round"><desc>optics:glass:{N1:.4f}</desc></path>\n')
    out += (f'    <path d="{path_str(pts_b)}" '
            f'style="fill:#b7c2dd;fill-opacity:0.7;stroke:#000;stroke-width:0.7px;'
            f'stroke-linejoin:round"><desc>optics:glass:{N1:.4f}</desc></path>\n')

    # Fuente (rojo)
    out += circulo(px(x_src), 0, px(1.6), fill="#dd2200")
    # Frente de onda plano indicativo a la derecha
    x_wave = px(x_exit - 5.0)
    out += linea(x_wave, -px(r_apb * 0.9), x_wave, px(r_apb * 0.9),
                 color="#2277cc", sw=0.8, dash="4,2", op=0.8)
    out += texto(x_wave, -px(r_apb * 0.9) - 4, "frente plano",
                 size=6, anchor="middle", color="#2277cc")
    out += texto(px(x_src), px(r_doc) - 6, "fuente puntual",
                 size=6, anchor="middle", color="#aa2200")
    out += texto(px((x_src + x_exit) / 2), -px(r_doc) + 10,
                 "LSOE 1         +         LSOE 2  →  haz colimado",
                 size=7, anchor="middle")

    # Rayos
    N_RAYOS = 11
    angs = np.linspace(-radians(8.0), radians(8.0), N_RAYOS)

    if not traced:
        for th in angs:
            xe = px(x_src) + px(40.0) * cos(th)
            ye = px(40.0) * sin(th)
            out += (f'    <path d="M {px(x_src):.3f},0 L {xe:.3f},{ye:.3f}" '
                    f'style="fill:none;stroke:#ff6600;stroke-width:1px"><desc>'
                    f'optics:beam</desc></path>\n')
    else:
        # Trazar los 4 dioptrios con el motor GOTS 3D (eje = z)
        sup0a = SuperficieCartesiana.desde_parametros_fisicos(N0, N1, Z0a, D0a, d1a)
        sup1a = SuperficieCartesiana.desde_parametros_fisicos(N1, N2, Z1a, d1a, D2a)
        sup0b = SuperficieCartesiana.desde_parametros_fisicos(N0, N1, Z0b, D0b, d1b)
        sup1b = SuperficieCartesiana.desde_parametros_fisicos(N1, N2, Z1b, d1b, D2b)
        sis = SistemaOptico()
        for s in (sup0a, sup1a, sup0b, sup1b):
            sis.agregar_superficie(s)
        for th in angs:
            r = Rayo(origen=np.array([0.0, 0.0, D0a]),
                     direccion=np.array([0.0, sin(th), cos(th)]))
            res = sis.trazar_rayo(r)
            pts_r = list(res.puntos)
            Pf, Df = res.puntos[-1], res.direcciones[-1]
            z_end = x_exit
            t_ext = (z_end - Pf[2]) / Df[2]
            Pend  = Pf + t_ext * Df
            pts_r.append(Pend)
            segs = " ".join(
                f"{'M' if i == 0 else 'L'} {px(P[2]):.3f},{-px(P[1]):.3f}"
                for i, P in enumerate(pts_r)
            )
            out += (f'    <path d="{segs}" '
                    f'style="fill:none;stroke:#ff6600;stroke-width:0.7px"/>\n')

    out += svg_footer()
    with open(out_path, "w") as f:
        f.write(out)
    print(f"{out_path}  (traced={traced})")


# ────────────────── Ejemplo 2 — Experimento multielemento ──────────────────
# Usa el motor canónico raytracing.World directamente.

def _segmento_bezier(p0, p1):
    """Representa el segmento recto P0-P1 como CubicBezier."""
    q0 = Vector(p0[0], p0[1])
    q3 = Vector(p1[0], p1[1])
    d  = q3 - q0
    q1 = q0 + Vector(d.x / 3.0, d.y / 3.0)
    q2 = q0 + Vector(2.0 * d.x / 3.0, 2.0 * d.y / 3.0)
    return CubicBezier(q0, q1, q2, q3)


def _oval_lsoe_bezier(p_front, p_back, r_ap, xshift, yshift=0.0):
    """Devuelve CompoundGeometricObject + lista de puntos SVG para dibujo."""
    r0, z0 = perfil_superficie(p_front, N=200, r_max=r_ap)
    r1, z1 = perfil_superficie(p_back,  N=200, r_max=r_ap)

    def pt(z, y): return (z + xshift, y + yshift)

    # Muestreamos el contorno completo (igual que _contorno_lsoe) y
    # unimos puntos consecutivos con CubicBezier rectos.
    contorno = [pt(z0[0], 0.0)]
    for z, r in zip(z0[1:], r0[1:]):
        contorno.append(pt(z, -r))
    contorno.append(pt(z1[-1], -r1[-1]))
    for z, r in zip(reversed(list(z1[:-1])), reversed(list(r1[:-1]))):
        contorno.append(pt(z, -r))
    for z, r in zip(z1[1:], r1[1:]):
        contorno.append(pt(z, r))
    contorno.append(pt(z0[-1], r0[-1]))
    for z, r in zip(reversed(list(z0[:-1])), reversed(list(r0[:-1]))):
        contorno.append(pt(z, r))

    beziers = [
        _segmento_bezier(contorno[i], contorno[(i + 1) % len(contorno)])
        for i in range(len(contorno))
    ]
    return CompoundGeometricObject(beziers), contorno


def generar_experimento(out_path, traced=False):
    """Experimento de un solo paso con todos los elementos de la extensión.

    Layout (convención SVG: y crece hacia abajo):

        fuente ─► LSOE colimadora ─► BS (45°) ─┬─► brazo transmitido
                                               │    (vidrio + beam dump)
                                               └─► brazo reflejado
                                                    (espejo pliegue + beam dump)

    Cada rayo atraviesa cada elemento a lo sumo una vez, evitando la
    recursión exponencial de un Michelson con recombinación.
    """
    # ── Lente colimadora (LSOE σ=0) ──
    N0, N1, N2 = 1.0, 1.5, 1.0
    Zf, Zb     = 0.0, 6.0
    D0,  D2    = -40.0, 1.0e5
    SIGMA      = 0.0
    d1         = calcular_d1_sigma(SIGMA, Zf, Zb, D0, D2, N0, N1, N2)
    pf         = calcular_gots(N0, N1, Zf, D0, d1)
    pb         = calcular_gots(N1, N2, Zb, d1, D2)
    r_ap       = encontrar_apertura(pf, pb) * 0.55

    lens_geo, lens_contorno = _oval_lsoe_bezier(pf, pb, r_ap, xshift=0.0)

    # ── BS a 45° en (X_BS, 0), pendiente +1, tangente (1,1)/√2 ──
    # Un rayo horizontal (+x) se refleja hacia +y (abajo en SVG).
    X_BS, L_BS = 70.0, 16.0
    d_bs = L_BS / 2.0 / np.sqrt(2.0)
    BS_A = (X_BS - d_bs, -d_bs)
    BS_B = (X_BS + d_bs, +d_bs)
    bs_geo = _segmento_bezier(BS_A, BS_B)

    # ── Bloque de vidrio en el brazo transmitido ──
    GX0, GX1 = 100.0, 120.0
    GY0, GY1 = -10.0, +10.0
    glass_corners = [(GX0, GY0), (GX1, GY0), (GX1, GY1), (GX0, GY1)]
    glass_geo = CompoundGeometricObject([
        _segmento_bezier(glass_corners[i],
                         glass_corners[(i + 1) % 4]) for i in range(4)
    ])

    # ── Beam dump transmitido (placa vertical al extremo derecho) ──
    X_DUMP_T = 160.0
    DT_A = (X_DUMP_T, -12.0)
    DT_B = (X_DUMP_T, +12.0)
    dump_t_geo = _segmento_bezier(DT_A, DT_B)

    # ── Espejo pliegue en el brazo reflejado (redirige +y → +x) ──
    # El rayo viene desde (X_BS, 0) en dirección (0, +1).  Lo cruza un
    # espejo con tangente (1, 1)/√2 en (X_BS, 40).  Reflexión:
    #     d = (0,1), n = (1,-1)/√2, d·n = −1/√2
    #     r = d − 2(d·n)n = (0,1) + (1,−1) = (1, 0) ✓
    Y_FOLD = 45.0
    L_FOLD = 14.0
    df = L_FOLD / 2.0 / np.sqrt(2.0)
    FOLD_A = (X_BS - df, Y_FOLD - df)
    FOLD_B = (X_BS + df, Y_FOLD + df)
    fold_geo = _segmento_bezier(FOLD_A, FOLD_B)

    # ── Beam dump reflejado (placa vertical a la derecha del pliegue) ──
    X_DUMP_R = 160.0
    DR_A = (X_DUMP_R, Y_FOLD - 12.0)
    DR_B = (X_DUMP_R, Y_FOLD + 12.0)
    dump_r_geo = _segmento_bezier(DR_A, DR_B)

    # ── Mundo ──
    world = World()
    world.add(OpticalObject(lens_geo,   Glass(N1)))
    world.add(OpticalObject(bs_geo,     BeamSplitter()))
    world.add(OpticalObject(glass_geo,  Glass(1.5)))
    world.add(OpticalObject(fold_geo,   Mirror()))
    world.add(OpticalObject(dump_t_geo, BeamDump()))
    world.add(OpticalObject(dump_r_geo, BeamDump()))

    # ── Rayos semilla ──
    N_RAYOS = 7
    angs = np.linspace(-radians(5.0), radians(5.0), N_RAYOS)
    fuente = Vector(D0, 0.0)

    # ── Encuadre ──
    xmin = px(D0 - 8.0);   xmax = px(X_DUMP_T + 10.0)
    ymin = px(-r_ap - 10.0); ymax = px(Y_FOLD + 20.0)

    out = svg_header(xmin, xmax, ymin, ymax)

    # Lente
    pts_svg_lens = [(px(x), px(y)) for (x, y) in lens_contorno]
    out += (f'    <path d="{path_str(pts_svg_lens)}" '
            f'style="fill:#b7c2dd;fill-opacity:0.7;stroke:#000;stroke-width:0.6px;'
            f'stroke-linejoin:round"><desc>optics:glass:{N1:.4f}</desc></path>\n')
    out += texto(px(Zf + 3.0), px(-r_ap - 4.0), "LSOE colimadora", size=6)

    # Beam splitter
    out += (f'    <path d="M {px(BS_A[0]):.3f},{px(BS_A[1]):.3f} '
            f'L {px(BS_B[0]):.3f},{px(BS_B[1]):.3f}" '
            f'style="fill:none;stroke:#4488aa;stroke-width:1.3px">'
            f'<desc>optics:beam_splitter</desc></path>\n')
    out += texto(px(X_BS + 8), px(-8), "BS 45°", size=6, anchor="start", color="#4488aa")

    # Vidrio rectangular
    gpts = [(px(x), px(y)) for (x, y) in glass_corners]
    out += (f'    <path d="M {gpts[0][0]:.3f},{gpts[0][1]:.3f} '
            f'L {gpts[1][0]:.3f},{gpts[1][1]:.3f} '
            f'L {gpts[2][0]:.3f},{gpts[2][1]:.3f} '
            f'L {gpts[3][0]:.3f},{gpts[3][1]:.3f} Z" '
            f'style="fill:#c8e0f0;fill-opacity:0.5;stroke:#3355aa;stroke-width:0.6px">'
            f'<desc>optics:glass:1.5</desc></path>\n')
    out += texto(px((GX0 + GX1) / 2), px(GY0 - 3), "vidrio n=1.5",
                 size=5, color="#3355aa")

    # Espejo pliegue
    out += (f'    <path d="M {px(FOLD_A[0]):.3f},{px(FOLD_A[1]):.3f} '
            f'L {px(FOLD_B[0]):.3f},{px(FOLD_B[1]):.3f}" '
            f'style="fill:none;stroke:#aa3377;stroke-width:1.4px">'
            f'<desc>optics:mirror</desc></path>\n')
    out += texto(px(X_BS - 10), px(Y_FOLD + 10), "M (pliegue)",
                 size=6, color="#aa3377", anchor="end")

    # Beam dumps
    for (A, B, label, tx, ty) in [
        (DT_A, DT_B, "dump T", px(X_DUMP_T + 3), px(0)),
        (DR_A, DR_B, "dump R", px(X_DUMP_R + 3), px(Y_FOLD)),
    ]:
        out += (f'    <path d="M {px(A[0]):.3f},{px(A[1]):.3f} '
                f'L {px(B[0]):.3f},{px(B[1]):.3f}" '
                f'style="fill:none;stroke:#222;stroke-width:2.0px">'
                f'<desc>optics:beam_dump</desc></path>\n')
        out += texto(tx, ty, label, size=5, color="#222", anchor="start")

    # Fuente
    out += circulo(px(D0), 0, px(1.5), fill="#dd2200")
    out += texto(px(D0), px(-5.0), "fuente", size=6, color="#aa2200")

    if not traced:
        for th in angs:
            xe = px(D0) + px(22.0) * cos(th)
            ye = px(22.0) * sin(th)
            out += (f'    <path d="M {px(D0):.3f},0 L {xe:.3f},{ye:.3f}" '
                    f'style="fill:none;stroke:#ff6600;stroke-width:1px">'
                    f'<desc>optics:beam</desc></path>\n')
    else:
        # Limitar recursión para proteger frente a lazos no previstos.
        world.max_recursion_depth = 40
        for th in angs:
            seed = Ray(fuente, UnitVector(cos(th), sin(th)))
            beams = world.propagate_beams(seed)
            for beam in beams:
                if not beam:
                    continue
                p0 = beam[0].origin
                d_path = f"M {px(p0.x):.3f},{px(p0.y):.3f}"
                for ray in beam:
                    travel = ray.travel if ray.travel > 0 else 30.0
                    p1 = ray.origin + travel * ray.direction
                    # Evita coordenadas absurdas (inf/nan)
                    if not (np.isfinite(p1.x) and np.isfinite(p1.y)):
                        break
                    d_path += f" L {px(p1.x):.3f},{px(p1.y):.3f}"
                out += (f'    <path d="{d_path}" '
                        f'style="fill:none;stroke:#ff6600;stroke-width:0.7px;'
                        f'opacity:0.85"/>\n')

    out += svg_footer()
    with open(out_path, "w") as f:
        f.write(out)
    print(f"{out_path}  (traced={traced})")


# ──────────────────────────────── main ──────────────────────────────────────

if __name__ == "__main__":
    generar_colimador(os.path.join(AQUI, "ejemplo_colimador.svg"),         traced=False)
    generar_colimador(os.path.join(AQUI, "ejemplo_colimador_traced.svg"),  traced=True)
    generar_experimento(os.path.join(AQUI, "ejemplo_experimento.svg"),         traced=False)
    generar_experimento(os.path.join(AQUI, "ejemplo_experimento_traced.svg"),  traced=True)
