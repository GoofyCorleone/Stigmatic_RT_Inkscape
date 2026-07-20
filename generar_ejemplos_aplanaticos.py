"""Ejemplos de lentes rigurosamente aplanáticas (cap. 4 de Silva-Lora 2024).

Reproduce las cuatro figuras de sistemas aplanáticos en todo rigor de la
tesis (Figs. 21, 24, 27 y 30): para cada tipo se traza el haz procedente
de dos puntos objeto, uno sobre el eje (H = 0) y otro fuera de eje
(H = 8 mm), usando el motor canónico ``raytracing.World`` sobre la MISMA
geometría Bézier que la extensión escribe en el SVG.

Al carecer de aberración esférica y de coma, ambos abanicos convergen en
un punto nítido; cuando la imagen es virtual se dibuja además la
prolongación de los rayos emergentes (línea punteada) hasta el punto
imagen, tal como hace el autor en la tesis.

Salida: ejemplo_aplanatica_tipo{0,1,2,3}.svg
"""

import os
import sys
from math import cos, sin, radians

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(AQUI, "inkscape-raytracing", "inkscape_raytracing"))

from gots_util import (
    disenar_aplanatica,
    perfil_superficie_aplanatica,
    puntos_a_bezier_path_str,
)

from raytracing import World, OpticalObject, Ray
from raytracing.vector import Vector, UnitVector
from raytracing.geometry import CubicBezier, CompoundGeometricObject
from raytracing.material import Glass, BeamDump


# ──────────────────────────── utilidades SVG ────────────────────────────────

SCALE = 2.0          # px / mm


def px(v):
    return v * SCALE


def svg_header(xmin, xmax, ymin, ymax):
    w, h = xmax - xmin, ymax - ymin
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" '
        f'version="1.1" width="{w:.2f}" height="{h:.2f}" '
        f'viewBox="0 0 {w:.2f} {h:.2f}">\n'
        f'  <g inkscape:groupmode="layer" inkscape:label="Capa 1" '
        f'id="layer1" transform="translate({-xmin:.3f},{-ymin:.3f})">\n'
    )


def svg_footer():
    return "  </g>\n</svg>\n"


def linea(x1, y1, x2, y2, color="#888", sw=0.5, dash=None, op=1.0):
    da = f";stroke-dasharray:{dash}" if dash else ""
    return (f'    <path d="M {x1:.3f},{y1:.3f} L {x2:.3f},{y2:.3f}" '
            f'style="fill:none;stroke:{color};stroke-width:{sw}px;'
            f'opacity:{op}{da}"/>\n')


def circulo(cx, cy, r, fill="#d22", stroke="none", sw=0.3, dash=None):
    da = f";stroke-dasharray:{dash}" if dash else ""
    return (f'    <circle cx="{cx:.3f}" cy="{cy:.3f}" r="{r:.3f}" '
            f'style="fill:{fill};stroke:{stroke};stroke-width:{sw}px{da}"/>\n')


def texto(x, y, s, size=7, anchor="middle", color="#333"):
    return (f'    <text x="{x:.2f}" y="{y:.2f}" '
            f'style="font:{size}px sans-serif;fill:{color};'
            f'text-anchor:{anchor}">{s}</text>\n')


# ────────────────────── geometría de la lente aplanática ────────────────────

def _contorno_aplanatica(diseno, r_ap, N=320):
    """Contorno cerrado (mm, y hacia abajo) con ζ₀ en el origen."""
    sup0, sup1 = diseno['superficies']
    zeta_ref = diseno['zetas'][0]
    r0, z0 = perfil_superficie_aplanatica(sup0, r_ap, N=N)
    r1, z1 = perfil_superficie_aplanatica(sup1, r_ap, N=N)
    z0 = z0 - zeta_ref
    z1 = z1 - zeta_ref

    pts = [(float(z0[0]), 0.0)]
    for z, r in zip(z0[1:], r0[1:]):
        pts.append((float(z), -float(r)))
    pts.append((float(z1[-1]), -float(r1[-1])))
    for z, r in zip(reversed(list(z1[:-1])), reversed(list(r1[:-1]))):
        pts.append((float(z), -float(r)))
    for z, r in zip(z1[1:], r1[1:]):
        pts.append((float(z), float(r)))
    pts.append((float(z0[-1]), float(r0[-1])))
    for z, r in zip(reversed(list(z0[:-1])), reversed(list(r0[:-1]))):
        pts.append((float(z), float(r)))
    return pts


def _apertura_maxima(diseno):
    """Radio máximo admisible: dominio de cada superficie, escala física
    del sistema y punto donde ambos perfiles se cruzarían."""
    topes = []
    for sup in diseno['superficies']:
        if sup.get('esfera'):
            topes.append(0.98 * abs(sup['R']))
        else:
            coef_rad = 2.0 * sup['S'] - sup['O'] * sup['OG']
            if coef_rad < -1e-15:
                topes.append(0.9 * float(np.sqrt(-1.0 / coef_rad)))
    # Escala física: no tiene sentido una apertura comparable al conjugado
    topes.append(0.35 * abs(diseno['d'][0] - diseno['zetas'][0]))
    r_tope = min(topes)

    # Recorte si las dos superficies se cruzan antes (lente sin espesor)
    sup0, sup1 = diseno['superficies']
    r0, z0 = perfil_superficie_aplanatica(sup0, r_tope, N=300)
    r1, z1 = perfil_superficie_aplanatica(sup1, r_tope, N=300)
    if len(r0) > 2 and len(r1) > 2:
        rp = np.linspace(r_tope * 1e-3, min(r0[-1], r1[-1]), 400)
        dz = np.interp(rp, r1, z1) - np.interp(rp, r0, z0)
        cruces = np.where(np.diff(np.sign(dz)) != 0)[0]
        if len(cruces) > 0:
            r_tope = min(r_tope, 0.92 * float(rp[cruces[0]]))
    return r_tope


def _path_a_beziers(d_str):
    """Convierte la cadena 'M .. C ..' en CubicBezier del motor canónico.

    Se traza exactamente la misma curva que se dibuja en el SVG, de modo
    que el trazado mostrado coincide con lo que vería la extensión.
    """
    toks = d_str.replace(",", " ").split()
    beziers, i, cur, inicio = [], 0, None, None
    while i < len(toks):
        if toks[i] == "M":
            cur = (float(toks[i + 1]), float(toks[i + 2]))
            inicio = cur
            i += 3
        elif toks[i] == "C":
            p1 = (float(toks[i + 1]), float(toks[i + 2]))
            p2 = (float(toks[i + 3]), float(toks[i + 4]))
            p3 = (float(toks[i + 5]), float(toks[i + 6]))
            beziers.append(CubicBezier(Vector(*cur), Vector(*p1),
                                       Vector(*p2), Vector(*p3)))
            cur = p3
            i += 7
        elif toks[i] in ("Z", "z"):
            if cur != inicio:
                d = (inicio[0] - cur[0], inicio[1] - cur[1])
                beziers.append(CubicBezier(
                    Vector(*cur),
                    Vector(cur[0] + d[0] / 3.0, cur[1] + d[1] / 3.0),
                    Vector(cur[0] + 2 * d[0] / 3.0, cur[1] + 2 * d[1] / 3.0),
                    Vector(*inicio)))
            i += 1
        else:
            i += 1
    return CompoundGeometricObject(beziers)


# ───────────────────────────── generación ───────────────────────────────────

DESCRIPCION = {
    0: ("Tipo-0 — esferas en los puntos de Young",
        "libre de esférica, coma y astigmatismo"),
    1: ("Tipo-1 — cónicas, interior colimado",
        "libre de esférica y coma"),
    2: ("Tipo-2 — ovoides de vértice plano (O = 0)",
        "libre de esférica y coma"),
    3: ("Tipo-3 — ovoides con G = 0",
        "libre de esférica y coma"),
}


def _semillas(z_obj, z_v1, H, n_rayos, semiangulo):
    """Rayos semilla de los dos abanicos (en eje y fuera de eje)."""
    for color, altura in [("#1144cc", 0.0), ("#dd8800", H)]:
        theta_c = np.arctan2(-altura, z_v1 * 0.5 - z_obj)
        for th in np.linspace(-radians(semiangulo), radians(semiangulo),
                              n_rayos):
            ang = theta_c + th
            yield color, altura, Ray(Vector(z_obj, altura),
                                     UnitVector(cos(ang), sin(ang)))


def _mundo(diseno, r_ap, x_fin, n_vidrio):
    """Mundo de trazado sobre la MISMA curva Bézier que se dibuja."""
    contorno = _contorno_aplanatica(diseno, r_ap)
    geo = _path_a_beziers(puntos_a_bezier_path_str(contorno, cerrar=True))
    world = World()
    world.add(OpticalObject(geo, Glass(n_vidrio)))
    world.add(OpticalObject(
        CompoundGeometricObject([CubicBezier(
            Vector(x_fin, -60.0), Vector(x_fin, -20.0),
            Vector(x_fin, 20.0), Vector(x_fin, 60.0))]),
        BeamDump()))
    world.max_recursion_depth = 20
    return world, contorno


def _radio_necesario(diseno, r_prueba, x_fin, n_vidrio, semillas):
    """Altura máxima que alcanzan los rayos sobre las dos superficies.

    Se traza con una apertura de prueba holgada para que ningún rayo
    tropiece con el borde; así se mide la apertura que realmente exige
    el haz (dentro del vidrio el cono se abre o se cierra según el tipo).
    """
    world, _ = _mundo(diseno, r_prueba, x_fin, n_vidrio)
    r_max = 0.0
    for _, _, seed in semillas:
        for beam in world.propagate_beams(seed):
            for rayo in beam[:-1]:      # descarta el tramo hacia la pantalla
                if rayo.travel > 0:
                    p = rayo.origin + rayo.travel * rayo.direction
                    if np.isfinite(p.y):
                        r_max = max(r_max, abs(p.y))
    return r_max


def generar_aplanatica(out_path, tipo, rama="real", r_ap=None,
                       H=8.0, n_rayos=11, semiangulo=7.0):
    """Figura de una lente rigurosamente aplanática con dos puntos objeto."""
    N0, N1 = 1.0, 1.8
    ZETA_0, ZETA_1, D_0 = 60.0, 70.0, 0.0

    diseno = disenar_aplanatica(tipo, N0, N1, ZETA_0, ZETA_1, D_0, rama=rama)
    d_0, d_1, d_2 = diseno['d']
    M_lat = diseno['magnificacion']

    # Coordenadas locales: vértice frontal en el origen
    z_obj = d_0 - ZETA_0
    z_img = d_2 - ZETA_0
    z_v1  = ZETA_1 - ZETA_0
    img_virtual = d_2 < ZETA_1

    # Pantalla absorbente al final del documento para acotar los rayos
    x_fin = max(z_img, z_v1) + 55.0
    x_ini = min(z_obj, z_img) - 12.0

    # ── Apertura: la exigida por el haz, no un valor fijo ────────────────
    r_tope = _apertura_maxima(diseno)
    if r_ap is None:
        semillas = list(_semillas(z_obj, z_v1, H, n_rayos, semiangulo))
        r_nec = _radio_necesario(diseno, 0.95 * r_tope, x_fin, N1, semillas)
        r_ap  = min(1.06 * r_nec, 0.95 * r_tope)

    world, contorno = _mundo(diseno, r_ap, x_fin, N1)

    # ── Encuadre ─────────────────────────────────────────────────────────
    r_doc = r_ap + max(abs(H), abs(M_lat * H)) + 6.0
    xmin, xmax = px(x_ini), px(x_fin + 4.0)
    ymin, ymax = -px(r_doc), px(r_doc + 6.0)

    out  = svg_header(xmin, xmax, ymin, ymax)
    out += linea(xmin, 0, xmax, 0, color="#aaaaaa", sw=0.4, dash="6,3")

    # ── Lente ────────────────────────────────────────────────────────────
    pts_svg = [(px(x), px(y)) for x, y in contorno]
    out += (f'    <path d="{puntos_a_bezier_path_str(pts_svg, cerrar=True)}" '
            f'style="fill:#c2ddb7;fill-opacity:0.7;stroke:#000;'
            f'stroke-width:0.6px;stroke-linejoin:round">'
            f'<desc>optics:glass:{N1:.4f}</desc></path>\n')

    # ── Trazado de los dos abanicos ──────────────────────────────────────
    for color, altura, seed in _semillas(z_obj, z_v1, H, n_rayos, semiangulo):
        for beam in world.propagate_beams(seed):
            if not beam:
                continue
            p0 = beam[0].origin
            d_path = f"M {px(p0.x):.3f},{px(p0.y):.3f}"
            for rayo in beam:
                travel = rayo.travel if rayo.travel > 0 else 40.0
                p1 = rayo.origin + travel * rayo.direction
                if not (np.isfinite(p1.x) and np.isfinite(p1.y)):
                    break
                d_path += f" L {px(p1.x):.3f},{px(p1.y):.3f}"
            out += (f'    <path d="{d_path}" style="fill:none;'
                    f'stroke:{color};stroke-width:0.55px;'
                    f'opacity:0.85"/>\n')

            # Prolongación hacia atrás hasta la imagen virtual
            if img_virtual and len(beam) >= 2:
                ult = beam[-1]
                a, v = ult.origin, ult.direction
                t = (z_img - a.x) / v.x if abs(v.x) > 1e-12 else 0.0
                if t < 0:
                    out += (f'    <path d="M {px(a.x):.3f},{px(a.y):.3f} '
                            f'L {px(a.x + t * v.x):.3f},'
                            f'{px(a.y + t * v.y):.3f}" '
                            f'style="fill:none;stroke:{color};'
                            f'stroke-width:0.35px;opacity:0.45;'
                            f'stroke-dasharray:2,2"/>\n')

    # ── Puntos objeto e imagen ───────────────────────────────────────────
    rad = px(1.3)
    for altura, color in [(0.0, "#1144cc"), (H, "#dd8800")]:
        out += circulo(px(z_obj), px(altura), rad, fill=color)
    for altura, color in [(0.0, "#1144cc"), (M_lat * H, "#dd8800")]:
        if img_virtual:
            out += circulo(px(z_img), px(altura), rad, fill="none",
                           stroke=color, sw=0.8, dash="1.5,1")
        else:
            out += circulo(px(z_img), px(altura), rad, fill=color)

    # ── Rótulos ──────────────────────────────────────────────────────────
    titulo, propiedad = DESCRIPCION[tipo]
    if tipo == 1:
        titulo += f"  ({'imagen real' if rama == 'real' else 'imagen virtual'})"
    out += texto((xmin + xmax) / 2, ymin + px(5.0), titulo, size=8)
    out += texto((xmin + xmax) / 2, ymin + px(10.0),
                 f"{propiedad}   ·   n₁ = {N1}   ·   aumento M = {M_lat:+.3f}",
                 size=6, color="#666")
    out += texto(px(z_obj), px(-1.5) - px(1.0), "objeto", size=6,
                 color="#1144cc")
    etiqueta_img = "imagen virtual" if img_virtual else "imagen"
    out += texto(px(z_img), ymax - px(3.0), etiqueta_img, size=6,
                 color="#008833")
    out += texto(px(z_v1 / 2.0), -px(r_ap + 2.0), "lente", size=6)

    out += svg_footer()
    with open(out_path, "w") as f:
        f.write(out)

    return diseno


def main():
    casos = [
        ("ejemplo_aplanatica_tipo0.svg", 0, "real"),
        ("ejemplo_aplanatica_tipo1.svg", 1, "real"),
        ("ejemplo_aplanatica_tipo2.svg", 2, "real"),
        ("ejemplo_aplanatica_tipo3.svg", 3, "real"),
    ]
    for nombre, tipo, rama in casos:
        ruta = os.path.join(AQUI, nombre)
        dis  = generar_aplanatica(ruta, tipo, rama=rama)
        d0, d1, d2 = dis['d']
        d1_txt = "∞" if np.isinf(d1) else f"{d1:.6f}"
        print(f"{ruta}\n    tipo-{tipo}: d₀={d0:.3f}  d₁={d1_txt}  "
              f"d₂={d2:.6f}  M={dis['magnificacion']:+.4f}")


if __name__ == "__main__":
    main()
