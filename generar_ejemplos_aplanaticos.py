"""Ejemplos de lentes rigurosamente aplanáticas (cap. 4 de Silva-Lora 2024).

Para cada uno de los cuatro tipos se trazan CUATRO fuentes puntuales
situadas a distintas alturas del plano objeto, limitadas por una pupila
física, y se localiza el foco meridional REAL de cada abanico.

Lo que hay que mirar en las figuras no es que los focos caigan sobre un
plano —no lo hacen—, sino que:

  · cada abanico converge en un punto nítido (spot de micras): eso es el
    aplanatismo, la ausencia simultánea de aberración esférica y coma,
    garantizada por la condición seno de Abbe que se anota en la figura;
  · el lugar geométrico de esos focos es una SUPERFICIE CURVA (curvatura
    de campo), dibujada punteada en morado y claramente separada del
    plano imagen paraxial, que se dibuja como línea gris de referencia.

La tesis muestra exactamente esto en la Fig. 17 («cuya imagen se forma
sobre una superficie curva») y lo discute en la p. 101.

El recuadro inferior amplía la zona de la imagen en el eje z para hacer
visible esa curvatura, que a escala real es de décimas de milímetro.

Salida: ejemplo_aplanatica_tipo{0,1,2,3}.svg
"""

import os
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(AQUI, "inkscape-raytracing", "inkscape_raytracing"))

from gots_util import (
    disenar_aplanatica,
    perfil_superficie_aplanatica,
    puntos_a_bezier_path_str,
    trazar_rayo,
    foco_meridional,
    superficie_imagen,
    seno_abbe,
    direcciones_pupila,
    condicion_aplanatismo,
    rms_aplanatismo,
)


# ──────────────────────────── utilidades SVG ────────────────────────────────

SCALE = 2.4          # px / mm

COLORES = ["#1144cc", "#00913f", "#dd8800", "#cc2222"]


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
        f'  <rect x="0" y="0" width="{w:.2f}" height="{h:.2f}" '
        f'fill="#ffffff"/>\n'
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


def polilinea(pts, color="#888", sw=0.5, dash=None, op=1.0):
    if len(pts) < 2:
        return ""
    da = f";stroke-dasharray:{dash}" if dash else ""
    d = "M " + " L ".join(f"{x:.3f},{y:.3f}" for x, y in pts)
    return (f'    <path d="{d}" style="fill:none;stroke:{color};'
            f'stroke-width:{sw}px;opacity:{op}{da}"/>\n')


def circulo(cx, cy, r, fill="#d22", stroke="none", sw=0.3):
    return (f'    <circle cx="{cx:.3f}" cy="{cy:.3f}" r="{r:.3f}" '
            f'style="fill:{fill};stroke:{stroke};stroke-width:{sw}px"/>\n')


def texto(x, y, s, size=7, anchor="middle", color="#333", peso="normal"):
    return (f'    <text x="{x:.2f}" y="{y:.2f}" '
            f'style="font:{peso} {size}px sans-serif;fill:{color};'
            f'text-anchor:{anchor}">{s}</text>\n')


# ────────────────────── geometría de la lente ───────────────────────────────

def _contorno(diseno, r_ap, N=320):
    """Contorno cerrado (mm, y hacia abajo) en coordenadas de diseño."""
    sup0, sup1 = diseno['superficies']
    r0, z0 = perfil_superficie_aplanatica(sup0, r_ap, N=N)
    r1, z1 = perfil_superficie_aplanatica(sup1, r_ap, N=N)

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
    """Radio máximo admisible: dominio, escala física y cruce de perfiles."""
    topes = []
    for sup in diseno['superficies']:
        if sup.get('esfera'):
            topes.append(0.95 * abs(sup['R']))
        else:
            coef_rad = 2.0 * sup['S'] - sup['O'] * sup['OG']
            if coef_rad < -1e-15:
                topes.append(0.9 * float(np.sqrt(-1.0 / coef_rad)))
    topes.append(0.35 * abs(diseno['d'][0] - diseno['zetas'][0]))
    r_tope = min(topes)

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


def _abanicos(superficies, d_obj, alturas, r_ap, z_pup, r_pup, n_rayos):
    """Traza los abanicos de todas las fuentes y localiza sus focos."""
    salida = []
    for h in alturas:
        origen = (d_obj, float(h))
        trayectorias, orig_sal, dir_sal = [], [], []
        for d0 in direcciones_pupila(d_obj, h, z_pup, r_pup, n_rayos):
            pts, dirs, _ = trazar_rayo(superficies, origen, d0, r_ap)
            if pts is None:
                continue
            trayectorias.append((pts, dirs))
            orig_sal.append(pts[-1])
            dir_sal.append(dirs[-1])
        if len(orig_sal) < 3:
            salida.append({'h': h, 'trayectorias': trayectorias,
                           'foco': None, 'spot': None})
            continue
        foco, spot = foco_meridional(orig_sal, dir_sal)
        salida.append({'h': h, 'trayectorias': trayectorias,
                       'foco': foco, 'spot': spot})
    return salida


def _radio_necesario(superficies, d_obj, alturas, r_prueba,
                     z_pup, r_pup, n_rayos):
    """Altura máxima que alcanzan los rayos sobre las superficies."""
    r_max = 0.0
    for h in alturas:
        for d0 in direcciones_pupila(d_obj, h, z_pup, r_pup, n_rayos):
            pts, _, _ = trazar_rayo(superficies, (d_obj, float(h)),
                                    d0, r_prueba)
            if pts is None:
                continue
            for _, r in pts[1:]:
                r_max = max(r_max, abs(r))
    return r_max


# ───────────────────────────── generación ───────────────────────────────────

DESCRIPCION = {
    0: ("Tipo-0 · esferas en los puntos de Young",
        "2Sk = Gk·Ok² — libre de esférica, coma y astigmatismo"),
    1: ("Tipo-1 · cónicas con interior colimado",
        "d1 → infinito — libre de esférica y coma"),
    2: ("Tipo-2 · ovoides de vértice plano",
        "Ok = 0 — libre de esférica y coma"),
    3: ("Tipo-3 · ovoides con G = 0",
        "Gk = 0 — libre de esférica y coma"),
}


def generar_aplanatica(out_path, tipo, rama="real",
                       alturas=(0.0, 3.0, 6.0, 9.0),
                       n_rayos=9, r_pupila=8.0):
    N0, N1 = 1.0, 1.8
    ZETA_0, ZETA_1, D_0 = 60.0, 70.0, 0.0
    Z_PUPILA = 55.0                     # pupila física, como en la tesis

    diseno = disenar_aplanatica(tipo, N0, N1, ZETA_0, ZETA_1, D_0, rama=rama)
    sups = diseno['superficies']
    d_0, d_1, d_2 = diseno['d']
    M_lat = diseno['magnificacion']

    # ── Apertura exigida por el haz que atraviesa la pupila ──────────────
    r_tope = _apertura_maxima(diseno)
    r_nec = _radio_necesario(sups, d_0, alturas, 0.95 * r_tope,
                             Z_PUPILA, r_pupila, n_rayos)
    r_ap = min(1.05 * r_nec, 0.95 * r_tope)

    # ── Trazado, superficie imagen y condición seno ──────────────────────
    fans = _abanicos(sups, d_0, alturas, r_ap, Z_PUPILA, r_pupila, n_rayos)
    curva = superficie_imagen(sups, d_0, np.linspace(0, max(alturas), 13),
                              r_ap, n_rayos=n_rayos,
                              z_pupila=Z_PUPILA, r_pupila=r_pupila)
    abbe = seno_abbe(sups, d_0, r_ap,
                     semiangulo_deg=np.degrees(
                         np.arctan(r_pupila / (Z_PUPILA - d_0))))

    focos = [f['foco'] for f in fans if f['foco'] is not None]
    dz_max = max(abs(e['foco'][0] - d_2) for e in curva) or 1e-6

    # ── Encuadre del panel principal ─────────────────────────────────────
    z_int = [d_0, ZETA_0, ZETA_1, d_2] + [f[0] for f in focos]
    x_ini = min(z_int) - 16.0
    x_fin = max(max(z_int), ZETA_1) + 42.0
    r_doc = max(r_ap, max(alturas), max(abs(f[1]) for f in focos)) + 8.0

    ALTO_TXT = 26.0                      # banda de cifras bajo el trazado
    xmin, xmax = px(x_ini), px(x_fin)
    ymin = -px(r_doc) - px(13.0)
    ymax = px(r_doc) + px(ALTO_TXT + 28.0)

    out = svg_header(xmin, xmax, ymin, ymax)
    out += linea(xmin, 0, xmax, 0, color="#bbbbbb", sw=0.4, dash="6,3")

    # ── Lente ────────────────────────────────────────────────────────────
    contorno = _contorno(diseno, r_ap)
    pts_svg = [(px(x), px(y)) for x, y in contorno]
    out += (f'    <path d="{puntos_a_bezier_path_str(pts_svg, cerrar=True)}" '
            f'style="fill:#c2ddb7;fill-opacity:0.75;stroke:#33691e;'
            f'stroke-width:0.7px;stroke-linejoin:round">'
            f'<desc>optics:glass:{N1:.4f}</desc></path>\n')

    # ── Pupila física ────────────────────────────────────────────────────
    # 2 % de holgura: si el borde coincide con el rayo marginal, el
    # trazador lo absorbe por ambigüedad numérica.
    for lado in (-1.0, +1.0):
        out += (f'    <path d="M {px(Z_PUPILA):.3f},'
                f'{lado * px(r_pupila * 1.02):.3f} '
                f'L {px(Z_PUPILA):.3f},{lado * px(r_doc * 0.92):.3f}" '
                f'style="fill:none;stroke:#333;stroke-width:1.6px">'
                f'<desc>optics:beam_dump</desc></path>\n')
    out += texto(px(Z_PUPILA), -px(r_doc * 0.92) - px(1.8), "pupila",
                 size=5.5, color="#333")

    # ── Plano imagen paraxial (referencia) ───────────────────────────────
    out += linea(px(d_2), -px(r_doc * 0.78), px(d_2), px(r_doc * 0.78),
                 color="#999999", sw=0.7, dash="3,2")
    out += texto(px(d_2), px(r_doc * 0.78) + px(3.5), "plano paraxial",
                 size=5.5, color="#888")

    # ── Abanicos ─────────────────────────────────────────────────────────
    for idx, fan in enumerate(fans):
        color = COLORES[idx % len(COLORES)]
        foco = fan['foco']
        for pts, dirs in fan['trayectorias']:
            camino = [(px(z), px(r)) for z, r in pts]
            zf, rf = pts[-1]
            dz, dr = dirs[-1]
            t = (x_fin - zf) / dz if abs(dz) > 1e-12 else 0.0
            if t > 0:
                camino.append((px(zf + t * dz), px(rf + t * dr)))
            out += polilinea(camino, color=color, sw=0.45, op=0.75)
            if foco is not None:
                t_f = (foco[0] - zf) * dz + (foco[1] - rf) * dr
                if t_f < 0:      # imagen virtual: prolongación hacia atrás
                    out += polilinea(
                        [(px(zf), px(rf)),
                         (px(zf + t_f * dz), px(rf + t_f * dr))],
                        color=color, sw=0.3, op=0.35, dash="2,2")

    # ── Fuentes y focos reales ───────────────────────────────────────────
    for idx, fan in enumerate(fans):
        color = COLORES[idx % len(COLORES)]
        out += circulo(px(d_0), px(fan['h']), px(1.2), fill=color)
        if fan['foco'] is not None:
            f = fan['foco']
            out += circulo(px(f[0]), px(f[1]), px(1.4), fill="none",
                           stroke=color, sw=1.3)

    # ── Superficie imagen curva (las dos ramas) ──────────────────────────
    # Se dibuja sobre un realce blanco para que destaque entre los rayos.
    rama_sup = [(px(e['foco'][0]), px(e['foco'][1])) for e in curva]
    rama_inf = [(x, -y) for x, y in rama_sup]
    for rama, op in ((rama_sup, 1.0), (rama_inf, 0.45)):
        out += polilinea(rama, color="#ffffff", sw=3.0, op=0.85 * op)
        out += polilinea(rama, color="#7b1fa2", sw=1.3, dash="4,2.5", op=op)

    # ── Rótulos ──────────────────────────────────────────────────────────
    titulo, subtitulo = DESCRIPCION[tipo]
    if tipo == 1:
        titulo += "  (imagen real)" if rama == "real" else "  (imagen virtual)"
    cx = (xmin + xmax) / 2
    out += texto(cx, ymin + px(5.0), titulo, size=9.5, peso="bold")
    out += texto(cx, ymin + px(10.0), subtitulo, size=6.5, color="#666")
    out += texto(px(d_0), -px(max(alturas)) - px(3.5),
                 "4 fuentes puntuales", size=6, color="#444")
    # Los focos y la superficie imagen se explican en el bloque de cifras
    # y en el recuadro de detalle; aquí sólo se marca el extremo de la rama
    # para que se entienda de qué curva se trata.
    x_lbl, y_lbl = rama_sup[-1]
    out += texto(x_lbl, y_lbl - px(3.5) if y_lbl <= 0 else y_lbl + px(5.5),
                 "focos reales", size=6, color="#7b1fa2")

    # ── Cifras de validación ─────────────────────────────────────────────
    x_txt = xmin + px(4.0)
    y_txt = px(r_doc) + px(8.0)
    peor_spot = max(f['spot'] for f in fans if f['spot'] is not None)
    filas = [
        (f"aumento M = {M_lat:+.3f}   ·   apertura {r_ap:.1f} mm   ·   "
         f"pupila r = {r_pupila:.0f} mm en z = {Z_PUPILA:.0f} mm", "#333"),
        (f"condición seno de Abbe:  sin u0 / sin uN = {abbe['media']:+.6f}, "
         f"constante a {abbe['desviacion_rel']:.0e}   ⇒  aplanático", "#333"),
        (f"spot meridional ≤ {peor_spot * 1000:.2f} µm en los cuatro focos"
         f"   ⇒  sin aberración esférica ni coma", "#1b5e20"),
        (f"los focos NO caen sobre el plano paraxial: se apartan hasta "
         f"{dz_max:.2f} mm   ⇒  curvatura de campo", "#7b1fa2"),
    ]
    for i, (s, col) in enumerate(filas):
        out += texto(x_txt, y_txt + px(5.2 * i), s, size=6.2,
                     anchor="start", color=col)

    # ── Recuadro de detalle: zona imagen ampliada en z ───────────────────
    zoom = 13.0 / dz_max
    ancho_ins = px(30.0)
    alto_ins = px(max(alturas)) * 0.62 + px(5.0)
    ins_x = cx
    ins_y = px(r_doc) + px(ALTO_TXT) + alto_ins + px(5.0)

    out += (f'    <rect x="{ins_x - ancho_ins:.2f}" '
            f'y="{ins_y - alto_ins:.2f}" width="{2 * ancho_ins:.2f}" '
            f'height="{2 * alto_ins:.2f}" rx="3" '
            f'style="fill:#fbfaff;stroke:#c9b8d6;stroke-width:0.6px"/>\n')

    def ins(z, r):
        return (ins_x + (z - d_2) * zoom * SCALE, ins_y + r * SCALE * 0.62)

    p_lo, p_hi = ins(d_2, -max(alturas)), ins(d_2, max(alturas))
    out += linea(p_lo[0], p_lo[1], p_hi[0], p_hi[1],
                 color="#999", sw=0.7, dash="3,2")
    pc = [ins(e['foco'][0], e['foco'][1]) for e in curva]
    out += polilinea(pc, color="#7b1fa2", sw=1.4, dash="4,2.5")
    out += polilinea([(x, 2 * ins_y - y) for x, y in pc],
                     color="#7b1fa2", sw=1.4, dash="4,2.5", op=0.4)
    for idx, fan in enumerate(fans):
        if fan['foco'] is None:
            continue
        x, y = ins(fan['foco'][0], fan['foco'][1])
        out += circulo(x, y, 2.4, fill=COLORES[idx % len(COLORES)])

    out += texto(ins_x, ins_y - alto_ins - px(2.0),
                 f"detalle de la zona imagen — eje z ampliado x{zoom:.0f}",
                 size=6, color="#555")
    out += texto(ins_x - ancho_ins + px(1.5), ins_y + alto_ins - px(1.5),
                 "gris: plano paraxial", size=5, anchor="start", color="#999")
    out += texto(ins_x + ancho_ins - px(1.5), ins_y - alto_ins + px(4.0),
                 "morado: superficie imagen", size=5, anchor="end",
                 color="#7b1fa2")

    out += svg_footer()
    with open(out_path, "w") as f:
        f.write(out)

    return diseno, fans, abbe, dz_max, r_ap


def main():
    for tipo in (0, 1, 2, 3):
        ruta = os.path.join(AQUI, f"ejemplo_aplanatica_tipo{tipo}.svg")
        dis, fans, abbe, dz_max, r_ap = generar_aplanatica(ruta, tipo)
        d0, d1, d2 = dis['d']
        d1_txt = "inf" if np.isinf(d1) else f"{d1:.6f}"
        peor = max(f['spot'] for f in fans if f['spot'] is not None)
        print(f"{os.path.basename(ruta)}")
        print(f"    d0={d0:.3f}  d1={d1_txt}  d2={d2:.6f}  "
              f"M={dis['magnificacion']:+.4f}  apertura={r_ap:.2f} mm")
        print(f"    seno de Abbe = {abbe['media']:+.9f} "
              f"(dispersion {abbe['desviacion_rel']:.1e})   "
              f"spot <= {peor * 1000:.3f} um   "
              f"curvatura de campo = {dz_max:.3f} mm")

    ruta = os.path.join(AQUI, "ejemplo_aplanatica_3lentes.svg")
    M_tot, abbe, rms, spot, dz_max, r_lente = generar_3_lentes(ruta)
    print(f"\n{os.path.basename(ruta)}")
    print(f"    M total = {M_tot:+.4f}   aperturas = "
          + ", ".join(f"{r:.1f}" for r in r_lente) + " mm")
    print(f"    seno de Abbe = {abbe['media']:+.9f} "
          f"(dispersion {abbe['desviacion_rel']:.1e})   "
          f"(M-1)_RMS = {rms:.2e}")
    print(f"    spot <= {spot * 1000:.3f} um   "
          f"curvatura de campo = {dz_max:.3f} mm")

    ruta_ink = os.path.join(AQUI, "ejemplo_aplanatica_3lentes_inkscape.svg")
    generar_3_lentes(ruta_ink, semillas=True)
    print(f"{os.path.basename(ruta_ink)}  (listo para "
          f"Extensions > Optics > Ray Tracing)")


# ═══════════════ Sistema aplanático de tres lentes ═════════════════════════
#
# La condición de aplanatismo M(N, rho_k) = 1 es un PRODUCTO sobre las
# superficies (Ec. 78), de modo que encadenar elementos rigurosamente
# aplanáticos da un sistema rigurosamente aplanático.  Aquí se encadenan
# tres lentes tipo-1 (biconvexas, interior colimado, M = -1 cada una)
# formando un relé 1:1: la imagen real de cada lente es el objeto de la
# siguiente.  El sistema completo tiene seis superficies cartesianas.

def _sistema_3_lentes():
    """Devuelve (superficies, disenos, conjugados) del relé de tres lentes."""
    N0, N1 = 1.0, 1.8
    ESPESOR, CONJUGADO = 10.0, 60.0

    disenos, superficies, conjugados = [], [], [0.0]
    z_vertice = CONJUGADO                       # primer vértice
    d_obj = 0.0
    for _ in range(3):
        dis = disenar_aplanatica(1, N0, N1, z_vertice, z_vertice + ESPESOR,
                                 d_obj, rama="real")
        disenos.append(dis)
        superficies.extend(dis['superficies'])
        d_obj = dis['d'][2]                     # la imagen es el objeto siguiente
        conjugados.append(d_obj)
        z_vertice = d_obj + CONJUGADO           # siguiente lente
    return superficies, disenos, conjugados


def generar_3_lentes(out_path, alturas=(0.0, 1.0, 2.0, 3.0),
                     n_rayos=9, r_pupila=4.0, semillas=False):
    """Relé aplanático de tres lentes tipo-1 con trazado de rayos.

    Con ``semillas=True`` no se dibujan las trayectorias ya calculadas
    sino sólo los segmentos ``optics:beam`` de partida, de modo que el
    archivo se puede abrir en Inkscape y trazar con
    Extensions → Optics → Ray Tracing.
    """
    superficies, disenos, conjugados = _sistema_3_lentes()
    Z_PUPILA = 55.0
    d_0 = 0.0
    d_final = conjugados[-1]
    M_total = float(np.prod([d['magnificacion'] for d in disenos]))

    # ── Apertura necesaria en cada lente ─────────────────────────────────
    # Se mide en UNA sola pasada con una apertura de prueba muy holgada.
    # Iterar sería un error: en cuanto la apertura de prueba se ajusta, los
    # rayos que la exceden quedan viñeteados y dejan de contar, con lo que
    # el bucle converge a un punto fijo falso y las lentes salen pequeñas.
    r_prueba = [120.0] * len(superficies)
    r_sup = [0.0] * len(superficies)
    for h in alturas:
        for d0 in direcciones_pupila(d_0, h, Z_PUPILA, r_pupila, n_rayos):
            pts, _, _ = trazar_rayo(superficies, (d_0, float(h)),
                                    d0, r_prueba)
            if pts is None:
                continue
            for i, (_, r) in enumerate(pts[1:]):
                r_sup[i] = max(r_sup[i], abs(r))
    # El tope físico es el radio al que las dos superficies de la lente se
    # cruzan: más allá el elemento no tiene espesor y el trazado se rompe.
    r_fisico = [_apertura_maxima(d) for d in disenos]
    r_lente = [min(1.12 * max(r_sup[2 * i], r_sup[2 * i + 1]), r_fisico[i])
               for i in range(3)]
    for i, (nec, tope) in enumerate(zip(r_lente, r_fisico)):
        if 1.12 * max(r_sup[2 * i], r_sup[2 * i + 1]) > tope + 1e-9:
            print(f"    aviso: L{i + 1} necesitaría "
                  f"{1.12 * max(r_sup[2 * i], r_sup[2 * i + 1]):.1f} mm pero "
                  f"sus superficies se cruzan en {tope:.1f} mm; habrá viñeteo")
    aperturas = [r_lente[i // 2] for i in range(6)]

    # ── Trazado de los cuatro abanicos ───────────────────────────────────
    fans = []
    for h in alturas:
        trayect, orig_sal, dir_sal = [], [], []
        for d0 in direcciones_pupila(d_0, h, Z_PUPILA, r_pupila, n_rayos):
            pts, dirs, _ = trazar_rayo(superficies, (d_0, float(h)),
                                       d0, aperturas)
            if pts is None:
                continue
            trayect.append((pts, dirs))
            orig_sal.append(pts[-1])
            dir_sal.append(dirs[-1])
        foco, spot = (foco_meridional(orig_sal, dir_sal)
                      if len(orig_sal) >= 3 else (None, None))
        fans.append({'h': h, 'trayectorias': trayect,
                     'foco': foco, 'spot': spot})

    curva = superficie_imagen(superficies, d_0,
                              np.linspace(0, max(alturas), 11), aperturas,
                              n_rayos=n_rayos, z_pupila=Z_PUPILA,
                              r_pupila=r_pupila)
    abbe = seno_abbe(superficies, d_0, aperturas,
                     semiangulo_deg=np.degrees(
                         np.arctan(r_pupila / (Z_PUPILA - d_0))))

    # ── Comprobación de la condición de aplanatismo del conjunto ─────────
    sistema = {'superficies': superficies}
    valores = []
    for d0 in direcciones_pupila(d_0, 0.0, Z_PUPILA, r_pupila, 11):
        _, _, rhos = trazar_rayo(superficies, (d_0, 0.0), d0, aperturas)
        if rhos is None:
            continue
        valores.append(condicion_aplanatismo(sistema, rhos))
    rms = rms_aplanatismo(valores)

    focos = [f['foco'] for f in fans if f['foco'] is not None]
    dz_max = max(abs(e['foco'][0] - d_final) for e in curva) or 1e-6

    # ── Encuadre ─────────────────────────────────────────────────────────
    x_ini = d_0 - 18.0
    x_fin = max(d_final, max(f[0] for f in focos)) + 30.0
    r_doc = max(max(r_lente), max(alturas)) + 8.0
    ALTO_TXT = 24.0
    xmin, xmax = px(x_ini), px(x_fin)
    ymin = -px(r_doc) - px(14.0)
    ymax = px(r_doc) + px(ALTO_TXT + 6.0)

    out = svg_header(xmin, xmax, ymin, ymax)
    out += linea(xmin, 0, xmax, 0, color="#bbbbbb", sw=0.4, dash="6,3")

    # ── Lentes ───────────────────────────────────────────────────────────
    for i, dis in enumerate(disenos):
        contorno = _contorno(dis, r_lente[i])
        pts_svg = [(px(x), px(y)) for x, y in contorno]
        out += (f'    <path d="{puntos_a_bezier_path_str(pts_svg, cerrar=True)}" '
                f'style="fill:#c2ddb7;fill-opacity:0.75;stroke:#33691e;'
                f'stroke-width:0.7px;stroke-linejoin:round">'
                f'<desc>optics:glass:1.8000</desc></path>\n')
        out += texto(px(dis['zetas'][0] + 5.0), -px(r_lente[i]) - px(3.0),
                     f"L{i + 1}", size=7, color="#33691e", peso="bold")

    # ── Pupila ───────────────────────────────────────────────────────────
    for lado in (-1.0, +1.0):
        out += (f'    <path d="M {px(Z_PUPILA):.3f},'
                f'{lado * px(r_pupila * 1.02):.3f} '
                f'L {px(Z_PUPILA):.3f},{lado * px(r_doc * 0.9):.3f}" '
                f'style="fill:none;stroke:#333;stroke-width:1.6px">'
                f'<desc>optics:beam_dump</desc></path>\n')
    out += texto(px(Z_PUPILA), -px(r_doc * 0.9) - px(1.8), "pupila",
                 size=5.5, color="#333")

    # ── Imágenes intermedias ─────────────────────────────────────────────
    for i, z in enumerate(conjugados[1:-1], start=1):
        out += linea(px(z), -px(r_doc * 0.55), px(z), px(r_doc * 0.55),
                     color="#bbb", sw=0.5, dash="2,2")
        out += texto(px(z), px(r_doc * 0.55) + px(3.0),
                     f"imagen intermedia {i}", size=5, color="#999")

    # ── Plano imagen paraxial final ──────────────────────────────────────
    out += linea(px(d_final), -px(r_doc * 0.8), px(d_final), px(r_doc * 0.8),
                 color="#999", sw=0.7, dash="3,2")

    # ── Abanicos ─────────────────────────────────────────────────────────
    for idx, fan in enumerate(fans):
        color = COLORES[idx % len(COLORES)]
        if semillas:
            # sólo el segmento inicial, etiquetado para el trazador
            for d0 in direcciones_pupila(d_0, fan['h'], Z_PUPILA,
                                         r_pupila, n_rayos):
                L = 0.85 * (Z_PUPILA - d_0)
                out += (f'    <path d="M {px(d_0):.3f},{px(fan["h"]):.3f} '
                        f'L {px(d_0 + L * d0[0]):.3f},'
                        f'{px(fan["h"] + L * d0[1]):.3f}" '
                        f'style="fill:none;stroke:{color};'
                        f'stroke-width:0.5px"><desc>optics:beam</desc>'
                        f'</path>\n')
            continue
        for pts, dirs in fan['trayectorias']:
            camino = [(px(z), px(r)) for z, r in pts]
            zf, rf = pts[-1]
            dz, dr = dirs[-1]
            t = (x_fin - zf) / dz if abs(dz) > 1e-12 else 0.0
            if t > 0:
                camino.append((px(zf + t * dz), px(rf + t * dr)))
            out += polilinea(camino, color=color, sw=0.4, op=0.75)

    # ── Fuentes y focos ──────────────────────────────────────────────────
    for idx, fan in enumerate(fans):
        color = COLORES[idx % len(COLORES)]
        out += circulo(px(d_0), px(fan['h']), px(1.1), fill=color)
        if fan['foco'] is not None:
            f = fan['foco']
            out += circulo(px(f[0]), px(f[1]), px(1.3), fill="none",
                           stroke=color, sw=1.2)

    rama = [(px(e['foco'][0]), px(e['foco'][1])) for e in curva]
    for r_, op in ((rama, 1.0), ([(x, -y) for x, y in rama], 0.45)):
        out += polilinea(r_, color="#ffffff", sw=3.0, op=0.85 * op)
        out += polilinea(r_, color="#7b1fa2", sw=1.3, dash="4,2.5", op=op)

    # ── Rótulos y cifras ─────────────────────────────────────────────────
    cx = (xmin + xmax) / 2
    out += texto(cx, ymin + px(5.0),
                 "Sistema aplanático de tres lentes — relé 1:1 de seis "
                 "superficies cartesianas", size=9.5, peso="bold")
    out += texto(cx, ymin + px(10.0),
                 "tres lentes tipo-1 encadenadas: la imagen real de cada una "
                 "es el objeto de la siguiente", size=6.5, color="#666")
    out += texto(px(d_0), -px(max(alturas)) - px(3.5),
                 "4 fuentes puntuales", size=6, color="#444")

    peor_spot = max(f['spot'] for f in fans if f['spot'] is not None)
    x_txt = xmin + px(4.0)
    y_txt = px(r_doc) + px(8.0)
    filas = [
        (f"aumento total M = {M_total:+.4f}   ·   conjugados en z = "
         + ", ".join(f"{c:.0f}" for c in conjugados) + " mm", "#333"),
        (f"condición seno de Abbe del conjunto:  sin u0 / sin uN = "
         f"{abbe['media']:+.6f}, constante a {abbe['desviacion_rel']:.0e}",
         "#333"),
        (f"(M − 1)_RMS sobre las seis superficies = {rms:.2e}   ⇒  el "
         f"encadenamiento conserva el aplanatismo riguroso", "#1b5e20"),
        (f"spot meridional ≤ {peor_spot * 1000:.2f} µm; los focos se apartan "
         f"del plano paraxial hasta {dz_max:.2f} mm (curvatura de campo)",
         "#7b1fa2"),
    ]
    for i, (s, col) in enumerate(filas):
        out += texto(x_txt, y_txt + px(5.2 * i), s, size=6.2,
                     anchor="start", color=col)

    out += svg_footer()
    with open(out_path, "w") as f:
        f.write(out)

    return M_total, abbe, rms, peor_spot, dz_max, r_lente


if __name__ == "__main__":
    main()
