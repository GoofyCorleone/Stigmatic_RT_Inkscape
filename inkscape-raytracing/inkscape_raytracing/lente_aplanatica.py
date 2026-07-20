"""Extensión para generar y analizar lentes singlete aplanáticas.

Dos modos de trabajo:

1. **Tipos rigurosamente aplanáticos** (§4.4 de la tesis de Silva-Lora,
   2024): los cuatro casos cerrados en los que la condición de
   aplanatismo generalizada M(ρₖ) = 1 (Ecs. 78/93) se cumple en todo
   rigor, de modo que la lente carece de aberración esférica Y de coma
   de manera exacta, no aproximada:

     tipo 0 — esferas en los puntos de Young (2Sₖ = GₖOₖ²); además libre
              de astigmatismo (anastigmática).
     tipo 1 — interior colimado (d₁ → ∞): superficies cónicas idénticas.
              Rama 'real' → biconvexa (M = −1); 'virtual' → menisco.
     tipo 2 — Oₖ = 0: ovoides de vértice plano; imagen virtual.
     tipo 3 — Gₖ = 0; menisco, imagen virtual.

2. **Superficies libres**: se introducen directamente los conjugados de
   cada superficie (ζₖ, dₖ, d_{k+1}) y la extensión deriva sus
   parámetros GOTS con ``calcular_gots``.  El singlete resultante es
   rigurosamente estigmático para el par axial, pero en general NO
   aplanático: la pestaña de análisis mide cuánto se aparta, con la
   métrica (M − 1)_RMS de la Ec. 94 y la condición seno de Abbe.  Es el
   modo para explorar diseños propios y reproducir barridos como el de
   la Fig. 16 de la tesis.

La pestaña **Análisis** puede además dibujar dónde se forma realmente la
imagen: en un sistema aplanático los puntos objeto fuera de eje enfocan
nítidamente, pero sobre una SUPERFICIE CURVA, no sobre el plano imagen
paraxial (Fig. 17 de la tesis).  La extensión traza los abanicos de
varios puntos de campo, marca sus focos meridionales reales y une esos
focos con una curva.
"""

import numpy as np
import inkex

from gots_util import (
    calcular_gots,
    disenar_aplanatica,
    perfil_superficie_aplanatica,
    puntos_a_bezier_path_str,
    superficie_imagen,
    seno_abbe,
    condicion_aplanatismo,
    rms_aplanatismo,
    trazar_rayo,
    alturas_estigmaticas,
)

# inkex serializa los números de un path con {:.6g}: para coordenadas de
# ~10 unidades eso cuantiza a ~1e-4 y corrompe las normales de los
# segmentos Bézier finos de la superficie (errores de refracción de
# hasta ~0.5°). Subimos la precisión global igual que en render.py.
try:
    from inkex.paths.interfaces import PathCommand as _PathCommand
    _PathCommand.number_template = "{:.10g}"
except ImportError:  # versiones antiguas de inkex
    pass


class LenteAplanatica(inkex.GenerateExtension):
    """Genera una lente aplanática con haz de prueba y análisis."""

    @property
    def style_lente(self):
        return {
            "stroke": "#33691e",
            "fill": "#c2ddb7",
            "fill-opacity": "0.75",
            "stroke-linejoin": "round",
            "stroke-width": "0.5pt",
        }

    @staticmethod
    def add_arguments(pars):
        pars.add_argument("--tab", type=str, default="diseno")

        # ── Modo y diseño por tipo ──────────────────────────────────────
        pars.add_argument("--modo",     type=str,   default="tipo")
        pars.add_argument("--tipo",     type=str,   default="1")
        pars.add_argument("--rama",     type=str,   default="real")
        pars.add_argument("--n0",       type=float, default=1.0)
        pars.add_argument("--n1",       type=float, default=1.8)
        pars.add_argument("--n2",       type=float, default=1.0)
        pars.add_argument("--zeta_0",   type=float, default=60.0)
        pars.add_argument("--zeta_1",   type=float, default=70.0)
        pars.add_argument("--d_objeto", type=float, default=0.0)
        pars.add_argument("--unidad",   type=str,   default="mm")

        # ── Superficies libres ──────────────────────────────────────────
        pars.add_argument("--d_inter",  type=float, default=599.172082)
        pars.add_argument("--d_imagen", type=float, default=200.0)

        # ── Apertura y haz ──────────────────────────────────────────────
        pars.add_argument("--r_apertura",     type=float, default=0.0)
        pars.add_argument("--n_rayos",        type=int,   default=9)
        pars.add_argument("--mostrar_eje",    type=inkex.Boolean, default=True)
        pars.add_argument("--mostrar_puntos", type=inkex.Boolean, default=True)

        # ── Pupila física ───────────────────────────────────────────────
        pars.add_argument("--mostrar_pupila", type=inkex.Boolean, default=True)
        pars.add_argument("--z_pupila",       type=float, default=55.0)
        pars.add_argument("--r_pupila",       type=float, default=8.0)

        # ── Análisis ────────────────────────────────────────────────────
        pars.add_argument("--n_campo",        type=int,   default=4)
        pars.add_argument("--h_campo",        type=float, default=9.0)
        pars.add_argument("--trazar_campo",   type=inkex.Boolean, default=True)
        pars.add_argument("--dibujar_imagen", type=inkex.Boolean, default=True)
        pars.add_argument("--anotar",         type=inkex.Boolean, default=True)

    # ── Helpers ─────────────────────────────────────────────────────────

    def _sv(self, val):
        return self.svg.viewport_to_unit(f"{val}{self.options.unidad}")

    def _X(self, z):
        """Coordenada SVG horizontal (vértice frontal en el origen)."""
        return self._sv(z - self.options.zeta_0)

    def _Y(self, r):
        """Coordenada SVG vertical (r hacia arriba ⇒ y hacia abajo)."""
        return self._sv(r)

    def _apertura_defecto(self, diseno):
        """Apertura automática: dominio de cada superficie, escala física
        del sistema y punto donde ambos perfiles se cruzarían."""
        sup0, sup1 = diseno['superficies']
        topes = []
        for sup in (sup0, sup1):
            if sup.get('esfera'):
                topes.append(0.75 * abs(sup['R']))
            else:
                coef_rad = 2.0 * sup['S'] - sup['O'] * sup['OG']
                if coef_rad < -1e-15:
                    topes.append(0.85 * float(np.sqrt(-1.0 / coef_rad)))
        topes.append(abs(diseno['d'][0] - diseno['zetas'][0]) / 6.0)
        r_ap = min(topes)

        try:
            r0, z0 = perfil_superficie_aplanatica(sup0, r_ap, N=400)
            r1, z1 = perfil_superficie_aplanatica(sup1, r_ap, N=400)
            rp = np.linspace(r_ap * 1e-3, min(r0[-1], r1[-1]), 400)
            dz = np.interp(rp, r1, z1) - np.interp(rp, r0, z0)
            cruces = np.where(np.diff(np.sign(dz)) != 0)[0]
            if len(cruces) > 0:
                r_ap = min(r_ap, 0.92 * float(rp[cruces[0]]))
        except Exception:
            pass
        return r_ap

    def _diseno_libre(self):
        """Singlete estigmático a partir de los conjugados introducidos."""
        o = self.options
        p0 = calcular_gots(o.n0, o.n1, o.zeta_0, o.d_objeto, o.d_inter)
        p1 = calcular_gots(o.n1, o.n2, o.zeta_1, o.d_inter, o.d_imagen)
        sup0 = dict(p0, n_in=o.n0, n_out=o.n1,
                    d_in=o.d_objeto, d_out=o.d_inter)
        sup1 = dict(p1, n_in=o.n1, n_out=o.n2,
                    d_in=o.d_inter, d_out=o.d_imagen)
        M_lat = ((o.n0 / o.n2) * (o.d_inter - o.zeta_0)
                 * (o.d_imagen - o.zeta_1)
                 / ((o.d_objeto - o.zeta_0) * (o.d_inter - o.zeta_1)))
        return {
            'tipo': 'libre',
            'n': (o.n0, o.n1, o.n2),
            'zetas': (o.zeta_0, o.zeta_1),
            'd': (o.d_objeto, o.d_inter, o.d_imagen),
            'magnificacion': M_lat,
            'superficies': [sup0, sup1],
        }

    # ── Generación ──────────────────────────────────────────────────────

    def generate(self):
        opts = self.options
        modo_libre = (opts.modo == "libre")
        tipo = int(opts.tipo)

        if abs(opts.n0 - 1.0) > 1e-9:
            inkex.utils.errormsg(
                "Aviso: el trazador de inkscape-raytracing supone que el\n"
                "medio exterior es aire (n = 1). Con n₀ ≠ 1 el trazado no\n"
                "coincidirá con el diseño."
            )
        if (not modo_libre) and tipo != 0 and abs(opts.n2 - opts.n0) > 1e-9:
            inkex.utils.errormsg(
                f"Aviso: el tipo {tipo} exige n₂ = n₀ para que M = 1\n"
                "(Ec. 112 de la tesis); se usará n₂ = n₀ y se ignorará el\n"
                "valor introducido."
            )

        # ── 1. Diseño ────────────────────────────────────────────────────
        try:
            if modo_libre:
                diseno = self._diseno_libre()
            else:
                diseno = disenar_aplanatica(
                    tipo, opts.n0, opts.n1,
                    opts.zeta_0, opts.zeta_1, opts.d_objeto,
                    n2=(opts.n2 if tipo == 0 else None),
                    rama=opts.rama,
                )
        except Exception as exc:
            inkex.utils.errormsg(f"Error en el diseño: {exc}")
            return

        sups = diseno['superficies']
        d_0, d_1, d_2 = diseno['d']
        M_lat = diseno['magnificacion']

        # ── 2. Apertura ──────────────────────────────────────────────────
        if opts.r_apertura > 1e-6:
            r_ap = opts.r_apertura
        else:
            r_ap = self._apertura_defecto(diseno)
        if r_ap < 1e-3:
            inkex.utils.errormsg(
                "La apertura calculada es demasiado pequeña.\n"
                "Introduzca r_apertura manualmente."
            )
            return

        # ── 3. Perfiles y contorno ───────────────────────────────────────
        r0_arr, z0_arr = perfil_superficie_aplanatica(sups[0], r_ap, N=800)
        r1_arr, z1_arr = perfil_superficie_aplanatica(sups[1], r_ap, N=800)
        if len(r0_arr) < 2 or len(r1_arr) < 2:
            inkex.utils.errormsg(
                "No se pudieron calcular los perfiles de las superficies."
            )
            return

        z0s = np.array([self._X(z) for z in z0_arr])
        z1s = np.array([self._X(z) for z in z1_arr])
        r0s = np.array([self._Y(r) for r in r0_arr])
        r1s = np.array([self._Y(r) for r in r1_arr])

        puntos = [(float(z0s[0]), 0.0)]
        for z, r in zip(z0s[1:], r0s[1:]):
            puntos.append((float(z), -float(r)))
        puntos.append((float(z1s[-1]), -float(r1s[-1])))
        for z, r in zip(reversed(list(z1s[:-1])), reversed(list(r1s[:-1]))):
            puntos.append((float(z), -float(r)))
        for z, r in zip(z1s[1:], r1s[1:]):
            puntos.append((float(z), float(r)))
        puntos.append((float(z0s[-1]), float(r0s[-1])))
        for z, r in zip(reversed(list(z0s[:-1])), reversed(list(r0s[:-1]))):
            puntos.append((float(z), float(r)))

        lente = inkex.PathElement()
        lente.style = self.style_lente
        lente.path = inkex.Path(puntos_a_bezier_path_str(puntos, cerrar=True))
        lente.desc = f"optics:glass:{opts.n1:.4f}"
        yield lente

        # ── 4. Eje óptico ────────────────────────────────────────────────
        x_obj = self._X(d_0)
        x_img = self._X(d_2)
        if opts.mostrar_eje:
            margen = self._sv(15.0)
            x_lo = min(x_obj, x_img, float(z0s[0])) - margen
            x_hi = max(x_obj, x_img, float(z1s[0])) + margen
            eje = inkex.PathElement()
            eje.style = {
                "stroke": "#aaaaaa", "stroke-width": "0.3pt", "fill": "none",
                "stroke-dasharray": f"{self._sv(3)},{self._sv(1.5)}",
            }
            eje.path = inkex.Path(f"M {x_lo:.4f},0 L {x_hi:.4f},0")
            yield eje

        # ── 5. Pupila física ─────────────────────────────────────────────
        if opts.mostrar_pupila:
            x_pup = self._X(opts.z_pupila)
            r_pup = self._Y(opts.r_pupila * 1.02)   # holgura: si el borde
            # coincide con el rayo marginal, el trazador lo absorbe
            largo = self._Y(opts.r_pupila * 0.8)
            for lado in (-1.0, +1.0):
                seg = inkex.PathElement()
                seg.style = {"stroke": "#222222", "stroke-width": "1.2pt",
                             "fill": "none"}
                seg.path = inkex.Path(
                    f"M {x_pup:.4f},{lado * r_pup:.4f} "
                    f"L {x_pup:.4f},{lado * (r_pup + largo):.4f}")
                seg.desc = "optics:beam_dump"
                yield seg

        # ── 6. Abanicos de los puntos de campo ───────────────────────────
        alturas = ([0.0] if opts.n_campo <= 1
                   else list(np.linspace(0.0, opts.h_campo, opts.n_campo)))
        paleta = ["#1144cc", "#00913f", "#dd8800", "#cc2222",
                  "#7b1fa2", "#00838f"]

        if opts.trazar_campo:
            for idx, h in enumerate(alturas):
                color = paleta[idx % len(paleta)]
                for y in np.linspace(-opts.r_pupila, opts.r_pupila,
                                     opts.n_rayos):
                    dz = opts.z_pupila - d_0
                    dr = y - h
                    norma = np.hypot(dz, dr)
                    # El trazador toma el endpoint del path como origen del
                    # rayo y −tangente final como dirección, así que basta
                    # dibujar el segmento fuente → pupila.
                    L = self._sv(min(abs(dz) * 0.9, 25.0))
                    haz = inkex.PathElement()
                    haz.style = {"stroke": color, "stroke-width": "0.4pt",
                                 "fill": "none"}
                    haz.path = inkex.Path(
                        f"M {x_obj:.4f},{self._Y(h):.4f} "
                        f"L {x_obj + L * dz / norma:.4f},"
                        f"{self._Y(h) + L * dr / norma:.4f}")
                    haz.desc = "optics:beam"
                    yield haz

        # ── 7. Marcadores objeto ─────────────────────────────────────────
        if opts.mostrar_puntos:
            for idx, h in enumerate(alturas):
                circ = inkex.Circle()
                circ.set('cx', f"{x_obj:.4f}")
                circ.set('cy', f"{self._Y(h):.4f}")
                circ.set('r', f"{self._sv(1.2):.4f}")
                circ.style = {"fill": paleta[idx % len(paleta)],
                              "stroke": "none"}
                yield circ

        # ── 8. Análisis: dónde se forma realmente la imagen ──────────────
        resultados = None
        if opts.dibujar_imagen or opts.anotar:
            try:
                resultados = superficie_imagen(
                    sups, d_0, np.linspace(0.0, opts.h_campo, 13), r_ap,
                    n_rayos=max(opts.n_rayos, 7),
                    z_pupila=opts.z_pupila, r_pupila=opts.r_pupila)
                abbe = seno_abbe(
                    sups, d_0, r_ap,
                    semiangulo_deg=np.degrees(
                        np.arctan(opts.r_pupila / abs(opts.z_pupila - d_0))))
            except Exception as exc:
                inkex.utils.errormsg(f"No se pudo analizar la imagen: {exc}")
                resultados = None

        if resultados and opts.dibujar_imagen:
            # plano imagen paraxial, como referencia
            plano = inkex.PathElement()
            plano.style = {"stroke": "#999999", "stroke-width": "0.4pt",
                           "fill": "none",
                           "stroke-dasharray": f"{self._sv(1.5)},{self._sv(1)}"}
            plano.path = inkex.Path(
                f"M {x_img:.4f},{-self._Y(opts.h_campo * 1.25):.4f} "
                f"L {x_img:.4f},{self._Y(opts.h_campo * 1.25):.4f}")
            yield plano

            # superficie imagen curva (ambas ramas)
            for signo in (+1.0, -1.0):
                pts = [(self._X(e['foco'][0]), signo * self._Y(e['foco'][1]))
                       for e in resultados]
                curva = inkex.PathElement()
                curva.style = {
                    "stroke": "#7b1fa2", "stroke-width": "0.6pt",
                    "fill": "none",
                    "stroke-dasharray": f"{self._sv(2)},{self._sv(1.2)}"}
                curva.path = inkex.Path(puntos_a_bezier_path_str(
                    pts, cerrar=False))
                yield curva

            # focos reales de cada punto de campo
            for idx, h in enumerate(alturas):
                cerc = min(resultados, key=lambda e: abs(e['h'] - h))
                circ = inkex.Circle()
                circ.set('cx', f"{self._X(cerc['foco'][0]):.4f}")
                circ.set('cy', f"{self._Y(cerc['foco'][1]):.4f}")
                circ.set('r', f"{self._sv(1.3):.4f}")
                circ.style = {"fill": "none",
                              "stroke": paleta[idx % len(paleta)],
                              "stroke-width": "0.5pt"}
                yield circ

        # ── 9. Anotación numérica de la verificación ─────────────────────
        if resultados and opts.anotar:
            rhos = [alturas_estigmaticas(diseno, rho)
                    for rho in np.linspace(r_ap * 0.15, r_ap * 0.95, 7)]
            valores = [condicion_aplanatismo(diseno, par) for par in rhos]
            rms = rms_aplanatismo(valores)
            dz_max = max(abs(e['foco'][0] - d_2) for e in resultados)
            spot_max = max(e['spot'] for e in resultados)
            aplanatico = rms < 1e-9

            lineas = [
                (f"{'lente aplanática tipo ' + str(tipo) if not modo_libre else 'singlete estigmático (superficies libres)'}"
                 f"   ·   n = {opts.n0:g}/{opts.n1:g}/{diseno['n'][2]:g}"),
                (f"conjugados: d0 = {d_0:g}, "
                 f"d1 = {'inf' if np.isinf(d_1) else format(d_1, '.6g')}, "
                 f"d2 = {d_2:.6g}   ·   aumento M = {M_lat:+.4f}"),
                (f"condición seno de Abbe: sin u0/sin uN = {abbe['media']:+.6f} "
                 f"(dispersión {abbe['desviacion_rel']:.1e})"),
                (f"(M − 1)_RMS = {rms:.3e}   ⇒   "
                 f"{'APLANÁTICO en todo rigor' if aplanatico else 'NO aplanático'}"),
                (f"spot meridional ≤ {spot_max * 1000:.2f} µm; "
                 f"los focos se apartan del plano paraxial hasta {dz_max:.3f}"),
            ]
            y0 = self._Y(max(r_ap, opts.h_campo)) + self._sv(9.0)
            for i, s in enumerate(lineas):
                t = inkex.TextElement()
                t.text = s
                t.set("x", f"{x_obj:.4f}")
                t.set("y", f"{y0 + i * self._sv(4.2):.4f}")
                t.style = {
                    "font-size": f"{self._sv(2.6):.4f}",
                    "font-family": "sans-serif",
                    "fill": "#7b1fa2" if i == 3 else "#333333",
                }
                yield t


if __name__ == "__main__":
    LenteAplanatica().run()
