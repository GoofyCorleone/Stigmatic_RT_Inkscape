"""Extensión para dibujar una superficie cartesiana (óvalo de Descartes).

Siempre dibuja el **óvalo cerrado completo** — la cuártica exacta — tanto
si los parámetros se dan en forma física (n₁, n₂, d₀, d₁, ζ) como en
forma GOTS directa (G, O, T, S, ζ).  La pestaña «Geometría» sólo fija
opciones de visualización y no actúa como modo.

El contorno se emite como curvas cúbicas Bézier con tangentes estimadas
por Catmull–Rom; esto reduce el error de aproximación poligonal de
O(h²) a O(h⁴) y elimina prácticamente el desplazamiento del foco que
producía un path de segmentos rectos.

El elemento se etiqueta ``optics:glass:{n₂}`` para ser procesado por la
extensión ray-tracing de Inkscape (damienBloch/inkscape-raytracing).
"""

import numpy as np
import inkex

from gots_util import (
    calcular_gots,
    perfil_ovoide_descartes,
    perfil_gots_oval_completo,
    puntos_a_bezier_path_str,
)


class SuperficieCartesiana(inkex.GenerateExtension):
    """Dibuja un óvalo de Descartes completo como lente refractante."""

    @property
    def style(self):
        return {
            "stroke": "#000000",
            "fill": "#b7c2dd",
            "fill-opacity": "0.75",
            "stroke-linejoin": "round",
            "stroke-width": "0.5pt",
        }

    @staticmethod
    def add_arguments(pars):
        pars.add_argument("--tab", type=str, default="fisico")

        # ── Modo «parámetros físicos» ──────────────────────────────────
        pars.add_argument("--n1",        type=float, default=1.0)
        pars.add_argument("--n2",        type=float, default=1.5)
        pars.add_argument("--d_objeto",  type=float, default=-100.0)
        pars.add_argument("--d_imagen",  type=float, default=200.0)
        pars.add_argument("--zeta",      type=float, default=0.0)

        # ── Modo «parámetros GOTS directos» ────────────────────────────
        pars.add_argument("--G_param",   type=float, default=1.5)
        pars.add_argument("--O_param",   type=float, default=0.015)
        pars.add_argument("--T_param",   type=float, default=1.5e-6)
        pars.add_argument("--S_param",   type=float, default=0.005)
        pars.add_argument("--zeta_gots", type=float, default=0.0)

        # ── Geometría / visualización (válido para ambos modos) ───────
        pars.add_argument("--unidad",         type=str,           default="mm")
        pars.add_argument("--mostrar_eje",    type=inkex.Boolean, default=True)
        pars.add_argument("--mostrar_puntos", type=inkex.Boolean, default=True)
        pars.add_argument("--generar_haz",    type=inkex.Boolean, default=True)
        pars.add_argument("--n_rayos",        type=int,           default=7)
        pars.add_argument("--angulo_max_deg", type=float,         default=10.0)

        # Compat: parámetro aceptado pero ignorado (el óvalo siempre se
        # dibuja completo).
        pars.add_argument("--r_apertura",     type=float, default=0.0)

    # ── helpers ───────────────────────────────────────────────────────

    def _sv(self, val):
        return self.svg.viewport_to_unit(f"{val}{self.options.unidad}")

    # ── generación ────────────────────────────────────────────────────

    def generate(self):
        opts = self.options

        # El «modo» lo determinan las pestañas Físico o GOTS; la pestaña
        # Geometría NO cambia el modo (es sólo visualización).  Para
        # cualquier otro valor de tab caemos al modo físico, que es el
        # caso de uso habitual.
        modo_gots = (opts.tab == "gots")

        # ── 1. Parámetros GOTS y focos ────────────────────────────────
        if modo_gots:
            OG = (opts.O_param * opts.G_param
                  if abs(opts.O_param) > 1e-30 else 0.0)
            params = {
                'G':  opts.G_param, 'O':  opts.O_param,
                'T':  opts.T_param, 'S':  opts.S_param,
                'OG': OG,           'zeta': opts.zeta_gots,
            }
            zeta_ref  = opts.zeta_gots
            n_vidrio  = opts.n2
            d_obj_fis = d_img_fis = None
            tiene_focos = False
        else:
            try:
                calcular_gots(
                    opts.n1, opts.n2, opts.zeta,
                    opts.d_objeto, opts.d_imagen,
                )
            except Exception as exc:
                inkex.utils.errormsg(f"Error al calcular GOTS: {exc}")
                return
            zeta_ref  = opts.zeta
            n_vidrio  = opts.n2
            d_obj_fis = opts.d_objeto
            d_img_fis = opts.d_imagen
            tiene_focos = True

        # ── 2. Perfil meridional — óvalo completo ─────────────────────
        if modo_gots:
            zs, rs = perfil_gots_oval_completo(params, N=400)
        else:
            zs, rs = perfil_ovoide_descartes(
                opts.n1, opts.n2, opts.zeta,
                opts.d_objeto, opts.d_imagen, N=400,
            )

        if zs is None or len(zs) < 4:
            inkex.utils.errormsg(
                "Los parámetros no definen un óvalo de Descartes cerrado.\n"
                "Revise que n₁ ≠ n₂, que ξ·η < 0 (objeto e imagen a lados\n"
                "opuestos del vértice) y que κ = n₂·η − n₁·ξ ≠ 0."
            )
            return

        # ── 3. Convertir a unidades SVG y construir el contorno cerrado
        z_svg = np.array([self._sv(z - zeta_ref) for z in zs])
        r_svg = np.array([self._sv(r)             for r in rs])

        puntos = [(float(z), -float(r)) for z, r in zip(z_svg, r_svg)]
        for z, r in zip(z_svg[-2:0:-1], r_svg[-2:0:-1]):
            puntos.append((float(z), float(r)))

        ovalo = inkex.PathElement()
        ovalo.style = self.style
        ovalo.path  = inkex.Path(
            puntos_a_bezier_path_str(puntos, cerrar=True)
        )
        ovalo.desc  = f"optics:glass:{n_vidrio:.4f}"
        yield ovalo

        # ── 4. Eje óptico ─────────────────────────────────────────────
        if opts.mostrar_eje:
            margen = self._sv(15.0)
            x_min  = float(np.min(z_svg)) - margen
            x_max  = float(np.max(z_svg)) + margen
            if tiene_focos:
                x_min = min(x_min, self._sv(d_obj_fis - zeta_ref) - margen)
                x_max = max(x_max, self._sv(d_img_fis - zeta_ref) + margen)
            eje = inkex.PathElement()
            eje.style = {
                "stroke": "#888888",
                "stroke-width": "0.3pt",
                "stroke-dasharray": f"{self._sv(3)},{self._sv(1.5)}",
                "fill": "none",
            }
            eje.path = inkex.Path(f"M {x_min:.4f},0 L {x_max:.4f},0")
            yield eje

        # ── 5. Marcadores objeto / imagen ─────────────────────────────
        if opts.mostrar_puntos and tiene_focos:
            radio_m = self._sv(1.2)
            for x_fis, color in [(d_obj_fis, "#dd2200"),
                                 (d_img_fis, "#00aa33")]:
                x = self._sv(x_fis - zeta_ref)
                circ = inkex.Circle()
                circ.set('cx', f"{x:.4f}")
                circ.set('cy', "0")
                circ.set('r',  f"{radio_m:.4f}")
                circ.style = {"fill": color, "stroke": "none"}
                yield circ

        # ── 6. Haz de rayos divergentes desde el punto objeto ─────────
        if opts.generar_haz and tiene_focos:
            angulo_max = opts.angulo_max_deg * np.pi / 180.0
            x_fuente   = self._sv(d_obj_fis - zeta_ref)
            margen_haz = self._sv(15.0)
            L_haz      = margen_haz * 0.85
            angulos    = np.linspace(-angulo_max, angulo_max, opts.n_rayos)

            for theta in angulos:
                # El ray-tracer toma como origen del rayo el endpoint del
                # path y la dirección −tangente final; así, para una
                # línea «M source L end», el rayo nace en ``end`` con
                # dirección (end − source).  La extrapolación inversa
                # pasa por la fuente, con lo que el resultado equivale
                # a un rayo originado en el punto objeto.
                x_end = float(x_fuente + L_haz * np.cos(theta))
                y_end = float(L_haz * np.sin(theta))
                haz   = inkex.PathElement()
                haz.style = {
                    "stroke": "#ff6600",
                    "stroke-width": "0.5pt",
                    "fill": "none",
                }
                haz.path = inkex.Path(
                    f"M {x_fuente:.4f},0 L {x_end:.4f},{y_end:.4f}"
                )
                haz.desc = "optics:beam"
                yield haz


if __name__ == "__main__":
    SuperficieCartesiana().run()
