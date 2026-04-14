"""Extensión para visualizar una superficie cartesiana (ovoide de Descartes).

Genera una lente plano-cartesiana en Inkscape: superficie cartesiana al
frente y plano plano al fondo.  El elemento se etiqueta como
optics:glass:{n2} para ser usado directamente con la extensión de
ray-tracing de Inkscape (damienBloch/inkscape-raytracing).

Opcionalmente dibuja el eje óptico, los puntos objeto e imagen y un
haz de rayos divergentes desde el punto objeto.
"""

import numpy as np
import inkex

from gots_util import (
    calcular_gots,
    perfil_superficie,
    perfil_ovoide_descartes,
    perfil_a_path_str,
)


class SuperficieCartesiana(inkex.GenerateExtension):
    """Dibuja una superficie cartesiana como lente plano-cartesiana."""

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
        # Selector de modo de entrada
        pars.add_argument("--tab", type=str, default="fisico")

        # ── Modo parámetros físicos ──────────────────────────────────────
        pars.add_argument("--n1", type=float, default=1.0)
        pars.add_argument("--n2", type=float, default=1.5)
        pars.add_argument("--d_objeto", type=float, default=-100.0)
        pars.add_argument("--d_imagen", type=float, default=200.0)
        pars.add_argument("--zeta", type=float, default=0.0)

        # ── Modo parámetros GOTS directos ───────────────────────────────
        pars.add_argument("--G_param", type=float, default=1.5)
        pars.add_argument("--O_param", type=float, default=0.015)
        pars.add_argument("--T_param", type=float, default=1.5e-6)
        pars.add_argument("--S_param", type=float, default=0.005)
        pars.add_argument("--zeta_gots", type=float, default=0.0)

        # ── Geometría de la lente ────────────────────────────────────────
        pars.add_argument("--r_apertura", type=float, default=0.0)
        # 0 = óvalo completo hasta su ecuador natural
        pars.add_argument("--unidad", type=str, default="mm")

        # ── Opciones de visualización ────────────────────────────────────
        pars.add_argument("--mostrar_eje",    type=inkex.Boolean, default=True)
        pars.add_argument("--mostrar_puntos", type=inkex.Boolean, default=True)
        pars.add_argument("--generar_haz",    type=inkex.Boolean, default=True)
        pars.add_argument("--n_rayos",        type=int,           default=7)
        pars.add_argument("--angulo_max_deg", type=float,         default=10.0)

    # ── Helpers ─────────────────────────────────────────────────────────────

    def _sv(self, val):
        """Convierte un valor en unidades 'u' a unidades de usuario SVG."""
        return self.svg.viewport_to_unit(f"{val}{self.options.unidad}")

    # ── Generación ──────────────────────────────────────────────────────────

    def generate(self):
        opts = self.options

        # ── 1. Calcular GOTS ─────────────────────────────────────────────
        if opts.tab in ("fisico", "fisica", "Físico"):
            try:
                params = calcular_gots(
                    opts.n1, opts.n2, opts.zeta, opts.d_objeto, opts.d_imagen
                )
            except Exception as exc:
                inkex.utils.errormsg(f"Error al calcular GOTS: {exc}")
                return
            n_vidrio   = opts.n2
            zeta_ref   = opts.zeta
            d_obj_fis  = opts.d_objeto
            d_img_fis  = opts.d_imagen
            tiene_puntos = True
        else:
            # Modo GOTS directo
            OG = (opts.O_param * opts.G_param
                  if abs(opts.O_param) > 1e-30 else 0.0)
            params = {
                'G': opts.G_param, 'O': opts.O_param,
                'T': opts.T_param, 'S': opts.S_param,
                'OG': OG, 'zeta': opts.zeta_gots,
            }
            n_vidrio     = opts.n2
            zeta_ref     = opts.zeta_gots
            d_obj_fis    = None
            d_img_fis    = None
            tiene_puntos = False

        # ── 2. Perfil del óvalo ──────────────────────────────────────────
        # Si r_apertura > 0 → rama cerca truncada con arista vertical (lente).
        # Si r_apertura = 0 y modo físico → óvalo cuártico de Descartes completo.
        # Si r_apertura = 0 y modo GOTS → rama cerca hasta ecuador natural.
        usa_apertura = opts.r_apertura > 1e-6

        if usa_apertura or not tiene_puntos:
            r_max_val = opts.r_apertura if usa_apertura else None
            r_arr, z_arr = perfil_superficie(params, N=400, r_max=r_max_val)

            if len(r_arr) < 2:
                inkex.utils.errormsg(
                    "No se pudo generar el perfil del óvalo de Descartes.\n"
                    "Compruebe que los parámetros son consistentes."
                )
                return

            z_svg = np.array([self._sv(z - zeta_ref) for z in z_arr])
            r_svg = np.array([self._sv(r)             for r in r_arr])

            puntos = []
            puntos.append((float(z_svg[0]), 0.0))
            for z, r in zip(z_svg[1:], r_svg[1:]):
                puntos.append((float(z), -float(r)))
            puntos.append((float(z_svg[-1]), float(r_svg[-1])))
            for z, r in zip(reversed(list(z_svg[:-1])), reversed(list(r_svg[:-1]))):
                puntos.append((float(z), float(r)))
        else:
            # Óvalo de Descartes completo (cuártica cerrada exacta)
            zs, rs = perfil_ovoide_descartes(
                opts.n1, opts.n2, opts.zeta, opts.d_objeto, opts.d_imagen,
                N=400,
            )
            if zs is None:
                inkex.utils.errormsg(
                    "El óvalo de Descartes no está cerrado para estos parámetros.\n"
                    "Pruebe con una apertura r_max > 0 para dibujar la lente truncada."
                )
                return

            z_svg = np.array([self._sv(z - zeta_ref) for z in zs])
            r_svg = np.array([self._sv(r)             for r in rs])

            # Recorrido cerrado: vértice frontal → mitad superior (r negativo en
            # SVG porque y↓) → vértice trasero → mitad inferior → cierre.
            puntos = []
            for z, r in zip(z_svg, r_svg):
                puntos.append((float(z), -float(r)))
            for z, r in zip(z_svg[-2:0:-1], r_svg[-2:0:-1]):
                puntos.append((float(z), float(r)))

        ovoide = inkex.PathElement()
        ovoide.style = self.style
        ovoide.path  = inkex.Path(perfil_a_path_str(puntos, cerrar=True))
        ovoide.desc  = f"optics:glass:{n_vidrio:.4f}"
        yield ovoide

        # ── 4. Eje óptico ────────────────────────────────────────────────
        if opts.mostrar_eje:
            margen = self._sv(15.0)
            x_izq  = (self._sv(d_obj_fis - zeta_ref) - margen
                      if tiene_puntos else -margen)
            x_der  = (self._sv(d_img_fis - zeta_ref) + margen
                      if tiene_puntos else float(np.max(z_svg)) + margen)
            eje = inkex.PathElement()
            eje.style = {
                "stroke": "#888888",
                "stroke-width": "0.3pt",
                "stroke-dasharray": f"{self._sv(3)},{self._sv(1.5)}",
                "fill": "none",
            }
            eje.path = inkex.Path(
                f"M {x_izq:.4f},0 L {x_der:.4f},0"
            )
            yield eje

        # ── 5. Marcadores objeto / imagen ────────────────────────────────
        if opts.mostrar_puntos and tiene_puntos:
            radio_m = self._sv(1.2)
            x_obj   = self._sv(d_obj_fis - zeta_ref)
            x_img   = self._sv(d_img_fis - zeta_ref)

            for x, color in [(x_obj, "#dd2200"), (x_img, "#00aa33")]:
                circ = inkex.Circle()
                circ.set('cx', f"{x:.4f}")
                circ.set('cy', "0")
                circ.set('r',  f"{radio_m:.4f}")
                circ.style = {"fill": color, "stroke": "none"}
                yield circ

        # ── 6. Haz de rayos desde el punto objeto ────────────────────────
        if opts.generar_haz and tiene_puntos:
            angulo_max = opts.angulo_max_deg * np.pi / 180.0
            x_fuente   = self._sv(d_obj_fis - zeta_ref)
            # Indicadores: desde la fuente al borde izquierdo (≈ margen)
            margen_haz = self._sv(15.0)
            L_haz      = margen_haz * 0.85
            angulos    = np.linspace(-angulo_max, angulo_max, opts.n_rayos)

            for theta in angulos:
                # Dibuja el indicador desde la fuente hacia la derecha (hacia
                # la lente). El ray-tracer toma como origen el punto final del
                # path y la dirección -tangente_final, por lo que para una
                # línea «M source L end» el rayo nace en ``end`` con dirección
                # (end − source), es decir, a la derecha. La extrapolación
                # inversa del rayo pasa por el punto fuente, así que es
                # equivalente a un rayo originado en la fuente.
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
