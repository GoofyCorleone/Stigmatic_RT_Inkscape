"""Extensión para generar una lente ovoide estigmática (LSOE) completa.

Genera el contorno cerrado de una lente singlete ovoide estigmática
(LSOE) a partir de sus parámetros de diseño, marcada como optics:glass
para la extensión inkscape-raytracing.  Además coloca un haz de rayos
divergentes desde el punto objeto, listos para el ray-tracing.

El factor de forma σ controla la distribución de curvatura:
  σ = −1 → plano-convexo (frente plano, atrás convexo)
  σ =  0 → biconvexo simétrico
  σ = +1 → convexo-plano (frente convexo, atrás plano)
"""

import numpy as np
import inkex

from gots_util import (
    calcular_gots,
    perfil_superficie,
    encontrar_apertura,
    calcular_d1_sigma,
    perfil_a_path_str,
)


class LenteOvoide(inkex.GenerateExtension):
    """Genera una lente LSOE con fuente de rayos en el plano objeto."""

    @property
    def style_lente(self):
        return {
            "stroke": "#000000",
            "fill": "#b7c2dd",
            "fill-opacity": "0.75",
            "stroke-linejoin": "round",
            "stroke-width": "0.5pt",
        }

    @staticmethod
    def add_arguments(pars):
        pars.add_argument("--tab", type=str, default="diseno")

        # ── Parámetros del sistema ───────────────────────────────────────
        pars.add_argument("--n0",      type=float, default=1.0)
        pars.add_argument("--n1",      type=float, default=1.6)
        pars.add_argument("--n2",      type=float, default=1.0)
        pars.add_argument("--zeta_0",  type=float, default=0.0)
        pars.add_argument("--zeta_1",  type=float, default=10.0)
        pars.add_argument("--d_objeto",type=float, default=-80.0)
        pars.add_argument("--d_imagen",type=float, default=110.0)
        pars.add_argument("--sigma",   type=float, default=0.0)
        pars.add_argument("--unidad",  type=str,   default="mm")

        # ── Apertura ────────────────────────────────────────────────────
        pars.add_argument("--r_apertura", type=float, default=0.0)
        # 0 = calcular automáticamente

        # ── Haz de rayos ────────────────────────────────────────────────
        pars.add_argument("--n_rayos",        type=int,           default=9)
        pars.add_argument("--angulo_max_deg",  type=float,         default=8.0)
        pars.add_argument("--mostrar_eje",     type=inkex.Boolean, default=True)
        pars.add_argument("--mostrar_puntos",  type=inkex.Boolean, default=True)

    # ── Helpers ─────────────────────────────────────────────────────────────

    def _sv(self, val):
        return self.svg.viewport_to_unit(f"{val}{self.options.unidad}")

    # ── Generación ──────────────────────────────────────────────────────────

    def generate(self):
        opts = self.options

        # ── 1. Calcular d₁ desde σ ───────────────────────────────────────
        try:
            d_1 = calcular_d1_sigma(
                opts.sigma,
                opts.zeta_0, opts.zeta_1,
                opts.d_objeto, opts.d_imagen,
                opts.n0, opts.n1, opts.n2,
            )
        except Exception as exc:
            inkex.utils.errormsg(f"Error calculando d₁ para σ={opts.sigma}: {exc}")
            return

        # ── 2. Parámetros GOTS de ambas superficies ──────────────────────
        try:
            params0 = calcular_gots(opts.n0, opts.n1, opts.zeta_0, opts.d_objeto, d_1)
            params1 = calcular_gots(opts.n1, opts.n2, opts.zeta_1, d_1, opts.d_imagen)
        except Exception as exc:
            inkex.utils.errormsg(f"Error calculando GOTS: {exc}")
            return

        # ── 3. Apertura ──────────────────────────────────────────────────
        if opts.r_apertura > 1e-6:
            r_ap = opts.r_apertura
        else:
            r_ap = encontrar_apertura(params0, params1)

        if r_ap < 1e-3:
            inkex.utils.errormsg(
                "La apertura calculada es demasiado pequeña.\n"
                "Introduzca r_apertura manualmente o ajuste los parámetros."
            )
            return

        # ── 4. Perfiles de ambas superficies ─────────────────────────────
        r0_arr, z0_arr = perfil_superficie(params0, N=400, r_max=r_ap)
        r1_arr, z1_arr = perfil_superficie(params1, N=400, r_max=r_ap)

        if len(r0_arr) < 2 or len(r1_arr) < 2:
            inkex.utils.errormsg(
                "No se pudieron calcular los perfiles de las superficies.\n"
                "Revise los parámetros del sistema."
            )
            return

        # ── 5. Convertir a unidades SVG (vértice frontal en el origen) ───
        zeta_ref = opts.zeta_0

        z0_svg = np.array([self._sv(z - zeta_ref) for z in z0_arr])
        z1_svg = np.array([self._sv(z - zeta_ref) for z in z1_arr])
        r0_svg = np.array([self._sv(r)             for r in r0_arr])
        r1_svg = np.array([self._sv(r)             for r in r1_arr])

        # ── 6. Construir el path cerrado de la lente ─────────────────────
        #
        # Recorrido (sentido en SVG, y hacia abajo):
        #
        #   vértice_frontal ──sup0_superior──> borde_superior
        #                                           │ rim
        #   vértice_trasero <──sup1_superior──  borde_superior
        #        │
        #   vértice_trasero ──sup1_inferior──> borde_inferior
        #                                           │ rim
        #   vértice_frontal <──sup0_inferior──  borde_inferior
        #
        puntos = []

        # vértice frontal
        puntos.append((float(z0_svg[0]), 0.0))

        # superficie frontal — mitad superior (r: 0 → r_ap)
        for z, r in zip(z0_svg[1:], r0_svg[1:]):
            puntos.append((float(z), -float(r)))

        # borde superior (conecta los dos bordes; puede ser un punto si se cruzan)
        puntos.append((float(z1_svg[-1]), -float(r1_svg[-1])))

        # superficie trasera — mitad superior invertida (r: r_ap → 0)
        for z, r in zip(reversed(list(z1_svg[:-1])), reversed(list(r1_svg[:-1]))):
            puntos.append((float(z), -float(r)))

        # vértice trasero (ya se llega aquí automáticamente)
        # superficie trasera — mitad inferior (r: 0 → r_ap)
        for z, r in zip(z1_svg[1:], r1_svg[1:]):
            puntos.append((float(z), float(r)))

        # borde inferior
        puntos.append((float(z0_svg[-1]), float(r0_svg[-1])))

        # superficie frontal — mitad inferior invertida (r: r_ap → 0)
        for z, r in zip(reversed(list(z0_svg[:-1])), reversed(list(r0_svg[:-1]))):
            puntos.append((float(z), float(r)))

        lente = inkex.PathElement()
        lente.style = self.style_lente
        lente.path  = inkex.Path(perfil_a_path_str(puntos, cerrar=True))
        lente.desc  = f"optics:glass:{opts.n1:.4f}"
        yield lente

        # ── 7. Eje óptico ────────────────────────────────────────────────
        if opts.mostrar_eje:
            margen  = self._sv(15.0)
            x_obj   = self._sv(opts.d_objeto - zeta_ref)
            x_img   = self._sv(opts.d_imagen - zeta_ref)
            eje = inkex.PathElement()
            eje.style = {
                "stroke": "#888888",
                "stroke-width": "0.3pt",
                "stroke-dasharray": f"{self._sv(3)},{self._sv(1.5)}",
                "fill": "none",
            }
            eje.path = inkex.Path(
                f"M {x_obj - margen:.4f},0 L {x_img + margen:.4f},0"
            )
            yield eje

        # ── 8. Marcadores objeto / imagen ────────────────────────────────
        if opts.mostrar_puntos:
            radio_m = self._sv(1.2)
            x_obj   = self._sv(opts.d_objeto - zeta_ref)
            x_img   = self._sv(opts.d_imagen - zeta_ref)

            for x, color in [(x_obj, "#dd2200"), (x_img, "#00aa33")]:
                circ = inkex.Circle()
                circ.set('cx', f"{x:.4f}")
                circ.set('cy', "0")
                circ.set('r',  f"{radio_m:.4f}")
                circ.style = {"fill": color, "stroke": "none"}
                yield circ

        # ── 9. Haz de rayos divergentes desde el punto objeto ────────────
        angulo_max  = opts.angulo_max_deg * np.pi / 180.0
        x_fuente    = self._sv(opts.d_objeto - zeta_ref)
        # Indicadores: distancia desde la fuente al borde izquierdo del doc
        # (aprox. = margen de 15 mm usado en el eje)
        margen_haz  = self._sv(15.0)
        L_haz       = margen_haz * 0.85
        angulos     = np.linspace(-angulo_max, angulo_max, opts.n_rayos)

        for theta in angulos:
            # Indicador dibujado desde la fuente hacia la lente (+x).
            # El ray-tracer usa el endpoint del path como origen y −tangente
            # final como dirección: para «M source L end» el rayo nace en
            # ``end`` con dirección (end−source). La línea extrapolada hacia
            # atrás pasa por la fuente, así que es equivalente a un rayo
            # originado en el punto objeto.
            x_end = float(x_fuente + L_haz * np.cos(theta))
            y_end = float(L_haz * np.sin(theta))   # SVG: y↓, óptica: r↑
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
    LenteOvoide().run()
