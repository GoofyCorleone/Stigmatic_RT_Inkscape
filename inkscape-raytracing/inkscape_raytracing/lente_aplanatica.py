"""Extensión para generar lentes singlete rigurosamente aplanáticas.

Implementa los cuatro tipos de sistemas aplanáticos en todo rigor del
capítulo 4 de la tesis de Silva-Lora (2024): lentes libres de aberración
esférica Y de coma de manera exacta (no aproximada), obtenidas como
soluciones cerradas de la condición de aplanatismo generalizada
M(ρₖ) = 1 (Ecs. 78/93):

  tipo 0 — esferas en los puntos de Young (2Sₖ = GₖOₖ²). Además libre
           de astigmatismo (anastigmática). Menisco, imagen virtual.
  tipo 1 — interior colimado (d₁ → ∞): superficies cónicas idénticas.
           Rama 'real' → biconvexa (imagen real, M = −1);
           rama 'virtual' → menisco (imagen virtual, M = +1).
  tipo 2 — Oₖ = 0: ovoides de vértice plano (generaliza la lámina
           plano-paralela). Imagen virtual, tipo lente oftálmica.
  tipo 3 — Gₖ = 0. Menisco, imagen virtual.

La lente se marca como optics:glass para la extensión de ray tracing, y
se emite un haz de prueba desde el punto objeto más, opcionalmente, una
pupila física (dos segmentos beam-dump) como en la tesis.
"""

import numpy as np
import inkex

from gots_util import (
    disenar_aplanatica,
    perfil_superficie_aplanatica,
    puntos_a_bezier_path_str,
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
    """Genera una lente rigurosamente aplanática con haz de prueba."""

    @property
    def style_lente(self):
        return {
            "stroke": "#000000",
            "fill": "#c2ddb7",
            "fill-opacity": "0.75",
            "stroke-linejoin": "round",
            "stroke-width": "0.5pt",
        }

    @staticmethod
    def add_arguments(pars):
        pars.add_argument("--tab", type=str, default="diseno")

        # ── Diseño ──────────────────────────────────────────────────────
        pars.add_argument("--tipo",     type=str,   default="1")
        pars.add_argument("--rama",     type=str,   default="real")
        pars.add_argument("--n0",       type=float, default=1.0)
        pars.add_argument("--n1",       type=float, default=1.8)
        pars.add_argument("--n2",       type=float, default=1.0)
        pars.add_argument("--zeta_0",   type=float, default=60.0)
        pars.add_argument("--zeta_1",   type=float, default=70.0)
        pars.add_argument("--d_objeto", type=float, default=0.0)
        pars.add_argument("--unidad",   type=str,   default="mm")

        # ── Apertura y haz ──────────────────────────────────────────────
        pars.add_argument("--r_apertura",     type=float, default=0.0)
        pars.add_argument("--n_rayos",        type=int,   default=9)
        pars.add_argument("--angulo_max_deg", type=float, default=8.0)
        pars.add_argument("--mostrar_eje",    type=inkex.Boolean, default=True)
        pars.add_argument("--mostrar_puntos", type=inkex.Boolean, default=True)

        # ── Pupila física ───────────────────────────────────────────────
        pars.add_argument("--mostrar_pupila", type=inkex.Boolean, default=False)
        pars.add_argument("--z_pupila",       type=float, default=55.0)
        pars.add_argument("--r_pupila",       type=float, default=10.0)

    # ── Helpers ─────────────────────────────────────────────────────────

    def _sv(self, val):
        return self.svg.viewport_to_unit(f"{val}{self.options.unidad}")

    def _apertura_defecto(self, diseno):
        """Apertura automática: donde ambos perfiles siguen siendo
        univaluados y la lente conserva un espesor de borde razonable."""
        sup0, sup1 = diseno['superficies']
        candidatos = []
        for sup in (sup0, sup1):
            if sup.get('esfera'):
                candidatos.append(0.75 * abs(sup['R']))
            else:
                coef_rad = 2.0 * sup['S'] - sup['O'] * sup['OG']
                if coef_rad < -1e-15:
                    # dominio limitado del radical → borde del óvalo
                    candidatos.append(0.85 * float(np.sqrt(-1.0 / coef_rad)))
        # Escala física del sistema como tope general
        d_0 = diseno['d'][0]
        z_0 = diseno['zetas'][0]
        candidatos.append(abs(d_0 - z_0) / 6.0)

        r_ap = min(candidatos)

        # Si las superficies se cruzan antes (lente biconvexa), recortar
        r_probe = np.linspace(r_ap * 1e-3, r_ap, 400)
        try:
            r0, z0 = perfil_superficie_aplanatica(sup0, r_ap, N=400)
            r1, z1 = perfil_superficie_aplanatica(sup1, r_ap, N=400)
            z0i = np.interp(r_probe, r0, z0)
            z1i = np.interp(r_probe, r1, z1)
            cruces = np.where(np.diff(np.sign(z1i - z0i)) != 0)[0]
            if len(cruces) > 0:
                r_ap = min(r_ap, 0.92 * float(r_probe[cruces[0]]))
        except Exception:
            pass
        return r_ap

    # ── Generación ──────────────────────────────────────────────────────

    def generate(self):
        opts = self.options
        tipo = int(opts.tipo)

        if abs(opts.n0 - 1.0) > 1e-9:
            inkex.utils.errormsg(
                "Aviso: el trazador de inkscape-raytracing supone que el\n"
                "medio exterior es aire (n = 1). Con n₀ ≠ 1 el trazado no\n"
                "coincidirá con el diseño aplanático."
            )
        if tipo != 0 and abs(opts.n2 - opts.n0) > 1e-9:
            inkex.utils.errormsg(
                f"Aviso: el tipo {tipo} exige n₂ = n₀ para que M = 1\n"
                "(Ec. 112 de la tesis); se usará n₂ = n₀ y se ignorará el\n"
                "valor introducido."
            )

        # ── 1. Diseño aplanático ─────────────────────────────────────────
        try:
            diseno = disenar_aplanatica(
                tipo, opts.n0, opts.n1,
                opts.zeta_0, opts.zeta_1, opts.d_objeto,
                n2=(opts.n2 if tipo == 0 else None),
                rama=opts.rama,
            )
        except Exception as exc:
            inkex.utils.errormsg(f"Error en el diseño aplanático: {exc}")
            return

        sup0, sup1 = diseno['superficies']
        d_0, d_1, d_2 = diseno['d']

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

        # ── 3. Perfiles ──────────────────────────────────────────────────
        r0_arr, z0_arr = perfil_superficie_aplanatica(sup0, r_ap, N=800)
        r1_arr, z1_arr = perfil_superficie_aplanatica(sup1, r_ap, N=800)
        if len(r0_arr) < 2 or len(r1_arr) < 2:
            inkex.utils.errormsg(
                "No se pudieron calcular los perfiles de las superficies."
            )
            return

        # ── 4. Unidades SVG (vértice frontal en el origen) ──────────────
        zeta_ref = opts.zeta_0
        z0_svg = np.array([self._sv(z - zeta_ref) for z in z0_arr])
        z1_svg = np.array([self._sv(z - zeta_ref) for z in z1_arr])
        r0_svg = np.array([self._sv(r)             for r in r0_arr])
        r1_svg = np.array([self._sv(r)             for r in r1_arr])

        # ── 5. Contorno cerrado (mismo recorrido que lente_ovoide) ──────
        puntos = [(float(z0_svg[0]), 0.0)]
        for z, r in zip(z0_svg[1:], r0_svg[1:]):
            puntos.append((float(z), -float(r)))
        puntos.append((float(z1_svg[-1]), -float(r1_svg[-1])))
        for z, r in zip(reversed(list(z1_svg[:-1])), reversed(list(r1_svg[:-1]))):
            puntos.append((float(z), -float(r)))
        for z, r in zip(z1_svg[1:], r1_svg[1:]):
            puntos.append((float(z), float(r)))
        puntos.append((float(z0_svg[-1]), float(r0_svg[-1])))
        for z, r in zip(reversed(list(z0_svg[:-1])), reversed(list(r0_svg[:-1]))):
            puntos.append((float(z), float(r)))

        lente = inkex.PathElement()
        lente.style = self.style_lente
        lente.path  = inkex.Path(puntos_a_bezier_path_str(puntos, cerrar=True))
        lente.desc  = f"optics:glass:{opts.n1:.4f}"
        yield lente

        # ── 6. Eje óptico ────────────────────────────────────────────────
        x_obj = self._sv(d_0 - zeta_ref)
        x_img = self._sv(d_2 - zeta_ref)
        if opts.mostrar_eje:
            margen = self._sv(15.0)
            x_lo = min(x_obj, x_img, float(z0_svg[0])) - margen
            x_hi = max(x_obj, x_img, float(z1_svg[0])) + margen
            eje = inkex.PathElement()
            eje.style = {
                "stroke": "#888888",
                "stroke-width": "0.3pt",
                "stroke-dasharray": f"{self._sv(3)},{self._sv(1.5)}",
                "fill": "none",
            }
            eje.path = inkex.Path(f"M {x_lo:.4f},0 L {x_hi:.4f},0")
            yield eje

        # ── 7. Marcadores objeto / imagen ────────────────────────────────
        if opts.mostrar_puntos:
            radio_m = self._sv(1.2)
            # imagen: relleno si es real (d₂ > ζ₁), hueco si es virtual
            img_real = d_2 > opts.zeta_1
            marcas = [(x_obj, "#dd2200", True), (x_img, "#00aa33", img_real)]
            for x, color, relleno in marcas:
                circ = inkex.Circle()
                circ.set('cx', f"{x:.4f}")
                circ.set('cy', "0")
                circ.set('r',  f"{radio_m:.4f}")
                if relleno:
                    circ.style = {"fill": color, "stroke": "none"}
                else:
                    circ.style = {"fill": "none", "stroke": color,
                                  "stroke-width": "0.6pt",
                                  "stroke-dasharray": f"{self._sv(1.2)},{self._sv(0.8)}"}
                yield circ

        # ── 8. Pupila física (dos segmentos beam-dump) ──────────────────
        if opts.mostrar_pupila:
            x_pup  = self._sv(opts.z_pupila - zeta_ref)
            r_pup  = self._sv(opts.r_pupila)
            altura = self._sv(opts.r_pupila * 0.8)
            for lado in (-1.0, +1.0):
                seg = inkex.PathElement()
                seg.style = {"stroke": "#222222", "stroke-width": "1.2pt",
                             "fill": "none"}
                seg.path = inkex.Path(
                    f"M {x_pup:.4f},{lado * r_pup:.4f} "
                    f"L {x_pup:.4f},{lado * (r_pup + altura):.4f}"
                )
                seg.desc = "optics:beam_dump"
                yield seg

        # ── 9. Haz de rayos desde el punto objeto ────────────────────────
        angulo_max = opts.angulo_max_deg * np.pi / 180.0
        L_haz      = self._sv(12.0)
        for theta in np.linspace(-angulo_max, angulo_max, opts.n_rayos):
            # El trazador toma el endpoint del path como origen del rayo y
            # −tangente final como dirección → rayo equivalente a uno
            # nacido en el punto objeto.
            x_end = float(x_obj + L_haz * np.cos(theta))
            y_end = float(L_haz * np.sin(theta))
            haz   = inkex.PathElement()
            haz.style = {"stroke": "#ff6600", "stroke-width": "0.5pt",
                         "fill": "none"}
            haz.path = inkex.Path(
                f"M {x_obj:.4f},0 L {x_end:.4f},{y_end:.4f}"
            )
            haz.desc = "optics:beam"
            yield haz


if __name__ == "__main__":
    LenteAplanatica().run()
