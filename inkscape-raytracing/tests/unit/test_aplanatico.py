"""Valida los diseños rigurosamente aplanáticos (cap. 4, Silva-Lora 2024)
contra los ejemplos numéricos de la tesis (Tablas 6-9) y la condición de
aplanatismo generalizada M(ρₖ) = 1 (Ecs. 78/93)."""

import numpy as np
from pytest import approx, raises

from inkscape_raytracing.gots_util import (
    disenar_aplanatica,
    condicion_aplanatismo,
    factor_aplanatismo,
    alturas_estigmaticas,
    rms_aplanatismo,
    perfil_superficie_aplanatica,
)

# Parámetros comunes de todos los ejemplos del §4.4 de la tesis:
# objeto en 0, vértices en 60 y 70, lente n₁ = 1.8 en aire.
ARGS = dict(n0=1.0, n1=1.8, zeta_0=60.0, zeta_1=70.0, d_0=0.0)


def _M_encadenado(diseno, rho_0):
    """M evaluado con las alturas físicas exactas del rayo (cadena
    estigmática), como exige la Definición 1 de la tesis."""
    rhos = alturas_estigmaticas(diseno, rho_0)
    return condicion_aplanatismo(diseno, rhos)


# ── Tipo-0: esferas en los puntos de Young (Tabla 6) ──────────────────────

def test_tipo0_conjugados_tabla6():
    d = disenar_aplanatica(0, **ARGS)
    assert d['d'][1] == approx(26.666667, abs=1e-5)
    assert d['d'][2] == approx(-8.0, abs=1e-9)
    assert d['superficies'][0]['R'] == approx(-60.0 / 2.8, abs=1e-9)
    assert d['superficies'][1]['R'] == approx(-27.857143, abs=1e-5)
    # O de la Tabla 6
    assert d['superficies'][0]['O'] == approx(-0.046667, abs=1e-6)
    assert d['superficies'][1]['O'] == approx(-0.035897, abs=1e-6)


def test_tipo0_factor_unidad_por_superficie():
    # En los puntos de Young cada superficie conserva el invariante de
    # Abbe por sí sola: su factor vale 1 para CUALQUIER ρ (Ec. 104).
    d = disenar_aplanatica(0, **ARGS)
    for sup in d['superficies']:
        for rho in (0.5, 3.0, 8.0, 12.0):
            assert factor_aplanatismo(sup, rho) == approx(1.0, abs=1e-12)


# ── Tipo-1: interior colimado, bicónicas (Tabla 7) ────────────────────────

def test_tipo1_real_tabla7():
    d = disenar_aplanatica(1, rama='real', **ARGS)
    assert np.isinf(d['d'][1])
    assert d['d'][2] == approx(130.0, abs=1e-9)
    s0, s1 = d['superficies']
    assert s0['G'] == approx(-3.24, abs=1e-9)
    assert s1['G'] == approx(-3.24, abs=1e-9)
    assert s0['O'] == approx(+0.020833, abs=1e-6)
    assert s1['O'] == approx(-0.020833, abs=1e-6)
    assert d['magnificacion'] == approx(-1.0, abs=1e-12)


def test_tipo1_virtual_menisco():
    d = disenar_aplanatica(1, rama='virtual', **ARGS)
    assert d['d'][2] == approx(10.0, abs=1e-9)      # imagen virtual
    assert d['magnificacion'] == approx(+1.0, abs=1e-12)


def test_tipo1_M_unidad():
    # El rayo interno viaja paralelo al eje y ambas superficies son
    # idénticas → ρ₁ = ρ₀ y el producto de factores vale 1 (Ec. 106).
    for rama in ('real', 'virtual'):
        d = disenar_aplanatica(1, rama=rama, **ARGS)
        for rho in (0.5, 4.0, 9.0, 14.0):
            assert _M_encadenado(d, rho) == approx(1.0, abs=1e-9)


# ── Tipo-2: Oₖ = 0, vértices planos (Tabla 8) ─────────────────────────────

def test_tipo2_tabla8():
    d = disenar_aplanatica(2, **ARGS)
    assert d['d'][1] == approx(-48.0, abs=1e-9)
    assert d['d'][2] == approx(4.444444, abs=1e-5)
    s0, s1 = d['superficies']
    assert s0['O'] == 0.0 and s1['O'] == 0.0
    assert s0['T'] == approx(-1.000228e-6, rel=1e-4)
    assert s0['S'] == approx(0.000259, abs=2e-6)
    assert s1['T'] == approx(-7.668749e-7, rel=1e-4)
    assert s1['S'] == approx(0.000217, abs=2e-6)


def test_tipo2_M_unidad():
    # Los factores individuales NO son constantes, pero el producto con
    # las alturas físicas del rayo vale 1 en todo rigor.
    d = disenar_aplanatica(2, **ARGS)
    valores = [_M_encadenado(d, rho) for rho in (0.5, 4.0, 9.0, 14.0)]
    assert rms_aplanatismo(valores) == approx(0.0, abs=1e-9)


# ── Tipo-3: Gₖ = 0 (Tabla 9) ──────────────────────────────────────────────

def test_tipo3_tabla9():
    d = disenar_aplanatica(3, **ARGS)
    assert d['d'][1] == approx(41.481481, abs=1e-5)
    assert d['d'][2] == approx(-22.4, abs=1e-6)
    s0, s1 = d['superficies']
    assert s0['O'] == approx(-0.100667, abs=1e-6)
    assert s1['O'] == approx(-0.065368, abs=1e-6)
    assert s0['T'] == approx(2.94e-5, rel=1e-2)
    assert s1['T'] == approx(8.049801e-6, rel=1e-4)
    assert abs(s0['G']) < 1e-9 and abs(s1['G']) < 1e-9


def test_tipo3_M_unidad():
    d = disenar_aplanatica(3, **ARGS)
    valores = [_M_encadenado(d, rho) for rho in (0.5, 4.0, 9.0, 14.0)]
    assert rms_aplanatismo(valores) == approx(0.0, abs=1e-9)


# ── Perfiles y validaciones ───────────────────────────────────────────────

def test_perfiles_alcanzan_apertura():
    for tipo, kw in [(0, {}), (1, {'rama': 'real'}),
                     (1, {'rama': 'virtual'}), (2, {}), (3, {})]:
        d = disenar_aplanatica(tipo, **ARGS, **kw)
        for sup in d['superficies']:
            r, z = perfil_superficie_aplanatica(sup, r_max=10.0, N=200)
            assert r[0] == approx(0.0, abs=1e-12)
            assert r[-1] == approx(10.0, rel=1e-3)
            assert np.all(np.diff(r) >= -1e-12)


def test_parametros_degenerados():
    with raises(ValueError):
        disenar_aplanatica(0, n0=1.0, n1=1.0, zeta_0=60, zeta_1=70, d_0=0)
    with raises(ValueError):
        disenar_aplanatica(2, n0=1.0, n1=1.8, zeta_0=60, zeta_1=70, d_0=60)
    with raises(ValueError):
        disenar_aplanatica(5, **ARGS)


# ── Trazador exacto y análisis de la imagen ───────────────────────────────

from inkscape_raytracing.gots_util import (          # noqa: E402
    trazar_rayo,
    foco_meridional,
    seno_abbe,
    superficie_imagen,
    normal_superficie,
    punto_superficie,
    condicion_aplanatismo,
    calcular_gots,
    _f_superficie,
)

D2_ESPERADO = {0: -8.0, 1: 130.0, 2: 4.4444444, 3: -22.4}


def _fan(diseno, r_max=9.0, grados=(-5, -2.5, 2.5, 5)):
    orig, dirs = [], []
    for g in grados:
        th = np.radians(g)
        pts, ds, _ = trazar_rayo(diseno['superficies'], (0.0, 0.0),
                                 (np.cos(th), np.sin(th)), r_max)
        assert pts is not None, f"rayo a {g}° perdido"
        orig.append(pts[-1])
        dirs.append(ds[-1])
    return orig, dirs


def test_normal_analitica_coincide_con_gradiente():
    # La normal se obtiene de la propiedad de V̄C̄ₖ (Ec. 90); debe coincidir
    # con el gradiente numérico del residuo de pertenencia a la superficie.
    for tipo in (0, 1, 2, 3):
        d = disenar_aplanatica(tipo, **ARGS)
        for sup in d['superficies']:
            for rho in (1.0, 5.0, 9.0):
                z, r = punto_superficie(sup, rho)
                n = np.array(normal_superficie(sup, z, r))
                h = 1e-6
                g = np.array([
                    (_f_superficie(sup, z + h, r) - _f_superficie(sup, z - h, r)),
                    (_f_superficie(sup, z, r + h) - _f_superficie(sup, z, r - h)),
                ])
                g /= np.hypot(*g)
                assert abs(abs(g @ n) - 1.0) < 1e-9


def test_foco_axial_exacto():
    # Estigmatismo riguroso: el abanico axial converge en un punto sin
    # extensión medible, en la posición de diseño.
    for tipo, d2 in D2_ESPERADO.items():
        d = disenar_aplanatica(tipo, **ARGS)
        foco, spot = foco_meridional(*_fan(d))
        assert foco[0] == approx(d2, abs=1e-6)
        assert foco[1] == approx(0.0, abs=1e-9)
        assert spot < 1e-9


def test_seno_de_abbe_igual_al_aumento():
    # Condición seno (Ec. 69/87): sin u₀/sin u_N = (n_N/n₀)·M, constante.
    for tipo in (0, 1, 2, 3):
        d = disenar_aplanatica(tipo, **ARGS)
        s = seno_abbe(d['superficies'], 0.0, 9.0)
        n0, _, nN = d['n']
        assert s['media'] == approx((nN / n0) * d['magnificacion'], abs=1e-9)
        assert s['desviacion_rel'] < 1e-12


def test_identidad_seno_aplanatismo():
    # Identidad de las Ecs. 77/86/87: sin u₀/sin u_N = (n_N/n₀)·M·M(ρₖ).
    # Se comprueba en un diseño aplanático y en otro que NO lo es, lo que
    # valida a la vez inversa_vc, condicion_aplanatismo y el trazador.
    # La apertura de prueba no puede exceder el dominio útil de cada
    # superficie, o el trazador saldría de la rama diseñada.
    casos = [(disenar_aplanatica(3, **ARGS), 9.0)]

    n0, n1, n2 = 1.0, 1.8, 1.0
    z0, z1, d0, d1, d2 = 60.0, 70.0, 0.0, 599.172082, 200.0
    p0 = calcular_gots(n0, n1, z0, d0, d1)
    p1 = calcular_gots(n1, n2, z1, d1, d2)
    casos.append(({
        'd': (d0, d1, d2), 'zetas': (z0, z1), 'n': (n0, n1, n2),
        'magnificacion': ((n0 / n2) * (d1 - z0) * (d2 - z1)
                          / ((d0 - z0) * (d1 - z1))),
        'superficies': [dict(p0, n_in=n0, n_out=n1, d_in=d0, d_out=d1),
                        dict(p1, n_in=n1, n_out=n2, d_in=d1, d_out=d2)],
    }, 20.0))

    comprobados = 0
    for d, r_max in casos:
        nn0, _, nnN = d['n']
        for grados in (1.5, 4.0, 6.0):
            th = np.radians(grados)
            pts, dirs, rhos = trazar_rayo(d['superficies'], (d['d'][0], 0.0),
                                          (np.cos(th), np.sin(th)), r_max)
            if pts is None:
                continue
            razon = np.sin(th) / dirs[-1][1]
            pred = (nnN / nn0) * d['magnificacion'] * condicion_aplanatismo(d, rhos)
            assert razon == approx(pred, rel=1e-11)
            comprobados += 1
    assert comprobados >= 5


def test_superficie_imagen_es_curva():
    # El aplanatismo NO pone los focos sobre un plano: los puntos de campo
    # enfocan nítidamente pero sobre una superficie curva.
    for tipo in (0, 1, 2, 3):
        d = disenar_aplanatica(tipo, **ARGS)
        res = superficie_imagen(d['superficies'], 0.0, [0.0, 3.0, 6.0], 9.0,
                                z_pupila=55.0, r_pupila=5.0)
        assert len(res) == 3
        # foco axial en la posición de diseño
        assert res[0]['foco'][0] == approx(D2_ESPERADO[tipo], abs=1e-5)
        # los focos fuera de eje se apartan del plano paraxial…
        desvios = [abs(e['foco'][0] - D2_ESPERADO[tipo]) for e in res]
        assert desvios[-1] > 10.0 * max(desvios[0], 1e-12)
        # …pero siguen siendo nítidos (sin coma)
        for e in res:
            assert e['spot'] < 0.05        # < 50 µm
