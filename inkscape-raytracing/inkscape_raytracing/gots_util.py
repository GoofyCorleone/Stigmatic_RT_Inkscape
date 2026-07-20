"""Utilidades GOTS auto-contenidas para las extensiones de Inkscape.

Implementa los cálculos de parámetros GOTS y perfil de superficie
cartesiana, adaptados de la tesis de Silva-Lora (2024) y del paquete
gots en Python/RayTracing/gots/.

Este módulo no tiene dependencias externas más allá de numpy.
"""

import numpy as np


# ── Cálculo de parámetros GOTS ─────────────────────────────────────────────


def calcular_gots(n_k, n_k1, zeta_k, d_k, d_k1):
    """Calcula los parámetros GOTS a partir de los parámetros físicos.

    Implementa las Ecs. 10-13 de Silva-Lora (2024).

    Args:
        n_k   : índice de refracción del medio objeto
        n_k1  : índice de refracción del medio imagen
        zeta_k: posición axial del vértice de la superficie
        d_k   : posición axial del punto objeto
        d_k1  : posición axial del punto imagen

    Returns:
        dict con claves: G, O, T, S, OG, zeta

    Raises:
        ValueError: si los parámetros son degenerados.
    """
    xi    = d_k  - zeta_k
    eta   = d_k1 - zeta_k
    kappa = n_k1 * eta - n_k * xi

    for nombre, val in [('xi',    xi),
                        ('eta',   eta),
                        ('n1-n0', n_k1 - n_k),
                        ('kappa', kappa)]:
        if abs(val) < 1e-30:
            raise ValueError(
                f"Parámetro degenerado: {nombre} = {val:.3e}  "
                "(verifique que n1≠n0, objeto≠vértice, imagen≠vértice y "
                "que κ = n1·η − n0·ξ ≠ 0)"
            )

    O_k  = (n_k1 * xi - n_k * eta) / (xi * eta * (n_k1 - n_k))
    T_k  = ((n_k1 - n_k) * (n_k1 + n_k) ** 2
            / (4.0 * n_k * n_k1 * xi * eta * kappa))
    S_k  = ((n_k1 + n_k) * (n_k1 ** 2 * eta - n_k ** 2 * xi)
            / (2.0 * n_k * n_k1 * xi * eta * kappa))
    OG_k = ((n_k1 ** 2 * eta - n_k ** 2 * xi) ** 2
            / (n_k * n_k1 * xi * eta * (n_k1 - n_k) * kappa))
    G_k  = OG_k / O_k if abs(O_k) > 1e-30 else float('inf')

    return {'G': G_k, 'O': O_k, 'T': T_k, 'S': S_k, 'OG': OG_k, 'zeta': zeta_k}


# ── Geometría de la superficie ─────────────────────────────────────────────


def _tau_de_rho(params, rho):
    """τ(ρ) = z(ρ) − ζ usando la Ec. 16 de Silva-Lora."""
    O, T, S, OG = params['O'], params['T'], params['S'], params['OG']
    rho2 = np.asarray(rho, dtype=float) ** 2
    coef_rad = 2.0 * S - O * OG
    radical  = np.sqrt(np.maximum(1.0 + coef_rad * rho2, 0.0))
    denom    = 1.0 + S * rho2 + radical
    numer    = (O + T * rho2) * rho2
    return np.where(np.abs(denom) < 1e-30, 0.0, numer / denom)


def perfil_ovoide_descartes(n0, n1, zeta, d_obj, d_img, N=400):
    """Contorno cerrado del óvalo de Descartes definido por

        m₀·|P−F₀| + m₁·|P−F₁| = k ,   m₀ = ±n₀,  m₁ = ±n₁,

    donde F₀=(d_obj,0), F₁=(d_img,0) y k se fija pasando por el vértice
    (ζ, 0). Esta es la cuártica exacta del óvalo refractante.

    Los signos de m₀ y m₁ codifican la naturaleza real o virtual de los
    focos (propagación hacia +z):
      objeto real   (d_obj < ζ) → m₀ = +n₀ ;  objeto virtual → m₀ = −n₀
      imagen real   (d_img > ζ) → m₁ = +n₁ ;  imagen virtual → m₁ = −n₁
    Para imagen virtual el invariante de Fermat es la DIFERENCIA de
    caminos ópticos (n₀·l₀ − n₁·l₁ = cte), no la suma; usar siempre la
    suma produce una superficie no estigmática.

    Args:
        n0, n1         : índices de refracción a cada lado de la superficie
        zeta           : posición axial del vértice frontal
        d_obj, d_img   : posiciones de los focos (puntos objeto e imagen)
        N              : número de puntos en el contorno superior

    Returns:
        (z_arr, r_arr): arrays con la rama superior del óvalo (r≥0),
                        desde el vértice frontal hasta el vértice trasero.
                        El contorno cerrado se obtiene reflejando r→−r.
        (None, None) si los parámetros no generan un óvalo cerrado.
    """
    n0 = float(n0); n1 = float(n1)
    zeta = float(zeta); d_obj = float(d_obj); d_img = float(d_img)

    # Índices con signo según foco real/virtual (propagación hacia +z)
    m0 = n0 if d_obj < zeta else -n0
    m1 = n1 if d_img > zeta else -n1

    k = m0 * abs(zeta - d_obj) + m1 * abs(zeta - d_img)

    # Raíces axiales: z tal que m₀|z−d_obj| + m₁|z−d_img| = k.
    # Resolución por regímenes de signos.
    raices = []
    for s0 in (+1.0, -1.0):
        for s1 in (+1.0, -1.0):
            denom = m0 * s0 + m1 * s1
            if abs(denom) < 1e-30:
                continue
            z_c = (k + m0 * s0 * d_obj + m1 * s1 * d_img) / denom
            # Verificar consistencia de signos
            if (z_c - d_obj) * s0 >= -1e-6 and (z_c - d_img) * s1 >= -1e-6:
                if not any(abs(z_c - r) < 1e-6 for r in raices):
                    raices.append(z_c)
    raices.sort()

    if len(raices) < 2:
        return None, None

    # Vértice frontal = raíz más cercana a ζ
    z_v0 = min(raices, key=lambda r: abs(r - zeta))
    otras = [r for r in raices if abs(r - z_v0) > 1e-6]
    if not otras:
        return None, None
    # Vértice trasero = la otra raíz más cercana a ζ (óvalo cerrado mínimo)
    z_v1 = min(otras, key=lambda r: abs(r - zeta))

    z_lo, z_hi = sorted([z_v0, z_v1])

    # Muestreo Chebyshev: puntos densos cerca de los vértices donde la
    # superficie es casi vertical, evitando errores de refracción al
    # aproximarla por segmentos rectos.
    i    = np.arange(N)
    cheb = -np.cos(np.pi * i / (N - 1))   # de −1 a +1
    zs   = 0.5 * (z_lo + z_hi) + 0.5 * (z_hi - z_lo) * cheb
    rs = np.zeros_like(zs)
    M  = n0 * n0 - n1 * n1

    for i, z in enumerate(zs):
        a = z - d_obj
        b = z - d_img
        # La cuadrática sólo depende de n₀², n₁² y k² (al elevar al
        # cuadrado se funden las dos ramas del óvalo); los signos de
        # m₀, m₁ actúan al seleccionar la raíz correcta más abajo.
        Nc = n0 * n0 * a * a - n1 * n1 * b * b - k * k
        Q2 = 4.0 * k * k * n1 * n1

        if abs(M) < 1e-30:
            # n₀ = n₁: degenerado (no hay refracción)
            rs[i] = 0.0
            continue

        # Cuadrática en u=r²: M² u² + (2MN − Q²) u + (N² − Q² b²) = 0
        Bq = 2.0 * M * Nc - Q2
        Cq = Nc * Nc - Q2 * b * b
        disc = Bq * Bq - 4.0 * M * M * Cq
        if disc < 0:
            rs[i] = 0.0
            continue
        sd = np.sqrt(disc)
        u1 = (-Bq + sd) / (2.0 * M * M)
        u2 = (-Bq - sd) / (2.0 * M * M)

        # Elegir la raíz ≥ 0 que satisface m₀√(a²+u) + m₁√(b²+u) = k
        mejor_u  = 0.0
        mejor_er = float('inf')
        for u in (u1, u2):
            if u < -1e-9:
                continue
            uu = max(u, 0.0)
            err = abs(m0 * np.sqrt(a * a + uu)
                      + m1 * np.sqrt(b * b + uu) - k)
            if err < mejor_er:
                mejor_er = err
                mejor_u  = uu
        rs[i] = np.sqrt(mejor_u)

    # Ordenar: queremos empezar en el vértice frontal
    if z_v0 > z_v1:
        zs = zs[::-1].copy()
        rs = rs[::-1].copy()
    # Asegurar r=0 exacto en los vértices
    rs[0]  = 0.0
    rs[-1] = 0.0
    return zs, rs


def perfil_gots_oval_completo(params, N=400):
    """Óvalo de Descartes COMPLETO a partir de los parámetros GOTS.

    Resuelve la cuadrática en τ (Eq. 9 de Silva-Lora) combinando las dos
    ramas τ₋ (cerca del vértice) y τ₊ (lejos del vértice) para obtener el
    contorno cerrado del óvalo en el plano meridiano.

    Devuelve:
        (z_arr, r_arr) con la rama superior (r ≥ 0) desde el vértice
        frontal hasta el vértice trasero del óvalo, ordenada por z
        creciente cuando el óvalo es convexo hacia +z.  El contorno
        cerrado se obtiene reflejando r → −r.
        (None, None) si los parámetros GOTS no definen un óvalo cerrado.
    """
    O  = float(params['O'])
    T  = float(params['T'])
    S  = float(params['S'])
    OG = float(params['OG'])
    zeta = float(params['zeta'])

    if abs(OG) < 1e-30:
        return None, None

    coef_rad = 2.0 * S - O * OG
    # Dominio de ρ: 1 + coef_rad·ρ² ≥ 0
    if coef_rad < -1e-15:
        rho_dom = np.sqrt(-1.0 / coef_rad)
    else:
        rho_dom = 1.0e4   # sin límite estricto

    # El óvalo real (rama «cercana al vértice», τ₋) se traza al variar ρ
    # desde 0 hasta ρ_max, donde r² = ρ² − τ₋² vuelve a anularse (vértice
    # trasero).  Primero hallamos ρ_max buscando el segundo cero de r(ρ).
    def _tau_menos(rho):
        rho2 = rho * rho
        radical = np.sqrt(max(1.0 + coef_rad * rho2, 0.0))
        denom = 1.0 + S * rho2 + radical
        if abs(denom) < 1e-30:
            return 0.0
        return (O + T * rho2) * rho2 / denom

    def _r2(rho):
        t = _tau_menos(rho)
        return rho * rho - t * t

    # Barrido grueso para encontrar ρ_max (segundo cero de r(ρ))
    M = 4000
    rho_busq = np.linspace(1e-6, rho_dom * 0.999999, M)
    r2_busq  = np.array([_r2(x) for x in rho_busq])
    sig      = np.sign(r2_busq)
    cambios  = np.where(np.diff(sig) != 0)[0]
    if len(cambios) == 0:
        # r² no cambia de signo: óvalo no cerrado con esta rama
        return None, None
    # primer cruce = vértice trasero del óvalo
    idx = cambios[0]
    # refinamos por bisección
    a, b = rho_busq[idx], rho_busq[idx + 1]
    for _ in range(60):
        c = 0.5 * (a + b)
        if _r2(c) > 0:
            a = c
        else:
            b = c
    rho_max = 0.5 * (a + b)

    # Muestreo Chebyshev denso en ambos extremos (vértice frontal y trasero)
    ks   = np.arange(N)
    cheb = 0.5 * (1.0 - np.cos(np.pi * ks / (N - 1)))
    rho  = rho_max * cheb

    rho2    = rho ** 2
    radical = np.sqrt(np.maximum(1.0 + coef_rad * rho2, 0.0))
    tau     = (O + T * rho2) * rho2 / (1.0 + S * rho2 + radical)
    r2      = np.maximum(rho2 - tau ** 2, 0.0)
    rs      = np.sqrt(r2)
    zs      = zeta + tau

    # forzar r=0 exacto en los vértices
    rs[0]  = 0.0
    rs[-1] = 0.0
    return zs, rs


def perfil_superficie(params, N=500, r_max=None):
    """Perfil meridional (r, z) de una superficie cartesiana.

    Devuelve únicamente la rama ascendente del óvalo (r creciente hasta
    el ecuador o hasta r_max, lo que sea menor).

    Args:
        params: dict GOTS (de ``calcular_gots`` o definido directamente)
        N     : número de puntos de muestreo en ρ
        r_max : apertura máxima (mismas unidades que los parámetros GOTS);
                None = sin límite explícito

    Returns:
        (r_arr, z_arr): arrays numpy 1D con r creciente desde 0
    """
    O, S, OG = params['O'], params['S'], params['OG']
    zeta     = params['zeta']

    # Dominio de ρ donde el radical es no negativo
    coef_rad = 2.0 * S - O * OG
    if coef_rad < -1e-15:
        rho_max_dominio = np.sqrt(-1.0 / coef_rad)
    else:
        rho_max_dominio = np.inf

    rho_lim = rho_max_dominio * 0.99 if np.isfinite(rho_max_dominio) else 500.0
    # Muestreo Chebyshev completo: nodos densos cerca de AMBOS extremos (ρ=0
    # en el vértice y ρ=ρ_lim cerca del ecuador/apertura), lo que mantiene
    # pequeña la desviación de las normales poligonales en todo el meridiano.
    i    = np.arange(N)
    cheb = 0.5 * (1.0 - np.cos(np.pi * i / (N - 1)))
    rho  = rho_lim * cheb

    tau  = _tau_de_rho(params, rho)
    z    = zeta + tau
    r    = np.sqrt(np.maximum(rho ** 2 - tau ** 2, 0.0))

    # Recortar a la rama ascendente (antes del ecuador del óvalo)
    dr       = np.diff(r)
    idx_desc = np.where(dr < -1e-12)[0]
    if len(idx_desc) > 0:
        corte = idx_desc[0] + 1
        r = r[:corte]
        z = z[:corte]

    # Recortar en r_max si se especifica
    if r_max is not None and r_max > 0:
        mask = r <= r_max + 1e-9
        r = r[mask]
        z = z[mask]

    return r, z


def encontrar_apertura(params0, params1, N=5000):
    """Radio de apertura donde los perfiles de dos superficies se interceptan.

    Encuentra el r donde z₀(r) = z₁(r) (los dos óvalos se cruzan en 3D).

    Args:
        params0, params1: dicts GOTS de las dos superficies
        N               : puntos de muestreo

    Returns:
        r_apertura (float, mismas unidades que los parámetros GOTS)
    """
    r0, z0 = perfil_superficie(params0, N=N)
    r1, z1 = perfil_superficie(params1, N=N)

    if len(r0) < 2 or len(r1) < 2:
        return 10.0

    r_max_comun = min(r0.max(), r1.max())
    if r_max_comun < 1e-6:
        return 10.0

    r_comun    = np.linspace(1e-3, r_max_comun * 0.999, 3000)
    z0_interp  = np.interp(r_comun, r0, z0)
    z1_interp  = np.interp(r_comun, r1, z1)
    diff       = z0_interp - z1_interp
    cruces     = np.where(np.diff(np.sign(diff)))[0]

    if len(cruces) > 0:
        idx   = cruces[0]
        denom = abs(diff[idx]) + abs(diff[idx + 1]) + 1e-30
        f     = abs(diff[idx]) / denom
        return float(r_comun[idx] + f * (r_comun[idx + 1] - r_comun[idx]))
    else:
        return float(r_max_comun * 0.9)


# ── Factor de forma σ → d₁ ────────────────────────────────────────────────


def calcular_d1_sigma(sigma, zeta_0, zeta_1, d_0, d_2, n_0, n_1, n_2):
    """Calcula d₁ a partir del factor de forma σ — Ecs. 39-43 de Silva-Lora.

    Args:
        sigma            : factor de forma (−1 = plano-convexo, 0 = biconvexo
                           simétrico, +1 = convexo-plano)
        zeta_0, zeta_1   : posiciones de los vértices de las dos superficies
        d_0              : posición del punto objeto
        d_2              : posición del punto imagen
        n_0, n_1, n_2    : índices de refracción

    Returns:
        d_1 (float): posición de la imagen intermedia dentro de la lente

    Raises:
        ValueError: si no existe solución real.
    """
    C2 = (n_2 * (d_0 - zeta_0) * (n_0 - n_1) * (sigma + 1)
          + n_0 * (d_2 - zeta_1) * (n_1 - n_2) * (sigma - 1))
    C1 = (n_1 * (d_0 - zeta_0) * (d_2 - zeta_1)
          * (n_2 * (sigma - 1) - n_0 * (sigma + 1) + 2.0 * n_1)
          + (zeta_0 + zeta_1)
          * (n_2 * (n_1 - n_0) * (sigma + 1) * (d_0 - zeta_0)
             + n_0 * (n_2 - n_1) * (sigma - 1) * (d_2 - zeta_1)))
    C0 = (n_1 * (d_0 - zeta_0) * (d_2 - zeta_1)
          * (zeta_0 * (n_0 - n_1) * (sigma + 1)
             + zeta_1 * (n_1 - n_2) * (sigma - 1))
          + zeta_0 * zeta_1
          * (n_2 * (d_0 - zeta_0) * (sigma + 1) * (n_0 - n_1)
             + n_0 * (d_2 - zeta_1) * (sigma - 1) * (n_1 - n_2)))

    disc = C1 ** 2 - 4.0 * C2 * C0
    if disc < 0:
        raise ValueError(
            f"No existe d₁ real para σ={sigma}. "
            "Pruebe con un valor de σ diferente o revise los parámetros del sistema."
        )

    sqrt_disc = np.sqrt(disc)
    d1_menos  = (-C1 - sqrt_disc) / (2.0 * C2)
    d1_mas    = (-C1 + sqrt_disc) / (2.0 * C2)

    lo, hi = min(zeta_0, zeta_1), max(zeta_0, zeta_1)
    dentro_mas   = lo < d1_mas   < hi
    dentro_menos = lo < d1_menos < hi

    # Seleccionar la raíz FUERA del intervalo [ζ₀, ζ₁]
    if not dentro_mas:
        return d1_mas
    elif not dentro_menos:
        return d1_menos
    else:
        # Ambas dentro (caso degenerado): usar la más alejada del intervalo
        dist_mas   = min(abs(d1_mas   - lo), abs(d1_mas   - hi))
        dist_menos = min(abs(d1_menos - lo), abs(d1_menos - hi))
        return d1_mas if dist_mas >= dist_menos else d1_menos


# ── Aplanatismo (capítulo 4 de Silva-Lora 2024) ───────────────────────────
#
# La condición de aplanatismo generalizada (Ec. 78/93) exige que
#
#     M(N, ρₖ) = ∏ₖ (n_{k+1}/n_k) · (1/(d_{k+1}−ζₖ) − 1/V̄C̄ₖ)
#                                   ─────────────────────────────  =  1
#                                   (1/(dₖ−ζₖ)     − 1/V̄C̄ₖ)
#
# para todos los rayos que acceden al sistema, donde V̄C̄ₖ es la distancia
# vértice → corte de la normal con el eje (Ec. 92), función de ρₖ.


def inversa_vc(sup, rho):
    """1/V̄C̄ₖ — Ec. 92 de Silva-Lora en forma numéricamente robusta.

    La forma publicada divide por G y por O, que degeneran en los tipos
    rigurosamente aplanáticos (tipo-2: O=0, G=∞; tipo-3: G=0, S=0).
    Usando la identidad T = S²/OG (válida idénticamente para las
    Ecs. 10–13) la expresión se reescribe sin esas divisiones:

        q(ρ) = [O + 2(2T − S·O)·ρ² / (1 + rad)] / rad ,
        rad  = √(1 + (2S − O·OG)·ρ²) ,

    que reproduce todos los límites: q(0) = O (curvatura paraxial),
    tipo-1 (S=T=0) → O/rad, tipo-2 (O=0) → 4Tρ²/((1+rad)·rad).

    Args:
        sup: dict GOTS (claves O, T, S, OG) o superficie esférica
             ({'esfera': True, 'O': 1/R})
        rho: distancia vértice-superficie del punto de incidencia
             (escalar o array)

    Returns:
        1/V̄C̄ (escalar o array según rho)
    """
    rho2 = np.asarray(rho, dtype=float) ** 2
    if sup.get('esfera'):
        # Esfera: la normal siempre pasa por el centro → V̄C̄ = R = cte.
        return np.broadcast_to(float(sup['O']), rho2.shape).copy() \
            if rho2.shape else float(sup['O'])
    O, T, S, OG = (float(sup['O']), float(sup['T']),
                   float(sup['S']), float(sup['OG']))
    rad = np.sqrt(np.maximum(1.0 + (2.0 * S - O * OG) * rho2, 0.0))
    return (O + 2.0 * (2.0 * T - S * O) * rho2 / (1.0 + rad)) / rad


def factor_aplanatismo(sup, rho):
    """Factor de la superficie k en la condición de aplanatismo (Ec. 78).

    Args:
        sup: dict de superficie con claves n_in, n_out, zeta, d_in, d_out
             y los parámetros GOTS (o 'esfera'/'O' para esferas);
             d_out = np.inf indica conjugado interno colimado (tipo-1).
        rho: distancia vértice-superficie del punto de incidencia

    Returns:
        factor (float o array)
    """
    q = inversa_vc(sup, rho)
    zeta = float(sup['zeta'])
    d_in, d_out = float(sup['d_in']), float(sup['d_out'])
    a = 0.0 if np.isinf(d_out) else 1.0 / (d_out - zeta)
    b = 0.0 if np.isinf(d_in)  else 1.0 / (d_in  - zeta)
    return (sup['n_out'] / sup['n_in']) * (a - q) / (b - q)


def condicion_aplanatismo(diseno, rhos):
    """M(N, ρₖ) — Ec. 78/93 — para un diseño de ``disenar_aplanatica``.

    Args:
        diseno: dict devuelto por ``disenar_aplanatica``
        rhos  : secuencia con el ρₖ del rayo en cada superficie (en el
                mismo orden que diseno['superficies'])

    Returns:
        M (float o array): el sistema es aplanático si M = 1 ∀ rayos.
    """
    M = 1.0
    for sup, rho in zip(diseno['superficies'], rhos):
        M = M * factor_aplanatismo(sup, rho)
    return M


def rms_aplanatismo(valores_M):
    """(M − 1)_RMS — Ec. 94 — sobre un conjunto de rayos."""
    v = np.asarray(valores_M, dtype=float)
    return float(np.sqrt(np.mean(np.abs(v - 1.0) ** 2)))


def _gots_limite_O_cero(n_k, n_k1, zeta_k, d_k):
    """Parámetros GOTS exactos en el límite O = 0 (superficie tipo-2).

    Con la condición (d_{k+1}−ζ)/n_{k+1} = (d_k−ζ)/n_k (Ec. 111) las
    Ecs. 10–13 degeneran (O = 0, G = ∞) pero T, S y OG tienen límites
    finitos que aquí se evalúan en forma cerrada.
    """
    xi = d_k - zeta_k
    T  = n_k * (n_k1 + n_k) / (4.0 * n_k1 ** 2 * xi ** 3)
    S  = (n_k1 ** 3 - n_k ** 3) / (2.0 * n_k1 ** 2 * xi ** 2 * (n_k1 - n_k))
    OG = ((n_k1 ** 3 - n_k ** 3) ** 2
          / (n_k * n_k1 ** 2 * xi * (n_k1 - n_k) * (n_k1 ** 2 - n_k ** 2)))
    return {'G': float('inf'), 'O': 0.0, 'T': T, 'S': S, 'OG': OG,
            'zeta': zeta_k}


def disenar_aplanatica(tipo, n0, n1, zeta_0, zeta_1, d_0,
                       n2=None, rama='real'):
    """Diseña una lente singlete rigurosamente aplanática (§4.4 de la tesis).

    Los cuatro tipos hacen M(ρₖ) ≡ 1 en todo rigor (sistemas libres de
    aberración esférica Y de coma simultáneamente):

      tipo 0: superficies esféricas con conjugados en los puntos de
              Young (2Sₖ = GₖOₖ², Ecs. 95-96). Además anastigmática.
      tipo 1: interior colimado (d₁ → ∞) → superficies cónicas idénticas
              (S=T=0, Ec. 108); exige n₂=n₀ y |d₂−ζ₁| = |d₀−ζ₀|
              (Ec. 107): rama='real' → biconvexa e imagen real;
              rama='virtual' → menisco e imagen virtual.
      tipo 2: Oₖ = 0 → ovoides de vértice plano (Ec. 113);
              d_{k+1} = ζₖ + (n_{k+1}/nₖ)(dₖ − ζₖ) (Ec. 111).
      tipo 3: Gₖ = 0 (Ec. 122);
              d_{k+1} = ζₖ + (nₖ/n_{k+1})²(dₖ − ζₖ) (Ec. 121).

    Args:
        tipo   : 0, 1, 2 o 3
        n0, n1 : índices del espacio objeto y de la lente
        zeta_0, zeta_1: vértices de las dos superficies (ζ₁ > ζ₀)
        d_0    : posición axial del punto objeto (d₀ < ζ₀ para objeto real)
        n2     : índice del espacio imagen. Los tipos 1-3 exigen n₂ = n₀
                 (se fuerza); el tipo 0 lo admite arbitrario (default n₀).
        rama   : sólo tipo 1: 'real' (biconvexa) o 'virtual' (menisco)

    Returns:
        dict con claves:
          tipo, n (n0,n1,n2), zetas (ζ₀,ζ₁), d (d₀,d₁,d₂),
          magnificacion (aumento lateral M = ∏gₖ, Ec. 89),
          superficies: [sup0, sup1] — cada sup con n_in, n_out, zeta,
          d_in, d_out y parámetros GOTS ('esfera'+'R' para tipo 0).

    Raises:
        ValueError: parámetros degenerados o tipo desconocido.
    """
    tipo = int(tipo)
    n0, n1 = float(n0), float(n1)
    zeta_0, zeta_1, d_0 = float(zeta_0), float(zeta_1), float(d_0)

    if abs(n1 - n0) < 1e-12:
        raise ValueError("Se requiere n1 ≠ n0 (sin salto de índice no hay "
                         "refracción).")
    if abs(d_0 - zeta_0) < 1e-9:
        raise ValueError("El objeto no puede coincidir con el vértice ζ₀.")
    if zeta_1 <= zeta_0:
        raise ValueError("Se requiere ζ₁ > ζ₀.")

    if tipo == 0:
        n2 = n0 if n2 is None else float(n2)
        # Puntos de Young (Ecs. 95-96): dₖ−ζ = R(nₖ+n_{k+1})/nₖ
        R0  = n0 * (d_0 - zeta_0) / (n0 + n1)
        d_1 = zeta_0 + (n0 / n1) * (d_0 - zeta_0)
        if abs(d_1 - zeta_1) < 1e-9:
            raise ValueError("Conjugado intermedio sobre el vértice ζ₁; "
                             "modifique ζ₁ o d₀.")
        R1  = n1 * (d_1 - zeta_1) / (n1 + n2)
        d_2 = zeta_1 + (n1 / n2) * (d_1 - zeta_1)
        sup0 = {'esfera': True, 'R': R0, 'O': 1.0 / R0, 'zeta': zeta_0,
                'n_in': n0, 'n_out': n1, 'd_in': d_0, 'd_out': d_1}
        sup1 = {'esfera': True, 'R': R1, 'O': 1.0 / R1, 'zeta': zeta_1,
                'n_in': n1, 'n_out': n2, 'd_in': d_1, 'd_out': d_2}

    elif tipo == 1:
        n2 = n0                     # exigido por M = 1 (Ec. 106)
        d_1 = float('inf')
        signo = -1.0 if rama == 'real' else +1.0     # Ec. 107
        d_2 = zeta_1 + signo * (d_0 - zeta_0)
        p0 = {'G': -(n1 / n0) ** 2,
              'O': -n0 / ((d_0 - zeta_0) * (n1 - n0)),
              'T': 0.0, 'S': 0.0, 'zeta': zeta_0}
        p0['OG'] = p0['O'] * p0['G']
        p1 = {'G': -(n1 / n2) ** 2,
              'O': -n2 / ((d_2 - zeta_1) * (n1 - n2)),
              'T': 0.0, 'S': 0.0, 'zeta': zeta_1}
        p1['OG'] = p1['O'] * p1['G']
        sup0 = dict(p0, n_in=n0, n_out=n1, d_in=d_0, d_out=d_1)
        sup1 = dict(p1, n_in=n1, n_out=n2, d_in=d_1, d_out=d_2)

    elif tipo == 2:
        n2 = n0                     # exigido por M = n_N/n₀ = 1 (Ec. 112)
        d_1 = zeta_0 + (n1 / n0) * (d_0 - zeta_0)    # Ec. 111/114
        d_2 = zeta_1 + (n2 / n1) * (d_1 - zeta_1)    # Ec. 111/115
        p0 = _gots_limite_O_cero(n0, n1, zeta_0, d_0)
        p1 = _gots_limite_O_cero(n1, n2, zeta_1, d_1)
        sup0 = dict(p0, n_in=n0, n_out=n1, d_in=d_0, d_out=d_1)
        sup1 = dict(p1, n_in=n1, n_out=n2, d_in=d_1, d_out=d_2)

    elif tipo == 3:
        n2 = n0                     # exigido por M = n_N/n₀ = 1 (Ec. 112)
        d_1 = zeta_0 + (n0 / n1) ** 2 * (d_0 - zeta_0)    # Ec. 121
        d_2 = zeta_1 + (n1 / n2) ** 2 * (d_1 - zeta_1)    # Ec. 121
        p0 = calcular_gots(n0, n1, zeta_0, d_0, d_1)
        p1 = calcular_gots(n1, n2, zeta_1, d_1, d_2)
        sup0 = dict(p0, n_in=n0, n_out=n1, d_in=d_0, d_out=d_1)
        sup1 = dict(p1, n_in=n1, n_out=n2, d_in=d_1, d_out=d_2)

    else:
        raise ValueError(f"Tipo aplanático desconocido: {tipo} (use 0-3).")

    # Aumento lateral M = ∏gₖ (Ec. 88-89), robusto para d₁ = ∞
    if np.isinf(d_1):
        M_lat = (n0 / n2) * (d_2 - zeta_1) / (d_0 - zeta_0)
    else:
        M_lat = ((n0 / n2) * (d_1 - zeta_0) * (d_2 - zeta_1)
                 / ((d_0 - zeta_0) * (d_1 - zeta_1)))

    return {
        'tipo': tipo,
        'n': (n0, n1, n2),
        'zetas': (zeta_0, zeta_1),
        'd': (d_0, d_1, d_2),
        'magnificacion': M_lat,
        'superficies': [sup0, sup1],
    }


def perfil_superficie_aplanatica(sup, r_max, N=500):
    """Perfil meridional (r, z) de una superficie aplanática hasta r_max.

    A diferencia de ``perfil_superficie``, el dominio de ρ se ajusta a la
    apertura pedida (las superficies tipo 1-3 tienen dominio ilimitado y
    un muestreo con tope fijo dejaría casi todos los nodos fuera de la
    apertura). Las esferas (tipo 0) se muestrean por ángulo (exacto).

    Args:
        sup  : superficie de un diseño de ``disenar_aplanatica``
        r_max: radio de apertura
        N    : número de puntos

    Returns:
        (r_arr, z_arr) con r creciente desde 0 hasta ≈ r_max
    """
    r_max = float(r_max)
    i    = np.arange(N)
    cheb = 0.5 * (1.0 - np.cos(np.pi * i / (N - 1)))    # denso en extremos

    if sup.get('esfera'):
        R = float(sup['R'])
        if r_max >= abs(R):
            r_max = 0.999 * abs(R)
        th_max = np.arcsin(r_max / abs(R))
        th = th_max * cheb
        r  = np.abs(R) * np.sin(th)
        z  = sup['zeta'] + R * (1.0 - np.cos(th))
        return r, z

    # Dominio admisible del radical
    coef_rad = 2.0 * sup['S'] - sup['O'] * sup['OG']
    rho_dom  = (np.sqrt(-1.0 / coef_rad) * 0.999999
                if coef_rad < -1e-15 else np.inf)

    # ρ_hi tal que r(ρ_hi) ≥ r_max: iteración de punto fijo
    # ρ² = r² + τ(ρ)² (converge en pocas iteraciones; τ varía suavemente)
    rho_hi = r_max
    for _ in range(80):
        rho_hi = min(rho_hi, rho_dom)
        tau    = float(_tau_de_rho(sup, rho_hi))
        nuevo  = np.sqrt(r_max ** 2 + tau ** 2)
        nuevo  = min(nuevo, rho_dom)
        if abs(nuevo - rho_hi) < 1e-12 * max(1.0, rho_hi):
            rho_hi = nuevo
            break
        rho_hi = nuevo

    rho  = rho_hi * cheb
    tau  = _tau_de_rho(sup, rho)
    z    = sup['zeta'] + tau
    r    = np.sqrt(np.maximum(rho ** 2 - tau ** 2, 0.0))

    # Rama ascendente y recorte en r_max
    dr       = np.diff(r)
    idx_desc = np.where(dr < -1e-12)[0]
    if len(idx_desc) > 0:
        corte = idx_desc[0] + 1
        r, z = r[:corte], z[:corte]
    mask = r <= r_max * (1.0 + 1e-9)
    return r[mask], z[mask]


def alturas_estigmaticas(diseno, rho_0):
    """ρ₁ del rayo que incide en la superficie 0 con parámetro ρ₀.

    Encadena la geometría estigmática exacta: el rayo refractado por la
    superficie 0 pasa (él o su prolongación) por el conjugado intermedio
    (d₁, 0), lo que determina el punto de incidencia sobre la
    superficie 1 sin necesidad de un trazador de rayos. Para d₁ = ∞
    (tipo-1) el rayo interno viaja paralelo al eje a la altura r₀.

    Args:
        diseno: dict de ``disenar_aplanatica``
        rho_0 : parámetro ρ del punto de incidencia sobre la superficie 0

    Returns:
        (rho_0, rho_1): parámetros de incidencia en ambas superficies
    """
    sup0, sup1 = diseno['superficies']
    d_1 = diseno['d'][1]

    def _punto(sup, rho):
        if sup.get('esfera'):
            R = float(sup['R'])
            # ρ = distancia vértice-punto (cuerda): ρ = 2|R|·sin(θ/2)
            th = 2.0 * np.arcsin(min(rho / (2.0 * abs(R)), 1.0))
            return sup['zeta'] + R * (1.0 - np.cos(th)), abs(R) * np.sin(th)
        tau = float(_tau_de_rho(sup, rho))
        return sup['zeta'] + tau, np.sqrt(max(rho * rho - tau * tau, 0.0))

    z0, r0 = _punto(sup0, rho_0)

    if np.isinf(d_1):
        # Rayo interno paralelo al eje: buscar ρ₁ con r₁(ρ₁) = r₀
        f = lambda rho: _punto(sup1, rho)[1] - r0
    else:
        # P₁ alineado con (d₁,0) y P₀: cruce del producto vectorial
        f = lambda rho: ((_punto(sup1, rho)[0] - d_1) * r0
                         - _punto(sup1, rho)[1] * (z0 - d_1))

    # Bracket por barrido geométrico + bisección
    lo, f_lo = 0.0, f(0.0)
    hi = None
    paso = max(rho_0, 1e-3)
    x = paso * 0.05
    for _ in range(400):
        f_x = f(x)
        if f_x == 0.0:
            return rho_0, x
        if np.sign(f_x) != np.sign(f_lo):
            hi = x
            break
        lo, f_lo = x, f_x
        x *= 1.15
    if hi is None:
        raise ValueError("No se encontró la intersección con la segunda "
                         "superficie (¿apertura excesiva?).")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if np.sign(f(mid)) == np.sign(f_lo):
            lo = mid
        else:
            hi = mid
    return rho_0, 0.5 * (lo + hi)


# ── Conversión de perfil a path SVG ───────────────────────────────────────


def puntos_a_bezier_path_str(puntos_xy, cerrar=True):
    """Convierte una secuencia de puntos (x, y) en un path SVG de curvas
    cúbicas Bézier (`M ... C ... C ... Z`) que interpola la curva suave
    que los une.

    Las tangentes se estiman con la fórmula de Bessel para nodos NO
    uniformes (mezcla de las pendientes de cuerda ponderada por las
    longitudes de cuerda adyacentes) y los handles Bézier del segmento
    [pᵢ, pⱼ] son  P₁ = pᵢ + mᵢ·hᵢ/3  y  P₂ = pⱼ − mⱼ·hᵢ/3, que es la
    conversión exacta spline→Bézier.  Esto da error O(h⁴) real incluso
    con muestreo Chebyshev (muy no uniforme en los extremos), a
    diferencia de la versión anterior con diferencias centradas y
    tensión fija, cuyos handles resultaban el doble de largos y
    producían un abombamiento O(h²) entre nodos (y con él un
    desplazamiento medible del foco).

    Args:
        puntos_xy: lista o array (N, 2)
        cerrar   : si True añade 'Z' y usa tangentes periódicas

    Returns:
        str con el path SVG
    """
    import math

    pts = [(float(x), float(y)) for x, y in puntos_xy]

    # En un contorno cerrado el punto inicial no debe repetirse al final:
    # la cuerda de cierre tendría longitud cero y su tangente de Bessel
    # degeneraría (pliegue microscópico en el vértice que desvía los
    # rayos axiales).
    if cerrar and len(pts) > 2:
        while len(pts) > 2 and math.hypot(pts[-1][0] - pts[0][0],
                                          pts[-1][1] - pts[0][1]) < 1e-12:
            pts.pop()

    # Decimación: eliminar puntos casi coincidentes.  El muestreo
    # Chebyshev produce segmentos de <1e-4 unidades en los extremos;
    # esos micro-segmentos no aportan geometría (la curvatura la
    # capturan los handles) y desestabilizan el solver de intersección
    # del trazador (tolerancias absolutas ~1e-8 sobre coeficientes que
    # escalan con h³).
    if len(pts) > 2:
        L_tot = sum(math.hypot(b[0] - a[0], b[1] - a[1])
                    for a, b in zip(pts[:-1], pts[1:]))
        umbral = max(L_tot * 1e-4, 1e-9)
        filtrados = [pts[0]]
        for p in pts[1:-1]:
            q = filtrados[-1]
            if math.hypot(p[0] - q[0], p[1] - q[1]) >= umbral:
                filtrados.append(p)
        # conservar el último punto siempre (extremo de la rama)
        q = filtrados[-1]
        if math.hypot(pts[-1][0] - q[0], pts[-1][1] - q[1]) < umbral:
            filtrados[-1] = pts[-1]
        else:
            filtrados.append(pts[-1])
        # en contornos cerrados, evitar que la cuerda de cierre quede
        # microscópica (desestabilizaría la tangente de Bessel del
        # punto inicial)
        if cerrar:
            while len(filtrados) > 2 and math.hypot(
                    filtrados[-1][0] - filtrados[0][0],
                    filtrados[-1][1] - filtrados[0][1]) < umbral:
                filtrados.pop()
        pts = filtrados

    n = len(pts)
    if n == 0:
        return ""
    if n == 1:
        return f"M {pts[0][0]:.9f},{pts[0][1]:.9f}"
    if n == 2:
        return (f"M {pts[0][0]:.9f},{pts[0][1]:.9f} "
                f"L {pts[1][0]:.9f},{pts[1][1]:.9f}"
                + (" Z" if cerrar else ""))

    # Cuerdas h_i = |p_{i+1} − p_i| (con envolvente si es cerrada) y
    # pendientes de cuerda d_i = (p_{i+1} − p_i)/h_i.
    n_seg = n if cerrar else n - 1
    h = [0.0] * n_seg
    d = [(0.0, 0.0)] * n_seg
    for i in range(n_seg):
        j = (i + 1) % n
        dx = pts[j][0] - pts[i][0]
        dy = pts[j][1] - pts[i][1]
        h[i] = math.hypot(dx, dy)
        if h[i] > 1e-30:
            d[i] = (dx / h[i], dy / h[i])

    # Tangentes de Bessel: m_i = (h_i·d_{i−1} + h_{i−1}·d_i)/(h_{i−1}+h_i)
    m = [(0.0, 0.0)] * n
    for i in range(n):
        if cerrar:
            i0 = (i - 1) % n_seg
            i1 = i % n_seg
        else:
            if i == 0 or i == n - 1:
                continue        # extremos: después
            i0 = i - 1
            i1 = i
        hs = h[i0] + h[i1]
        if hs < 1e-30:
            continue
        m[i] = ((h[i1] * d[i0][0] + h[i0] * d[i1][0]) / hs,
                (h[i1] * d[i0][1] + h[i0] * d[i1][1]) / hs)

    if not cerrar:
        # Extremos abiertos: extrapolación parabólica m₀ = 2d₀ − m₁
        m[0] = (2.0 * d[0][0] - m[1][0], 2.0 * d[0][1] - m[1][1])
        m[n - 1] = (2.0 * d[n_seg - 1][0] - m[n - 2][0],
                    2.0 * d[n_seg - 1][1] - m[n - 2][1])

    # 9 decimales: con muestreo denso los segmentos pueden medir ~1e-4
    # unidades y el redondeo de las coordenadas contaminaría las
    # tangentes (y con ellas las normales de refracción).
    partes = [f"M {pts[0][0]:.9f},{pts[0][1]:.9f}"]
    for i in range(n_seg):
        j = (i + 1) % n
        p0x, p0y = pts[i]
        p3x, p3y = pts[j]
        hi3 = h[i] / 3.0
        p1x = p0x + m[i][0] * hi3
        p1y = p0y + m[i][1] * hi3
        p2x = p3x - m[j][0] * hi3
        p2y = p3y - m[j][1] * hi3
        partes.append(
            f"C {p1x:.9f},{p1y:.9f} {p2x:.9f},{p2y:.9f} "
            f"{p3x:.9f},{p3y:.9f}"
        )
    if cerrar:
        partes.append("Z")
    return " ".join(partes)


def perfil_a_path_str(puntos_xy, cerrar=True):
    """Convierte lista de (x, y) a cadena de path SVG (M ... L ... Z).

    Args:
        puntos_xy: lista o array de forma (N, 2)
        cerrar   : si True, añade 'Z' al final

    Returns:
        str con el path SVG
    """
    if len(puntos_xy) == 0:
        return ""
    x0, y0 = puntos_xy[0]
    partes  = [f"M {x0:.5f},{y0:.5f}"]
    for x, y in puntos_xy[1:]:
        partes.append(f"L {x:.5f},{y:.5f}")
    if cerrar:
        partes.append("Z")
    return " ".join(partes)
