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

        n₀·|P−F₀| + n₁·|P−F₁| = k ,

    donde F₀=(d_obj,0), F₁=(d_img,0) y k se fija pasando por el vértice
    (ζ, 0). Esta es la cuártica exacta del óvalo refractante.

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

    k = n0 * abs(zeta - d_obj) + n1 * abs(zeta - d_img)

    # Raíces axiales: z tal que n₀|z−d_obj| + n₁|z−d_img| = k.
    # Resolución por regímenes de signos.
    raices = []
    for s0 in (+1.0, -1.0):
        for s1 in (+1.0, -1.0):
            denom = n0 * s0 + n1 * s1
            if abs(denom) < 1e-30:
                continue
            z_c = (k + n0 * s0 * d_obj + n1 * s1 * d_img) / denom
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

        # Elegir la raíz ≥ 0 que satisface n₀√(a²+u) + n₁√(b²+u) = k
        mejor_u  = 0.0
        mejor_er = float('inf')
        for u in (u1, u2):
            if u < -1e-9:
                continue
            uu = max(u, 0.0)
            err = abs(n0 * np.sqrt(a * a + uu)
                      + n1 * np.sqrt(b * b + uu) - k)
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


# ── Conversión de perfil a path SVG ───────────────────────────────────────


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
