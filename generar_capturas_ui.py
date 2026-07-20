"""Genera imágenes de referencia de los diálogos de las extensiones.

Las imágenes NO son capturas de pantalla: se dibujan a partir de los
propios archivos ``.inx``, que son la definición canónica de la interfaz
que Inkscape construye.  De ese modo la documentación no puede quedar
desincronizada con la extensión: si se añade un parámetro, basta con
volver a ejecutar este script.

Salida: docs/img/ui-<extension>-<pestaña>.svg
"""

import os
import xml.etree.ElementTree as ET

AQUI = os.path.dirname(os.path.abspath(__file__))
INX_DIR = os.path.join(AQUI, "inkscape-raytracing", "inkscape_raytracing")
IMG_DIR = os.path.join(AQUI, "docs", "img")

NS = {"inx": "http://www.inkscape.org/namespace/inkscape/extension"}

# ── Métricas del diálogo ────────────────────────────────────────────────────
ANCHO      = 700
MARGEN     = 16
ALTO_TIT   = 30
ALTO_TAB   = 28
ALTO_FILA  = 26
ALTO_CAB   = 26
ALTO_SEP   = 12
ALTO_TXT   = 15
ALTO_PIE   = 44
CTRL_ANCHO = 215

C_FONDO   = "#f6f5f4"
C_BORDE   = "#c8c5c2"
C_TIT     = "#e1dedb"
C_TEXTO   = "#2e3436"
C_SUAVE   = "#6b6b6b"
C_CAMPO   = "#ffffff"
C_ACENTO  = "#3584e4"


def esc(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def ancho_texto(s, size):
    return len(s) * size * 0.52


def envolver(texto, size, ancho_max):
    """Parte un texto en líneas que quepan en ancho_max."""
    lineas = []
    for parrafo in texto.split("\n"):
        palabras = parrafo.split()
        if not palabras:
            lineas.append("")
            continue
        actual = palabras[0]
        for w in palabras[1:]:
            if ancho_texto(actual + " " + w, size) <= ancho_max:
                actual += " " + w
            else:
                lineas.append(actual)
                actual = w
        lineas.append(actual)
    return lineas


# ── Lectura del .inx ────────────────────────────────────────────────────────

def leer_inx(ruta):
    """Devuelve (nombre, [(pestaña, [elementos])]).

    Cada elemento es (clase, datos) con clase en:
    'cabecera', 'texto', 'separador', 'campo', 'casilla', 'lista'.
    """
    arbol = ET.parse(ruta)
    raiz = arbol.getroot()
    nombre = raiz.find("inx:name", NS).text

    def elementos_de(nodo):
        salida = []
        for hijo in nodo:
            etiqueta = hijo.tag.split("}")[-1]
            if etiqueta == "label":
                texto = (hijo.text or "").strip()
                if hijo.get("appearance") == "header":
                    salida.append(("cabecera", texto))
                elif texto:
                    salida.append(("texto", texto))
            elif etiqueta == "separator":
                salida.append(("separador", None))
            elif etiqueta == "param":
                tipo = hijo.get("type")
                gui = (hijo.get("gui-text") or hijo.get("name") or "").strip()
                if tipo in ("float", "int", "string"):
                    salida.append(("campo", (gui, (hijo.text or "").strip())))
                elif tipo == "boolean":
                    val = (hijo.text or "").strip().lower() == "true"
                    salida.append(("casilla", (gui, val)))
                elif tipo in ("enum", "optiongroup"):
                    opciones = [(o.text or "").strip() for o in hijo]
                    salida.append(("lista", (gui, opciones)))
        return salida

    notebook = raiz.find(".//inx:param[@type='notebook']", NS)
    if notebook is None:
        return nombre, [(None, elementos_de(raiz))]

    paginas = []
    for pagina in notebook.findall("inx:page", NS):
        paginas.append((pagina.get("gui-text"), elementos_de(pagina)))
    return nombre, paginas


# ── Dibujo ──────────────────────────────────────────────────────────────────

def _texto(x, y, s, size=12, color=C_TEXTO, peso="normal", anchor="start"):
    return (f'  <text x="{x:.1f}" y="{y:.1f}" font-family="Cantarell, '
            f'DejaVu Sans, sans-serif" font-size="{size}" fill="{color}" '
            f'font-weight="{peso}" text-anchor="{anchor}">{esc(s)}</text>\n')


def _rect(x, y, w, h, fill, stroke=None, rx=0, sw=1):
    s = f' stroke="{stroke}" stroke-width="{sw}"' if stroke else ""
    return (f'  <rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" '
            f'height="{h:.1f}" rx="{rx}" fill="{fill}"{s}/>\n')


def alto_elementos(elementos):
    alto = 0
    for clase, datos in elementos:
        if clase == "cabecera":
            alto += ALTO_CAB
        elif clase == "texto":
            alto += ALTO_TXT * len(envolver(datos, 10.5, ANCHO - 2 * MARGEN)) + 6
        elif clase == "separador":
            alto += ALTO_SEP
        else:
            alto += ALTO_FILA
    return alto


def dibujar(nombre, pestanas, indice_activa, salida):
    pestana_activa, elementos = pestanas[indice_activa]
    alto_cont = alto_elementos(elementos) + 2 * MARGEN
    alto_tabs = ALTO_TAB if pestana_activa else 0
    alto = ALTO_TIT + alto_tabs + alto_cont + ALTO_PIE

    out = (f'<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
           f'<svg xmlns="http://www.w3.org/2000/svg" width="{ANCHO}" '
           f'height="{alto}" viewBox="0 0 {ANCHO} {alto}">\n')

    # marco y barra de título
    out += _rect(0.5, 0.5, ANCHO - 1, alto - 1, C_FONDO, C_BORDE, rx=8)
    out += _rect(0.5, 0.5, ANCHO - 1, ALTO_TIT, C_TIT, C_BORDE, rx=8)
    out += _rect(0.5, ALTO_TIT - 8, ANCHO - 1, 8, C_TIT)
    out += _texto(ANCHO / 2, ALTO_TIT / 2 + 4.5, nombre, 12.5, C_TEXTO,
                  "bold", "middle")
    for i, cx in enumerate((ANCHO - 22,)):
        out += (f'  <path d="M {cx - 5},{ALTO_TIT / 2 - 5} l 10,10 M '
                f'{cx + 5},{ALTO_TIT / 2 - 5} l -10,10" stroke="{C_SUAVE}" '
                f'stroke-width="1.4" fill="none"/>\n')

    y = ALTO_TIT

    # pestañas
    if pestana_activa:
        x = MARGEN
        for i, (nom, _) in enumerate(pestanas):
            w = ancho_texto(nom, 11.5) + 26
            activa = (i == indice_activa)
            out += _rect(x, y + 3, w, ALTO_TAB - 3,
                         "#ffffff" if activa else C_TIT,
                         C_BORDE, rx=5)
            if activa:
                out += _rect(x + 6, y + 3, w - 12, 2.5, C_ACENTO)
            out += _texto(x + w / 2, y + ALTO_TAB / 2 + 5, nom, 11.5,
                          C_TEXTO if activa else C_SUAVE,
                          "bold" if activa else "normal", "middle")
            x += w + 4
        out += (f'  <line x1="0.5" y1="{y + ALTO_TAB:.1f}" x2="{ANCHO - 0.5}" '
                f'y2="{y + ALTO_TAB:.1f}" stroke="{C_BORDE}"/>\n')
        y += ALTO_TAB

    # contenido
    y += MARGEN
    for clase, datos in elementos:
        if clase == "cabecera":
            out += _texto(MARGEN, y + 15, datos, 11.5, C_TEXTO, "bold")
            y += ALTO_CAB
        elif clase == "texto":
            for linea in envolver(datos, 10.5, ANCHO - 2 * MARGEN):
                out += _texto(MARGEN, y + 11, linea, 10.5, C_SUAVE)
                y += ALTO_TXT
            y += 6
        elif clase == "separador":
            out += (f'  <line x1="{MARGEN}" y1="{y + ALTO_SEP / 2:.1f}" '
                    f'x2="{ANCHO - MARGEN}" y2="{y + ALTO_SEP / 2:.1f}" '
                    f'stroke="{C_BORDE}"/>\n')
            y += ALTO_SEP
        elif clase == "campo":
            gui, valor = datos
            out += _texto(MARGEN, y + 16, gui, 11)
            cx = ANCHO - MARGEN - CTRL_ANCHO
            out += _rect(cx, y + 3, CTRL_ANCHO, 20, C_CAMPO, C_BORDE, rx=4)
            out += _texto(cx + 8, y + 17, valor, 11)
            out += (f'  <path d="M {ANCHO - MARGEN - 16},{y + 10} l 4,-4 l 4,4 '
                    f'M {ANCHO - MARGEN - 16},{y + 16} l 4,4 l 4,-4" '
                    f'stroke="{C_SUAVE}" stroke-width="1.2" fill="none"/>\n')
            y += ALTO_FILA
        elif clase == "casilla":
            gui, val = datos
            out += _rect(MARGEN, y + 5, 15, 15,
                         C_ACENTO if val else C_CAMPO, C_BORDE, rx=3)
            if val:
                out += (f'  <path d="M {MARGEN + 3.5},{y + 12.5} l 3,3 l 5,-6" '
                        f'stroke="#ffffff" stroke-width="2" fill="none"/>\n')
            out += _texto(MARGEN + 24, y + 17, gui, 11)
            y += ALTO_FILA
        elif clase == "lista":
            gui, opciones = datos
            out += _texto(MARGEN, y + 16, gui, 11)
            cx = ANCHO - MARGEN - CTRL_ANCHO
            out += _rect(cx, y + 3, CTRL_ANCHO, 20, C_CAMPO, C_BORDE, rx=4)
            elegido = opciones[0] if opciones else ""
            if ancho_texto(elegido, 11) > CTRL_ANCHO - 30:
                while elegido and ancho_texto(elegido + "…", 11) > CTRL_ANCHO - 30:
                    elegido = elegido[:-1]
                elegido += "…"
            out += _texto(cx + 8, y + 17, elegido, 11)
            out += (f'  <path d="M {ANCHO - MARGEN - 18},{y + 11} l 5,5 l 5,-5" '
                    f'stroke="{C_SUAVE}" stroke-width="1.3" fill="none"/>\n')
            y += ALTO_FILA

    # pie: vista previa y botones
    y_pie = alto - ALTO_PIE
    out += (f'  <line x1="0.5" y1="{y_pie:.1f}" x2="{ANCHO - 0.5}" '
            f'y2="{y_pie:.1f}" stroke="{C_BORDE}"/>\n')
    out += _rect(MARGEN, y_pie + 14, 15, 15, C_CAMPO, C_BORDE, rx=3)
    out += _texto(MARGEN + 24, y_pie + 26, "Live preview", 11, C_SUAVE)
    for etiqueta, dx, principal in (("Apply", 168, True), ("Close", 84, False)):
        bx = ANCHO - MARGEN - dx
        out += _rect(bx, y_pie + 11, 76, 22,
                     C_ACENTO if principal else C_CAMPO, C_BORDE, rx=5)
        out += _texto(bx + 38, y_pie + 26, etiqueta, 11,
                      "#ffffff" if principal else C_TEXTO, "bold", "middle")

    out += "</svg>\n"
    with open(salida, "w") as f:
        f.write(out)
    return salida


EXTENSIONES = [
    ("render.inx",               "raytracing"),
    ("set_material.inx",         "set-material"),
    ("superficie_cartesiana.inx", "superficie"),
    ("lente_ovoide.inx",         "lsoe"),
    ("lente_aplanatica.inx",     "aplanatica"),
]


def main():
    os.makedirs(IMG_DIR, exist_ok=True)
    for archivo, mote in EXTENSIONES:
        nombre, pestanas = leer_inx(os.path.join(INX_DIR, archivo))
        for i, (pest, _) in enumerate(pestanas):
            sufijo = f"-{i + 1}" if pest else ""
            ruta = os.path.join(IMG_DIR, f"ui-{mote}{sufijo}.svg")
            dibujar(nombre, pestanas, i, ruta)
            print(f"{os.path.relpath(ruta, AQUI)}"
                  + (f"   [{pest}]" if pest else ""))


if __name__ == "__main__":
    main()
