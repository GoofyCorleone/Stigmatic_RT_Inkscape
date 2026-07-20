#!/usr/bin/env bash
#
# instalar.sh — copia las extensiones de este repo a la carpeta de extensiones
# de usuario de Inkscape y así Inkscape las ejecuta.
#
# Recuerda: editar los archivos del repo NO actualiza Inkscape. Inkscape solo
# lee las extensiones desde su carpeta de usuario, y solo al arrancar. Ejecuta
# este script cada vez que cambies el código y REINICIA Inkscape después.
#
# Uso:
#   ./instalar.sh                      # detecta el sistema y copia
#   INKSCAPE_EXT=/ruta ./instalar.sh   # fuerza una carpeta de destino concreta
#
set -euo pipefail

# Carpeta de este script (raíz del repo), para poder ejecutarlo desde cualquier sitio.
DIR_REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ORIGEN="$DIR_REPO/inkscape-raytracing/inkscape_raytracing"

if [[ ! -d "$ORIGEN" ]]; then
    echo "ERROR: no encuentro el código de la extensión en:" >&2
    echo "  $ORIGEN" >&2
    echo "¿Estás ejecutando el script desde dentro del repositorio?" >&2
    exit 1
fi

# --- Determinar la carpeta de extensiones de Inkscape ---
if [[ -n "${INKSCAPE_EXT:-}" ]]; then
    DESTINO="$INKSCAPE_EXT"
else
    case "$(uname -s)" in
        Darwin)
            DESTINO="$HOME/Library/Application Support/org.inkscape.Inkscape/config/inkscape/extensions"
            ;;
        Linux)
            DESTINO="${XDG_CONFIG_HOME:-$HOME/.config}/inkscape/extensions"
            ;;
        *)
            echo "ERROR: sistema no reconocido ($(uname -s))." >&2
            echo "Indica la carpeta a mano:  INKSCAPE_EXT=/ruta ./instalar.sh" >&2
            echo "En Windows usa la instrucción de PowerShell del README." >&2
            exit 1
            ;;
    esac
fi

echo "Origen:  $ORIGEN"
echo "Destino: $DESTINO"
echo

mkdir -p "$DESTINO"

# --- Copiar (excluyendo cachés de Python) ---
if command -v rsync >/dev/null 2>&1; then
    rsync -a --exclude='__pycache__' --exclude='*.pyc' "$ORIGEN"/ "$DESTINO"/
else
    # Alternativa sin rsync.
    cp -R "$ORIGEN"/. "$DESTINO"/
    find "$DESTINO" -name '__pycache__' -type d -prune -exec rm -rf {} + 2>/dev/null || true
fi

echo "Instalado. Extensiones .inx en la carpeta de Inkscape:"
for f in superficie_cartesiana lente_ovoide lente_aplanatica; do
    if [[ -f "$DESTINO/$f.inx" ]]; then
        echo "  ✓ $f.inx"
    else
        echo "  ✗ $f.inx  (¡no se copió!)"
    fi
done

echo
echo "Listo. Ahora CIERRA Inkscape por completo y vuelve a abrirlo."
echo "Las verás en:  Extensions → Optics"
