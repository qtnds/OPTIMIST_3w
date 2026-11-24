"""Lecture des arguments de type liste et des longueurs avec unité."""

from __future__ import annotations

import re

_LENGTH_UNITS = {
    "m": 1.0,
    "mm": 1e-3,
    "um": 1e-6,
    "µm": 1e-6,
    "μm": 1e-6,
    "nm": 1e-9,
    "pm": 1e-12,
}

_NUM_RE = re.compile(r"^\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)\s*([a-zA-Zµμ]*)\s*$")


def parse_list(value) -> list[str] | None:
    """'[a, b,c]' ou 'a,b,c' -> ['a', 'b', 'c'].  None -> None.

    Une liste déjà construite (lue dans un fichier JSON) est renvoyée en chaînes.
    """
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        return [str(v).strip() for v in value]
    s = str(value).strip()
    if s.startswith("[") and s.endswith("]"):
        s = s[1:-1]
    if not s:
        return []
    return [item.strip() for item in s.split(",")]


def is_default(token: str | None) -> bool:
    """Un élément '-', '' ou 'none' signifie : garder la valeur de la base de matériaux."""
    return token is None or str(token).strip().lower() in {"", "-", "none", "null", "default"}


def parse_length(token, default_unit: str = "nm") -> float:
    """'50nm' -> 5e-8 ; '300um' -> 3e-4 ; '20' -> 20 nm (unité par défaut)."""
    if isinstance(token, (int, float)):
        return float(token) * _LENGTH_UNITS[default_unit]
    m = _NUM_RE.match(str(token))
    if not m:
        raise ValueError(f"Longueur illisible : {token!r} (exemples valides : 50nm, 0.3um, 300um, 1e-6m)")
    value, unit = float(m.group(1)), m.group(2) or default_unit
    if unit not in _LENGTH_UNITS:
        raise ValueError(f"Unité inconnue {unit!r} dans {token!r} (unités : {', '.join(_LENGTH_UNITS)})")
    return value * _LENGTH_UNITS[unit]


def parse_float_list(value, n: int, name: str) -> list[float | None]:
    """Liste de flottants de longueur n ; les éléments '-' deviennent None (valeur par défaut)."""
    items = parse_list(value)
    if items is None:
        return [None] * n
    if len(items) != n:
        raise ValueError(f"--{name} attend {n} valeurs (une par couche), {len(items)} reçues : {items}")
    return [None if is_default(t) else float(t) for t in items]


def fmt_length(x: float) -> str:
    """Affichage lisible d'une longueur en mètres."""
    ax = abs(x)
    if ax >= 1e-3:
        return f"{x * 1e3:.3g} mm"
    if ax >= 1e-6:
        return f"{x * 1e6:.3g} µm"
    return f"{x * 1e9:.3g} nm"
