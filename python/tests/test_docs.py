"""The Fortran and C references cover every public name in src/."""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FORTRAN_DOCS = ROOT / "docs" / "reference" / "fortran"
C_DOCS = ROOT / "docs" / "reference" / "c_api.rst"
CAPI = {"statespace_capi", "statespace_capi_models", "statespace_capi_extras"}


def _join(text):
    """Fortran source with continuation lines joined."""
    return re.sub(r"&[ \t]*\n\s*&?", " ", text)


def public_names(text):
    text = _join(text)
    names = set()
    for m in re.finditer(r"^\s*public\s*::\s*(.*)$", text, re.M):
        names |= {n.strip() for n in m.group(1).split("!")[0].split(",") if n.strip()}
    for m in re.finditer(
        r"^\s*(?:integer|real\(dp\)|logical)(?:\(c_int\))?,\s*parameter,\s*"
        r"public\s*::\s*(.*)$",
        text,
        re.M,
    ):
        names |= {
            p.split("=")[0].strip()
            for p in m.group(1).split("!")[0].split(",")
            if "=" in p
        }
    for m in re.finditer(r"^\s*type\s*,[^:]*public[^:]*::\s*(\w+)", text, re.M):
        names.add(m.group(1))
    return names


def test_fortran_reference_is_complete():
    missing = []
    for src in sorted((ROOT / "src").glob("*.f90")):
        mod = src.stem
        if mod in CAPI or mod == "statespace":
            continue
        page = FORTRAN_DOCS / f"{mod}.rst"
        assert page.exists(), f"no reference page for {mod}"
        doc = page.read_text()
        missing += [
            f"{mod}.{n}"
            for n in public_names(src.read_text())
            if not re.search(r"\b" + re.escape(n) + r"\b", doc)
        ]
    assert not missing, f"undocumented: {missing}"


def test_c_reference_is_complete():
    doc = C_DOCS.read_text()
    missing = []
    for mod in CAPI:
        text = (ROOT / "src" / f"{mod}.f90").read_text()
        for name in re.findall(r'bind\(C, name="(ss_\w+)"\)', text):
            if not re.search(r"\b" + name + r"\b", doc):
                missing.append(name)
    assert not missing, f"undocumented: {missing}"


def _signatures(text):
    """{name: argument list} for subroutine and function statements."""
    text = _join(text)
    sigs = {}
    for m in re.finditer(
        r"\b(?:subroutine|function)\s+(\w+)\s*\(([^)]*)\)", text, re.I
    ):
        sigs.setdefault(m.group(1).lower(), []).append(
            [a.strip().lower() for a in m.group(2).split(",") if a.strip()]
        )
    return sigs


def test_fortran_signatures_match_source():
    source = {}
    for src in (ROOT / "src").glob("*.f90"):
        for name, arglists in _signatures(src.read_text()).items():
            source.setdefault(name, []).extend(arglists)
    wrong = []
    for page in FORTRAN_DOCS.glob("*.rst"):
        for name, arglists in _signatures(page.read_text()).items():
            # A type-bound name may be implemented as <prefix>_<name>
            # (e.g. loglike => model_loglike), and a deferred one by its
            # abstract interface <name>_iface.
            known = [
                a
                for n, lists in source.items()
                if n in (name, name + "_iface") or n.endswith("_" + name)
                for a in lists
            ]
            for args in arglists:
                if known and args not in known:
                    wrong.append(f"{page.stem}: {name}({', '.join(args)})")
    assert not wrong, f"signatures differ from src/: {wrong}"
