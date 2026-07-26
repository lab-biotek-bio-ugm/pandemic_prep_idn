"""The QC decision registry must not drift from the code it documents.

qc/decision_criteria.yaml is the cited source of truth the QC dossier renders.
If a threshold there disagrees with the constant the pipeline actually applies,
the dossier would justify a gate the code does not enforce -- an audit lie, which
for a national-surveillance assay is the worst failure mode. This test pins every
numeric value in the registry to the code constant it names, and checks that every
citation resolves to a real bib entry.

Only pyyaml plus the two pure-stdlib tool modules are needed, so this runs in the
contract job. iupac_utils imports primer3 lazily (inside functions), so importing
it at module scope here is safe.
"""
import re

import pytest
import yaml

import candidate_selection
import iupac_utils
from conftest import REPO

pytestmark = pytest.mark.contract

REGISTRY = yaml.safe_load(
    (REPO / "qc" / "decision_criteria.yaml").read_text(encoding="utf-8")
)
GATES = REGISTRY["gates"]
_MODULES = {"iupac_utils": iupac_utils, "candidate_selection": candidate_selection}
_BIBKEYS = set(re.findall(r"@\w+\{([^,]+),",
                          (REPO / "report" / "references.bib").read_text(encoding="utf-8")))


def _resolve(code_constant):
    """Resolve a code_constant reference to its live value.

    Forms handled:
      module.CONST
      module.LISTCONST[*].field   -> [item[field] for item in LISTCONST]
    """
    m = re.fullmatch(r"(\w+)\.(\w+)\[\*\]\.(\w+)", code_constant)
    if m:
        mod, listname, field = m.groups()
        return [item.get(field) for item in getattr(_MODULES[mod], listname)]
    mod, name = code_constant.split(".", 1)
    return getattr(_MODULES[mod], name)


_WITH_CONSTANT = [g for g in GATES if g.get("code_constant") and "value" in g]


@pytest.mark.parametrize("gate", _WITH_CONSTANT, ids=lambda g: g["id"])
def test_registry_value_matches_the_code_constant(gate):
    actual = _resolve(gate["code_constant"])
    expected = gate["value"]
    assert actual == expected, (
        f"gate {gate['id']}: qc/decision_criteria.yaml says {expected!r} but "
        f"{gate['code_constant']} is {actual!r}. Update the registry (and its "
        f"rationale/citation if the science changed), or revert the code change."
    )


@pytest.mark.parametrize("gate", GATES, ids=lambda g: g["id"])
def test_every_citation_resolves_to_a_bib_entry(gate):
    missing = [k for k in (gate.get("references") or []) if k not in _BIBKEYS]
    assert not missing, (
        f"gate {gate['id']} cites {missing}, which are not in report/references.bib. "
        f"Add the entry (verify the DOI first)."
    )


def test_no_unresolved_pending_references():
    pending = (REGISTRY.get("meta") or {}).get("pending_references") or []
    unresolved = [k for k in pending if k not in _BIBKEYS]
    assert not unresolved, (
        f"meta.pending_references still lists {unresolved}; add them to "
        f"references.bib and clear the list."
    )


def test_every_gate_has_rationale_and_class():
    for g in GATES:
        assert g.get("class") in (1, 2), f"gate {g.get('id')} has no valid class"
        assert g.get("rationale", "").strip(), f"gate {g.get('id')} has no rationale"
        assert g.get("references"), f"gate {g.get('id')} has no citation"
