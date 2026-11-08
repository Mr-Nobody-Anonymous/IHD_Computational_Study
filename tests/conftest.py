"""Shared fixtures. Session-scoped registries keep the suite fast."""

from __future__ import annotations

import numpy as np
import pytest

from ihd_solutions.catalog import Regimen, load_catalog, validate_combination
from ihd_solutions.evidence import load_anchors
from ihd_solutions.population import background_regimen_ids, load_populations
from ihd_solutions.units import time_grid


@pytest.fixture(scope="session")
def catalog():
    return load_catalog()


@pytest.fixture(scope="session")
def anchors():
    return load_anchors()


@pytest.fixture(scope="session")
def populations():
    return load_populations()


@pytest.fixture(scope="session")
def population(populations):
    return populations["stable_ccs"]


@pytest.fixture(scope="session")
def days() -> np.ndarray:
    return time_grid(5.0, 30.0)


@pytest.fixture(scope="session")
def short_days() -> np.ndarray:
    return time_grid(2.0, 60.0)


@pytest.fixture(scope="session")
def background(catalog) -> Regimen:
    return validate_combination(catalog, background_regimen_ids())


@pytest.fixture(scope="session")
def make_regimen(catalog):
    def _make(*therapy_ids: str) -> Regimen:
        return validate_combination(catalog, therapy_ids)

    return _make


@pytest.fixture(scope="session")
def with_background(catalog):
    def _make(*therapy_ids: str) -> Regimen:
        return validate_combination(catalog, background_regimen_ids() + therapy_ids)

    return _make
