# -*- coding: utf-8 -*-

"""Conversion utilities for BEL."""

from typing import Mapping

from pybel import BELGraph
from pybel.dsl import BaseEntity


def get_relationship(graph: BELGraph, u: BaseEntity, v: BaseEntity, key: str) -> Mapping[str, str]:
    """Convert a graph."""
    raise NotImplementedError
