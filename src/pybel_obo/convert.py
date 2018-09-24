# -*- coding: utf-8 -*-

"""Conversion utilities for BEL."""

from pybel import BELGraph
from pybel.dsl import BaseEntity


def convert(graph: BELGraph, u: BaseEntity, v: BaseEntity, key: str):
    """Convert a graph."""
    raise NotImplementedError
