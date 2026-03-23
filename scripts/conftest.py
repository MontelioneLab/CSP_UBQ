"""
Pytest configuration and fixtures for CSP pipeline tests.
"""

import pytest


def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line(
        "markers", "integration: marks tests as integration (deselect with '-m \"not integration\"')"
    )
