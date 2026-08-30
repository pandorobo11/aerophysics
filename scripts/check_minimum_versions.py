"""Confirm the public direct dependencies resolve to their declared minima."""

from importlib.metadata import version

EXPECTED_MINIMUMS = {
    "numpy": "2.0.0",
    "plotly": "6.9.0",
    "scipy": "1.14.0",
    "streamlit": "1.59.0",
}


def main() -> None:
    """Fail when lowest-direct resolution no longer selects a declared minimum."""
    mismatches = {
        package: (expected, version(package))
        for package, expected in EXPECTED_MINIMUMS.items()
        if version(package) != expected
    }
    if mismatches:
        detail = ", ".join(
            f"{package}: expected {expected}, resolved {actual}"
            for package, (expected, actual) in sorted(mismatches.items())
        )
        raise SystemExit(f"direct dependency minima were not selected: {detail}")

    for package, expected in sorted(EXPECTED_MINIMUMS.items()):
        print(f"{package}=={expected}")


if __name__ == "__main__":
    main()
