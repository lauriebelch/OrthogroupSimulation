#!/usr/bin/env python3
'''
"Script to run the unit testing
'''

import glob
import importlib
import os
import sys
import traceback

# make sure path to orthosim repo is stored and accessible
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

# folder
TESTS_DIR = os.path.join(REPO_ROOT, "tests")


## function to find the individual unit test scripts, which are in format test_*.py
def discover_test_modules():
    paths = sorted(glob.glob(os.path.join(TESTS_DIR, "test_*.py")))
    module_names = []
    for path in paths:
        file_name = os.path.basename(path)
        module_names.append("tests." + file_name[:-3])
    return module_names

## function to run a test module
def run_module(module_name, results):
    # first check that its all importable and runnable
    try:
        module = importlib.import_module(module_name)
    except Exception as exc:
        print(f"\n{module_name}: IMPORT ERROR -> {exc!r}")
        traceback.print_exc()
        results["import_errors"].append(module_name)
        return
    # then find some test functions to run
    # all these must begin with 'test_'
    test_function_names = sorted(
        name for name in dir(module)
        if name.startswith("test_") and callable(getattr(module, name))
    )
    if not test_function_names:
        print(f"\n{module_name}: no test_* functions found")
        return
    # this bit will run the tests, and capture if the test passes or fails
    print(f"\n{module_name}:")
    for name in test_function_names:
        fn = getattr(module, name)
        try:
            fn()
            print(f"  PASS  {name}")
            results["passed"].append(f"{module_name}.{name}")
        except Exception as exc:
            print(f"  FAIL  {name} -> {exc!r}")
            results["failed"].append((f"{module_name}.{name}", exc))

# main will do the tests and print the results
def main():
    results = {"passed": [], "failed": [], "import_errors": []}
    for module_name in discover_test_modules():
        run_module(module_name, results)
    total = len(results["passed"]) + len(results["failed"])
    print("\n" + "=" * 60)
    print(f"Modules with import errors: {len(results['import_errors'])}")
    for module_name in results["import_errors"]:
        print(f"  - {module_name}")
    print(f"Tests run:     {total}")
    print(f"Tests passed:  {len(results['passed'])}")
    print(f"Tests failed:  {len(results['failed'])}")
    for test_name, exc in results["failed"]:
        print(f"  - {test_name}: {exc!r}")
    print("=" * 60)
    if results["failed"] or results["import_errors"]:
        sys.exit(1)
    sys.exit(0)


if __name__ == "__main__":
    main()
