# README

## About

The `src/` directory contains the `aiidalab_qe_workchain` package that must be installed as a dependency for the QE app as declared in the [setup.cfg](setup.cfg) file.

The app distributes its own AiiDA workchain implementation which constitutes an app dependency and must therefore be globally installed so that the AiiDA daemon is able to import it.
The workchain package is automatically built as a wheel and attached to every app release via the [GitHub actions `Release` workflow](https://github.com/aiidalab/aiidalab-qe/blob/master/.github/workflows/release.yml).
In this way one can install the wheel directly from the GitHub release via a PEP 508 compliant URL, example: `aiidalab-qe-workchain @ https://github.com/aiidalab/aiidalab-qe/releases/download/v22.01.0/aiidalab_qe_workchain-1.0-py3-none-any.whl`

## How to update the `aiidalab_qe_workchain` package

Any updates to the workchain that have not been released yet must be installed manually.

```console
cd src
pip install  .
```

Additional notes:
- Consider to use the editable mode (`pip install -e .`) while actively developing the workchain.
- CI tests will likely fail for a branch with not yet released workchain updates.
- A prerelease should be created shortly after merging workchain updates into the default branch to simplify development.

## Note on alternative approaches for distributing the workchain package

The following alternatives approaches for the distribution of the workchain wheel could be considered (in rough order of preference at the time of writing):

1. Install the package directly from the app directory (something like: `aiidalab-qe-workchain@file://./src/dist/aiidalab_qe_workchain-1.0-py3-none-any.whl`).
   However this is currently not possible, because it would be difficult to reliably determine the absolute location of the package and non-local URIs are not universally supported
2. Distribute the wheel [via GitHub packages](https://github.com/orgs/aiidalab/packages) once [Python packages are supported](https://github.com/github/roadmap/issues/94).
3. Distribute the workchain as part of its own dedicated package or as part of the [aiida-quantumespresso](https://github.com/aiidateam/aiida-quantumespresso) package.
4. Modify aiidalab to support the installation of additional packages not declared as part of [setup.cfg](setup.cfg).
