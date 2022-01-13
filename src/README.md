# README

## About

The `src/` directory contains the `aiidalab_qe_workchain` package that must be installed as a dependency for the QE app as declared in the [setup.cfg](setup.cfg) file.

The app distributes its own AiiDA workchain implementation which constitutes an app dependency and must therefore be globally installed so that the AiiDA daemon is able to import it.
The workchain package is automatically built as a wheel and attached to every app release via the [GitHub actions `Release` workflow](https://github.com/aiidalab/aiidalab-qe/blob/master/.github/workflows/release.yml).
In this way one can install the wheel directly from the GitHub release via a PEP 508 compliant URL, example: `aiidalab-qe-workchain @ https://github.com/aiidalab/aiidalab-qe/releases/download/v22.01.0/aiidalab_qe_workchain-1.0-py3-none-any.whl`

## How to update the `aiidalab_qe_workchain` package

If changes to the workchain package cannot be avoided, this is how you can update it:

1. Develop any fixes or features in conjunction with necessary updates to the workchain package on a dedicated branch.*
2. Once development is finished merge the branch.
3. Create a _prerelease_ (e.g. 22.01.0rc0) that contains the update including the revised workchain with: `bumpver --tag rc`.
4. Update the workchain dependency in [setup.cfg](setup.cfg) to install the revised workchain from the prerelease.
5. Create the final release (e.g. 22.01.0) with: `bumpver --tag final`.
6. Update the workchain dependency in [setup.cfg](setup.cfg) to install the revised workchain from the latest final release.†

_*) This will require a manual installation of the package during the development process and CI tests might fail._

_†) This step is technically optional and will only take effect on the next release._

## Note on alternative approaches for distributing the workchain package

The rather involved approach for updating the package is considered sub-optimal, however alternative (simpler) approaches for distributing and installing the wheel are currently not available or not known.

The following alternatives approaches could be considered (in rough order of preference at the time of writing):

1. Install the package directly from the app directory (something like: `aiidalab-qe-workchain@file://./src/dist/aiidalab_qe_workchain-1.0-py3-none-any.whl`).
   However this is currently not possible, because it would be difficult to reliably determine the absolute location of the package and non-local URIs are not universally supported
2. Distribute the wheel [via GitHub packages](https://github.com/orgs/aiidalab/packages) once [Python packages are supported](https://github.com/github/roadmap/issues/94).
3. Distribute the workchain as part of its own dedicated package or as part of the [aiida-quantumespresso](https://github.com/aiidateam/aiida-quantumespresso) package.
4. Modify aiidalab to support the installation of additional packages not declared as part of [setup.cfg](setup.cfg).
