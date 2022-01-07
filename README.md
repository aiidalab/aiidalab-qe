# AiiDAlab Quantum ESPRESSO application

## About

This is a early-development implementation of an AiiDAlab application for Quantum ESPRESSO workflow.
The app allows the execution of a workflow with Quantum ESPRESSO that includes the selection of an input structure, its relaxation, and the bands structure calculation.

**The app is currently in an early development stage!**

## For developers

The package uses pre-commit hooks to check the style consistency of all commits.
To use those you need to first install the pre-commit package itself, e.g. with:
```
pip install .[dev]
```
and then install the pre-commit hooks with
```
pre-commit install
```
The pre-commit checks should now be automatically executed prior to each commit.

## For maintainers

To create a new release, clone the repository, install development dependencies with `pip install -e '.[dev]'`, and then execute `bumpver update`.
This will:

  1. Create a tagged release with bumped version and push it to the repository.
  2. Trigger a GitHub actions workflow that creates a GitHub release.

Additional notes:

  - Use the `--dry` option to preview the release change.
  - The release tag (e.g. a/b/rc) is determined from the last release.
    Use the `--tag` option to switch the release tag.

In case that the [`aiidalab_qe_workchain`](src/) is changed, the corresponding dependency entry in [setup.cfg](setup.cfg) must also be changed.
Unfortunately this is a bit of a chicken and egg problem: the work chain package must be released first (e.g. via a release candidate), and the dependency is then updated for the next proper app release.

## Acknowledgements

This work is supported by the
[MARVEL National Centre for Competency in Research](<http://nccr-marvel.ch>) funded by the [Swiss National Science Foundation](<http://www.snf.ch/en>),
the MARKETPLACE project funded by [Horizon 2020](https://ec.europa.eu/programmes/horizon2020/) under the H2020-NMBP-25-2017 call (Grant No. 760173),
as well as by the [MaX
European Centre of Excellence](<http://www.max-centre.eu/>) funded by the Horizon 2020 EINFRA-5 program,
Grant No. 676598.

<div style="text-align:center">
 <img src="miscellaneous/logos/MARVEL.png" alt="MARVEL" height="75px">
 <img src="miscellaneous/logos/MaX.png" alt="MaX" height="75px">
 <img src="miscellaneous/logos/MarketPlace.png" alt="MarketPlace" height="75px">
</div>
