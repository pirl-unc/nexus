# Contributing to Nexus

We would love your help in making `Nexus` a useful resource for the community. 
No contribution is too small, and we especially appreciate usability improvements 
like better documentation, tutorials, tests, or code cleanup.

## Project scope
We hope `Nexus` will grow to include **state-of-the-art bioinformatics pipelines 
for cancer immunogenomics**. This includes DNA/RNA alignment, *de novo* assembly, 
DNA/RNA variant calling, quantification of *de novo* (variant) and canonical 
isoforms as well as peptides, prediction of immunogenic peptides, and any useful downstream 
analysis. The pipelines implemented in `Nexus` support both long-read and short-read 
sequencing data to increase the sensitivity of immunologically targetable alterations 
in cancer.

All committed code to `Nexus` should be suitable for regular research use by practioners.

If you are contemplating a large contribution, such as the addition of a new workflow, it probably makes sense to reach out on the Github issue tracker (or email us at ajslee@unc.edu) to discuss and coordinate the work.

## Making a contribution
All contributions can be made as pull requests on Github. One of the core developers will review your contribution. As needed the core contributors will also make releases and submit to PyPI.

A few other guidelines:

 * `Nexus` is written for Python3 on Linux and OS X. We can't guarantee support for Windows.
 * All workflows should be documented using [numpy-style docstrings](https://numpydoc.readthedocs.io/en/latest/format.html) and associated with unit tests.
 * Bugfixes should be accompanied with test that illustrates the bug when feasible.
 * Contributions are licensed under Apache 2.0
 * Please adhere to our [code of conduct](https://github.com/pirl-unc/nexus/code-of-conduct.md).
