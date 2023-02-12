---
layout: default
nav_order: 2
title: About
---

## About *gpaw-tools*

`gpaw-tools` is a tool developed by the [LRG](https://www.lrgresearch.org) research group at Gazi University for conducting Density Functional Theory (DFT) calculations. The LRG group conducts both experimental and computational studies, and prior to 2020, had been using the DFT tool ATK (later named QuantumATK) for their research on low-dimensional structures.

However, due to licensing issues with the Europractice program at Gazi University, the LRG group began searching for a new DFT tool at the beginning of 2020. After careful consideration, they decided to use ASE and GPAW because they could be integrated and were flexible to work with Python. ASE is a powerful software and the ability to use GPAW and various modes such as PW, FD, LCAO, and Exact-exchange, provides significant flexibility. Additionally, by using ASE+ASAP+OpenKIM interatomic potentials for geometry pre-optimization, it is very similar to the desired features of their previous workflows.

One of the important considerations for the LRG group when choosing a DFT tool was ease of use. While ASE and GPAW are powerful, they do require knowledge of Python. Therefore, the `gpaw-tools` project was started to create a user interface for ASE and GPAW, making it more accessible for users who are not proficient in Python but still want to conduct materials research.

The `gpaw-tools` project is still in its early stages and the LRG group hopes to continue to improve and expand the capabilities of the tool. They welcome feedback and contributions from the community, and encourage users to use/extend the software and create issues and requests on the [Github page](https://github.com/lrgresearch/gpaw-tools).

Thank you for choosing `gpaw-tools` for your DFT and/or MD calculations!

[back](./)
