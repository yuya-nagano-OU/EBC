# EBC — Effective Beam Convolver for CMB Science

Julia tools for effective beam convolution on the sphere, with a focus on CMB temperature/polarization analysis.

This repository is an early-stage research codebase for handling:
- spherical-harmonic coefficients of sky and beam maps
- spin-weighted combinations for polarization
- Wigner-D based rotations
- effective/local beam convolution on HEALPix pixels
- map-making related local rotation handling

## Overview

EBC is designed for CMB beam convolution studies in harmonic space.

The current implementation includes:
- structures for sky/beam/convolution settings
- slicing of harmonic coefficients by multipole range
- conversion between polarization-related harmonic components
- global and local effective Wigner-D based rotation utilities
- one-pixel convolution routines
- HEALPix ring/pixel helper functions

This repository is currently closer to a **research prototype** than a polished public package.  
The API may change as the implementation is cleaned up and extended.

