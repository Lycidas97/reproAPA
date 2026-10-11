# Third-party sources

## Code

- Module 03 keeps the captured pyclustree implementation and its license.
- The upstream scripts keep their original notices.

## APADETECT

| Version | Location | Core SHA-256 | Use |
| --- | --- | --- | --- |
| 0.3.0 (legacy) | `code/vendor/apadetect_annotation/postprocessing/` | `3e19278cb018644add4707da57380f69fbe0867c8f3977431182b78498180615` | Comparison with the legacy annotation |
| 0.4.1 (current) | External; not in this repository | `868c8cc8f6fac250605bc40d5150b8bc6bdb46635a3ec9b025440f9df654c223` | Current common annotation |

Version 0.4.1 uses a different dispatch for batches and bedtools. The [evidence record](../environments/evidence/pas-runtime-and-structural-cache.json) gives the shared helper semantics, the constants and the actual cache identities.

Limits:

- The full installed 0.4.1 package was not exported.
- The captured local copy of APADETECT has no license file.

## Resources and license

The authors' code, configurations and documentation are released under the MIT License (`LICENSE`). The third-party components above keep their own licenses and notices, and the scientific reference resources keep their own terms of use.
