# Paper Manifest

Licensed PDFs and downloaded paper images in this directory are intentionally
ignored by git. Keep bibliographic metadata, access notes, and checksum records
in this manifest; do not commit publisher PDFs.

Use checksums to identify local copies, not to imply redistribution rights:

```bash
sha256sum docs/papers/<file>
```

## Required Source Inventory

| Local status | Expected local filename | Citation | DOI / URL | Access notes | SHA-256 |
| --- | --- | --- | --- | --- | --- |
| Missing as of 2026-06-01 | `DecentKing2008_IMAJAM_73_1_37-68_hxm043.pdf` | S. P. Decent and A. C. King, "Surface-tension-driven flow in a slender cone," *IMA Journal of Applied Mathematics* **73**(1), 37--68, 2008. | DOI [`10.1093/imamat/hxm043`](https://doi.org/10.1093/imamat/hxm043) | Primary target for source-fidelity work. Published online in 2007 and in the February 2008 issue. Requires publisher/institutional access unless open via OUP. | TODO once acquired locally |
| Present locally; ignored | `Billingham1999_JFM_397_45.pdf` | J. Billingham, "Surface-tension-driven flow in fat fluid wedges and cones," *Journal of Fluid Mechanics* **397**, 45--71, 1999. | DOI [`10.1017/S0022112099006011`](https://doi.org/10.1017/S0022112099006011) | Local copy should be treated as licensed publisher content. | `a4c99bd1dff44abcccd2c8cbce1f6bfa0cf6817ffbf53cb62b607838f61e88ef` |
| Present locally; ignored | `DecentKing2001_IUTAM.pdf` | S. P. Decent and A. C. King, "The Recoil of A Broken Liquid Bridge," in *IUTAM Symposium on Free Surface Flows*, Springer, 2001. | DOI [`10.1007/978-94-010-0796-2_4`](https://doi.org/10.1007/978-94-010-0796-2_4) | Local copy is a large Springer proceedings PDF; verify page/chapter extraction before using for source-backed tests. | `efd95dddbb6571ce324a3115a2812df37c3265d296d4818f96d8c804fa7ab524` |
| Present locally; ignored | `Eggers1997_RMP_69_865.pdf` | J. Eggers, "Nonlinear dynamics and breakup of free-surface flows," *Reviews of Modern Physics* **69**, 865--929, 1997. | DOI [`10.1103/RevModPhys.69.865`](https://doi.org/10.1103/RevModPhys.69.865) | Local APS review copy; use for contextual literature, not Decent--King benchmarks. | `7e5afdd865742d778c7955b359ae53e3a4dd48b77820cc4c91ca774d94349685` |
| Present locally; ignored | `KellerMiksis1983_SIAMJAM_43_268.pdf` | J. B. Keller and M. J. Miksis, "Surface tension driven flows," *SIAM Journal on Applied Mathematics* **43**(2), 268--277, 1983. | DOI [`10.1137/0143018`](https://doi.org/10.1137/0143018) | Local SIAM copy; source for the `t^(2/3)` similarity scaling. | `2942ea6bdddd6f8a2bb2c5735363f4208de1850f3e111b329f72e9cf12f89a21` |
| PDF missing; page image present | `KellerKingTing1995_PoF_7_226.pdf` | J. B. Keller, A. C. King, and L. Ting, "Blob formation," *Physics of Fluids* **7**(1), 226--228, 1995. | DOI [`10.1063/1.868513`](https://doi.org/10.1063/1.868513) | Only `KellerKingTing1995_page.png` is currently present. Acquire and checksum the PDF before citing quantitative details. | TODO once acquired locally |

## Quarantined Mislabeled Artifact

The file previously named `DecentKing2008_QJMAM_61_1.pdf` was not the
Decent--King cone paper. It was renamed locally with:

```bash
mv -f docs/papers/DecentKing2008_QJMAM_61_1.pdf docs/papers/NOT_DecentKing_PrestonJensenRichardson2008_QJMAM_61_1_hbm021.pdf
```

| Local filename | Actual metadata | DOI / URL | Why quarantined | SHA-256 |
| --- | --- | --- | --- | --- |
| `NOT_DecentKing_PrestonJensenRichardson2008_QJMAM_61_1_hbm021.pdf` | S. P. Preston, O. E. Jensen, and G. Richardson, "Buckling of an axisymmetric vesicle under compression: the effects of resistance to shear," *The Quarterly Journal of Mechanics and Applied Mathematics* **61**(1), 1--24, 2008. | DOI [`10.1093/qjmam/hbm021`](https://doi.org/10.1093/qjmam/hbm021) | It is a vesicle-compression paper, not a slender-cone recoil paper. Keep it quarantined or remove the local copy. | `8ebdf98aa7593bd04a16d078e7aaf94477451474f34730f222d2fc71792fa52e` |

The previously cited DOI `10.1093/qjmam/hbm028` is also not the cone paper.
Crossref/OUP metadata identify it as M. H. B. M. Shariff, "Nonlinear
transversely isotropic elastic solids: an alternative representation," *The
Quarterly Journal of Mechanics and Applied Mathematics* **61**(2), 129--149,
2008.

## Page Images And Fetch Diagnostics

These ignored PNGs are local acquisition aids rather than source artifacts.
Do not cite them as references.

| Local file | Status | SHA-256 |
| --- | --- | --- |
| `KellerKingTing1995_page.png` | Single page image for the missing AIP paper. | `07a8a63bd9be4e3baf063bb8e903ea472e973a9b28141a81eb50fe1f6e8de7a8` |
| `KellerMiksis1983_directfail.png` | Fetch diagnostic screenshot. | `4a2cb1e2b1b99c89b9f07e400e9d24cfbc1c7f147ee8f1c72659402cb67f33e4` |
| `KellerMiksis1983_page.png` | Fetch diagnostic screenshot. | `51857e1bbe9b7cd097f4b701f55bff89f69dd7432476329dcf493765bbc5d1aa` |
| `doi_test.png` | Fetch diagnostic screenshot. | `4862ce112e959576907def262b634e770223438cf9a7be7feccdc4be8ba821f0` |
| `scholar_test.png` | Fetch diagnostic screenshot. | `495255093c53ecd3b13bef494c2c8dcd7f97f771ca52800e23fc9df01564e72c` |
| `tib_search.png` | Fetch diagnostic screenshot. | `91d02e0af16a459116ff5c9584c2a46a2e6337d46242a0dea058d7051fe8060b` |

## Fetching Policy

Paper fetching is best effort because publisher flows, VPN entitlements, and
Cloudflare challenges change outside this repository. The Node/Playwright tool
chain is locked by `package-lock.json`, which should be tracked; use `npm ci`
before running fetch tooling so browser automation versions are reproducible.
Downloaded PDFs, browser profiles, and diagnostic images remain ignored by git.
