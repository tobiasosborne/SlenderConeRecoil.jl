# Paper Manifest

Licensed PDFs and downloaded paper images in this directory are intentionally ignored by git. Keep bibliographic metadata in this manifest, and do not commit publisher PDFs.

## Primary Target

| Status | Expected local filename | Metadata | Notes |
| --- | --- | --- | --- |
| Not present locally as of 2026-06-01 | `DecentKing2008_IMAJAM_73_1_37-68_hxm043.pdf` | S. P. Decent and A. C. King, "Surface-tension-driven flow in a slender cone," *IMA Journal of Applied Mathematics* **73**(1), 37--68, 2008. DOI [`10.1093/imamat/hxm043`](https://doi.org/10.1093/imamat/hxm043). | Verified from DOI content negotiation, Oxford/OUP metadata, and the University of Dundee record. The article was published online in 2007 and appears in the February 2008 issue. |

## Quarantined Mislabeled Artifact

The file previously named `DecentKing2008_QJMAM_61_1.pdf` was not the Decent--King cone paper. It was renamed locally with:

```bash
mv -f docs/papers/DecentKing2008_QJMAM_61_1.pdf docs/papers/NOT_DecentKing_PrestonJensenRichardson2008_QJMAM_61_1_hbm021.pdf
```

| Local filename | Actual metadata | Why quarantined |
| --- | --- | --- |
| `NOT_DecentKing_PrestonJensenRichardson2008_QJMAM_61_1_hbm021.pdf` | S. P. Preston, O. E. Jensen, and G. Richardson, "Buckling of an axisymmetric vesicle under compression: the effects of resistance to shear," *The Quarterly Journal of Mechanics and Applied Mathematics* **61**(1), 1--24, 2008. DOI `10.1093/qjmam/hbm021`. | It is a vesicle-compression paper, not a slender-cone recoil paper. |

The previously cited DOI `10.1093/qjmam/hbm028` is also not the cone paper. Crossref/OUP metadata identify it as M. H. B. M. Shariff, "Nonlinear transversely isotropic elastic solids: an alternative representation," *The Quarterly Journal of Mechanics and Applied Mathematics* **61**(2), 129--149, 2008.

## Other Local Inventory

These ignored files were present locally during the 2026-06-01 source-fidelity cleanup. Their detailed metadata was not revalidated in this issue.

| Local file | Status |
| --- | --- |
| `Billingham1999_JFM_397_45.pdf` | Present locally; ignored by git. |
| `DecentKing2001_IUTAM.pdf` | Present locally; ignored by git. |
| `Eggers1997_RMP_69_865.pdf` | Present locally; ignored by git. |
| `KellerMiksis1983_SIAMJAM_43_268.pdf` | Present locally; ignored by git. |
| `KellerKingTing1995_page.png` | Page image only; full AIP paper not present locally. |
