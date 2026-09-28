# L-BFGS-B

The Fortran sources of [jacobwilliams/lbfgsb](https://github.com/jacobwilliams/lbfgsb)
at commit `ce1a9322016d6f31b652393e000a914b63af6fef`, unmodified, with their
`License.txt` (BSD-3-Clause). A modernization of L-BFGS-B 3.0 by Zhu, Byrd, Lu,
Nocedal and Morales; see `docs/references.md` for the papers.

They are copied here so that the library builds without network access (for
example when pip builds it from the source distribution). Both builds use this
copy: `fpm.toml` as a path dependency and `CMakeLists.txt` directly. To update,
replace `src/` with the sources of a newer commit and record it above.
