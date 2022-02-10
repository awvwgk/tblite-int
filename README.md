# Integral calculator for xTB Hamiltonian

[![Apache-2.0](https://img.shields.io/github/license/awvwgk/tblite-int)](LICENSE)


## Building with fpm

Invoke [fpm](https://fpm.fortran-lang.org) in the project root with

```
fpm build
```

You can access the ``tblite-int`` program using the run subcommand

```
fpm run -- --help
```

To calculate the overlap integrals for a pair of two hyrogen atoms at 0.7 Ångström distance

```
fpm run -- H H 0,0,0.7
```


## License

Licensed under the Apache License, Version 2.0 (the “License”);
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an *“as is” basis*,
*without warranties or conditions of any kind*, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in this project by you, as defined in the
Apache-2.0 license, shall be licensed as above, without any additional
terms or conditions.
