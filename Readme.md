# ptope

Library to handle computations using polytopes and their associated Gram
matrices.

### Build

The build system used is `make`. Running `make` will compile the dynamic library
`libptope.so`. `make test` will compile and run the test suite.

### Dependencies

`ptope` uses:

* [armadillo linear algebra library][arma] to handle matrix
operations and decompositions.
* [GoogleTest] to provide unit tests.

### License

`ptope` is provided under the Apache 2.0 license. See License for more details.

```
   Copyright 2015 John Lawson

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
```

[arma]: http://arma.sourceforge.net
[GoogleTest]: https://github.com/google/googletest

