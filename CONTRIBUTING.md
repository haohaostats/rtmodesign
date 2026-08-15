# Contributing to rtmodesign

Bug reports and focused feature requests are welcome through GitHub Issues.
For code contributions, please open an issue before a large change, add or
update tests, and confirm that `R CMD check --no-manual` completes with
`Status: OK`.

The public interface is intentionally limited to ordinary-user workflows.
New exported arguments should be proposed only when they represent a clinical
input that cannot be inferred safely by the package.

